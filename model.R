library(epiworldR)

#' Build a measles model with community and school layers
#'
#' Latent, prodromal and rash durations are rounded lognormal draws. The
#' incubation standard deviation is informed by Lessler et al. (2009), the
#' prodromal standard deviation by the Hagelloch data, and the rash standard
#' deviation is an assumption. Vaccination and contact assignments are given
#' in `setup`; the surveillance event records each run in `state$runs`.
#' @param p Disease and surveillance parameters.
#' @param setup Population, contacts, vaccination and shedding inputs.
#' @param state Environment that receives one record per simulation.
#' @param scenario Baseline, wastewater, or cases.
#' @return An epiworldR model ready for [epiworldR::run_multiple()].
#' @references Lessler et al. (2009), doi:10.1016/S1473-3099(09)70069-6.
#'   Meyer et al. (2017), doi:10.18637/jss.v077.i11.
#' @examples
#' parameters <- read.csv("data/parameters.csv", stringsAsFactors = FALSE)
#' p <- lapply(parameters$default, function(x) eval(parse(text = x)))
#' names(p) <- parameters$r_argument
#' p$n <- 100
#' p$ndays <- 30
#' setup <- list(
#'   persons = list(age_group = rep(0L, p$n), school = rep(0L, p$n)),
#'   community = matrix(2, 1, 1),
#'   school_probability = 0.05,
#'   vaccinated = integer(),
#'   kernel = rep(1 / 19, 19)
#' )
#' model <- build_model(p, setup, new.env(parent = emptyenv()))
#' epiworldR::get_stages(model)
build_model <- function(
    p,
    setup,
    state,
    scenario = "baseline"
) {
  population         <- setup$persons
  community_matrix   <- setup$community
  school_probability <- setup$school_probability
  vaccinated         <- setup$vaccinated
  stopifnot(
    length(population$age_group) == p$n,
    length(population$school) == p$n
  )

  # Define disease stages and their progression using the model builder.
  virus <- virus(name = "Measles") |>
    add_stage(label = "Susceptible") |>
    add_stage(label = "Latent", carrier = TRUE) |>
    add_stage(label = "Prodromal", infectious = TRUE) |>
    add_stage(label = "Rash", infectious = TRUE) |>
    add_stage(label = "Rash isolated", infectious = TRUE) |>
    add_stage(label = "Hospitalized", carrier = TRUE) |>
    add_stage(label = "Recovered") |>
    set_contact(stage = "Rash isolated", weight = 0) |>
    set_contact(stage = "Hospitalized", weight = 0) |>
    add_infection(
      from = "Susceptible",
      to   = "Latent",
      rate = p$transmission_rate
    ) |>
    add_after(
      from     = "Latent",
      to       = "Prodromal",
      duration = duration_lognormal(
        mean = p$incubation_period,
        sd   = p$incubation_period_sd
      )
    ) |>
    add_after(
      from     = "Prodromal",
      to       = "Rash",
      duration = duration_lognormal(
        mean = p$prodromal_period,
        sd   = p$prodromal_period_sd
      )
    )

  # Both rash stages share the recovery and hospitalization mechanisms.
  for (stage in c("Rash", "Rash isolated")) {
    virus <- virus |>
      add_after(
        from     = stage,
        to       = "Recovered",
        duration = duration_lognormal(
          mean = p$rash_period,
          sd   = p$rash_period_sd
        )
      )
  }
  virus <- virus |>
    add_external(
      from            = "Rash",
      to              = "Rash isolated",
      transition_name = "Isolate rash"
    ) |>
    add_rate(
      from = c("Rash", "Rash isolated"),
      to   = "Hospitalized",
      rate = rash_hospitalization_rate(
        probability = p$prob_hospitalization,
        mean        = p$rash_period,
        sd          = p$rash_period_sd
      )
    ) |>
    add_rate(
      from = "Hospitalized",
      to   = "Recovered",
      rate = 1 / p$hospitalization_period
    )

  # Add the two contact layers and seed infections.
  model <- Model(n = p$n) |>
    add_layer(
      layer = layer_mixing(
        name           = "community",
        groups         = as.integer(population$age_group),
        contact_matrix = community_matrix
      )
    ) |>
    add_layer(
      layer = layer_mixing(
        name         = "schools",
        groups       = as.integer(population$school),
        contact_prob = school_probability
      )
    ) |>
    add_virus(virus = virus, assign = assign_count(n = p$n_seed))

  # Apply MMR protection to the supplied agents.
  if (length(vaccinated)) {
    vaccine <- tool(
      name                     = "MMR",
      susceptibility_reduction = p$vax_efficacy,
      transmission_reduction   = 0,
      recovery_enhancer        = 0,
      death_reduction          = 0
    )
    model <- model |>
      add_tool(
        tool   = vaccine,
        assign = assign_agents(ids = as.integer(vaccinated))
      )
  }
  # Attach surveillance without capturing the model in the callback.
  model <- model |>
    add_globalevent(
      event = surveillance_event(
        p        = p,
        setup    = setup,
        scenario = scenario,
        state    = state
      )
    )
  model
}

#' Convert a target hospitalization probability to a daily rash hazard
#'
#' Assumes the rounded lognormal rash duration is independent of the daily
#' hospitalization draw. Numerical summation includes the far tail.
#' @param probability Probability of hospitalization during rash.
#' @param mean,sd Mean and standard deviation of rash duration in days.
#' @return Daily hospitalization probability.
#' @examples
#' rash_hospitalization_rate(probability = 0.2, mean = 4, sd = 1)
rash_hospitalization_rate <- function(probability, mean, sd) {
  stopifnot(probability >= 0, probability < 1, mean > 0, sd >= 0)

  # Handle degenerate distributions before numerical integration.
  if (probability == 0) return(0)
  if (sd == 0) return(1 - (1 - probability)^(1 / max(1, round(mean))))

  # Integrate the rounded lognormal duration to recover a daily hazard.
  s2 <- log1p((sd / mean)^2)
  days <- seq_len(1000L)
  mass <- diff(c(0, stats::plnorm(days + .5, log(mean) - s2 / 2, sqrt(s2))))
  stats::uniroot(
    f        = function(h) 1 - sum(mass * (1 - h)^days) - probability,
    interval = c(0, 1)
  )$root
}

#' Create the wastewater and case-surveillance event
#'
#' Sampling, reporting and interventions share one R event. The event records
#' daily series and trigger day separately for every run. R events require
#' `run_multiple()` to run on a single thread.
#' @param p Disease and surveillance parameters.
#' @param setup Contact and shedding inputs.
#' @param scenario Baseline, wastewater, or cases.
#' @param state Environment that receives the replicate records.
#' @return An epiworldR global event.
#' @examples
#' p <- list(ndays = 30, ww_log10_g_mean = 10, ww_log10_g_sd = 1)
#' state <- new.env(parent = emptyenv())
#' event <- surveillance_event(p, list(), "baseline", state)
#' class(event)
surveillance_event <- function(p, setup, scenario, state) {

  state$runs <- list()
  n_days <- p$ndays
  meanlog <- log(10) * p$ww_log10_g_mean
  sdlog   <- log(10) * p$ww_log10_g_sd

  reset <- function() {
    variables <- c("cohort_g", "cum_inf", "cum_prod", "cum_rash")
    for (v in variables)
      assign(v, numeric(n_days), envir = state)
    state$reported    <- numeric(n_days + p$case_report_delay)
    state$ww_obs      <- rep(NA_real_, n_days)
    state$trigger_day <- NA_integer_
    state$rash_seen <- rep(FALSE, p$n)
    # Day 0: seeded infections start in the latent state
    state$prev <- c(inf = p$n_seed, prod = 0, rash = 0)
  }

  fun <- function(m) {

    t <- epiworldR::today(m)
    if (t == 1) {
      reset()
      state$run <- length(state$runs) + 1L
    }

    # The final callback saves the replicate even when no trigger fires.
    on.exit(
      if (t == n_days) {
        state$runs[[state$run]] <- list(
          trigger_day = state$trigger_day,
          ww_obs      = state$ww_obs,
          reported    = state$reported[seq_len(n_days)],
          cum_inf     = state$cum_inf,
          cum_prod    = state$cum_prod,
          cum_rash    = state$cum_rash
        )
      }
    )

    # The model has one status, so stage totals describe the whole population.
    totals <- epiworldR::get_today_total(m)
    tot <- tapply(totals$counts, totals$stage, sum)

    # Scan individuals only on days when someone has a rash.
    if (tot[["Rash"]] > 0) {
      rash <- get_agents_stages(m) == "Rash"
      new_rash <- which(rash & !state$rash_seen)
      state$rash_seen[new_rash] <- TRUE
      selected <- new_rash[
        runif(length(new_rash)) < p$rash_reduction_contact_rate
      ]
      for (i in selected) fire(m, "Isolate rash", i - 1L)
    }

    cum <- c(inf = p$n - tot[["Susceptible"]])
    cum["prod"] <- cum[["inf"]] - tot[["Latent"]]
    cum["rash"] <- cum[["prod"]] - tot[["Prodromal"]]

    new <- cum - state$prev
    state$prev <- cum

    state$cum_inf[t]  <- cum[["inf"]]
    state$cum_prod[t] <- cum[["prod"]]
    state$cum_rash[t] <- cum[["rash"]]

    # Draw detection for each rash case and record it after the reporting delay.
    r <- t + p$case_report_delay
    state$reported[r] <- state$reported[r] +
      sum(stats::runif(new[["rash"]]) < p$case_report_prob)

    # Shedding: individual total genome copies are log-normal
    state$cohort_g[t] <- sum(stats::rlnorm(new[["prod"]], meanlog, sdlog))
    lag <- seq_len(min(t, length(setup$kernel))) - 1L
    conc <- sum(setup$kernel[lag + 1L] * state$cohort_g[t - lag]) /
      (p$ww_alpha * p$n)
    # Draw measurement error daily and use it on scheduled sampling days.
    eps <- stats::rnorm(1, 0, p$ww_sigma)
    if ((t %% 7) %in% p$ww_sampling_days) {
      obs <- conc * exp(eps)
      state$ww_obs[t] <- if (obs >= p$ww_lod) obs else 0
    }

    if (!is.na(state$trigger_day) || scenario == "baseline")
      return(invisible(NULL))

    fire <- if (scenario == "cases") {
      sum(state$reported[1:t]) >= p$trigger_cases
    } else {
      s <- t - p$ww_lab_delay
      s >= 1 && !is.na(state$ww_obs[s]) && state$ww_obs[s] >= p$trigger_ww
    }

    if (fire) {
      state$trigger_day <- t
      epiworldR::set_contact_prob(
        layer        = epiworldR::get_layer(m, "schools"),
        contact_prob = setup$school_probability *
          (1 - p$school_contact_reduction)
      )
    }

    invisible(NULL)

  }

  epiworldR::globalevent_fun(fun, name = "Surveillance", day = -99)
}


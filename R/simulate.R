# Simulation of measles outbreaks with wastewater (WW) or case-based triggers.
#
# Scenarios:
#  - "baseline":   no surveillance-triggered intervention.
#  - "wastewater": WW samples are collected on sampling days and become
#                  available `ww_lab_delay` days later. The intervention starts
#                  the first time a result is at or above `trigger_ww`.
#  - "cases":      each case is detected during its rash period with
#                  probability `case_report_prob` (default 70%) and reported
#                  `case_report_delay` days after rash onset. The intervention
#                  starts when cumulative reports reach `trigger_cases`.
#
# The intervention reduces contacts in the school layer for the rest of the
# simulation; community contacts remain unchanged.
#
# Scenarios run with the same seed share the same trajectory up to the
# trigger day (both the model and R's RNG are seeded identically, and the
# surveillance draws are the same in all scenarios), so differences between
# them are paired.

#' Default parameters
#'
#' All defaults can be overridden through `...`. See README.qmd for sources.
#' @param ... Named values replacing the defaults.
measles_ww_params <- function(...) {

  defaults <- list(
    # Population and transmission
    n                           = 1e4,
    n_seed                      = 1,
    R0                          = 7,
    prop_vaccinated             = 0,
    vax_efficacy                = 0.97,
    transmission_rate           = 0.9,
    incubation_period           = 10,
    incubation_period_sd        = 2.5,
    prodromal_period            = 4,
    prodromal_period_sd         = 1.7,
    rash_period                 = 4,
    rash_period_sd              = 1,
    rash_reduction_contact_rate = 0.5,
    prob_hospitalization        = 0.2,
    hospitalization_period      = 7,
    ndays                       = 365,
    # Wastewater shedding (days counted from prodrome onset)
    ww_t_peak                   = 4,
    ww_t_end                    = 18,
    ww_log10_vpeak              = 3,
    ww_log10_g_mean             = 10,
    ww_log10_g_sd               = 1,
    ww_alpha                    = 378500,
    ww_sigma                    = 0.5,
    ww_lod                      = 1,
    ww_sampling_days            = c(1, 3, 5),
    ww_lab_delay                = 3,
    # Clinical surveillance: probability that a case is detected (and
    # reported) during its rash period. Varied in a sensitivity analysis.
    case_report_prob            = 0.7,
    case_report_delay           = 3,
    # Triggers and intervention
    trigger_cases               = 5,
    trigger_ww                  = 1,
    school_contact_reduction    = 0.8
  )

  utils::modifyList(defaults, list(...))

}

#' Prepare the application population and calibrate transmission
#' @param p Model parameters.
#' @param inputs Application inputs for [layered_setup()].
#' @return Calibrated layered setup.
build_setup <- function(p, inputs) {
  setup <- do.call(layered_setup, c(list(p = p), inputs))
  calibrate_r0(p, setup)
}

#' Calibrate both layers against the simulated reproduction number
#' Estimates include index cases with zero secondary infections.
calibrate_r0 <- function(p, setup, reps = 50, iter = 3, n_seed = 10) {
  p$n_seed <- n_seed
  unvaccinated <- setup
  unvaccinated$vaccinated <- integer()
  estimate <- function(x) {
    # Estimate each seed's secondary infections in a susceptible population.
    values <- vapply(
      X         = seq_len(reps),
      FUN       = function(seed) {
        m <- build_model(p, x)
        epiworldR::run(m, ndays = 60, seed = seed)
        rt <- as.data.frame(epiworldR::get_reproductive_number(m))
        sum(rt$rt[rt$source_exposure_date == 0]) / n_seed
      },
      FUN.VALUE = numeric(1)
    )
    mean(values)
  }
  for (i in seq_len(iter)) {
    r <- estimate(unvaccinated)
    if (!is.finite(r) || r <= 0) stop("R0 calibration produced no transmission")
    factor <- p$R0 / r
    unvaccinated$community <- unvaccinated$community * factor
    unvaccinated$school_probability <- pmin(
      1,
      unvaccinated$school_probability * factor
    )
  }
  setup$community <- unvaccinated$community
  setup$school_probability <- unvaccinated$school_probability
  setup$R0_sim <- estimate(unvaccinated)
  unvaccinated$school_probability <- unvaccinated$school_probability *
    (1 - p$school_contact_reduction)
  setup$R_int <- estimate(unvaccinated)
  setup
}

#' Build the two-layer measles model used by all scenarios
build_model <- function(p, setup) {
  model <- layered_measles(
    p                  = p,
    population         = setup$persons,
    community_matrix   = setup$community,
    school_probability = setup$school_probability,
    vaccinated         = setup$vaccinated
  )
  epiworldR::verbose_off(model)
  model
}

#' Adds the surveillance process as a global event
#'
#' Every day the event computes new infections, prodrome onsets, and rash
#' onsets from the state counts, draws case reports and shedding, and
#' (on sampling days) a WW measurement. Results are stored in `state`.
#' @param m Model.
#' @param p Parameters.
#' @param setup Output of `build_setup()`.
#' @param scenario One of "baseline", "wastewater", or "cases".
#' @param state Environment where the daily series are stored.
add_surveillance <- function(m, p, setup, scenario, state) {

  n_days <- p$ndays
  meanlog <- log(10) * p$ww_log10_g_mean
  sdlog   <- log(10) * p$ww_log10_g_sd

  reset <- function() {
    variables <- c(
      "new_inf", "new_prod", "new_rash", "cohort_g", "ww_true",
      "cum_inf", "cum_prod", "cum_rash"
    )
    for (v in variables)
      assign(v, numeric(n_days), envir = state)
    state$reported    <- numeric(n_days + p$case_report_delay)
    state$ww_obs      <- rep(NA_real_, n_days)
    state$trigger_day <- NA_integer_
    # Day 0: seeded infections start in the latent state
    state$prev <- c(inf = p$n_seed, prod = 0, rash = 0)
  }

  fun <- function(m) {

    t <- epiworldR::today(m)
    if (t == 1)
      reset()

    tot <- epiworldR::get_today_total(m)
    if (is.data.frame(tot)) {
      counts <- tapply(tot$counts, tot$stage, sum)
      tot <- c(
        Susceptible               = 0,
        `Quarantined Susceptible` = 0,
        Latent                    = 0,
        `Quarantined Latent`      = 0,
        Prodromal                 = 0,
        `Quarantined Prodromal`   = 0
      )
      tot[names(counts)] <- counts
    }
    cum <- c(
      inf  = p$n - tot[["Susceptible"]] - tot[["Quarantined Susceptible"]]
    )
    cum["prod"] <- cum[["inf"]] - tot[["Latent"]] - tot[["Quarantined Latent"]]
    cum["rash"] <- cum[["prod"]] - tot[["Prodromal"]] -
      tot[["Quarantined Prodromal"]]

    new <- cum - state$prev
    state$prev <- cum

    state$new_inf[t]  <- new[["inf"]]
    state$new_prod[t] <- new[["prod"]]
    state$new_rash[t] <- new[["rash"]]
    state$cum_inf[t]  <- cum[["inf"]]
    state$cum_prod[t] <- cum[["prod"]]
    state$cum_rash[t] <- cum[["rash"]]

    # Clinical surveillance: each new rash case is detected with probability
    # case_report_prob, and reports arrive after the reporting delay. One
    # uniform is drawn per case (instead of a single binomial draw) so that
    # the number of random draws does not depend on case_report_prob. This
    # keeps the WW draws identical across detection probabilities, and the
    # cases detected at a lower probability are a subset of those detected
    # at a higher one (common random numbers for the sensitivity analysis).
    r <- t + p$case_report_delay
    state$reported[r] <- state$reported[r] +
      sum(stats::runif(new[["rash"]]) < p$case_report_prob)

    # Shedding: individual total genome copies are log-normal
    state$cohort_g[t] <- sum(stats::rlnorm(new[["prod"]], meanlog, sdlog))
    conc <- expected_concentration(
      state$cohort_g,
      t,
      kernel = setup$kernel,
      alpha  = p$ww_alpha,
      n      = p$n
    )
    state$ww_true[t] <- conc

    # The measurement error is drawn every day, and used only on sampling
    # days, so the number of random draws does not depend on the sampling
    # schedule (common random numbers across schedules)
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
      if (!is.null(setup$school_probability)) {
        epiworldR::set_contact_prob(
          layer        = epiworldR::get_layer(m, "schools"),
          contact_prob = setup$school_probability *
            (1 - p$school_contact_reduction)
        )
      } else {
        epiworldR::set_contact_matrix(m, setup$cm_int)
      }
    }

    invisible(NULL)

  }

  epiworldR::add_globalevent(
    model = m,
    event = epiworldR::globalevent_fun(fun, name = "Surveillance", day = -99)
  )

  # The model holds `fun`, whose environment is this frame, so the frame must
  # not reference the model; otherwise the model can never be freed
  rm(m)
  invisible(NULL)

}

#' Runs a single simulation
#' @param p Parameters.
#' @param setup Output of `build_setup()`.
#' @param scenario One of "baseline", "wastewater", or "cases".
#' @param seed Seed for the C++ model. The caller controls R's RNG.
#' @return A list with the daily series (`daily`) and the trigger day.
simulate_once <- function(p, setup, scenario, seed) {

  m     <- build_model(p, setup)
  # Empty parent: the default (this frame) would let the surveillance event
  # reach `m` and keep the model alive (see add_surveillance())
  state <- new.env(parent = emptyenv())
  add_surveillance(m, p, setup, scenario, state)

  epiworldR::run(m, ndays = p$ndays, seed = seed)

  days <- seq_len(p$ndays)
  list(
    daily = data.frame(
      day      = days,
      new_inf  = state$new_inf,
      new_prod = state$new_prod,
      new_rash = state$new_rash,
      reported = state$reported[days],
      ww_true  = state$ww_true,
      ww_obs   = state$ww_obs,
      cum_inf  = state$cum_inf,
      cum_prod = state$cum_prod,
      cum_rash = state$cum_rash
    ),
    trigger_day = state$trigger_day
  )

}

#' Runs paired simulations for several scenarios
#' @param p Parameters.
#' @param setup Output of `build_setup()`.
#' @param nsims Number of simulations per scenario.
#' @param scenarios Scenarios to run.
#' @param seed0 Seeds are `seed0 + 1:nsims`, shared across scenarios.
#' @param runner Execution callback; scripts can initialize R RNG per replicate.
#' @param cores Number of cores passed to `parallel::mclapply()`.
#' @return A data.table with one row per simulation and scenario, with the
#' daily series in the `daily` list column.
run_scenarios <- function(
    p,
    setup,
    nsims     = 100,
    scenarios = c("baseline", "wastewater", "cases"),
    seed0     = 1000,
    cores     = 1,
    runner    = simulate_once
) {

  grid <- expand.grid(
    sim              = seq_len(nsims),
    scenario         = scenarios,
    stringsAsFactors = FALSE
  )

  res <- parallel::mclapply(
    X        = seq_len(nrow(grid)),
    FUN      = function(i) {
      runner(p, setup, grid$scenario[i], seed0 + grid$sim[i])
    },
    mc.cores = cores
  )

  ans <- data.table::as.data.table(grid)
  ans[, trigger_day := vapply(res, `[[`, integer(1), "trigger_day")]
  ans[, daily := lapply(res, `[[`, "daily")]
  ans[, final_size := vapply(daily, function(d) d$cum_inf[nrow(d)], numeric(1))]
  ans[, reported := vapply(daily, function(d) sum(d$reported), numeric(1))]
  ans[, peak_day := vapply(
    X         = daily,
    FUN       = function(d) d$day[which.max(d$new_inf)],
    FUN.VALUE = numeric(1)
  )]
  ans[]

}

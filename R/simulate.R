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
# The intervention reduces non-household contacts of school-age groups (see
# `intervention_matrix()`) for the rest of the simulation.
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
    prodromal_period            = 4,
    rash_period                 = 4,
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

#' Pre-computes the population, contact matrices, and shedding kernel
#'
#' The contact matrix is scaled so the basic reproduction number equals
#' `R0`: R0 = transmission_rate * D * rho(M), with D the effective infectious
#' period (prodromal + rash period weighted by the remaining contact rate).
#' @param p Parameters from `measles_ww_params()`.
build_setup <- function(p) {

  pop <- city_population(p$n)
  cm  <- contact_matrices(pop)

  d_eff <- p$prodromal_period +
    (1 - p$rash_reduction_contact_rate) * p$rash_period

  rho   <- max(Re(eigen(cm$Total, only.values = TRUE)$values))
  scale <- p$R0 / (p$transmission_rate * d_eff * rho)

  cm_int <- intervention_matrix(cm, reduction = p$school_contact_reduction)

  setup <- list(
    pop     = pop,
    cm      = cm,
    cm_base = cm$Total * scale,
    cm_int  = cm_int * scale,
    R_int   = p$R0 * max(Re(eigen(cm_int, only.values = TRUE)$values)) / rho,
    kernel  = shedding_kernel(
      t_peak      = p$ww_t_peak,
      t_end       = p$ww_t_end,
      log10_vpeak = p$ww_log10_vpeak
    )
  )

  calibrate_r0(p, setup)

}

#' Rescales the contact matrices so the simulated R0 matches the target
#'
#' The analytic scaling in `build_setup()` overstates the R0 realized by the
#' agent-based model, so we estimate R0 by simulation (mean number of
#' secondary infections of `n_seed` index cases, without vaccination) and
#' rescale the contact matrices proportionally, a few times.
#' @param p Parameters.
#' @param setup Setup list with `cm_base` and `cm_int`.
#' @param reps Simulations per iteration.
#' @param iter Number of iterations.
#' @param n_seed Number of index cases per simulation.
#' @return The setup with rescaled matrices and the estimate in `R0_sim`.
calibrate_r0 <- function(
    p,
    setup,
    reps   = 100,
    iter   = 3,
    n_seed = 10
) {

  # R0 is defined for a fully susceptible population
  p$n_seed          <- n_seed
  p$prop_vaccinated <- 0

  estimate <- function() {
    mean(vapply(seq_len(reps), function(s) {
      m <- build_model(p, setup)
      epiworldR::run(m, ndays = 60, seed = s)
      rt <- as.data.frame(epiworldR::get_reproductive_number(m))
      mean(rt$rt[rt$source_exposure_date == 0])
    }, numeric(1)))
  }

  for (i in seq_len(iter)) {
    r <- estimate()
    setup$cm_base <- setup$cm_base * p$R0 / r
    setup$cm_int  <- setup$cm_int * p$R0 / r
  }

  setup$R0_sim <- estimate()
  setup

}

#' Builds the measles mixing model with one entity per age group
#' @param p Parameters.
#' @param setup Output of `build_setup()`.
build_model <- function(p, setup) {

  # Hospitalization rate matching the target probability (see
  # ?measles::ModelMeaslesMixing)
  h_rate <- p$prob_hospitalization * (1 / p$rash_period) /
    (1 - p$prob_hospitalization)

  m <- measles::ModelMeaslesMixing(
    n                            = p$n,
    prevalence                   = p$n_seed / p$n,
    contact_matrix               = setup$cm_base,
    transmission_rate            = p$transmission_rate,
    prop_vaccinated              = p$prop_vaccinated,
    vax_efficacy                 = p$vax_efficacy,
    incubation_period            = p$incubation_period,
    prodromal_period             = p$prodromal_period,
    rash_period                  = p$rash_period,
    rash_reduction_contact_rate  = p$rash_reduction_contact_rate,
    hospitalization_rate         = h_rate,
    hospitalization_period       = p$hospitalization_period,
    # Individual-level public health responses are off in all scenarios
    quarantine_willingness       = 0,
    isolation_willingness        = 0,
    contact_tracing_success_rate = 0
  )

  for (g in names(setup$pop)) {
    epiworldR::add_entity(
      m,
      epiworldR::entity(g, setup$pop[[g]], as_proportion = FALSE)
    )
  }

  epiworldR::verbose_off(m)
  m

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
    for (v in c("new_inf", "new_prod", "new_rash", "cohort_g", "ww_true",
                "cum_inf", "cum_prod", "cum_rash"))
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

    if ((t %% 7) %in% p$ww_sampling_days) {
      obs <- conc * exp(stats::rnorm(1, 0, p$ww_sigma))
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
      epiworldR::set_contact_matrix(m, setup$cm_int)
    }

    invisible(NULL)

  }

  epiworldR::add_globalevent(
    m,
    epiworldR::globalevent_fun(fun, name = "Surveillance", day = -99)
  )

  invisible(m)

}

#' Runs a single simulation
#' @param p Parameters.
#' @param setup Output of `build_setup()`.
#' @param scenario One of "baseline", "wastewater", or "cases".
#' @param seed Seed for both the model and R's RNG.
#' @return A list with the daily series (`daily`) and the trigger day.
simulate_once <- function(p, setup, scenario, seed) {

  m     <- build_model(p, setup)
  state <- new.env()
  add_surveillance(m, p, setup, scenario, state)

  set.seed(seed)
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
#' @param cores Number of cores passed to `parallel::mclapply()`.
#' @return A data.table with one row per simulation and scenario, with the
#' daily series in the `daily` list column.
run_scenarios <- function(
    p,
    setup,
    nsims     = 100,
    scenarios = c("baseline", "wastewater", "cases"),
    seed0     = 1000,
    cores     = 1
) {

  grid <- expand.grid(
    sim              = seq_len(nsims),
    scenario         = scenarios,
    stringsAsFactors = FALSE
  )

  res <- parallel::mclapply(
    seq_len(nrow(grid)),
    function(i) simulate_once(p, setup, grid$scenario[i], seed0 + grid$sim[i]),
    mc.cores = cores
  )

  ans <- data.table::as.data.table(grid)
  ans[, trigger_day := vapply(res, `[[`, integer(1), "trigger_day")]
  ans[, daily := lapply(res, `[[`, "daily")]
  ans[, final_size := vapply(daily, function(d) d$cum_inf[nrow(d)], numeric(1))]
  ans[, reported := vapply(daily, function(d) sum(d$reported), numeric(1))]
  ans[, peak_day := vapply(daily, function(d) d$day[which.max(d$new_inf)], numeric(1))]
  ans[]

}

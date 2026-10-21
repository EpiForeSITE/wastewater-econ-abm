# Runs the scenarios and stores aggregate results in analysis/results.rds.
# Usage (from the project root): Rscript analysis/run-scenarios.R

library(data.table)
for (f in list.files("R", full.names = TRUE)) source(f)

nsims        <- 400
cores        <- max(1, parallel::detectCores() - 1)
coverages    <- c(0, 0.85)
lab_delays   <- c(0, 1, 3, 5, 7)
outbreak_min <- 50 # baseline infections for a run to count as an outbreak

summarise_runs <- function(runs, baseline, p) {

  base <- baseline[, .(sim, base_final = final_size, base_daily = daily)]
  runs <- merge(runs, base, by = "sim")

  # Paired lead time: case-based trigger day minus WW trigger day, both
  # computed on the shared baseline trajectory
  runs[, lead := vapply(base_daily, function(d) {
    tr <- offline_triggers(d, p)
    tr$t_cases - tr$t_ww
  }, numeric(1))]

  at_trigger <- function(d, t, f) if (is.na(t)) NA_real_ else f(d, t)

  runs[, infected_at_trigger := mapply(
    at_trigger,
    daily,
    trigger_day,
    MoreArgs = list(f = function(d, t) d$cum_inf[t])
  )]
  runs[, reported_at_trigger := mapply(
    at_trigger,
    daily,
    trigger_day,
    MoreArgs = list(f = function(d, t) sum(d$reported[seq_len(t)]))
  )]

  # Summaries are conditional on runs where the baseline had an outbreak
  pct_outbreak <- 100 * mean(baseline$final_size >= outbreak_min)
  runs <- runs[base_final >= outbreak_min]

  runs[, .(
    pct_outbreak        = pct_outbreak,
    nsims               = .N,
    pct_trigger         = 100 * mean(!is.na(trigger_day)),
    trigger_day         = median(as.numeric(trigger_day), na.rm = TRUE),
    final_size          = median(final_size),
    final_q25           = quantile(final_size, .25),
    final_q75           = quantile(final_size, .75),
    averted             = mean(base_final - final_size),
    averted_q25         = quantile(base_final - final_size, .25),
    averted_q75         = quantile(base_final - final_size, .75),
    lead_median         = median(lead, na.rm = TRUE),
    infected_at_trigger = median(infected_at_trigger, na.rm = TRUE),
    reported_at_trigger = median(reported_at_trigger, na.rm = TRUE)
  ), by = scenario]

}

mean_curve <- function(runs, label) {
  x <- runs[, .(
    day     = daily[[1]]$day,
    new_inf = daily[[1]]$new_inf,
    ww_true = daily[[1]]$ww_true
  ), by = .(sim, scenario)]
  x[, .(new_inf = mean(new_inf), ww_true = mean(ww_true)), by = .(scenario, day)][
    , coverage := label][]
}

results <- list()

for (cov in coverages) {

  message("Coverage: ", cov)
  p     <- measles_ww_params(prop_vaccinated = cov)
  setup <- build_setup(p)

  baseline <- run_scenarios(p, setup, nsims, "baseline", cores = cores)

  cal <- calibrate_ww_threshold(baseline, p)
  p$trigger_ww <- cal$threshold

  triggered <- run_scenarios(
    p,
    setup,
    nsims,
    c("wastewater", "cases"),
    cores = cores
  )

  # Sensitivity: WW lab turnaround time (threshold kept fixed)
  lab <- rbindlist(lapply(lab_delays, function(ld) {
    pl  <- utils::modifyList(p, list(ww_lab_delay = ld))
    out <- run_scenarios(pl, setup, nsims, "wastewater", cores = cores)
    s   <- summarise_runs(out, baseline, pl)
    s[, ww_lab_delay := ld][]
  }))

  # Trigger-time diagnostics on the baseline trajectories
  diag <- rbindlist(lapply(baseline$daily, function(d) {
    tr <- offline_triggers(d, p)
    data.table(
      t_cases         = tr$t_cases,
      t_ww            = tr$t_ww,
      seen_cases      = tr$seen_cases,
      seen_ww         = tr$seen_ww,
      reported_at_ww  = if (is.na(tr$t_ww)) NA_real_ else
        sum(d$reported[seq_len(tr$t_ww)]),
      infected_at_ww  = if (is.na(tr$t_ww)) NA_real_ else d$cum_inf[tr$t_ww],
      infected_at_cas = if (is.na(tr$t_cases)) NA_real_ else
        d$cum_inf[tr$t_cases]
    )
  }))

  results[[as.character(cov)]] <- list(
    params      = p,
    R0_sim      = setup$R0_sim,
    R_int       = setup$R_int,
    pop         = setup$pop,
    calibration = cal,
    summary     = summarise_runs(
      rbind(baseline, triggered),
      baseline,
      p
    ),
    lab_delay   = lab,
    diagnostics = diag,
    curves      = mean_curve(rbind(baseline, triggered), cov)
  )

}

saveRDS(results, "analysis/results.rds")

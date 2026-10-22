# Runs the scenarios and stores aggregate results in analysis/results.rds.
# Usage (from the project root): Rscript analysis/run-scenarios.R

library(data.table)
for (f in list.files("R", full.names = TRUE)) source(f)

nsims        <- 400
cores        <- max(1, parallel::detectCores() - 1)
coverages    <- c(0, 0.85)
lab_delays   <- c(0, 1, 3, 5, 7, 10)
# WW sampling schedules (days of the week, t %% 7)
sampling_schedules <- list(
  "1" = 1,
  "2" = c(1, 4),
  "3" = c(1, 3, 5),
  "5" = 1:5
)
# Probability that a case is detected during its rash period. The reference
# value is used for the main scenarios and to calibrate the WW threshold;
# the sensitivity analysis varies it with the WW threshold held fixed.
report_prob_ref <- 0.7
report_probs    <- c(0.7, 0.5, 0.3, 0.1)
outbreak_min <- 50 # baseline infections for a run to count as an outbreak

# `ww_days` (optional): data.table with `sim` and `t_ww`, the WW trigger day
# of each run. When given, the lead is the run's own trigger day minus
# `t_ww`; this is used for case-based runs whose detection probability
# differs from the one used to generate the baseline reports.
# `case_days` (optional): data.table with `sim` and `t_cases`, the case-based
# trigger day of each run. When given, the lead is `t_cases` minus the run's
# own trigger day; this is used for WW runs whose sampling schedule differs
# from the one used to generate the baseline WW series.
summarise_runs <- function(
    runs,
    baseline,
    p,
    ww_days   = NULL,
    case_days = NULL
) {

  base <- baseline[, .(sim, base_final = final_size, base_daily = daily)]
  runs <- merge(runs, base, by = "sim")

  # Paired lead time: case-based trigger day minus WW trigger day, both
  # computed on the shared baseline trajectory
  if (!is.null(ww_days)) {
    runs <- merge(runs, ww_days, by = "sim")
    runs[, lead := as.numeric(trigger_day) - as.numeric(t_ww)]
  } else if (!is.null(case_days)) {
    runs <- merge(runs, case_days, by = "sim")
    runs[, lead := as.numeric(t_cases) - as.numeric(trigger_day)]
  } else {
    runs[, lead := vapply(base_daily, function(d) {
      tr <- offline_triggers(d, p)
      tr$t_cases - tr$t_ww
    }, numeric(1))]
  }

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
  p     <- measles_ww_params(
    prop_vaccinated  = cov,
    case_report_prob = report_prob_ref
  )
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

  # Sensitivity: WW sampling frequency (threshold kept fixed). The WW
  # measurement error is drawn every day (see add_surveillance()), so the
  # schedules share their draws and differ only in which days are sampled.
  # The baseline WW series uses the default schedule, so the lead is paired
  # with the case-based trigger of the main runs instead.
  case_days <- triggered[scenario == "cases", .(sim, t_cases = trigger_day)]
  sampling <- rbindlist(lapply(names(sampling_schedules), function(k) {
    ps  <- utils::modifyList(
      p,
      list(ww_sampling_days = sampling_schedules[[k]])
    )
    out <- run_scenarios(ps, setup, nsims, "wastewater", cores = cores)
    s   <- summarise_runs(out, baseline, ps, case_days = case_days)
    s[, samples_per_week := as.integer(k)][]
  }))

  # Sensitivity: probability of detecting a case during the rash period.
  # The WW threshold stays at its calibrated value, so only the case-based
  # trigger changes. Case detection uses common random numbers (see
  # add_surveillance()), so WW trigger days do not depend on the detection
  # probability and the lead is paired with the main WW runs.
  ww_days <- triggered[scenario == "wastewater", .(sim, t_ww = trigger_day)]
  report <- rbindlist(lapply(report_probs, function(rp) {
    pr  <- utils::modifyList(p, list(case_report_prob = rp))
    out <- run_scenarios(pr, setup, nsims, "cases", cores = cores)
    s   <- summarise_runs(out, baseline, pr, ww_days = ww_days)
    s[, case_report_prob := rp][]
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
    sampling    = sampling,
    report_prob = report,
    diagnostics = diag,
    curves      = mean_curve(rbind(baseline, triggered), cov)
  )

}

saveRDS(results, "analysis/results.rds")

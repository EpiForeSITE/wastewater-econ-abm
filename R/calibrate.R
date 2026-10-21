# Calibration of the WW threshold against the case-based threshold.
#
# Both triggers should fire when the same number of cases has been "seen" by
# each surveillance system. The clinical system sees a case at rash onset;
# WW sees it from the onset of shedding (prodrome). So, using baseline runs
# (no intervention), we compare:
#
#  - cases: cumulative rash onsets on the rash day of the report that crosses
#    `trigger_cases` (trigger day - case_report_delay);
#  - WW:    cumulative prodrome onsets on the collection day of the first
#    sample at or above the WW threshold (trigger day - ww_lab_delay);
#
# and pick the WW threshold whose median matches the case-based median.
# Because the comparison uses the collection day, the calibrated threshold
# does not depend on `ww_lab_delay`.

#' Trigger days computed from a baseline daily series
#' @param d Daily data frame from `simulate_once()`.
#' @param p Parameters.
#' @param trigger_ww WW threshold (copies/mL).
#' @return A list with the trigger days and the cases seen at each trigger.
offline_triggers <- function(d, p, trigger_ww = p$trigger_ww) {

  ndays <- nrow(d)

  t_cases <- which(cumsum(d$reported) >= p$trigger_cases)[1]
  seen_cases <- if (is.na(t_cases)) NA_real_ else
    d$cum_rash[max(1, t_cases - p$case_report_delay)]

  s_ww <- which(!is.na(d$ww_obs) & d$ww_obs >= trigger_ww &
    d$day + p$ww_lab_delay <= ndays)[1]
  seen_ww <- if (is.na(s_ww)) NA_real_ else d$cum_prod[s_ww]

  list(
    t_cases    = t_cases,
    t_ww       = s_ww + p$ww_lab_delay,
    seen_cases = seen_cases,
    seen_ww    = seen_ww
  )

}

#' Calibrates the WW threshold
#' @param baseline Output of `run_scenarios()` for the baseline scenario.
#' @param p Parameters.
#' @param grid Candidate thresholds (copies/mL); values below the LOD are
#' dropped.
#' @return A list with the selected threshold and the calibration table.
calibrate_ww_threshold <- function(
    baseline,
    p,
    grid = 10^seq(0, 3, by = 0.05)
) {

  grid <- grid[grid >= p$ww_lod]

  seen_cases <- vapply(
    baseline$daily,
    function(d) offline_triggers(d, p)$seen_cases,
    numeric(1)
  )
  target <- stats::median(seen_cases, na.rm = TRUE)

  tab <- data.frame(
    threshold = grid,
    seen_ww   = vapply(grid, function(thr) {
      stats::median(vapply(
        baseline$daily,
        function(d) offline_triggers(d, p, thr)$seen_ww,
        numeric(1)
      ), na.rm = TRUE)
    }, numeric(1))
  )

  best <- which.min(abs(tab$seen_ww - target))

  list(
    threshold = tab$threshold[best],
    target    = target,
    table     = tab
  )

}

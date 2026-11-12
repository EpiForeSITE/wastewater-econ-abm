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
# Both medians are taken over the same runs: those where the case trigger
# fires. A run where WW never crosses the threshold counts as WW never
# seeing the outbreak (`Inf`), so thresholds that miss outbreaks are
# penalized instead of being judged only on the runs they detect.
# Because the comparison uses the collection day, the calibrated threshold
# does not depend on `ww_lab_delay`.

#' Trigger days computed from a baseline daily series
#' @param d Daily data frame from `simulate_once()`.
#' @param p Parameters.
#' @param trigger_ww WW threshold (copies/mL).
#' @return A list with the trigger days and the cases seen at each trigger.
#' @examples
#' p <- measles_ww_params(
#'   trigger_cases = 1,
#'   case_report_delay = 0,
#'   ww_lab_delay = 0
#' )
#' d <- data.frame(
#'   day      = 1:3,
#'   reported = c(0, 1, 0),
#'   ww_obs   = c(0, 2, NA),
#'   cum_rash = c(0, 1, 1),
#'   cum_prod = c(0, 1, 1)
#' )
#' offline_triggers(d = d, p = p, trigger_ww = 1)
offline_triggers <- function(d, p, trigger_ww = p$trigger_ww) {

  ndays <- nrow(d)

  t_cases <- which(cumsum(d$reported) >= p$trigger_cases)[1]
  seen_cases <- if (is.na(t_cases)) NA_real_ else
    d$cum_rash[max(1, t_cases - p$case_report_delay)]

  s_ww <- which(
    !is.na(d$ww_obs) & d$ww_obs >= trigger_ww &
      d$day + p$ww_lab_delay <= ndays
  )[1]
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
#' @return A list with the selected threshold, the target, the number of runs
#' used (those where the case trigger fires), and the calibration table, which
#' also gives the percentage of those runs where WW crosses each threshold.
#' @examples
#' p <- measles_ww_params(
#'   trigger_cases = 1,
#'   case_report_delay = 0,
#'   ww_lab_delay = 0
#' )
#' d <- data.frame(
#'   day      = 1:3,
#'   reported = c(0, 1, 0),
#'   ww_obs   = c(0, 2, NA),
#'   cum_rash = c(0, 1, 1),
#'   cum_prod = c(0, 1, 1)
#' )
#' baseline <- data.frame(sim = 1)
#' baseline$daily <- list(d)
#' calibrate_ww_threshold(baseline = baseline, p = p, grid = c(1, 2))$table
calibrate_ww_threshold <- function(
    baseline,
    p,
    grid = 10^seq(0, 3, by = 0.05)
) {

  grid <- grid[grid >= p$ww_lod]

  seen_cases <- vapply(
    X         = baseline$daily,
    FUN       = function(d) offline_triggers(d, p)$seen_cases,
    FUN.VALUE = numeric(1)
  )
  eligible <- !is.na(seen_cases)
  if (!any(eligible))
    stop("The case trigger never fires in the baseline runs.")
  target <- stats::median(seen_cases[eligible])

  seen_ww <- vapply(
    X         = grid,
    FUN       = function(thr) {
      vapply(
        X         = baseline$daily[eligible],
        FUN       = function(d) offline_triggers(d, p, thr)$seen_ww,
        FUN.VALUE = numeric(1)
      )
    },
    FUN.VALUE = numeric(sum(eligible))
  )
  seen_ww <- matrix(seen_ww, ncol = length(grid))
  seen_ww[is.na(seen_ww)] <- Inf

  tab <- data.frame(
    threshold   = grid,
    seen_ww     = apply(seen_ww, 2, stats::median),
    pct_trigger = 100 * colMeans(is.finite(seen_ww))
  )

  best <- which.min(abs(tab$seen_ww - target))

  list(
    threshold = tab$threshold[best],
    target    = target,
    n_runs    = sum(eligible),
    table     = tab
  )

}

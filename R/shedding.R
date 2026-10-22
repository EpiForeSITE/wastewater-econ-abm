# Wastewater shedding and measurement, following the structure of the CDC
# wwinference model (https://github.com/CDCgov/ww-inference-model):
#
#   C(t) = 1 / (alpha * N) * sum_tau s(tau) * G(t - tau)
#
# where s(tau) is a normalized shedding profile that is piecewise linear in
# log10, G(d) is the total number of genome copies shed by the cohort that
# started shedding on day d, alpha is the wastewater volume (mL) per person
# per day, and N is the population size. Unlike the COVID model, where tau is
# counted from infection, here tau is counted from the onset of the
# prodromal period (measles RNA appears in urine during the prodrome).
# Observations are log-normal around C(t) and censored below the LOD.

#' Normalized shedding profile
#'
#' log10 s(tau) rises linearly from 0 at tau = 0 to `log10_vpeak` at
#' `t_peak`, then declines linearly back to 0 at `t_end`. The profile is
#' normalized to sum to one, so `G` is the total shed per infection.
#' @param t_peak Days from prodrome onset to peak shedding.
#' @param t_end Days from prodrome onset to the end of shedding.
#' @param log10_vpeak Peak (log10) relative to the start and end of shedding.
#' @return A numeric vector for tau = 0, ..., t_end.
shedding_kernel <- function(
    t_peak      = 4,
    t_end       = 18,
    log10_vpeak = 3
) {

  tau <- 0:t_end
  lv  <- ifelse(
    tau <= t_peak,
    log10_vpeak * tau / t_peak,
    log10_vpeak * (1 - (tau - t_peak) / (t_end - t_peak))
  )

  s <- 10^lv
  s / sum(s)

}

#' Expected concentration (copies/mL) on day t
#' @param cohort_g Vector of genome copies by cohort day (index 1 = day 1).
#' @param t Day.
#' @param kernel Output of `shedding_kernel()`.
#' @param alpha Wastewater volume per person per day (mL).
#' @param n Population size.
expected_concentration <- function(cohort_g, t, kernel, alpha, n) {

  tau <- seq_along(kernel) - 1
  tau <- tau[tau < t]
  sum(kernel[tau + 1] * cohort_g[t - tau]) / (alpha * n)

}

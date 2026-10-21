# End-to-end tests of the measles wastewater simulation.
# Run from the project root with: tinytest::run_test_dir("tests")

library(data.table)
for (f in list.files(file.path("..", "R"), full.names = TRUE)) source(f)
old <- setwd("..")

# Shedding profile is a proper distribution peaking at t_peak
k <- shedding_kernel(t_peak = 4, t_end = 18)
expect_equal(sum(k), 1)
expect_equal(which.max(k) - 1, 4)

p     <- measles_ww_params(ndays = 150, trigger_ww = 5)
setup <- build_setup(p)

expect_equal(sum(setup$pop), p$n)
expect_true(abs(setup$R0_sim - p$R0) / p$R0 < 0.1)
expect_true(setup$R_int < p$R0)

runs <- run_scenarios(p, setup, nsims = 20)

# Baseline never triggers; intervention scenarios trigger in outbreaks
expect_true(all(is.na(runs[scenario == "baseline", trigger_day])))
expect_true(any(!is.na(runs[scenario == "wastewater", trigger_day])))

# Paired runs share the trajectory up to the trigger day
b <- runs[scenario == "baseline"]
w <- runs[scenario == "wastewater"]
for (i in which(!is.na(w$trigger_day))) {
  td <- w$trigger_day[i]
  expect_equal(b$daily[[i]]$cum_inf[1:td], w$daily[[i]]$cum_inf[1:td])
}

# The in-model triggers match those recomputed from the baseline series
tr <- lapply(b$daily, offline_triggers, p = p)
expect_equal(
  w$trigger_day,
  vapply(tr, function(x) as.integer(x$t_ww), integer(1))
)
expect_equal(
  runs[scenario == "cases", trigger_day],
  vapply(tr, function(x) as.integer(x$t_cases), integer(1))
)

# Intervening never increases infections much on average
expect_true(
  mean(runs[scenario == "wastewater", final_size]) <=
    mean(b$final_size)
)

# Calibration returns a threshold at or above the LOD
cal <- calibrate_ww_threshold(b, p)
expect_true(cal$threshold >= p$ww_lod)

setwd(old)

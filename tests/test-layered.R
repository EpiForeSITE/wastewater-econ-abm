# Run from repository root: Rscript tests/test-layered.R
if (!dir.exists("R")) setwd("..")
for (f in list.files("R", full.names = TRUE)) source(f)
p <- measles_ww_params(n = 4000, ndays = 80, n_seed = 20,
  prop_vaccinated = 0, trigger_cases = 1, case_report_prob = 1,
  case_report_delay = 0)
setup <- layered_setup(p)
model <- build_model(p, setup)
stopifnot(inherits(epiworldR::get_layer(model, "community"), "epiworld_layer_mixing"),
          inherits(epiworldR::get_layer(model, "schools"), "epiworld_layer_mixing"))
stopifnot(sum(setup$pop) == p$n, length(setup$persons$school) == p$n,
          all(setup$school_probability >= 0 & setup$school_probability <= 1),
          all(setup$persons$coverage[!is.na(setup$persons$school)] >= 0),
          abs(sum(school_weights()) - 1) < 1e-10)
base <- simulate_once(p, setup, "baseline", 42)
cases <- simulate_once(p, setup, "cases", 42)
ww <- simulate_once(p, setup, "wastewater", 42)
stopifnot(is.na(base$trigger_day), !is.na(cases$trigger_day),
          all(diff(base$daily$cum_inf) >= 0),
          identical(base$daily$cum_inf[seq_len(cases$trigger_day)],
                    cases$daily$cum_inf[seq_len(cases$trigger_day)]),
          identical(base$daily$ww_true, ww$daily$ww_true) || !is.na(ww$trigger_day))
cat("Layered model checks passed.\n")

# A seeded cohort with contacts and hospitalization disabled exercises the
# actual duration transitions. A fixed duration must not have the early and
# long tail produced by a geometric waiting time.
q <- measles_ww_params(n = 1000, n_seed = 1000, incubation_period_sd = 0,
  prodromal_period_sd = 0, rash_period_sd = 0, prob_hospitalization = 0,
  rash_reduction_contact_rate = 0)
population <- list(age_group = rep(0L, q$n), school = rep(0L, q$n))
fixed <- layered_measles(q, population, matrix(0, 1, 1), 0)
epiworldR::verbose_off(fixed)
epiworldR::run(fixed, ndays = 40, seed = 9)
hist <- epiworldR::get_hist_total(fixed)
latent <- hist[hist$stage == "Latent", ]
stopifnot(all(latent$counts[latent$date < q$incubation_period] == q$n),
          all(latent$counts[latent$date > q$incubation_period + 1] == 0))
q$incubation_period_sd <- 2.5
varying <- layered_measles(q, population, matrix(0, 1, 1), 0)
epiworldR::verbose_off(varying)
epiworldR::run(varying, ndays = 40, seed = 9)
hist_varying <- epiworldR::get_hist_total(varying)
latent_varying <- hist_varying[hist_varying$stage == "Latent", ]
stopifnot(any(latent_varying$counts[latent_varying$date < q$incubation_period] < q$n),
          any(latent_varying$counts[latent_varying$date > q$incubation_period + 1] > 0))

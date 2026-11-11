# Run from the repository root with epiworldR v1-layers installed.
for (f in list.files("R", full.names = TRUE)) source(f)

p <- measles_ww_params(n = 10000, ndays = 120, n_seed = 10,
                       prop_vaccinated = 0.85)
setup <- layered_setup(p)
baseline <- simulate_layered_once(p, setup, "baseline", seed = 42)
cases <- simulate_layered_once(p, setup, "cases", seed = 42)
wastewater <- simulate_layered_once(p, setup, "wastewater", seed = 42)
print(data.frame(
  scenario = c("baseline", "cases", "wastewater"),
  trigger_day = c(baseline$trigger_day, cases$trigger_day, wastewater$trigger_day),
  infections = c(tail(baseline$daily$cum_inf, 1),
                 tail(cases$daily$cum_inf, 1), tail(wastewater$daily$cum_inf, 1))
))

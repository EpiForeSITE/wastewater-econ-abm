# Test the main model builder and the native multi-run workflow.
source("../model.R")
set.seed(137)
parameters <- read.csv("../data/parameters.csv", stringsAsFactors = FALSE)
p <- lapply(parameters$default, function(x) eval(parse(text = x)))
names(p) <- parameters$r_argument
p$n <- 1000
p$n_seed <- 5
p$ndays <- 30
p$ww_sampling_days <- 1:7
p$ww_lab_delay <- 0
p$ww_lod <- 0

# A small synthetic population lets the tests run without external data.
setup <- list(
  persons = list(
    age_group = rep(0L, p$n),
    school    = rep(0:9, each = 100)
  ),
  community          = matrix(2, 1, 1),
  school_probability = rep(0.04, 10),
  vaccinated         = integer(),
  kernel             = rep(1 / 19, 19)
)

# Baseline runs must retain one surveillance record per replicate.
state <- new.env(parent = emptyenv())
model <- build_model(p, setup, state, scenario = "baseline")
expect_true(inherits(get_layer(model, "community"), "epiworld_layer_mixing"))
expect_true(inherits(get_layer(model, "schools"), "epiworld_layer_mixing"))
run_multiple(
  m        = model,
  ndays    = p$ndays,
  nsims    = 3,
  seed     = 42,
  saver    = make_saver("total_hist"),
  nthreads = 1,
  verbose  = FALSE
)
history <- run_multiple_get_results(model, nthreads = 1)$total_hist
expect_equal(length(state$runs), 3L)
expect_equal(length(unique(history$sim_num)), 3L)
expect_true(all(vapply(state$runs, function(x) is.na(x$trigger_day), TRUE)))
valid_runs <- vapply(
  X         = state$runs,
  FUN       = function(x) {
    length(x$cum_inf) == p$ndays && all(diff(x$cum_inf) >= 0)
  },
  FUN.VALUE = TRUE
)
expect_true(all(valid_runs))
counts <- aggregate(counts ~ sim_num + date, data = history, FUN = sum)
expect_true(all(counts$counts == p$n))
susceptible <- history[
  history$stage == "Susceptible" & history$date == p$ndays,
]
susceptible <- susceptible[order(susceptible$sim_num), ]
recorded <- vapply(state$runs, function(x) tail(x$cum_inf, 1), numeric(1))
expect_equal(recorded, as.numeric(p$n - susceptible$counts))

# Forced triggers exercise the school-contact intervention in both scenarios.
p$trigger_ww <- 0
p$trigger_cases <- 0
for (scenario in c("wastewater", "cases")) {
  state <- new.env(parent = emptyenv())
  model <- build_model(p, setup, state, scenario = scenario)
  run_multiple(
    m        = model,
    ndays    = p$ndays,
    nsims    = 2,
    seed     = 42,
    saver    = make_saver("total_hist"),
    nthreads = 1,
    verbose  = FALSE
  )
  expect_equal(length(state$runs), 2L)
  expect_true(all(vapply(state$runs, function(x) x$trigger_day == 1L, TRUE)))
  expect_equal(
    get_contact_prob(get_layer(model, "schools")),
    setup$school_probability * (1 - p$school_contact_reduction)
  )
}

# Hospitalization hazard increases with the target hospitalization probability.
expect_equal(rash_hospitalization_rate(0, 4, 1), 0)
expect_true(rash_hospitalization_rate(0.2, 4, 1) > 0)

# Fixed and variable latent durations must use the configured distribution.
p$n_seed <- p$n
p$ndays <- 40
p$incubation_period_sd <- 0
p$rash_reduction_contact_rate <- 0
p$prob_hospitalization <- 0
setup$community[,] <- 0
setup$school_probability[] <- 0
state <- new.env(parent = emptyenv())
model <- build_model(p, setup, state)
run(model, ndays = p$ndays, seed = 9)
history <- get_hist_total(model)
latent <- history[history$stage == "Latent", ]
expect_true(all(latent$counts[latent$date < p$incubation_period] == p$n))
expect_true(all(latent$counts[latent$date > p$incubation_period + 1] == 0))

p$incubation_period_sd <- 2.5
state <- new.env(parent = emptyenv())
model <- build_model(p, setup, state)
run(model, ndays = p$ndays, seed = 9)
history <- get_hist_total(model)
latent <- history[history$stage == "Latent", ]
expect_true(any(latent$counts[latent$date < p$incubation_period] < p$n))
expect_true(any(latent$counts[latent$date > p$incubation_period + 1] > 0))

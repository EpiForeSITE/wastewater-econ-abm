# Test the main model builder, the response event, and the population helpers.
source("../model.R")
set.seed(137)
defaults <- read_parameters("../data/parameters.csv")
p <- defaults
p$detection_rate <- 0.5
p$alert_detection_rate <- 0.6
p$n <- 2000
p$n_seed <- 5
p$ndays <- 60

# Four schools of 250 students, half vaccinated, and 1,000 non-students.
school <- c(rep(0:3, each = 250), rep(NA_integer_, 1000))
setup <- list(
  persons = list(
    age_group = rep(0:1, each = 1000),
    school    = school
  ),
  community          = matrix(c(3, 1, 1, 3), 2),
  school_probability = rep(0.02, 4),
  vaccinated         = seq(0L, 998L, by = 2L),
  kernel             = rep(1 / 19, 19)
)
p$transmission_rate <- transmission_for_r0(p, setup$community)
expect_true(p$transmission_rate > 0 && p$transmission_rate < 1)

# Baseline runs keep one consistent record per replicate.
state <- new.env(parent = emptyenv())
model <- build_model(p, setup, state, response = "none")
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
counts <- aggregate(counts ~ sim_num + date, data = history, FUN = sum)
expect_true(all(counts$counts == p$n))
susceptible <- history[
  history$stage == "Susceptible" & history$date == p$ndays,
]
susceptible <- susceptible[order(susceptible$sim_num), ]
recorded <- vapply(
  X         = state$runs,
  FUN       = function(x) p$n_seed + sum(x$daily$new_infections),
  FUN.VALUE = numeric(1)
)
expect_equal(recorded, as.numeric(p$n - susceptible$counts))
expect_equal(
  vapply(state$runs, function(x) sum(x$infections_by_age), numeric(1)),
  recorded
)
expect_true(all(vapply(state$runs, function(x) sum(x$daily$excluded), 0) == 0))

# Paired runs share every draw until the first response acts.
run_paired <- function(scenario) {
  set.seed(11)
  state <- new.env(parent = emptyenv())
  model <- build_model(p, setup, state, scenario, seeds = c(1L, 3L, 5L))
  run(model, ndays = p$ndays, seed = 11)
  state$runs[[1]]
}
none <- run_paired("none")
cases <- run_paired("cases")
expect_false(is.na(cases$first_exclusion))
before <- seq_len(cases$first_exclusion)
signals <- setdiff(names(none$daily), "excluded")
expect_equal(none$daily[before, signals], cases$daily[before, signals])

# One detected student excludes the school's other unvaccinated students
# for exactly the response period, and nobody else.
q <- p
q$transmission_rate <- 0
q$detection_rate <- 1
q$case_report_delay <- 0
q$prob_hospitalization <- 0
q$ndays <- 100
state <- new.env(parent = emptyenv())
model <- build_model(q, setup, state, response = "cases", seeds = 1L)
run(model, ndays = q$ndays, seed = 3)
run1 <- state$runs[[1]]
unvaccinated_mates <- sum(school[seq(2, 250, by = 2)] == 0) - 1
expect_equal(run1$episodes, 1L)
expect_equal(max(run1$daily$excluded), unvaccinated_mates)
expect_equal(sum(run1$daily$excluded > 0), q$response_days)
expect_equal(
  which(run1$daily$excluded > 0)[1],
  which(run1$daily$reported > 0)[1]
)

# Wastewater alerts change contacts and detection only in their scenario.
q <- p
q$ww_false_positive <- 1
q$ww_sampling_days <- 0:6
q$ndays <- 30
for (scenario in c("cases", "wastewater")) {
  state <- new.env(parent = emptyenv())
  model <- build_model(q, setup, state, response = scenario)
  run(model, ndays = q$ndays, seed = 5)
  run1 <- state$runs[[1]]
  acting <- scenario == "wastewater"
  expect_equal(run1$first_alert, 1L + q$ww_lab_delay)
  expect_true(all(run1$daily$alert_on[-seq_len(q$ww_lab_delay)]))
  expect_equal(run1$alerts, 1L)
  expect_equal(
    get_param(model, "Community contact scale"),
    if (acting) 1 - q$community_contact_reduction else 1
  )
  expect_equal(
    get_param(model, "Detection rate"),
    if (acting) q$alert_detection_rate else q$detection_rate
  )
}

# Students are dealt into schools of about the target size.
level <- rep(c(1, 2, NA), c(2600, 1400, 1000))
assigned <- assign_schools(
  level        = level,
  source_level = c(1, 1, 1, 1, 1, 1, 2, 2),
  weight       = c(1, 2, 3, 1, 1, 1, 1, 1),
  size         = 500
)
sizes <- table(assigned$school)
expect_equal(length(sizes), 5L + 3L)
expect_true(all(abs(sizes - 500) <= 70))
expect_true(all(is.na(assigned$school[is.na(level)])))
expect_equal(length(assigned$source), 8L)

# Hospitalization hazard increases with the target hospitalization probability.
expect_equal(rash_hospitalization_rate(0, 4, 1), 0)
expect_true(rash_hospitalization_rate(0.2, 4, 1) > 0)

# Fixed and variable latent durations must use the configured distribution.
q <- p
q$n_seed <- q$n
q$ndays <- 40
q$incubation_period_sd <- 0
q$prob_hospitalization <- 0
q$transmission_rate <- 0
state <- new.env(parent = emptyenv())
model <- build_model(q, setup, state)
run(model, ndays = q$ndays, seed = 9)
history <- get_hist_total(model)
latent <- history[history$stage == "Latent", ]
expect_true(all(latent$counts[latent$date < q$incubation_period] == q$n))
expect_true(all(latent$counts[latent$date > q$incubation_period + 1] == 0))

q$incubation_period_sd <- 2.5
state <- new.env(parent = emptyenv())
model <- build_model(q, setup, state)
run(model, ndays = q$ndays, seed = 9)
history <- get_hist_total(model)
latent <- history[history$stage == "Latent", ]
expect_true(any(latent$counts[latent$date < q$incubation_period] < q$n))
expect_true(any(latent$counts[latent$date > q$incubation_period + 1] > 0))

# Splitting contacts into community and school layers keeps age totals and R0.
p <- defaults
p$detection_rate <- 0.5
p$n <- 20000
p$n_seed <- 20
p$ndays <- 40
p$R0 <- 3
age_group <- rep(0:1, c(8000, 12000))
total <- matrix(c(10, 2, 3, 8), 2)
school_matrix <- matrix(c(6, 0, 0, 0), 2)
school_id <- rep(NA_integer_, p$n)
school_id[1:6000] <- rep(0:29, each = 200)
attendance <- as.numeric(tapply(!is.na(school_id), age_group, mean))
contact_split <- split_school_contacts(total, school_matrix, attendance)
probability <- school_contact_prob(
  school_id, age_group, contact_split$per_student
)
school_contacts <- probability[school_id[1:6000] + 1L] * 199
expect_equal(
  rowSums(contact_split$community) +
    attendance * c(mean(school_contacts), 0),
  rowSums(total)
)
expect_equal(
  contact_split$community * c(8000, 12000),
  t(contact_split$community * c(8000, 12000))
)

p$transmission_rate <- transmission_for_r0(p, total)
seed_r0 <- function(setup) {
  model <- build_model(p, setup, new.env(parent = emptyenv()))
  run_multiple(
    m        = model,
    ndays    = p$ndays,
    nsims    = 40,
    seed     = 1,
    saver    = make_saver("reproductive"),
    nthreads = 1,
    verbose  = FALSE
  )
  r <- run_multiple_get_results(model, nthreads = 1)$reproductive
  sum(r$rt[r$source >= 0 & r$source_exposure_date == 0]) / (40 * p$n_seed)
}
layered <- seed_r0(list(
  persons            = list(age_group = age_group, school = school_id),
  community          = contact_split$community,
  school_probability = probability,
  vaccinated         = integer(),
  kernel             = rep(1 / 19, 19)
))
pooled <- seed_r0(list(
  persons            = list(
    age_group = age_group,
    school    = rep(NA_integer_, p$n)
  ),
  community          = total,
  school_probability = rep(0, 30),
  vaccinated         = integer(),
  kernel             = rep(1 / 19, 19)
))
expect_true(abs(layered / pooled - 1) < 0.15)
expect_true(abs(layered / p$R0 - 1) < 0.2)

# Non-student coverage follows the age rules, and older agents count as
# vaccinated.
coverage <- community_coverage(
  age         = c(0, 3, 20, 30, 50, 69, 70),
  elementary  = 0.9,
  high_school = 0.8,
  adult       = 0.95,
  immune_age  = 70
)
expect_equal(coverage, c(0, 0.9, 0.8, 0.875, 0.95, 0.95, 1))

# Every scenario in scenarios.yml, and in the Makefile, loads with known
# parameters; misspelled parameters stop.
scenario_names <- names(yaml::read_yaml("../scenarios.yml"))
for (name in scenario_names) {
  scenario <- read_scenario("../scenarios.yml", name, defaults)
  expect_true(scenario$response %in% c("none", "cases", "wastewater"))
  expect_equal(names(scenario$p), names(defaults))
}
makefile <- readLines("../Makefile")
targets <- grep("^SCENARIOS := ", makefile, value = TRUE)
targets <- sub("^SCENARIOS := ", "", targets)
expect_equal(strsplit(targets, " +")[[1]], scenario_names)
bad <- tempfile(fileext = ".yml")
writeLines(c("typo:", "  ww_lab_dealy: 2"), bad)
expect_error(read_scenario(bad, "typo", defaults), "ww_lab_dealy")
expect_error(read_scenario(bad, "missing", defaults), "not defined")

# Detection plus hospitalization reports the target share of cases.
q <- p
q$n <- 4000
q$n_seed <- q$n
q$ndays <- 40
q$transmission_rate <- 0
setup_r <- setup
setup_r$persons <- list(
  age_group = rep(0L, q$n),
  school    = rep(NA_integer_, q$n)
)
setup_r$community <- matrix(1, 1, 1)
q$detection_rate <- rash_detection_rate(
  reporting       = 0.25,
  hospitalization = rash_hospitalization_rate(
    probability = q$prob_hospitalization,
    mean        = q$rash_period,
    sd          = q$rash_period_sd
  ),
  mean            = q$rash_period,
  sd              = q$rash_period_sd
)
state <- new.env(parent = emptyenv())
model <- build_model(q, setup_r, state)
run(model, ndays = q$ndays, seed = 21)
detected <- sum(state$runs[[1]]$daily$new_detected)
expect_true(abs(detected / q$n - 0.25) < 0.03)
expect_error(
  rash_detection_rate(0.25, hospitalization = 0.2, mean = 4, sd = 1),
  "Hospitalizations alone"
)

# Shedding that starts before prodrome onset moves the signal earlier.
q <- p
q$n <- length(setup$persons$age_group)
q$n_seed <- 50
q$ndays <- 40
q$transmission_rate <- 0
q$ww_sampling_days <- 0:6
q$ww_sigma <- 0
q$ww_false_positive <- 0
first_positive <- function(start) {
  q$ww_t_start <- start
  setup$kernel <- rep(1 / (start + 19), start + 19)
  state <- new.env(parent = emptyenv())
  run(build_model(q, setup, state), ndays = q$ndays, seed = 4)
  which(state$runs[[1]]$daily$ww_positive)[1]
}
expect_true(first_positive(3) < first_positive(0))
q$ww_t_start <- q$ww_lab_delay + 1
expect_error(build_model(q, setup, new.env(parent = emptyenv())))

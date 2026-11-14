library(epiworldR)

#' Build a measles model with community and school layers
#'
#' Latent, prodromal and rash durations are rounded lognormal draws. The
#' incubation standard deviation is informed by Lessler et al. (2009), the
#' prodromal standard deviation by the Hagelloch data, and the rash standard
#' deviation is an assumption. Each day of rash, a case is detected with
#' probability "Detection rate" and isolated; isolation keeps the rash clock
#' that started at rash onset. Vaccination and contact assignments are given
#' in `setup`; the response event records each run in `state$runs`.
#' @param p Disease, surveillance and response parameters.
#' @param setup Population, contacts, vaccination and shedding inputs.
#' @param state Environment that receives one record per simulation.
#' @param response `"none"` (record only), `"cases"` (exclude unvaccinated
#'   students from schools with reported cases), or `"wastewater"` (the case
#'   response plus a county-wide alert after positive samples).
#' @param seeds Agent ids (from 0) infected on day 0; `NULL` draws `p$n_seed`
#'   agents at random in every run.
#' @return An epiworldR model ready for [epiworldR::run()] or
#'   [epiworldR::run_multiple()].
#' @references Lessler et al. (2009), doi:10.1016/S1473-3099(09)70069-6.
#'   Meyer et al. (2017), doi:10.18637/jss.v077.i11.
#' @examples
#' parameters <- read.csv("data/parameters.csv", stringsAsFactors = FALSE)
#' p <- lapply(parameters$default, function(x) eval(parse(text = x)))
#' names(p) <- parameters$r_argument
#' p$n <- 100
#' p$ndays <- 30
#' p$transmission_rate <- 0.1
#' setup <- list(
#'   persons = list(age_group = rep(0L, p$n), school = rep(0L, p$n)),
#'   community = matrix(2, 1, 1),
#'   school_probability = 0.05,
#'   vaccinated = integer(),
#'   kernel = rep(1 / 19, 19)
#' )
#' model <- build_model(p, setup, new.env(parent = emptyenv()))
#' get_stages(model)
build_model <- function(
    p,
    setup,
    state,
    response = "none",
    seeds    = NULL
) {
  population <- setup$persons
  stopifnot(
    response %in% c("none", "cases", "wastewater"),
    length(population$age_group) == p$n,
    length(population$school) == p$n
  )
  rash_duration <- duration_lognormal(
    mean = p$rash_period,
    sd   = p$rash_period_sd
  )

  # Define disease stages and their progression using the model builder.
  virus <- virus(name = "Measles") |>
    add_stage(label = "Susceptible") |>
    add_stage(label = "Latent", carrier = TRUE) |>
    add_stage(label = "Prodromal", infectious = TRUE) |>
    add_stage(label = "Rash", infectious = TRUE) |>
    add_stage(label = "Rash isolated", infectious = TRUE) |>
    add_stage(label = "Hospitalized", carrier = TRUE) |>
    add_stage(label = "Recovered") |>
    set_contact(stage = "Rash isolated", weight = 0) |>
    set_contact(stage = "Hospitalized", weight = 0) |>
    add_infection(
      from = "Susceptible",
      to   = "Latent",
      rate = p$transmission_rate
    ) |>
    add_after(
      from     = "Latent",
      to       = "Prodromal",
      duration = duration_lognormal(
        mean = p$incubation_period,
        sd   = p$incubation_period_sd
      )
    ) |>
    add_after(
      from     = "Prodromal",
      to       = "Rash",
      duration = duration_lognormal(
        mean = p$prodromal_period,
        sd   = p$prodromal_period_sd
      )
    )

  # Detection isolates a case without restarting its rash period.
  virus <- virus |>
    add_after(from = "Rash", to = "Recovered", duration = rash_duration) |>
    add_after(
      from     = "Rash isolated",
      to       = "Recovered",
      duration = rash_duration,
      since    = "Rash"
    ) |>
    add_rate(from = "Rash", to = "Rash isolated", rate = "Detection rate") |>
    add_rate(
      from = c("Rash", "Rash isolated"),
      to   = "Hospitalized",
      rate = rash_hospitalization_rate(
        probability = p$prob_hospitalization,
        mean        = p$rash_period,
        sd          = p$rash_period_sd
      )
    ) |>
    add_rate(
      from = "Hospitalized",
      to   = "Recovered",
      rate = 1 / p$hospitalization_period
    )

  # Add the two contact layers; the response event edits both during a run.
  community <- layer_mixing(
    name           = "community",
    groups         = as.integer(population$age_group),
    contact_matrix = setup$community
  )
  schools <- layer_mixing(
    name         = "schools",
    groups       = as.integer(population$school),
    contact_prob = setup$school_probability
  )
  placement <- if (is.null(seeds)) {
    assign_count(n = p$n_seed)
  } else {
    assign_agents(ids = as.integer(seeds))
  }
  model <- Model(n = p$n) |>
    add_param("Detection rate", p$detection_rate) |>
    add_layer(layer = community) |>
    add_layer(layer = schools) |>
    add_virus(virus = virus, assign = placement)

  # Apply MMR protection to the supplied agents.
  if (length(setup$vaccinated)) {
    vaccine <- tool(
      name                     = "MMR",
      susceptibility_reduction = p$vax_efficacy,
      transmission_reduction   = 0,
      recovery_enhancer        = 0,
      death_reduction          = 0
    )
    model <- model |>
      add_tool(
        tool   = vaccine,
        assign = assign_agents(ids = as.integer(setup$vaccinated))
      )
  }

  # Attach surveillance and response without capturing the model.
  model |>
    add_globalevent(
      event = response_event(
        p        = p,
        setup    = setup,
        response = response,
        state    = state,
        n_seed   = if (is.null(seeds)) p$n_seed else length(seeds)
      )
    )
}

#' Group students into schools of a target size
#'
#' Students of each level are split at random into schools of about `size`.
#' Each simulated school takes the attributes (e.g., MMR coverage) of a source
#' school of the same level, drawn with probability `weight`, without
#' replacement unless the level has fewer source schools than needed.
#' @param level School level of each agent; `NA` for non-students.
#' @param source_level Level of each source school.
#' @param weight Sampling weight of each source school.
#' @param size Target number of students per school.
#' @return A list with `school`, the simulated school of each agent (from 0;
#'   `NA` for non-students), and `source`, the source school of each
#'   simulated school.
#' @examples
#' schools <- assign_schools(
#'   level        = rep(c(1, 2, NA), c(1000, 600, 400)),
#'   source_level = c(1, 1, 1, 2, 2),
#'   weight       = c(1, 2, 1, 1, 1),
#'   size         = 500
#' )
#' table(schools$school)
#' schools$source
assign_schools <- function(
    level,
    source_level,
    weight,
    size = 500
) {
  school <- rep(NA_integer_, length(level))
  source <- integer()
  for (j in sort(unique(level[!is.na(level)]))) {
    students <- which(level == j)
    choices <- which(source_level == j & !is.na(weight) & weight > 0)
    if (!length(choices)) stop("No source school for level ", j)
    k <- max(1L, round(length(students) / size))

    # Draw source schools, then deal the students evenly among them.
    picked <- choices[sample.int(
      n       = length(choices),
      size    = k,
      replace = k > length(choices),
      prob    = weight[choices]
    )]
    school[students] <- length(source) +
      sample(rep_len(seq_len(k), length(students))) - 1L
    source <- c(source, picked)
  }
  list(school = school, source = source)
}

#' Split age-group contacts between the community and school layers
#'
#' The school layer takes the within-age school contacts of every age group
#' that attends school. All other contacts stay in the community matrix,
#' including school contacts between age groups (e.g., with teachers).
#' Averaged over an age group, community plus school contacts equal the row
#' sums of `total`. Removing a diagonal keeps the matrix reciprocal.
#' @param total Daily contacts across all settings, by age group.
#' @param school Daily school-setting contacts, already included in `total`.
#' @param attendance Share of each age group that attends school.
#' @return A list with the `community` matrix and the expected daily school
#'   contacts `per_student` of a student in each age group.
#' @examples
#' total <- matrix(c(8, 2, 2, 6), 2)
#' school <- matrix(c(4, 0.5, 0.5, 0), 2)
#' split_school_contacts(total, school, attendance = c(0.8, 0))
split_school_contacts <- function(total, school, attendance) {
  stopifnot(
    dim(total) == dim(school),
    length(attendance) == nrow(total),
    attendance >= 0,
    attendance <= 1
  )

  # Only groups with students move their within-age school contacts.
  within <- ifelse(attendance > 0, diag(school), 0)
  list(
    community   = total - diag(within, nrow = length(within)),
    per_student = ifelse(attendance > 0, within / attendance, 0)
  )
}

#' Daily contact probability between members of each school
#'
#' Every pair of schoolmates meets with the same daily probability, chosen so
#' that a student's expected school contacts equal the school's mean of
#' `per_student`. Schools with one student have no school contacts.
#' @param school School of each agent, starting at 0; `NA` for non-students.
#' @param age_group Age group of each agent, starting at 0.
#' @param per_student Expected daily school contacts by age group.
#' @return One contact probability per school.
#' @examples
#' school_contact_prob(
#'   school      = c(0, 0, 0, 1, 1, NA),
#'   age_group   = c(1, 1, 1, 2, 2, 3),
#'   per_student = c(0, 1, 0.5, 0)
#' )
school_contact_prob <- function(school, age_group, per_student) {
  stopifnot(length(school) == length(age_group))

  # Average the target over each school's students.
  students <- which(!is.na(school))
  groups <- factor(
    school[students],
    levels = seq_len(max(school[students]) + 1L) - 1L
  )
  size <- tabulate(groups, nbins = nlevels(groups))
  target <- tapply(per_student[age_group[students] + 1L], groups, mean)
  target[is.na(target)] <- 0

  pmin(1, ifelse(size > 1, as.numeric(target) / (size - 1), 0))
}

#' MMR coverage of agents outside school, by single year of age
#'
#' Follows the adult immunization assumptions of the multigroup.vaccine
#' vignette "Utah/Arizona border": children under 1 are unvaccinated, ages
#' 1-4 have the elementary-school coverage, ages 5-24 the high-school
#' coverage, ages 25-44 the mean of the high-school and adult coverage, and
#' ages 45 and older the adult coverage. Agents aged `immune_age` and older
#' count as vaccinated, a stand-in for presumptive immunity (in the United
#' States, birth before 1957). The model applies the vaccine efficacy.
#' @param age Age in years.
#' @param elementary,high_school MMR coverage in elementary and high schools.
#' @param adult Coverage of adults aged 45 to `immune_age`.
#' @param immune_age Age from which agents count as vaccinated.
#' @return Share vaccinated for each `age`.
#' @references multigroup.vaccine vignette "utaz_age_school", commit
#'   ddd9201f28b9908c033542f6a2445b84619af1c7.
#' @examples
#' community_coverage(c(0, 3, 20, 30, 50, 80), elementary = 0.9,
#'   high_school = 0.86)
community_coverage <- function(
    age,
    elementary,
    high_school,
    adult      = 0.95,
    immune_age = 70
) {
  stopifnot(
    age >= 0,
    elementary >= 0,
    elementary <= 1,
    high_school >= 0,
    high_school <= 1,
    adult >= 0,
    adult <= 1
  )

  # Assign coverage from the oldest band down to the youngest.
  coverage <- rep(adult, length(age))
  coverage[age < 45] <- (high_school + adult) / 2
  coverage[age < 25] <- high_school
  coverage[age < 5] <- elementary
  coverage[age < 1] <- 0
  coverage[age >= immune_age] <- 1
  coverage
}

#' Read default parameters
#'
#' @param file CSV with columns `r_argument` and `default`, where `default`
#'   holds an R expression.
#' @return A named list of parameter values.
#' @examples
#' p <- read_parameters("data/parameters.csv")
#' p$R0
read_parameters <- function(file) {
  parameters <- utils::read.csv(file, stringsAsFactors = FALSE)
  p <- lapply(parameters$default, function(x) eval(parse(text = x)))
  names(p) <- parameters$r_argument
  p
}

#' Read one scenario and apply it to the default parameters
#'
#' A scenario file is a YAML mapping from scenario names to their settings:
#' a `description`, a `response` (`"none"`, `"cases"` or `"wastewater"`), and
#' any parameters that differ from `defaults`. Unknown parameter names stop
#' with an error, so a misspelled parameter cannot be silently ignored.
#' @param file Path to the YAML file.
#' @param name Scenario name.
#' @param defaults Named list of default parameters.
#' @return A list with the scenario `name`, `description`, `response`, and
#'   parameters `p`.
#' @examples
#' scenario <- read_scenario(
#'   file     = "scenarios.yml",
#'   name     = "baseline",
#'   defaults = read_parameters("data/parameters.csv")
#' )
#' scenario$response
read_scenario <- function(file, name, defaults) {
  scenarios <- yaml::read_yaml(file)
  if (!name %in% names(scenarios)) {
    stop("Scenario '", name, "' is not defined in ", file)
  }
  spec <- scenarios[[name]]
  settings <- c("description", "response")
  unknown <- setdiff(names(spec), c(settings, names(defaults)))
  if (length(unknown)) {
    stop("Unknown parameters in '", name, "': ", toString(unknown))
  }
  response <- if (is.null(spec$response)) "none" else spec$response
  stopifnot(response %in% c("none", "cases", "wastewater"))

  # Override the defaults with the scenario's parameters.
  p <- defaults
  for (key in setdiff(names(spec), settings)) p[[key]] <- unlist(spec[[key]])
  list(
    name        = name,
    description = if (is.null(spec$description)) "" else spec$description,
    response    = response,
    p           = p
  )
}

#' Per-contact transmission probability for a target R0
#'
#' Uses the next-generation approximation `R0 = p * D * rho(C)`, where `C` is
#' the population-average daily contact matrix and `D` the expected days an
#' agent transmits: the prodromal period plus the expected rash days before
#' detection, `sum((1 - h)^(0:(L - 1)))` for daily detection probability `h`
#' and rash period `L`. It ignores hospitalization, clustering within schools,
#' and susceptible depletion, so simulations should refine the result.
#' @param p Disease parameters, including `R0` and `detection_rate`.
#' @param contacts Daily contacts by age group across all layers.
#' @return Probability of transmission per contact.
#' @examples
#' p <- list(
#'   R0               = 7,
#'   prodromal_period = 4,
#'   rash_period      = 4,
#'   detection_rate   = 0.5
#' )
#' transmission_for_r0(p, matrix(10, 1, 1))
transmission_for_r0 <- function(p, contacts) {
  rash_days <- seq_len(max(1, round(p$rash_period))) - 1
  days <- p$prodromal_period + sum((1 - p$detection_rate)^rash_days)
  r1 <- compute_reproduction_number(
    contacts,
    transmission_prob      = 1,
    infectious_period_days = days
  )$R
  p$R0 / r1
}

#' Convert a target hospitalization probability to a daily rash hazard
#'
#' Assumes the rounded lognormal rash duration is independent of the daily
#' hospitalization draw. Numerical summation includes the far tail.
#' @param probability Probability of hospitalization during rash.
#' @param mean,sd Mean and standard deviation of rash duration in days.
#' @return Daily hospitalization probability.
#' @examples
#' rash_hospitalization_rate(probability = 0.2, mean = 4, sd = 1)
rash_hospitalization_rate <- function(probability, mean, sd) {
  stopifnot(probability >= 0, probability < 1, mean > 0, sd >= 0)

  # Handle degenerate distributions before numerical integration.
  if (probability == 0) return(0)
  if (sd == 0) return(1 - (1 - probability)^(1 / max(1, round(mean))))

  # Integrate the rounded lognormal duration to recover a daily hazard.
  s2 <- log1p((sd / mean)^2)
  days <- seq_len(1000L)
  mass <- diff(c(0, stats::plnorm(days + .5, log(mean) - s2 / 2, sqrt(s2))))
  stats::uniroot(
    f        = function(h) 1 - sum(mass * (1 - h)^days) - probability,
    interval = c(0, 1)
  )$root
}

#' Daily rash detection probability for a target reporting fraction
#'
#' A case is reported if it is detected or hospitalized during its rash,
#' so the daily detection probability `h` solves
#' `1 - E[((1 - h) (1 - r))^D] = reporting`, where `r` is the daily
#' hospitalization probability and `D` the rounded lognormal rash duration.
#' Competing daily events are treated as independent.
#' @param reporting Target share of cases that are reported.
#' @param hospitalization Daily hospitalization probability during rash.
#' @param mean,sd Mean and standard deviation of rash duration in days.
#' @return Daily detection probability.
#' @examples
#' r <- rash_hospitalization_rate(probability = 0.02, mean = 4, sd = 1)
#' rash_detection_rate(reporting = 0.25, hospitalization = r, mean = 4,
#'   sd = 1)
rash_detection_rate <- function(reporting, hospitalization, mean, sd) {
  stopifnot(reporting > 0, reporting < 1, mean > 0, sd > 0)

  # Rounded lognormal rash duration, as in rash_hospitalization_rate().
  s2 <- log1p((sd / mean)^2)
  days <- seq_len(1000L)
  mass <- diff(c(0, stats::plnorm(days + .5, log(mean) - s2 / 2, sqrt(s2))))
  reported <- function(h) {
    1 - sum(mass * ((1 - h) * (1 - hospitalization))^days)
  }
  if (reported(0) >= reporting) {
    stop("Hospitalizations alone exceed the reporting fraction")
  }
  stats::uniroot(
    f        = function(h) reported(h) - reporting,
    interval = c(0, 1)
  )$root
}

#' Create the surveillance and response event
#'
#' Every response records the same surveillance signals, so paired runs stay
#' identical until a response acts. Detected cases (entries into isolation or
#' hospital) are reported after `case_report_delay` days. In the `"cases"`
#' and `"wastewater"` responses, a reported student case excludes the
#' school's unvaccinated students who have not been detected as cases;
#' they return `response_days` after the school's last reported case. Any
#' reported case also starts or extends a county-wide public health alert of
#' `response_days` days, which raises the daily detection probability to
#' `alert_detection_rate`.
#' Wastewater is sampled on `ww_sampling_days` of each week. Agents shed
#' from `ww_t_start` days before prodrome onset, following `setup$kernel`
#' (one entry per day from that start); because prodrome onsets are known
#' only when they happen, a sample is measured `ww_t_start` days after it is
#' taken, which requires `ww_t_start <= ww_lab_delay`. A sample is
#' positive when its concentration reaches `trigger_ww` or, with probability
#' `ww_false_positive`, falsely. In the `"wastewater"` response, a positive
#' result, available after `ww_lab_delay` days, also starts or extends the
#' alert; every response records the samples. R events require
#' single-threaded runs.
#' @param p Disease, surveillance and response parameters.
#' @param setup Population, vaccination and shedding inputs.
#' @param response `"none"`, `"cases"`, or `"wastewater"`.
#' @param state Environment that receives the replicate records.
#' @param n_seed Number of agents infected on day 0.
#' @return An epiworldR global event.
#' @examples
#' p <- list(
#'   ndays           = 30,
#'   ww_log10_g_mean = 10,
#'   ww_log10_g_sd   = 1,
#'   ww_t_start      = 0,
#'   ww_lab_delay    = 3
#' )
#' setup <- list(persons = list(school = NA))
#' event <- response_event(p, setup, "none", new.env(), n_seed = 1)
#' class(event)
response_event <- function(p, setup, response, state, n_seed) {
  stopifnot(
    p$ww_t_start >= 0,
    p$ww_t_start <= p$ww_lab_delay
  )

  state$runs <- list()
  n_days  <- p$ndays
  meanlog <- log(10) * p$ww_log10_g_mean
  sdlog   <- log(10) * p$ww_log10_g_sd
  school  <- setup$persons$school
  respond <- response != "none"

  # Unvaccinated students of each school, by agent id (from 0).
  n_schools <- max(c(school, -1L), na.rm = TRUE) + 1L
  keep <- !is.na(school) & !(seq_along(school) - 1L) %in% setup$vaccinated
  pupils <- split(
    which(keep) - 1L,
    factor(school[keep], levels = seq_len(n_schools) - 1L)
  )

  reset <- function() {
    daily <- c(
      "new_infections", "new_rash", "new_detected", "reported",
      "admissions", "hospitalized", "isolated", "excluded", "cohort_g"
    )
    for (v in daily) assign(v, numeric(n_days), envir = state)
    state$ever_hospital <- rep(FALSE, length(school))
    state$alert_on    <- logical(n_days)
    state$ww_obs      <- rep(NA_real_, n_days)
    state$ww_positive <- rep(NA, n_days)
    state$ww_false    <- rep(NA, n_days)
    state$queue       <- vector("list", n_days + p$case_report_delay)
    state$detected    <- rep(FALSE, length(school))
    state$excluded_ids <- list()
    state$exclusion_end <- numeric()
    state$episodes    <- 0L
    state$first_exclusion <- NA_integer_
    state$alert_end   <- 0
    state$acting      <- FALSE
    state$alerts      <- 0L
    state$first_alert <- NA_integer_
    state$first_trigger <- NA_character_
    state$prev <- c(inf = n_seed, prod = 0, rash = 0)
  }

  fun <- function(m) {

    t <- today(m)
    if (t == 1) reset()

    # Daily stage totals; agents are scanned only while cases are detected.
    totals <- get_today_total(m)
    tot <- tapply(totals$counts, totals$stage, sum)
    stages <- NULL
    if (tot[["Rash isolated"]] + tot[["Hospitalized"]] > 0) {
      stages <- get_agents_stages(m)
      found <- which(
        stages %in% c("Rash isolated", "Hospitalized") & !state$detected
      )
      state$detected[found] <- TRUE
      state$new_detected[t] <- length(found)
      admitted <- which(stages == "Hospitalized" & !state$ever_hospital)
      state$ever_hospital[admitted] <- TRUE
      state$admissions[t] <- length(admitted)
      r <- t + p$case_report_delay
      state$queue[[r]] <- c(state$queue[[r]], found - 1L)
    }

    # Cumulative infections, prodromal onsets, and rash onsets.
    cum <- c(inf = p$n - tot[["Susceptible"]])
    cum["prod"] <- cum[["inf"]] - tot[["Latent"]]
    cum["rash"] <- cum[["prod"]] - tot[["Prodromal"]]
    new <- cum - state$prev
    state$prev <- cum
    state$new_infections[t] <- new[["inf"]]
    state$new_rash[t]       <- new[["rash"]]
    state$hospitalized[t]   <- tot[["Hospitalized"]]
    state$isolated[t]       <- tot[["Rash isolated"]]

    # Wastewater: shedding cohorts, sampling error, and false positives.
    # Shedding starts ww_t_start days before prodrome onset, so the sample
    # of day u = t - ww_t_start is measured on day t, once every cohort
    # shedding into it is known; kernel entry k + 1 weights onset day t - k.
    state$cohort_g[t] <- sum(stats::rlnorm(new[["prod"]], meanlog, sdlog))
    lag <- seq_len(min(t, length(setup$kernel))) - 1L
    conc <- sum(setup$kernel[lag + 1L] * state$cohort_g[t - lag]) /
      (p$ww_alpha * p$n)
    eps <- stats::rnorm(1, 0, p$ww_sigma)
    false_draw <- stats::runif(1) < p$ww_false_positive
    u <- t - p$ww_t_start
    if (u >= 1 && (u %% 7) %in% p$ww_sampling_days) {
      obs <- conc * exp(eps)
      obs <- if (obs >= p$ww_lod) obs else 0
      true_positive <- obs >= p$trigger_ww
      state$ww_obs[u]      <- obs
      state$ww_positive[u] <- true_positive || false_draw
      state$ww_false[u]    <- !true_positive && false_draw
    }

    # A reported case, or in the wastewater response a positive result,
    # starts or extends the alert. Actions at the end of day t apply from
    # day t + 1.
    reported <- state$queue[[t]]
    state$reported[t] <- length(reported)
    s <- t - p$ww_lab_delay
    ww_trigger <- response == "wastewater" && s >= 1 &&
      isTRUE(state$ww_positive[s])
    case_trigger <- respond && length(reported) > 0
    if (case_trigger || ww_trigger) {
      if (t >= state$alert_end) state$alerts <- state$alerts + 1L
      if (is.na(state$first_alert)) {
        state$first_alert <- t
        state$first_trigger <- if (case_trigger) "case" else "wastewater"
      }
      state$alert_end <- t + p$response_days
    }
    on <- t < state$alert_end
    state$alert_on[t] <- on
    if (on != state$acting) {
      set_param(
        m,
        "Detection rate",
        if (on) p$alert_detection_rate else p$detection_rate
      )
      state$acting <- on
    }

    # Reported student cases exclude their school's unvaccinated students.
    if (respond) {
      layer <- get_layer(m, "schools")
      for (s_id in unique(stats::na.omit(school[reported + 1L]))) {
        key <- as.character(s_id)
        if (is.null(state$excluded_ids[[key]])) {
          ids <- pupils[[s_id + 1L]]
          ids <- ids[!state$detected[ids + 1L]]
          for (i in ids) set_group(layer, i, NA)
          state$excluded_ids[[key]] <- ids
          state$episodes <- state$episodes + 1L
          if (is.na(state$first_exclusion)) state$first_exclusion <- t
        }
        state$exclusion_end[key] <- t + p$response_days
      }

      # Students return once their school's exclusion period ends.
      for (key in names(state$excluded_ids)) {
        if (t >= state$exclusion_end[[key]]) {
          for (i in state$excluded_ids[[key]]) {
            set_group(layer, i, as.integer(key))
          }
          state$excluded_ids[[key]] <- NULL
        }
      }
      state$excluded[t] <- sum(lengths(state$excluded_ids))
    }

    # The final day saves the replicate, including age-specific infections.
    if (t == n_days) {
      if (is.null(stages)) stages <- get_agents_stages(m)
      ever <- stages != "Susceptible"
      active <- setdiff(names(tot), c("Susceptible", "Recovered"))
      state$runs[[length(state$runs) + 1L]] <- list(
        daily = data.frame(
          day            = seq_len(n_days),
          new_infections = state$new_infections,
          new_rash       = state$new_rash,
          new_detected   = state$new_detected,
          reported       = state$reported,
          admissions     = state$admissions,
          hospitalized   = state$hospitalized,
          isolated       = state$isolated,
          excluded       = state$excluded,
          alert_on       = state$alert_on,
          ww_obs         = state$ww_obs,
          ww_positive    = state$ww_positive,
          ww_false       = state$ww_false
        ),
        infections_by_age = tabulate(
          setup$persons$age_group[ever] + 1L,
          nbins = max(setup$persons$age_group) + 1L
        ),
        active_at_end   = sum(tot[active]) > 0,
        episodes        = state$episodes,
        first_exclusion = state$first_exclusion,
        alerts          = state$alerts,
        first_alert     = state$first_alert,
        first_trigger   = state$first_trigger
      )
    }

    invisible(NULL)

  }

  globalevent_fun(fun, name = "Surveillance and response", day = -99)
}

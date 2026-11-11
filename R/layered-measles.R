#' Construct the two-layer measles model
#'
#' Latent and prodromal periods are lognormal draws on entry, parameterized
#' by the mean and standard deviation in days. Lessler et al. (2009) describe
#' exposure to first symptoms, not every disease stage. The prodromal standard
#' deviation comes from the Hagelloch data fit; the rash value is an assumption.
#' The engine rounds durations to whole days, with a minimum of one.
#' @param p Parameters from [measles_ww_params()].
#' @param population Agent age groups and school memberships.
#' @param community_matrix Community contacts by age group.
#' @param school_probability Within-school contact probabilities.
#' @param vaccinated Zero-based IDs receiving the MMR tool.
#' @return An epiworldR model with community and school layers.
#' @references Lessler et al. (2009), doi:10.1016/S1473-3099(09)70069-6.
#'   Meyer et al. (2017), doi:10.18637/jss.v077.i11.
layered_measles <- function(p, population, community_matrix, school_probability,
                            vaccinated = integer()) {
  stopifnot(length(population$age_group) == p$n,
            length(population$school) == p$n)
  virus <- epiworldR::virus("Measles")
  epiworldR::add_stage(virus, "Susceptible")
  epiworldR::add_stage(virus, "Latent", carrier = TRUE)
  epiworldR::add_stage(virus, "Prodromal", infectious = TRUE)
  epiworldR::add_stage(virus, "Rash", infectious = TRUE)
  epiworldR::add_stage(virus, "Rash isolated", infectious = TRUE)
  epiworldR::add_stage(virus, "Hospitalized", carrier = TRUE)
  epiworldR::add_stage(virus, "Recovered")
  # A random fraction of rash cases leave all layers. This matches the
  # requested average reduction but is not identical to fractional contacts
  # for every case.
  epiworldR::set_contact(virus, "Rash isolated", 0)
  epiworldR::set_contact(virus, "Hospitalized", 0)
  epiworldR::add_infection(virus, "Susceptible", "Latent", p$transmission_rate)
  epiworldR::add_after(virus, "Latent", "Prodromal",
    epiworldR::duration_lognormal(p$incubation_period, p$incubation_period_sd))
  epiworldR::add_after(virus, "Prodromal", "Rash",
    epiworldR::duration_lognormal(p$prodromal_period, p$prodromal_period_sd))
  for (stage in c("Rash", "Rash isolated"))
    epiworldR::add_after(virus, stage, "Recovered",
      epiworldR::duration_lognormal(p$rash_period, p$rash_period_sd))
  epiworldR::add_external(virus, "Rash", "Rash isolated", "Isolate rash")
  epiworldR::add_rate(virus, c("Rash", "Rash isolated"), "Hospitalized",
    rash_hospitalization_rate(p$prob_hospitalization, p$rash_period, p$rash_period_sd))
  epiworldR::add_rate(virus, "Hospitalized", "Recovered", 1 / p$hospitalization_period)

  model <- epiworldR::Model(p$n)
  epiworldR::add_layer(model, epiworldR::layer_mixing("community",
    groups = as.integer(population$age_group), contact_matrix = community_matrix))
  epiworldR::add_layer(model, epiworldR::layer_mixing("schools",
    groups = as.integer(population$school), contact_prob = school_probability))
  epiworldR::add_virus(model, virus, epiworldR::assign_count(p$n_seed))
  epiworldR::add_globalevent(model, rash_isolation_event(p$n, p$rash_reduction_contact_rate))
  if (length(vaccinated)) {
    vaccine <- epiworldR::tool("MMR", p$vax_efficacy, 0, 0, 0)
    epiworldR::add_tool(model, vaccine, epiworldR::assign_agents(as.integer(vaccinated)))
  }
  model
}

#' Create the rash-isolation event
#'
#' An eligible case is selected once. The event applies its transition at the
#' end of the day and the isolated stage draws a new rash duration; this is an
#' approximation to fractional contact reduction throughout rash.
#' @param n Number of agents.
#' @param probability Probability of selecting a rash case for isolation.
#' @return An epiworldR global event.
rash_isolation_event <- function(n, probability) {
  was_rash <- rep(FALSE, n)
  isolate <- function(m) {
    if (epiworldR::today(m) == 1) was_rash[] <<- FALSE
    current <- epiworldR::get_agents_stages(m) == "Rash"
    new <- which(current & !was_rash)
    was_rash[new] <<- TRUE
    chosen <- new[stats::runif(length(new)) < probability]
    for (i in chosen) epiworldR::fire(m, "Isolate rash", i - 1L)
  }
  epiworldR::globalevent_fun(isolate, name = "Rash isolation", day = -99)
}

#' Convert a target hospitalization probability to a daily rash hazard
#'
#' Assumes the rounded lognormal rash duration is independent of the daily
#' hospitalization draw. Numerical summation includes the far tail.
#' @param probability Probability of hospitalization during rash.
#' @param mean,sd Mean and standard deviation of rash duration in days.
#' @return Daily hospitalization probability.
rash_hospitalization_rate <- function(probability, mean, sd) {
  stopifnot(probability >= 0, probability < 1, mean > 0, sd >= 0)
  if (probability == 0) return(0)
  if (sd == 0) return(1 - (1 - probability)^(1 / max(1, round(mean))))
  s2 <- log1p((sd / mean)^2)
  days <- seq_len(1000L)
  mass <- diff(c(0, stats::plnorm(days + .5, log(mean) - s2 / 2, sqrt(s2))))
  stats::uniroot(function(h) 1 - sum(mass * (1 - h)^days) - probability,
                 c(0, 1))$root
}

#' Prepare agents and contacts for the layered model
#'
#' School contact probabilities target `school_contacts` per student where
#' group size permits; they are capped at one contact per pair. Student
#' vaccination uses school coverage and other agents use `prop_vaccinated`.
#' @param p Parameters from [measles_ww_params()].
#' @param school_contacts Target daily contacts per student.
#' @param community_scale Multiplier for the initial community matrix.
#' @param seed Seed for population and vaccination assignment.
#' @return A list of population, contact and shedding inputs.
layered_setup <- function(p, school_contacts = 4, community_scale = 0.2,
                          seed = 137L) {
  pop <- city_population(p$n)
  persons <- assign_schools(pop, seed)
  cm <- contact_matrices(pop)
  community <- (cm$Household + cm$Work + cm$Other) * community_scale
  sizes <- tabulate(persons$school + 1L, nbins = max(persons$school, na.rm = TRUE) + 1L)
  school_prob <- pmin(1, school_contacts / pmax(1, sizes - 1L))
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, envir = .GlobalEnv))
  set.seed(seed + 1L)
  coverage <- persons$coverage
  coverage[is.na(coverage)] <- p$prop_vaccinated
  vaccinated <- which(stats::runif(p$n) < coverage) - 1L
  list(pop = pop, persons = persons, community = community,
       school_probability = school_prob, vaccinated = vaccinated,
       kernel = shedding_kernel(p$ww_t_peak, p$ww_t_end, p$ww_log10_vpeak))
}

#' Run one layered scenario
#' @param p Parameters from [measles_ww_params()].
#' @param setup A layered setup.
#' @param scenario Baseline, wastewater, or cases.
#' @param seed Random seed.
#' @return The output of [simulate_once()].
simulate_layered_once <- function(p, setup, scenario = "baseline", seed = 1L) {
  simulate_once(p, setup, scenario, seed)
}

# Measles natural history for epiworldR's v1-layers model builder.
# The caller supplies the population assignments and contact intensities.
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
  epiworldR::add_rate(virus, "Latent", "Prodromal", 1 / p$incubation_period)
  epiworldR::add_rate(virus, "Prodromal", "Rash",
    (1 - p$rash_reduction_contact_rate) / p$prodromal_period)
  epiworldR::add_rate(virus, "Prodromal", "Rash isolated",
    p$rash_reduction_contact_rate / p$prodromal_period)
  epiworldR::add_rate(virus, c("Rash", "Rash isolated"), "Recovered", 1 / p$rash_period)
  epiworldR::add_rate(virus, c("Rash", "Rash isolated"), "Hospitalized",
    p$prob_hospitalization / (p$rash_period * (1 - p$prob_hospitalization)))
  epiworldR::add_rate(virus, "Hospitalized", "Recovered", 1 / p$hospitalization_period)

  model <- epiworldR::Model(p$n)
  epiworldR::add_layer(model, epiworldR::layer_mixing("community",
    groups = as.integer(population$age_group), contact_matrix = community_matrix))
  epiworldR::add_layer(model, epiworldR::layer_mixing("schools",
    groups = as.integer(population$school), contact_prob = school_probability))
  epiworldR::add_virus(model, virus, epiworldR::assign_count(p$n_seed))
  if (length(vaccinated)) {
    vaccine <- epiworldR::tool("MMR", p$vax_efficacy, 0, 0, 0)
    epiworldR::add_tool(model, vaccine, epiworldR::assign_agents(as.integer(vaccinated)))
  }
  model
}

# Build a reproducible county population. School mixing probabilities target
# `school_contacts` per student where group size permits; small simulated
# schools cannot attain that target and are capped at one contact per pair.
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

# Run the existing surveillance process against the two-layer model.
simulate_layered_once <- function(p, setup, scenario = "baseline", seed = 1L) {
  stopifnot(scenario %in% c("baseline", "wastewater", "cases"))
  m <- layered_measles(p, setup$persons, setup$community,
                       setup$school_probability, setup$vaccinated)
  state <- new.env(parent = emptyenv())
  add_surveillance(m, p, setup, scenario, state)
  set.seed(seed)
  epiworldR::run(m, ndays = p$ndays, seed = seed)
  days <- seq_len(p$ndays)
  list(daily = data.frame(day = days, new_inf = state$new_inf,
    new_prod = state$new_prod, new_rash = state$new_rash,
    reported = state$reported[days], ww_true = state$ww_true,
    ww_obs = state$ww_obs, cum_inf = state$cum_inf,
    cum_prod = state$cum_prod, cum_rash = state$cum_rash),
    trigger_day = state$trigger_day)
}

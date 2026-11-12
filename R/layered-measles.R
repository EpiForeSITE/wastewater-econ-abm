library(epiworldR)

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
#' @examples
#' p <- measles_ww_params(n = 1000, ndays = 30, n_seed = 5)
#' inputs <- salt_lake_inputs(n = p$n)
#' setup <- do.call(layered_setup, c(list(p = p), inputs))
#' model <- layered_measles(
#'   p                  = p,
#'   population         = setup$persons,
#'   community_matrix   = setup$community,
#'   school_probability = setup$school_probability,
#'   vaccinated         = setup$vaccinated
#' )
#' epiworldR::get_stages(model)
layered_measles <- function(
    p,
    population,
    community_matrix,
    school_probability,
    vaccinated = integer()
) {
  stopifnot(
    length(population$age_group) == p$n,
    length(population$school) == p$n
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

  # Both rash stages share the recovery and hospitalization mechanisms.
  for (stage in c("Rash", "Rash isolated")) {
    virus <- virus |>
      add_after(
        from     = stage,
        to       = "Recovered",
        duration = duration_lognormal(
          mean = p$rash_period,
          sd   = p$rash_period_sd
        )
      )
  }
  virus <- virus |>
    add_external(
      from            = "Rash",
      to              = "Rash isolated",
      transition_name = "Isolate rash"
    ) |>
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

  # Add the two contact layers, seed infections, and register rash isolation.
  model <- Model(n = p$n) |>
    add_layer(
      layer = layer_mixing(
        name           = "community",
        groups         = as.integer(population$age_group),
        contact_matrix = community_matrix
      )
    ) |>
    add_layer(
      layer = layer_mixing(
        name         = "schools",
        groups       = as.integer(population$school),
        contact_prob = school_probability
      )
    ) |>
    add_virus(virus = virus, assign = assign_count(n = p$n_seed)) |>
    add_globalevent(
      event = rash_isolation_event(
        n           = p$n,
        probability = p$rash_reduction_contact_rate
      )
    )

  # Apply MMR protection to the supplied agents.
  if (length(vaccinated)) {
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
        assign = assign_agents(ids = as.integer(vaccinated))
      )
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
#' @examples
#' event <- rash_isolation_event(n = 1000, probability = 0.5)
#' class(event)
rash_isolation_event <- function(n, probability) {
  # Track rash onsets separately for each simulation run.
  was_rash <- rep(FALSE, n)
  isolate <- function(m) {
    if (today(m) == 1) was_rash[] <<- FALSE

    # Draw isolation only for agents newly observed with rash.
    current <- get_agents_stages(m) == "Rash"
    new <- which(current & !was_rash)
    was_rash[new] <<- TRUE
    chosen <- new[stats::runif(length(new)) < probability]
    for (i in chosen) fire(m, "Isolate rash", i - 1L)
  }

  globalevent_fun(isolate, name = "Rash isolation", day = -99)
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

#' Prepare agents and contacts for the layered model
#'
#' School contact probabilities target `school_contacts` per student where
#' group size permits; they are capped at one contact per pair. Student
#' vaccination uses school coverage and other agents use `prop_vaccinated`.
#' @param p Parameters from [measles_ww_params()].
#' @param school_contacts Target daily contacts per student.
#' @param community_scale Multiplier for the initial community matrix.
#' @param pop Age-group population counts.
#' @param schools School records with id, level, and vaxrate columns.
#' @param weights School sampling weights, in school-record order.
#' @param cm Contact matrices by setting.
#' @return A list of population, contact and shedding inputs.
#' @examples
#' p <- measles_ww_params(n = 1000)
#' inputs <- salt_lake_inputs(n = p$n)
#' setup <- do.call(layered_setup, c(list(p = p), inputs))
#' length(setup$persons$school)
#' head(setup$school_probability)
layered_setup <- function(
    p,
    pop,
    schools,
    weights,
    cm,
    school_contacts = 4,
    community_scale = 0.2
) {
  # Construct agent memberships from the application inputs.
  persons <- assign_schools(
    pop     = pop,
    schools = schools,
    weights = weights
  )

  # Scale community contacts and cap within-school probabilities.
  community <- (cm$Household + cm$Work + cm$Other) * community_scale
  sizes <- tabulate(
    bin   = persons$school + 1L,
    nbins = length(unique(schools$id))
  )
  school_prob <- pmin(1, school_contacts / pmax(1, sizes - 1L))

  # Use the caller's RNG stream to assign vaccination.
  coverage <- persons$coverage
  coverage[is.na(coverage)] <- p$prop_vaccinated
  vaccinated <- which(stats::runif(p$n) < coverage) - 1L

  # Return inputs shared by calibration and the scenario runner.
  list(
    pop                = pop,
    persons            = persons,
    community          = community,
    school_probability = school_prob,
    vaccinated         = vaccinated,
    kernel             = shedding_kernel(
      t_peak      = p$ww_t_peak,
      t_end       = p$ww_t_end,
      log10_vpeak = p$ww_log10_vpeak
    )
  )
}

#' Run one layered scenario
#' @param p Parameters from [measles_ww_params()].
#' @param setup A layered setup.
#' @param scenario Baseline, wastewater, or cases.
#' @param seed Random seed.
#' @return The output of [simulate_once()].
#' @examples
#' p <- measles_ww_params(n = 1000, ndays = 30, n_seed = 5)
#' inputs <- salt_lake_inputs(n = p$n)
#' setup <- do.call(layered_setup, c(list(p = p), inputs))
#' result <- simulate_layered_once(p, setup, scenario = "baseline", seed = 42)
#' head(result$daily)
simulate_layered_once <- function(p, setup, scenario = "baseline", seed = 1L) {
  simulate_once(p, setup, scenario, seed)
}

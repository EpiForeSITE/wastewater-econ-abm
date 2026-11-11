#' Load the current Salt Lake County application inputs
#'
#' This adapter contains location-specific data choices. Other applications
#' supply the same pop, schools, weights, and cm inputs to [layered_setup()].
#' @param n Number of simulated agents.
#' @return A list of inputs for the generic layered setup.
salt_lake_inputs <- function(n) {

  # Read the current application's demographic and school data.
  pop <- city_population(
    n          = n,
    county     = "Salt Lake County",
    state_fips = "49",
    year       = 2025
  )
  schools <- utah_schools()

  list(
    pop     = pop,
    schools = schools,
    weights = utah_school_weights(schools = schools),
    cm      = contact_matrices(pop = pop)
  )
}

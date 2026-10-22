# Population and contact structure for a synthetic city.
#
# Age groups follow the Epistorm-Mix 5-year brackets (0-4, ..., 70-74, 75+).
# The age distribution comes from the U.S. Census county estimates shipped with
# the multigroup.vaccine package, and the contact matrices from Epistorm-Mix
# (data/epistorm). The national matrices are re-weighted to the city's age
# distribution and made reciprocal (symmetric in total contacts).

age_limits <- seq(0, 75, by = 5)

settings <- c("Household", "School", "Work", "Other")

#' Read an Epistorm-Mix setting-specific M-matrix
#' @param setting One of `settings`.
#' @param path Directory with the CSV files.
read_epistorm_matrix <- function(setting, path = "data/epistorm") {

  x <- utils::read.csv(
    file.path(path, sprintf("%s-M-by5_75-matrix.csv", setting)),
    check.names = FALSE
  )

  ans <- as.matrix(x[, -1])
  rownames(ans) <- x[[1]]
  ans

}

#' Group sizes for a synthetic city of size `n`
#'
#' Uses the age distribution of a U.S. county and rounds the group sizes so
#' they add up to `n` (largest remainder).
#' @param n Population size.
#' @param county County name as in the census file.
#' @param state_fips State FIPS code.
#' @param year Census estimate year.
#' @return A named integer vector with one entry per age group.
city_population <- function(
    n          = 1e4,
    county     = "Salt Lake County",
    state_fips = "49",
    year       = 2025
) {

  census <- multigroup.vaccine::getCensusData(
    state_fips  = state_fips,
    county_name = county,
    year        = year,
    age_groups  = age_limits,
    csv_path    = system.file(
      "extdata",
      sprintf("cc-est%d-syasex-%s.csv", year, state_fips),
      package = "multigroup.vaccine"
    )
  )

  prop <- census$age_pops / sum(census$age_pops)
  ans  <- floor(prop * n)
  left <- n - sum(ans)
  if (left > 0) {
    extra      <- order(prop * n - ans, decreasing = TRUE)[seq_len(left)]
    ans[extra] <- ans[extra] + 1
  }

  stats::setNames(as.integer(ans), rownames(read_epistorm_matrix("Household")))

}

#' U.S. age distribution (proportions) over the same age groups
#' @param year Census estimate year.
us_age_distribution <- function(year = 2025) {

  us <- utils::read.csv(system.file(
    "extdata",
    sprintf("nc-est%d-agesex-res.csv", year),
    package = "multigroup.vaccine"
  ))

  # SEX == 0 is both sexes; AGE 999 is the total row
  us  <- us[us$SEX == 0 & us$AGE < 999, ]
  pop <- tapply(
    us[[sprintf("POPESTIMATE%d", year)]],
    findInterval(us$AGE, age_limits),
    sum
  )

  as.numeric(pop / sum(pop))

}

#' City-level contact matrices by setting
#'
#' Each Epistorm-Mix matrix is re-weighted so contacts with group `j` scale
#' with the city's share of `j` relative to the U.S. share, and then made
#' reciprocal with `measles::make_cmat_symmetric()`.
#' @param pop Group sizes (from `city_population()`).
#' @return A list with one matrix per setting plus `Total`.
contact_matrices <- function(pop) {

  w <- (pop / sum(pop)) / us_age_distribution()

  ans <- lapply(settings, function(s) {
    m <- t(t(read_epistorm_matrix(s)) * w)
    m <- measles::make_cmat_symmetric(m, as.numeric(pop))
    dimnames(m) <- list(names(pop), names(pop))
    m
  })
  names(ans) <- settings

  ans$Total <- Reduce(`+`, ans)
  ans

}

#' Contact matrix under the school-age intervention
#'
#' Non-household contacts (school, work, other) that involve at least one
#' school-age group are reduced by `reduction`. Household contacts are kept,
#' since children kept out of school still mix at home.
#' @param cm List from `contact_matrices()`.
#' @param reduction Proportional reduction in [0, 1].
#' @param school_groups Names of the school-age groups.
intervention_matrix <- function(
    cm,
    reduction     = 0.8,
    school_groups = c("5-9", "10-14", "15-19")
) {

  involved <- outer(
    rownames(cm$Total) %in% school_groups,
    rownames(cm$Total) %in% school_groups,
    `|`
  )

  outside <- cm$School + cm$Work + cm$Other
  cm$Household + outside * (1 - reduction * involved)

}

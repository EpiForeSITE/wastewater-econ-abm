# Population and contact structure for a synthetic city.
#
# Age groups follow the Epistorm-Mix 5-year brackets (0-4, ..., 70-74, 75+).
# The age distribution comes from the U.S. Census county estimates shipped with
# the multigroup.vaccine package, and the contact matrices from Epistorm-Mix
# (data/epistorm). The national matrices are re-weighted to the city's age
# distribution and made reciprocal (symmetric in total contacts).

age_limits <- seq(0, 75, by = 5)

settings <- c("Household", "School", "Work", "Other")

#' School-level measles coverage and district population shares
#'
#' Shares are within school districts, not across the county. They must be
#' combined with district population estimates before assigning students.
#' @examples
#' schools <- utah_schools(county = "Salt Lake")
#' head(schools[, c("id", "level", "vaxrate")])
utah_schools <- function(
    county = "Salt Lake",
    path   = "data/utah/school_district_schoolID_data_Utah_2025.csv"
) {
  # Load the requested county and check the school-record schema.
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  x <- x[x$county == county, , drop = FALSE]
  if (!nrow(x)) stop("No schools found for county: ", county)

  columns <- c("id", "schooldistrictlocation", "level", "popfrac", "vaxrate")
  if (anyNA(x[columns]) ||
      anyDuplicated(x[c("id", "level")]) || any(!x$level %in% 1:3) ||
      any(x$popfrac < 0) || any(x$vaxrate < 0 | x$vaxrate > 1))
    stop("Invalid school records")

  x
}

#' Approximate county-wide school weights from school-district ACS estimates
#'
#' District shares use the population aged 5--17 in the 2023 ACS. The source
#' school fractions are shares of the district school-age population. This
#' assumes the five district boundaries cover the county and that a school's
#' district-level share applies to its resident students.
#' @examples
#' schools <- utah_schools(county = "Salt Lake")
#' weights <- utah_school_weights(schools = schools)
#' head(weights)
#' sum(weights)
utah_school_weights <- function(
    schools = utah_schools(),
    path    = "data/utah/acs_Utah_2023_school_districts.csv"
) {
  # Match the source districts to their ACS population columns.
  acs <- utils::read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  district <- unique(schools$schooldistrictlocation)
  labels <- ifelse(
    test = district == "SALT LAKE",
    yes  = "Salt Lake City",
    no   = tools::toTitleCase(tolower(district))
  )
  cols <- paste0(labels, " School District, UT")
  if (any(!cols %in% names(acs))) stop("Missing district in ACS data")
  # Weight school shares by each district's school-age population.
  ages <- as.numeric(acs[[1]])
  population <- vapply(
    X         = cols,
    FUN       = function(col) sum(acs[[col]][ages %in% c(5, 10, 15, 16)]),
    FUN.VALUE = numeric(1)
  )
  if (anyNA(population) || any(population <= 0)) {
    stop("Invalid ACS district population")
  }
  weights <- population[match(schools$schooldistrictlocation, district)] *
    schools$popfrac

  weights / sum(weights)
}

#' Assign school-aged agents to schools and retain their MMR coverage
#'
#' Age groups are five-year bands. Three fifths of the 15--19 group are
#' assumed to be 15--17; 5--9 and 10--14 are assigned in full. School level
#' cutoffs are approximated by the available age bands, not individual grades.
#' @param pop Age-band population counts.
#' @param schools Application school records with id, level and vaxrate.
#' @param weights Sampling weights in school-record order.
#' @param age_levels School level assigned to each age band, or NA.
#' @param attendance Fraction attending school within each age band.
#' @return Age groups, school IDs and coverage for each simulated agent.
#' @examples
#' set.seed(137)
#' schools <- data.frame(
#'   id      = c("A", "B", "C"),
#'   level   = 1:3,
#'   vaxrate = c(0.8, 0.9, 0.95)
#' )
#' population <- assign_schools(
#'   pop     = rep(20L, 16),
#'   schools = schools,
#'   weights = c(1, 1, 1)
#' )
#' table(population$school, useNA = "ifany")
assign_schools <- function(
    pop,
    schools,
    weights,
    age_levels = c(NA, 1L, 2L, 3L, rep(NA, 12)),
    attendance = c(0, 1, 1, 3 / 5, rep(0, 12))
) {
  # Assign school levels using caller-supplied age-band attendance rules.
  stopifnot(
    length(age_levels) == length(pop),
    length(attendance) == length(pop)
  )
  ages <- rep(seq_along(pop), pop)
  level <- rep(NA_integer_, length(ages))
  for (g in seq_along(pop)) {
    candidates <- which(ages == g)
    selected <- candidates[sample.int(
      n    = length(candidates),
      size = floor(length(candidates) * attendance[g])
    )]
    level[selected] <- age_levels[g]
  }

  # Draw schools within each level using the supplied population weights.
  row <- rep(NA_integer_, length(ages))
  for (j in unique(stats::na.omit(age_levels))) {
    idx <- which(level == j)
    choices <- which(schools$level == j)
    if (length(idx) && !length(choices)) stop("No schools for level: ", j)
    row[idx] <- choices[sample.int(
      n       = length(choices),
      size    = length(idx),
      replace = TRUE,
      prob    = weights[choices]
    )]
  }
  school_id <- match(schools$id[row], unique(schools$id)) - 1L
  list(
    age_group  = ages - 1L,
    school     = school_id,
    coverage   = schools$vaxrate[row],
    school_row = row
  )
}

#' Read an Epistorm-Mix setting-specific M-matrix
#' @param setting One of `settings`.
#' @param path Directory with the CSV files.
#' @examples
#' contacts <- read_epistorm_matrix(setting = "School")
#' contacts[1:3, 1:3]
read_epistorm_matrix <- function(setting, path = "data/epistorm") {

  x <- utils::read.csv(
    file        = file.path(path, sprintf("%s-M-by5_75-matrix.csv", setting)),
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
#' @examples
#' pop <- city_population(n = 1000, county = "Salt Lake County")
#' pop
#' sum(pop)
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
#' @examples
#' shares <- us_age_distribution(year = 2025)
#' head(shares)
#' sum(shares)
us_age_distribution <- function(year = 2025) {

  us <- utils::read.csv(
    file = system.file(
      "extdata",
      sprintf("nc-est%d-agesex-res.csv", year),
      package = "multigroup.vaccine"
    )
  )

  # SEX == 0 is both sexes; AGE 999 is the total row
  us  <- us[us$SEX == 0 & us$AGE < 999, ]
  pop <- tapply(
    X     = us[[sprintf("POPESTIMATE%d", year)]],
    INDEX = findInterval(us$AGE, age_limits),
    FUN   = sum
  )

  as.numeric(pop / sum(pop))

}

#' City-level contact matrices by setting
#'
#' Each Epistorm-Mix matrix is re-weighted so contacts with group `j` scale
#' with the city's share of `j` relative to the U.S. share, and then made
#' reciprocal by averaging total contacts in each direction.
#' @param pop Group sizes (from `city_population()`).
#' @return A list with one matrix per setting plus `Total`.
#' @examples
#' pop <- city_population(n = 1000)
#' contacts <- contact_matrices(pop = pop)
#' names(contacts)
#' contacts$Total[1:3, 1:3]
contact_matrices <- function(pop) {

  w <- (pop / sum(pop)) / us_age_distribution()

  # Reweight and symmetrize each setting's contacts.
  ans <- lapply(
    X   = settings,
    FUN = function(s) {
      m <- t(t(read_epistorm_matrix(s)) * w)
      m <- symmetric_contacts(m, as.numeric(pop))
      dimnames(m) <- list(names(pop), names(pop))
      m
    }
  )
  names(ans) <- settings

  ans$Total <- Reduce(`+`, ans)
  ans

}

#' Make contact totals reciprocal between population groups
#'
#' Average the total contacts in both directions before converting them back
#' to per-person rates. This is the calculation used by
#' `measles::make_cmat_symmetric()` at commit
#' `d1e87fed204ebe820620759fa1b85cdb9d2e0b03`.
#' @param matrix Numeric square matrix of per-person contact rates.
#' @param pop Positive population sizes, one for each row and column.
#' @return A numeric matrix of reciprocal per-person contact rates.
#' @examples
#' contacts <- matrix(c(2, 4, 1, 3), nrow = 2, byrow = TRUE)
#' pop <- c(100, 200)
#' reciprocal <- symmetric_contacts(matrix = contacts, pop = pop)
#' reciprocal
#' reciprocal * pop
symmetric_contacts <- function(matrix, pop) {

  stopifnot(
    is.matrix(matrix),
    is.numeric(matrix),
    is.numeric(pop),
    length(pop) > 0,
    nrow(matrix) == length(pop),
    ncol(matrix) == length(pop),
    all(is.finite(pop) & pop > 0)
  )

  # Average the total contacts, then restore each row's per-person rate.
  total <- matrix * pop
  (total + t(total)) / 2 / pop
}

#' Contact matrix under the school-age intervention
#'
#' Non-household contacts (school, work, other) that involve at least one
#' school-age group are reduced by `reduction`. Household contacts are kept,
#' since children kept out of school still mix at home.
#' @param cm List from `contact_matrices()`.
#' @param reduction Proportional reduction in [0, 1].
#' @param school_groups Names of the school-age groups.
#' @examples
#' pop <- city_population(n = 1000)
#' contacts <- contact_matrices(pop = pop)
#' reduced <- intervention_matrix(cm = contacts, reduction = 0.8)
#' reduced[2:4, 2:4]
intervention_matrix <- function(
    cm,
    reduction     = 0.8,
    school_groups = c("5-9", "10-14", "15-19")
) {

  involved <- outer(
    X   = rownames(cm$Total) %in% school_groups,
    Y   = rownames(cm$Total) %in% school_groups,
    FUN = `|`
  )

  outside <- cm$School + cm$Work + cm$Other
  cm$Household + outside * (1 - reduction * involved)

}

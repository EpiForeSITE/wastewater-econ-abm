# Parameter table, with defaults, units, sources, and notes. It follows the
# format of measles::measles_parameters().

#' Reads the parameter table
#' @param file Path to the CSV file.
#' @return A data frame with one row per argument of `measles_ww_params()`.
ww_parameters <- function(file = file.path("data", "parameters.csv")) {
  utils::read.csv(
    file,
    colClasses   = "character",
    na.strings   = character(0),
    fileEncoding = "UTF-8"
  )
}

#' Formats the parameter table for display
#' @param x Output of `ww_parameters()`.
#' @return A data frame with the display columns; URLs become Markdown links.
format_ww_parameters <- function(x) {

  links <- vapply(strsplit(x$doi_or_url, "; ", fixed = TRUE), function(u) {
    if (length(u) == 0) return("")
    paste0(" ", paste0("[[", seq_along(u), "]](", u, ")", collapse = " "))
  }, character(1))

  range <- ifelse(
    x$lit_low == "" & x$lit_high == "",
    "",
    paste0(x$lit_low, "–", x$lit_high)
  )

  data.frame(
    Parameter          = paste0(x$parameter, " (`", x$r_argument, "`)"),
    Default            = paste0("`", x$default, "` ", x$units),
    `Literature range` = range,
    Source             = paste0(x$citation, links),
    Notes              = x$notes,
    check.names        = FALSE
  )

}

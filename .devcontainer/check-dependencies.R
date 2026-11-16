# Verify the direct R dependencies used to run, test, and render the analysis.
# CRAN versions are from the 2026-10-20 Posit Package Manager snapshot.
stopifnot(getRversion() == package_version("4.6.1"))

versions <- c(
  epiworldR          = "1.0.0-0",
  multigroup.vaccine = "0.1.1",
  cpp11              = "0.5.5",
  data.table         = "1.18.6.1",
  ggplot2            = "4.0.3",
  tinytest           = "1.4.3",
  remotes            = "2.5.0",
  knitr              = "1.52",
  rmarkdown          = "2.32",
  quarto             = "1.5.1",
  socialmixr         = "0.7.0",
  Rcpp               = "1.1.2",
  deSolve            = "1.42",
  shiny              = "1.14.0",
  bslib              = "0.12.0",
  htmltools          = "0.5.9"
)

for (package in names(versions)) {
  actual <- as.character(packageVersion(package))
  if (actual != as.character(package_version(versions[[package]]))) {
    stop(package, ": expected ", versions[[package]], ", found ", actual)
  }
}

# Package versions alone cannot distinguish development commits.
commits <- c(
  epiworldR = "2f1a5acc7835c44908fab53261d879afc4335e51",
  multigroup.vaccine = "ddd9201f28b9908c033542f6a2445b84619af1c7"
)
for (package in names(commits)) {
  actual <- packageDescription(package)$RemoteSha
  if (is.null(actual) || actual != commits[[package]]) {
    stop(package, ": installed Git commit does not match ", commits[[package]])
  }
}

# Check the layered model API and the county Census data used by the example.
api <- c(
  "Model", "layer_mixing", "add_layer", "add_stage", "add_after",
  "duration_lognormal", "assign_agents", "get_agents_stages"
)
stopifnot(all(api %in% getNamespaceExports("epiworldR")))
stopifnot(nzchar(system.file(
  "extdata", "cc-est2025-syasex-49.csv", package = "multigroup.vaccine"
)))

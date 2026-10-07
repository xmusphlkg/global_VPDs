# Shared definitions for the active revised analysis.
CORE_DISEASES <- c("Measles", "Pertussis", "Diphtheria", "Rubella", "Mumps")

DISEASE_PARAMETERS <- tibble::tribble(
  ~disease, ~vaccine, ~window_years,
  "Measles", "MCV1", 10L,
  "Pertussis", "DTP3", 5L,
  "Diphtheria", "DTP3", 10L,
  "Rubella", "RCV1", 10L,
  "Mumps", "MCV1", 10L
)

assert_unique <- function(d, keys) {
  if (anyDuplicated(d[keys])) stop("Duplicate keys: ", paste(keys, collapse = ", "))
}
calendar_summary <- function(year, value, window, fun = mean) {
  stopifnot(!anyDuplicated(year), length(year) == length(value), length(window) == 1)
  vapply(year, function(t) {
    idx <- match(seq.int(t - window, t - 1), year)
    if (anyNA(idx) || any(!is.finite(value[idx]))) return(NA_real_)
    fun(value[idx])
  }, numeric(1))
}
metric_values <- function(actual, prediction) {
  stopifnot(length(actual) == length(prediction), all(is.finite(actual)), all(is.finite(prediction)))
  sst <- sum((actual - mean(actual))^2)
  tibble::tibble(MAE = mean(abs(actual - prediction)),
    RMSE = sqrt(mean((actual - prediction)^2)),
    R2 = if (sst > 0) 1 - sum((actual - prediction)^2) / sst else NA_real_,
    correlation_squared = if (sd(actual) > 0 && sd(prediction) > 0) cor(actual, prediction)^2 else NA_real_)
}
# Fixed physical units avoid train/test standardisation leakage.
prepare_core <- function(d) {
  dplyr::mutate(d,
    coverage_10 = coverage / 10, vulnerability_10 = historical_vulnerability_index / 10,
    prior_burden = log1p(cases_lag1), log_gdp = log(gdp_per_capita_constant_ppp),
    uhc_10 = uhc_service_coverage_index / 10, urban_10 = urban_population_pct / 10,
    health_exp = health_expenditure_pct_gdp, time = (year - 2010) / 10,
    country_f = factor(country), disease_f = factor(disease, levels = CORE_DISEASES),
    country_disease = interaction(country, disease, drop = TRUE))
}
required_core <- c("cases", "population", "coverage_10", "vulnerability_10", "prior_burden",
  "log_gdp", "uhc_10", "urban_10", "health_exp")
complete_core <- function(d) dplyr::filter(d, year >= 2002, population > 0, cases >= 0,
  dplyr::if_all(dplyr::all_of(required_core), is.finite))

#!/usr/bin/env Rscript
# Sensitivity analysis: pooled EU/EEA monthly hurdle models. `mcv2_coverage_lag`
# is the preceding-year annual WUENIC estimate (not a cohort-weighted average).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(MASS)
  library(sandwich)
  library(glmmTMB)
})

out_dir <- "output/revised/measles_resurgence_full"
panel <- read_csv(file.path(out_dir, "measles_monthly_model_panel.csv"), show_col_types = FALSE)
panel$country <- factor(panel$country)
panel$month <- factor(panel$month)
panel$year <- factor(panel$year)

extract_fit <- function(fit, vcov_mat, region, outcome, term, n, n_country, family_name, aic) {
  b <- coef(fit)[term]
  se <- sqrt(diag(vcov_mat))[term]
  scale <- if (grepl("gap_|coverage", term)) 0.1 else 1
  data.frame(
    region = region, outcome = outcome, term = term, family = family_name,
    estimate = unname(b), se = unname(se),
    effect_per_0.10 = exp(unname(b) * scale),
    lower = exp((unname(b) - 1.96 * unname(se)) * scale),
    upper = exp((unname(b) + 1.96 * unname(se)) * scale),
    n = n, countries = n_country, AIC = aic,
    row.names = NULL
  )
}

fit_one <- function(dat, region, outcome, term, positive_only = FALSE) {
  d <- dat %>% filter(is.finite(.data[[term]]),
                      is.finite(mcv2_coverage_lag), is.finite(external_pressure_lag),
                      is.finite(offset_population))
  if (positive_only) d <- d %>% filter(cases > 0)
  if (nrow(d) < 40 || n_distinct(d$country) < 5) return(NULL)
  rhs <- paste0(term, " + mcv2_coverage_lag + external_pressure_lag + month + year + (1|country)")
  if (outcome == "cases") {
    form <- as.formula(paste0("cases ~ ", rhs, " + offset(log(offset_population))"))
    fit <- tryCatch(glmmTMB::glmmTMB(form, data = d, family = glmmTMB::nbinom2()), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    fam <- "negative-binomial"
  } else {
    form <- as.formula(paste0(outcome, " ~ ", rhs))
    fit <- tryCatch(glmmTMB::glmmTMB(form, data = d, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    fam <- "binomial"
  }
  sm <- summary(fit)$coefficients$cond
  b <- sm[term, "Estimate"]
  se <- sm[term, "Std. Error"]
  scale <- if (grepl("gap_|coverage", term)) 0.1 else 1
  data.frame(region = region, outcome = outcome, term = term, family = fam,
             estimate = b, se = se, effect_per_0.10 = exp(b * scale),
             lower = exp((b - 1.96 * se) * scale), upper = exp((b + 1.96 * se) * scale),
             n = nrow(d), countries = n_distinct(d$country), AIC = AIC(fit), row.names = NULL)
}

primary <- bind_rows(lapply("EU/EEA", function(region) {
  bind_rows(
    fit_one(panel, region, "outbreak_month", "gap_1_14_lag"),
    fit_one(panel, region, "sustained_2month", "gap_1_14_lag"),
    fit_one(panel, region, "cases", "gap_1_14_lag", positive_only = TRUE)
  )
}))
write_csv(primary, file.path(out_dir, "measles_hurdle_model_results.csv"))

message("Wrote pooled hurdle model results")

#!/usr/bin/env Rscript
# Within-country age alignment of reported measles cases with reconstructed
# cohort gaps, EU/EEA (ECDC TESSy case-based age distributions, 2023–2025).
#
# cases[k, a] ~ NB(mu), log mu = log pop[k, a] + alpha_k + gamma_a + beta * gap[k, a]
#
# k indexes country-years, so alpha_k absorbs outbreak size and reporting
# completeness; gamma_a are shared age-band effects. beta is therefore
# identified only from within-country-year deviations of the gap across ages.
# The primary fit uses closed 1–29-year bands and same-birth-cohort exposure
# alignment; the open 30+ definition is retained as a sensitivity.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(MASS)
  library(jsonlite)
})

out_dir <- "output/revised/measles_resurgence_full"
top10 <- read_csv(file.path(out_dir, "measles_hotspot_country_summary.csv"), show_col_types = FALSE) %>%
  filter(hotspot) %>% pull(country)
panel_raw <- read_csv(file.path(out_dir, "measles_case_age_panel.csv"), show_col_types = FALSE) %>%
  filter(is.finite(population), population > 0) %>%
  mutate(cases = round(cases), cy = paste(country, outcome_year), top10 = country %in% top10)

# Primary analysis: closed, exactly matched age bands 1–29. The open-ended
# 30+ case numerator cannot be paired with a 30–39 denominator and exposure.
panel <- panel_raw %>% filter(!open_ended, is.finite(gap)) %>%
  mutate(gap10 = gap * 10) %>%
  group_by(cy) %>% filter(n() >= 2, sum(cases) >= 10) %>% ungroup() %>%
  mutate(cy = factor(cy), band = factor(band))

# Open-ended 30+ rows are retained only as a sensitivity result with the
# documented 30–39 denominator/exposure truncation.
panel_open <- panel_raw %>% filter(is.finite(gap)) %>%
  mutate(gap10 = gap * 10) %>%
  group_by(cy) %>% filter(n() >= 2, sum(cases) >= 10) %>% ungroup() %>%
  mutate(cy = factor(cy), band = factor(band))

irr <- function(fit, term = "gap10") {
  s <- summary(fit)$coefficients[term, ]
  c(irr = exp(s[[1]]), lower = exp(s[[1]] - 1.96 * s[[2]]), upper = exp(s[[1]] + 1.96 * s[[2]]))
}
fit <- function(d, rhs = "gap10 + band") glm.nb(as.formula(paste("cases ~", rhs, "+ cy + offset(log(population))")), data = droplevels(d))
describe <- function(d) data.frame(n_rows = nrow(d), n_country_years = n_distinct(d$cy), n_countries = n_distinct(d$country),
  n_cases = sum(d$cases),
  # Identifying variation: within-country-year spread of the gap across bands.
  median_within_gap_range_pp = median(tapply(d$gap, as.character(d$cy), function(g) diff(range(g)))) * 100)
row <- function(model, d) cbind(data.frame(model = model), t(irr(fit(d))), describe(d))

rows <- list(row("eu_eea_band_fe", panel), row("top10_band_fe", filter(panel, top10)),
             row("other_band_fe", filter(panel, !top10)),
             row("open_band_sensitivity", panel_open))
for (y in sort(unique(panel$outcome_year))) rows[[length(rows) + 1]] <- row(paste0("year_", y), filter(panel, outcome_year == y))
loo <- sapply(unique(panel$country), function(cc) irr(fit(filter(panel, country != cc)))[["irr"]])
results <- bind_rows(rows) %>% mutate(
  loo_min = if_else(model == "eu_eea_band_fe", min(loo), NA_real_),
  loo_max = if_else(model == "eu_eea_band_fe", max(loo), NA_real_),
  loo_min_country = if_else(model == "eu_eea_band_fe", names(loo)[which.min(loo)], NA_character_),
  loo_max_country = if_else(model == "eu_eea_band_fe", names(loo)[which.max(loo)], NA_character_))
write_csv(results, file.path(out_dir, "measles_age_alignment_fe.csv"))

# Added-variable coordinates: observed/expected from the model without the
# gap, against the gap after removing country-year and age-band effects.
null_fit <- fit(panel, "band")
panel$expected <- fitted(null_fit)
panel$log2_oe <- log2((panel$cases + .5) / (panel$expected + .5))
panel$gap_resid <- resid(lm(gap ~ cy + band, data = panel))
write_csv(panel %>% dplyr::select(country, outcome_year, top10, band, age_lo, age_hi_used, open_ended, cases,
                                  expected, population, gap, exposure_alignment,
                                  gap_resid, log2_oe),
          file.path(out_dir, "measles_age_alignment_fe_points.csv"))
print(as.data.frame(results), digits = 3)

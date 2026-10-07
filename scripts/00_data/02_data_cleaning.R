# Canonical four-disease panel. Rebuild from source files on every run.
suppressPackageStartupMessages(library(tidyverse))
source("scripts/00_data/00_core_analysis.R")
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("output/revised/tables", recursive = TRUE, showWarnings = FALSE)
coverage <- read_csv("data/raw/coverage_wuenic_all.csv", show_col_types = FALSE) %>%
  filter(group == "COUNTRIES", coverage_category == "WUENIC", year %in% 2000:2025) %>%
  transmute(country, year, vaccine = antigen, coverage, coverage_source = coverage_category)
assert_unique(coverage, c("country", "year", "vaccine"))
stopifnot(all(is.na(coverage$coverage) | between(coverage$coverage, 0, 100)))
cases <- read_csv("data/raw/case_data.csv", show_col_types = FALSE) %>%
  mutate(disease = str_to_title(disease)) %>%
  filter(disease %in% CORE_DISEASES, year %in% 2000:2025)
assert_unique(cases, c("country", "disease", "year"))
stopifnot(all(is.na(cases$cases) | cases$cases >= 0))
pop <- read_csv("data/raw/population_data.csv", show_col_types = FALSE)
meta <- read_csv("data/raw/wpp_country_metadata.csv", show_col_types = FALSE)
covars <- read_csv("data/raw/country_covariates.csv", show_col_types = FALSE)
assert_unique(pop, c("country", "year"))
assert_unique(covars, c("country", "year"))
universe <- Reduce(intersect, list(coverage$country, cases$country, pop$country, meta$country))
# Use a full calendar grid: absent reports remain NA; reported zeros remain zero.
panel <- expand_grid(country = sort(universe), disease = CORE_DISEASES, year = 2000:2025) %>%
  left_join(DISEASE_PARAMETERS, by = "disease", relationship = "many-to-one") %>%
  left_join(coverage, by = c("country", "year", "vaccine"), relationship = "many-to-one") %>%
  left_join(cases, by = c("country", "year", "disease"), relationship = "one-to-one") %>%
  left_join(pop, by = c("country", "year"), relationship = "many-to-one") %>%
  left_join(covars, by = c("country", "year"), relationship = "many-to-one") %>%
  mutate(incidence_per_100k = if_else(population > 0, cases / population * 1e5, NA_real_))
assert_unique(panel, c("country", "disease", "year"))
write_csv(panel, "data/processed/analysis_data.csv")
write_csv(DISEASE_PARAMETERS, "data/processed/disease_parameters.csv")
write_csv(tibble(country = sort(universe)), "output/revised/tables/country_universe.csv")
cat("WUENIC panel:", nrow(panel), "rows;", length(universe), "countries/areas.\n")

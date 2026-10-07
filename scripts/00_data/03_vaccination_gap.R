# Exposure history uses coverage alone, independently of case reporting.
suppressPackageStartupMessages(library(tidyverse))
source("scripts/00_data/00_core_analysis.R")
d <- read_csv("data/processed/analysis_data.csv", show_col_types = FALSE)
assert_unique(d, c("country", "disease", "year"))
d <- d %>% group_by(country, disease) %>% group_modify(~ {
  x <- arrange(.x, year)
  x$historical_vulnerability_index <- calendar_summary(x$year, 100 - x$coverage, x$window_years[1])
  for (w in c(3, 5, 7, 10, 15)) x[[paste0("gap_window_", w)]] <- calendar_summary(x$year, 100-x$coverage, w)
  idx <- match(x$year - 1, x$year)
  x$cases_lag1 <- x$cases[idx]
  x$coverage_change_1yr <- x$coverage - x$coverage[idx]
  x$vaccination_gap_mean <- x$historical_vulnerability_index
  x
}) %>% ungroup()
write_csv(d, "data/processed/data_with_vaccination_gaps.csv")
m <- complete_core(prepare_core(d))
write_csv(m, "data/processed/core_model_data.csv")
write_csv(m, "data/processed/model_ready_data.csv")
write_csv(m, "data/processed/complete_case_data.csv")
write_csv(tibble(variable = c("coverage", "historical_vulnerability_index", "vaccination_gap_mean", "cases_lag1", "population"),
  definition = c("WUENIC dose coverage among the target group (%)", "Mean complement of coverage in all preceding window calendar years; not population susceptibility",
    "Compatibility alias for historical_vulnerability_index", "Reported cases at calendar year t-1; missing stays missing", "WPP midyear population in persons; 2024 onward projections")), "data/processed/variable_description.csv")
cat("Main complete cases:", nrow(m), "from", n_distinct(m$country), "countries/areas; years", min(m$year), "to", max(m$year), "\n")

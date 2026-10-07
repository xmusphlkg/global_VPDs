# ==============================================================================
# Prepare Measles Dose-Specific Coverage Data for Cohort Reconstruction
# ==============================================================================
# Purpose: Extract MCV1 and MCV2 coverage from WUENIC data and prepare for
#          cohort-based immunity gap reconstruction
#
# Input: data/raw/coverage_wuenic_all.csv
# Output: data/processed/cohort_immunity_profiles/measles_dose_coverage.csv
# ==============================================================================

library(tidyverse)

# Create output directory
dir.create("data/processed/cohort_immunity_profiles", recursive = TRUE, showWarnings = FALSE)

message("=============================================================")
message("PHASE 2A: MEASLES DOSE-SPECIFIC COVERAGE DATA PREPARATION")
message("=============================================================\n")

# ------------------------------------------------------------------------------
# Load WUENIC Data
# ------------------------------------------------------------------------------

message("Loading WUENIC coverage data...")
wuenic_raw <- read_csv("data/raw/coverage_wuenic_all.csv",
                       show_col_types = FALSE)

message("Raw WUENIC data: ", nrow(wuenic_raw), " rows")
message("Columns: ", paste(names(wuenic_raw), collapse = ", "))

# ------------------------------------------------------------------------------
# Extract MCV1 and MCV2 Data
# ------------------------------------------------------------------------------

message("\nExtracting MCV1 and MCV2 coverage...")

# Filter for measles-containing vaccines
mcv_data <- wuenic_raw %>%
  filter(antigen %in% c("MCV1", "MCV2")) %>%
  filter(coverage_category == "WUENIC") %>%  # WHO/UNICEF estimates only
  select(
    country = country,
    country_name = country_name,
    year = year,
    dose = antigen,
    coverage = coverage
  ) %>%
  mutate(
    coverage = as.numeric(coverage),
    year = as.integer(year)
  ) %>%
  filter(!is.na(coverage), !is.na(year))

message("Extracted MCV data: ", nrow(mcv_data), " rows")
message("Countries: ", n_distinct(mcv_data$country))
message("Year range: ", min(mcv_data$year), " to ", max(mcv_data$year))

# Check coverage by dose
dose_summary <- mcv_data %>%
  group_by(dose) %>%
  summarise(
    n_obs = n(),
    n_countries = n_distinct(country),
    mean_coverage = mean(coverage, na.rm = TRUE),
    min_year = min(year),
    max_year = max(year),
    .groups = "drop"
  )

message("\nDose-specific summary:")
print(dose_summary)

# ------------------------------------------------------------------------------
# Data Quality Checks
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("DATA QUALITY CHECKS")
message(rep("-", 60))

# Check 1: Are there country-years with MCV2 but no MCV1?
mcv_wide <- mcv_data %>%
  pivot_wider(
    names_from = dose,
    values_from = coverage
  )

invalid_mcv2 <- mcv_wide %>%
  filter(is.na(MCV1) & !is.na(MCV2))

if (nrow(invalid_mcv2) > 0) {
  message("\nWARNING: ", nrow(invalid_mcv2), " country-years have MCV2 without MCV1:")
  print(invalid_mcv2 %>% select(country_name, year, MCV2) %>% head(10))
} else {
  message("\n✓ No country-years with MCV2 but no MCV1")
}

# Check 2: Are there MCV2 > MCV1 cases?
invalid_sequence <- mcv_wide %>%
  filter(!is.na(MCV1) & !is.na(MCV2) & MCV2 > MCV1)

if (nrow(invalid_sequence) > 0) {
  message("\nWARNING: ", nrow(invalid_sequence), " country-years have MCV2 > MCV1:")
  print(invalid_sequence %>%
          select(country_name, year, MCV1, MCV2) %>%
          arrange(desc(MCV2 - MCV1)) %>%
          head(10))
  message("\nNote: This can happen due to:")
  message("  1. Different denominators (birth cohort size estimates)")
  message("  2. Catch-up vaccination reaching older children")
  message("  3. Data quality issues")
} else {
  message("\n✓ All MCV2 values ≤ MCV1")
}

# Check 3: Coverage gaps by year
coverage_gaps <- mcv_wide %>%
  group_by(year) %>%
  summarise(
    n_countries_total = n(),
    n_with_mcv1 = sum(!is.na(MCV1)),
    n_with_mcv2 = sum(!is.na(MCV2)),
    pct_with_mcv1 = 100 * n_with_mcv1 / n_countries_total,
    pct_with_mcv2 = 100 * n_with_mcv2 / n_countries_total,
    .groups = "drop"
  )

message("\nCoverage data availability by year:")
message("(Showing every 5 years)")
print(coverage_gaps %>% filter(year %% 5 == 0))

# Plot coverage availability over time
p_availability <- coverage_gaps %>%
  select(year, MCV1 = pct_with_mcv1, MCV2 = pct_with_mcv2) %>%
  pivot_longer(cols = c(MCV1, MCV2), names_to = "dose", values_to = "pct_countries") %>%
  ggplot(aes(x = year, y = pct_countries, color = dose)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "MCV1/MCV2 Data Availability Over Time",
    x = "Year",
    y = "% of Countries with Data",
    color = "Dose"
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 100))

ggsave("output/diagnostics/mcv_data_availability.png",
       p_availability, width = 10, height = 6)

# ------------------------------------------------------------------------------
# Historical Coverage Length by Country
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("HISTORICAL DATA LENGTH BY COUNTRY")
message(rep("-", 60))

country_history <- mcv_data %>%
  filter(dose == "MCV1") %>%  # Use MCV1 as reference
  group_by(country, country_name) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    n_years = n(),
    mean_coverage = mean(coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(first_year, country_name)

message("\nCountries by data length:")
length_dist <- country_history %>%
  mutate(years_group = cut(n_years,
                           breaks = c(0, 5, 10, 15, 20, 30, 50),
                           labels = c("1-5", "6-10", "11-15", "16-20", "21-30", "31+"))) %>%
  count(years_group)

print(length_dist)

# Countries with longest history (for cohort reconstruction)
message("\nCountries with longest MCV1 history (top 20):")
print(country_history %>%
        arrange(desc(n_years)) %>%
        select(country_name, first_year, last_year, n_years, mean_coverage) %>%
        head(20))

# Countries with shortest history
message("\nCountries with shortest MCV1 history (bottom 10):")
print(country_history %>%
        arrange(n_years) %>%
        select(country_name, first_year, last_year, n_years, mean_coverage) %>%
        head(10))

# ------------------------------------------------------------------------------
# Save Prepared Data
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("SAVING PREPARED DATA")
message(rep("-", 60))

# Save long format (easier for cohort reconstruction)
write_csv(mcv_data, "data/processed/cohort_immunity_profiles/measles_dose_coverage_long.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/measles_dose_coverage_long.csv")

# Save wide format (easier for analysis)
write_csv(mcv_wide, "data/processed/cohort_immunity_profiles/measles_dose_coverage_wide.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/measles_dose_coverage_wide.csv")

# Save country history metadata
write_csv(country_history, "data/processed/cohort_immunity_profiles/measles_country_history.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/measles_country_history.csv")

# ------------------------------------------------------------------------------
# Summary for Next Steps
# ------------------------------------------------------------------------------

message("\n", rep("=", 60))
message("SUMMARY AND NEXT STEPS")
message(rep("=", 60))

message("\nData prepared:")
message("  • MCV1 coverage: ", sum(mcv_data$dose == "MCV1"), " country-years")
message("  • MCV2 coverage: ", sum(mcv_data$dose == "MCV2"), " country-years")
message("  • Countries: ", n_distinct(mcv_data$country))
message("  • Year range: ", min(mcv_data$year), "-", max(mcv_data$year))

message("\nFor cohort reconstruction, we need:")
message("  1. ✓ MCV1 and MCV2 coverage by country-year (COMPLETE)")
message("  2. ⏳ Immunization schedules (target ages for each dose)")
message("  3. ⏳ Birth cohort populations")
message("  4. ⏳ SIA data (supplementary immunization activities)")

message("\nNext script: 02b_prepare_immunization_schedules.R")
message("\n", rep("=", 60))

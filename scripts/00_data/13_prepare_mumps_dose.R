# ==============================================================================
# Prepare Mumps Dose-Specific Coverage Data for Cohort Reconstruction
# ==============================================================================
# Purpose: Extract MMR coverage for mumps analysis
#          Mumps uses the SAME MMR vaccine as measles, so we use MCV1/MCV2
#          as direct proxies
#
# Input: data/raw/coverage_wuenic_all.csv
# Output: data/processed/cohort_immunity_profiles/mumps_dose_coverage_*.csv
# ==============================================================================

library(tidyverse)

# Create output directory
dir.create("data/processed/cohort_immunity_profiles", recursive = TRUE, showWarnings = FALSE)

message("=============================================================")
message("PHASE 2A: MUMPS DOSE-SPECIFIC COVERAGE DATA PREPARATION")
message("=============================================================\n")

message("NOTE: Mumps uses MMR vaccine (same as measles)")
message("      Using MCV1/MCV2 data as direct proxy for mumps coverage\n")

# ------------------------------------------------------------------------------
# Load WUENIC Data
# ------------------------------------------------------------------------------

message("Loading WUENIC coverage data...")
wuenic_raw <- read_csv("data/raw/coverage_wuenic_all.csv",
                       show_col_types = FALSE)

message("Raw WUENIC data: ", nrow(wuenic_raw), " rows")

# ------------------------------------------------------------------------------
# Extract MCV1 and MCV2 Data (for Mumps)
# ------------------------------------------------------------------------------

message("\nExtracting MCV1 and MCV2 coverage (MMR vaccine)...")

# Filter for measles-containing vaccines (which include mumps)
mmr_data <- wuenic_raw %>%
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

message("Extracted MMR data: ", nrow(mmr_data), " rows")
message("Countries: ", n_distinct(mmr_data$country))
message("Year range: ", min(mmr_data$year), " to ", max(mmr_data$year))

# Check coverage by dose
dose_summary <- mmr_data %>%
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
mmr_wide <- mmr_data %>%
  pivot_wider(
    names_from = dose,
    values_from = coverage
  )

invalid_mcv2 <- mmr_wide %>%
  filter(is.na(MCV1) & !is.na(MCV2))

if (nrow(invalid_mcv2) > 0) {
  message("\nWARNING: ", nrow(invalid_mcv2), " country-years have MCV2 without MCV1")
} else {
  message("\n✓ No country-years with MCV2 but no MCV1")
}

# Check 2: Are there MCV2 > MCV1 cases?
invalid_sequence <- mmr_wide %>%
  filter(!is.na(MCV1) & !is.na(MCV2) & MCV2 > MCV1)

if (nrow(invalid_sequence) > 0) {
  message("\nWARNING: ", nrow(invalid_sequence), " country-years have MCV2 > MCV1")
  message("\nNote: This can happen due to:")
  message("  1. Different denominators (birth cohort size estimates)")
  message("  2. Catch-up vaccination reaching older children")
  message("  3. Data quality issues")
} else {
  message("\n✓ All MCV2 values ≤ MCV1")
}

# ------------------------------------------------------------------------------
# Historical Coverage Length by Country
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("HISTORICAL DATA LENGTH BY COUNTRY")
message(rep("-", 60))

country_history <- mmr_data %>%
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

# ------------------------------------------------------------------------------
# Save Prepared Data
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("SAVING PREPARED DATA")
message(rep("-", 60))

# Save long format (easier for cohort reconstruction)
write_csv(mmr_data, "data/processed/cohort_immunity_profiles/mumps_dose_coverage_long.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/mumps_dose_coverage_long.csv")

# Save wide format (easier for analysis)
write_csv(mmr_wide, "data/processed/cohort_immunity_profiles/mumps_dose_coverage_wide.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/mumps_dose_coverage_wide.csv")

# Save country history metadata
write_csv(country_history, "data/processed/cohort_immunity_profiles/mumps_country_history.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/mumps_country_history.csv")

# ------------------------------------------------------------------------------
# Summary for Next Steps
# ------------------------------------------------------------------------------

message("\n", rep("=", 60))
message("SUMMARY AND NEXT STEPS")
message(rep("=", 60))

message("\nData prepared:")
message("  • MCV1 coverage (MMR dose 1): ", sum(mmr_data$dose == "MCV1"), " country-years")
message("  • MCV2 coverage (MMR dose 2): ", sum(mmr_data$dose == "MCV2"), " country-years")
message("  • Countries: ", n_distinct(mmr_data$country))
message("  • Year range: ", min(mmr_data$year), "-", max(mmr_data$year))

message("\nIMPORTANT NOTE:")
message("  Mumps uses the SAME MMR vaccine as measles")
message("  → Can reuse measles immunization schedule")
message("  → Main difference: LOWER vaccine efficacy for mumps")
message("     - MMR dose 1: ~88% for mumps (vs 85% for measles)")
message("     - MMR dose 2: ~95% for mumps (vs 95% for measles)")
message("  → Consider waning immunity (outbreaks in high-coverage populations)")

message("\nNext: Copy measles schedule for mumps")
message("\n", rep("=", 60))

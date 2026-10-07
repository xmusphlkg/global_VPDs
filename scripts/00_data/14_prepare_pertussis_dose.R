# ==============================================================================
# Prepare Pertussis (DTP) Dose-Specific Coverage Data for Cohort Reconstruction
# ==============================================================================
# Purpose: Extract DTP1, DTP3, DTP4 coverage from WUENIC data and prepare for
#          cohort-based immunity gap reconstruction with WANING IMMUNITY
#
# Input: data/raw/coverage_wuenic_all.csv
# Output: data/processed/cohort_immunity_profiles/pertussis_dose_coverage_long.csv
#
# KEY DIFFERENCE FROM MEASLES: Pertussis protection wanes over time
#   - DTP3 efficacy ~90% initially, drops to ~50% after 5 years
# ==============================================================================

library(tidyverse)

message("=============================================================")
message("PERTUSSIS DOSE-SPECIFIC COVERAGE DATA PREPARATION")
message("=============================================================\n")

# ------------------------------------------------------------------------------
# Load WUENIC Data
# ------------------------------------------------------------------------------

message("Loading WUENIC coverage data...")
wuenic_raw <- read_csv("data/raw/coverage_wuenic_all.csv",
                       show_col_types = FALSE)

message("Raw WUENIC data: ", nrow(wuenic_raw), " rows")

# ------------------------------------------------------------------------------
# Extract DTP Coverage Data
# ------------------------------------------------------------------------------

message("\nExtracting DTP coverage (pertussis is in DTP vaccine)...")

# Filter for DTP-containing vaccines
# DTP1 and DTP3 = primary series (DTP4 not available in WUENIC data)
dtp_data <- wuenic_raw %>%
  filter(antigen %in% c("DTP1", "DTP3")) %>%
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

message("Extracted DTP data: ", nrow(dtp_data), " rows")
message("Countries: ", n_distinct(dtp_data$country))
message("Year range: ", min(dtp_data$year), " to ", max(dtp_data$year))

# Check coverage by dose
dose_summary <- dtp_data %>%
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

# Check 1: Are there country-years with DTP3 but no DTP1?
dtp_wide <- dtp_data %>%
  pivot_wider(
    names_from = dose,
    values_from = coverage
  )

invalid_dtp3 <- dtp_wide %>%
  filter(is.na(DTP1) & !is.na(DTP3))

if (nrow(invalid_dtp3) > 0) {
  message("\nWARNING: ", nrow(invalid_dtp3), " country-years have DTP3 without DTP1:")
  print(invalid_dtp3 %>% select(country_name, year, DTP3) %>% head(10))
} else {
  message("\n✓ No country-years with DTP3 but no DTP1")
}

# Check 2: Are there DTP3 > DTP1 cases?
invalid_sequence <- dtp_wide %>%
  filter(!is.na(DTP1) & !is.na(DTP3) & DTP3 > DTP1)

if (nrow(invalid_sequence) > 0) {
  message("\nWARNING: ", nrow(invalid_sequence), " country-years have DTP3 > DTP1:")
  print(invalid_sequence %>%
          select(country_name, year, DTP1, DTP3) %>%
          arrange(desc(DTP3 - DTP1)) %>%
          head(10))
  message("\nNote: This can happen due to different denominators")
} else {
  message("\n✓ All DTP3 values ≤ DTP1")
}

# Check 3: Note about boosters
message("\nNote on DTP4 (booster doses):")
message("  DTP4 data not available in WUENIC dataset")
message("  Analysis focuses on DTP3 (primary series completion)")
message("  This provides a conservative baseline estimate")
message("  Booster doses would provide additional protection not captured here")

# ------------------------------------------------------------------------------
# Historical Coverage Length by Country
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("HISTORICAL DATA LENGTH BY COUNTRY")
message(rep("-", 60))

country_history <- dtp_data %>%
  filter(dose == "DTP3") %>%  # Use DTP3 as reference
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

# Countries with longest history
message("\nCountries with longest DTP3 history (top 20):")
print(country_history %>%
        arrange(desc(n_years)) %>%
        select(country_name, first_year, last_year, n_years, mean_coverage) %>%
        head(20))

# ------------------------------------------------------------------------------
# Save Prepared Data
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("SAVING PREPARED DATA")
message(rep("-", 60))

# Save long format (easier for cohort reconstruction)
write_csv(dtp_data, "data/processed/cohort_immunity_profiles/pertussis_dose_coverage_long.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/pertussis_dose_coverage_long.csv")

# Save wide format (easier for analysis)
write_csv(dtp_wide, "data/processed/cohort_immunity_profiles/pertussis_dose_coverage_wide.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/pertussis_dose_coverage_wide.csv")

# Save country history metadata
write_csv(country_history, "data/processed/cohort_immunity_profiles/pertussis_country_history.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/pertussis_country_history.csv")

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

message("\n", rep("=", 60))
message("SUMMARY AND NEXT STEPS")
message(rep("=", 60))

message("\nData prepared:")
message("  • DTP1 coverage: ", sum(dtp_data$dose == "DTP1"), " country-years")
message("  • DTP3 coverage: ", sum(dtp_data$dose == "DTP3"), " country-years")
message("  • Countries: ", n_distinct(dtp_data$country))
message("  • Year range: ", min(dtp_data$year), "-", max(dtp_data$year))

message("\n⚠️  KEY DIFFERENCE FROM MEASLES:")
message("  Pertussis protection WANES over time!")
message("  DTP3 efficacy: ~90% (year 1) → ~50% (year 5+)")
message("  This requires special handling in cohort reconstruction")

message("\nNext steps:")
message("  1. Create pertussis immunization schedule (DTP1-3 timing)")
message("  2. Implement cohort reconstruction WITH waning immunity")
message("  3. Compare with measles to show impact of waning")

message("\n", rep("=", 60))

# ==============================================================================
# Prepare Rubella Dose-Specific Coverage Data for Cohort Reconstruction
# ==============================================================================
# Purpose: Extract RCV1 and RCV2 coverage from WUENIC data and prepare for
#          cohort-based immunity gap reconstruction
#
# Input: data/raw/coverage_wuenic_all.csv
# Output: data/processed/cohort_immunity_profiles/rubella_dose_coverage_*.csv
# ==============================================================================

library(tidyverse)

# Create output directory
dir.create("data/processed/cohort_immunity_profiles", recursive = TRUE, showWarnings = FALSE)

message("=============================================================")
message("RUBELLA DOSE-SPECIFIC COVERAGE DATA PREPARATION")
message("=============================================================\n")

# ------------------------------------------------------------------------------
# Load WUENIC Data
# ------------------------------------------------------------------------------

message("Loading WUENIC coverage data...")
wuenic_raw <- read_csv("data/raw/coverage_wuenic_all.csv",
                       show_col_types = FALSE)

message("Raw WUENIC data: ", nrow(wuenic_raw), " rows")

# ------------------------------------------------------------------------------
# Extract RCV1 and RCV2 Data
# ------------------------------------------------------------------------------

message("\nExtracting RCV1 and RCV2 coverage...")

# Check what rubella antigens are available
rubella_antigens <- wuenic_raw %>%
  filter(str_detect(antigen, "^RCV")) %>%
  distinct(antigen, antigen_description)

message("\nAvailable rubella antigens:")
print(rubella_antigens)

# Filter for rubella-containing vaccines
rcv_data <- wuenic_raw %>%
  filter(antigen %in% c("RCV1", "RCV2")) %>%
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

message("Extracted RCV data: ", nrow(rcv_data), " rows")
message("Countries: ", n_distinct(rcv_data$country))
message("Year range: ", min(rcv_data$year), " to ", max(rcv_data$year))

# Check coverage by dose
dose_summary <- rcv_data %>%
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
# Handle Missing RCV2 Data
# ------------------------------------------------------------------------------

# Check if RCV2 exists
has_rcv2 <- "RCV2" %in% unique(rcv_data$dose)

if (!has_rcv2) {
  message("\nNOTE: RCV2 data not found in WUENIC dataset.")
  message("This is expected as rubella is often delivered via MMR vaccine.")
  message("Creating placeholder RCV2 data structure for consistency.")

  # Create empty RCV2 structure
  rcv2_placeholder <- rcv_data %>%
    filter(dose == "RCV1") %>%
    mutate(
      dose = "RCV2",
      coverage = NA_real_
    ) %>%
    head(0)  # Empty dataframe with correct structure

  rcv_data <- bind_rows(rcv_data, rcv2_placeholder)
}

# ------------------------------------------------------------------------------
# Data Quality Checks
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("DATA QUALITY CHECKS")
message(rep("-", 60))

# Create wide format for checks
rcv_wide <- rcv_data %>%
  pivot_wider(
    names_from = dose,
    values_from = coverage
  )

# Check 1: Are there country-years with RCV2 but no RCV1?
invalid_rcv2 <- rcv_wide %>%
  filter(is.na(RCV1) & !is.na(RCV2))

if (nrow(invalid_rcv2) > 0) {
  message("\nWARNING: ", nrow(invalid_rcv2), " country-years have RCV2 without RCV1:")
  print(invalid_rcv2 %>% select(country_name, year, RCV2) %>% head(10))
} else {
  message("\n✓ No country-years with RCV2 but no RCV1")
}

# Check 2: Are there RCV2 > RCV1 cases?
invalid_sequence <- rcv_wide %>%
  filter(!is.na(RCV1) & !is.na(RCV2) & RCV2 > RCV1)

if (nrow(invalid_sequence) > 0) {
  message("\nWARNING: ", nrow(invalid_sequence), " country-years have RCV2 > RCV1:")
  print(invalid_sequence %>%
          select(country_name, year, RCV1, RCV2) %>%
          arrange(desc(RCV2 - RCV1)) %>%
          head(10))
} else {
  message("\n✓ All RCV2 values ≤ RCV1 (or RCV2 not available)")
}

# Check 3: Coverage gaps by year
coverage_gaps <- rcv_wide %>%
  group_by(year) %>%
  summarise(
    n_countries_total = n(),
    n_with_rcv1 = sum(!is.na(RCV1)),
    n_with_rcv2 = sum(!is.na(RCV2)),
    pct_with_rcv1 = 100 * n_with_rcv1 / n_countries_total,
    pct_with_rcv2 = 100 * n_with_rcv2 / n_countries_total,
    .groups = "drop"
  )

message("\nCoverage data availability by year (every 5 years):")
print(coverage_gaps %>% filter(year %% 5 == 0))

# ------------------------------------------------------------------------------
# Historical Coverage Length by Country
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("HISTORICAL DATA LENGTH BY COUNTRY")
message(rep("-", 60))

country_history <- rcv_data %>%
  filter(dose == "RCV1") %>%  # Use RCV1 as reference
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

# Top countries with longest history
message("\nCountries with longest RCV1 history (top 20):")
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

# Save long format
write_csv(rcv_data, "data/processed/cohort_immunity_profiles/rubella_dose_coverage_long.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/rubella_dose_coverage_long.csv")

# Save wide format
write_csv(rcv_wide, "data/processed/cohort_immunity_profiles/rubella_dose_coverage_wide.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/rubella_dose_coverage_wide.csv")

# Save country history metadata
write_csv(country_history, "data/processed/cohort_immunity_profiles/rubella_country_history.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/rubella_country_history.csv")

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

message("\n", rep("=", 60))
message("SUMMARY")
message(rep("=", 60))

message("\nData prepared:")
message("  • RCV1 coverage: ", sum(rcv_data$dose == "RCV1"), " country-years")
message("  • RCV2 coverage: ", sum(rcv_data$dose == "RCV2" & !is.na(rcv_data$coverage)), " country-years")
message("  • Countries: ", n_distinct(rcv_data$country))
message("  • Year range: ", min(rcv_data$year), "-", max(rcv_data$year))

message("\nNext: 02b_prepare_rubella_schedules.R")
message(rep("=", 60), "\n")

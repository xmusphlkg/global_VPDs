# ==============================================================================
# Prepare Birth Cohort Population Data
# ==============================================================================
# Purpose: Extract annual births by country to serve as denominators for
#          cohort-based vaccination coverage and immunity gap calculations
#
# Input: data/raw/births_data.csv (from UN WPP)
# Output: data/processed/cohort_immunity_profiles/birth_cohorts.csv
# ==============================================================================

library(tidyverse)

message("=============================================================")
message("PHASE 2C: BIRTH COHORT POPULATION DATA PREPARATION")
message("=============================================================\n")

# ------------------------------------------------------------------------------
# Load Birth Data
# ------------------------------------------------------------------------------

message("Loading birth cohort data...")
births_raw <- read_csv("data/raw/births_data.csv",
                       show_col_types = FALSE)

message("Raw births data: ", nrow(births_raw), " rows")
message("Columns: ", paste(names(births_raw), collapse = ", "))

# Check structure
glimpse(births_raw)

# ------------------------------------------------------------------------------
# Check Country Coverage
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("COUNTRY COVERAGE CHECK")
message(rep("-", 60))

# Load country list from coverage data
coverage_countries <- read_csv("data/processed/cohort_immunity_profiles/measles_dose_coverage_wide.csv",
                               show_col_types = FALSE) %>%
  distinct(country, country_name)

message("\nCountries with coverage data: ", nrow(coverage_countries))

# Check which countries have birth data
countries_with_births <- births_raw %>%
  distinct(country)

message("Countries with birth data: ", nrow(countries_with_births))

# Find missing
missing_births <- coverage_countries %>%
  anti_join(countries_with_births, by = "country")

if (nrow(missing_births) > 0) {
  message("\nWARNING: ", nrow(missing_births), " countries have coverage but no birth data:")
  print(missing_births)
} else {
  message("\n✓ All countries with coverage data have birth data")
}

# ------------------------------------------------------------------------------
# Prepare Birth Cohort Dataset
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("PREPARING BIRTH COHORT DATASET")
message(rep("-", 60))

# Filter and clean birth data
birth_cohorts <- births_raw %>%
  # Keep only countries with coverage data
  semi_join(coverage_countries, by = "country") %>%
  # Add country_name from coverage data
  left_join(coverage_countries, by = "country") %>%
  select(
    country,
    country_name,
    year = year,
    births = births
  ) %>%
  filter(!is.na(births), births > 0) %>%
  arrange(country, year)

message("\nBirth cohort data prepared:")
message("  Countries: ", n_distinct(birth_cohorts$country))
message("  Year range: ", min(birth_cohorts$year), " to ", max(birth_cohorts$year))
message("  Total country-years: ", nrow(birth_cohorts))

# ------------------------------------------------------------------------------
# Data Quality Checks
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("DATA QUALITY CHECKS")
message(rep("-", 60))

# Check 1: Are there gaps in time series?
year_gaps <- birth_cohorts %>%
  group_by(country) %>%
  arrange(year) %>%
  mutate(
    year_diff = year - lag(year),
    has_gap = !is.na(year_diff) & year_diff > 1
  ) %>%
  filter(has_gap) %>%
  ungroup()

if (nrow(year_gaps) > 0) {
  message("\nWARNING: ", n_distinct(year_gaps$country), " countries have gaps in time series:")
  year_gaps_summary <- year_gaps %>%
    group_by(country, country_name) %>%
    summarise(
      n_gaps = n(),
      max_gap = max(year_diff),
      .groups = "drop"
    ) %>%
    arrange(desc(max_gap))
  print(head(year_gaps_summary, 10))
} else {
  message("\n✓ No gaps in birth time series")
}

# Check 2: Distribution of births
birth_stats <- birth_cohorts %>%
  group_by(country, country_name) %>%
  summarise(
    n_years = n(),
    first_year = min(year),
    last_year = max(year),
    mean_births = mean(births),
    median_births = median(births),
    min_births = min(births),
    max_births = max(births),
    cv_births = sd(births) / mean(births),
    .groups = "drop"
  )

message("\nBirth statistics by country:")
message("Summary of mean annual births:")
print(summary(birth_stats$mean_births))

# Countries with most births
message("\nCountries with highest mean annual births:")
print(birth_stats %>%
        arrange(desc(mean_births)) %>%
        select(country_name, mean_births, n_years) %>%
        head(10))

# Countries with most variation
message("\nCountries with highest birth variability (CV):")
print(birth_stats %>%
        filter(!is.na(cv_births), cv_births < Inf) %>%
        arrange(desc(cv_births)) %>%
        select(country_name, cv_births, mean_births) %>%
        head(10))

# ------------------------------------------------------------------------------
# Check Alignment with Coverage Years
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("ALIGNMENT WITH COVERAGE DATA")
message(rep("-", 60))

# Load coverage data to check year alignment
coverage_years <- read_csv("data/processed/cohort_immunity_profiles/measles_dose_coverage_wide.csv",
                           show_col_types = FALSE) %>%
  filter(!is.na(MCV1)) %>%
  select(country, country_name, coverage_year = year)

# For each coverage year, we need births from cohorts that could be vaccinated
# MCV1 typically at 9-12 months, so coverage in year Y targets births from year Y-1 mostly
# MCV2 typically at 15-18 months, so coverage in year Y targets births from year Y-2 mostly

alignment_check <- coverage_years %>%
  mutate(
    # Need births from 3 years before (for MCV2 in older children)
    # to year of coverage (for late vaccinations)
    needed_birth_year_start = coverage_year - 3,
    needed_birth_year_end = coverage_year
  ) %>%
  rowwise() %>%
  mutate(
    has_all_births = all(
      needed_birth_year_start:needed_birth_year_end %in%
        birth_cohorts$year[birth_cohorts$country == country]
    )
  ) %>%
  ungroup()

missing_alignment <- alignment_check %>%
  filter(!has_all_births)

if (nrow(missing_alignment) > 0) {
  message("\nWARNING: ", nrow(missing_alignment),
          " country-years have coverage but incomplete birth cohort data:")
  missing_summary <- missing_alignment %>%
    group_by(country_name) %>%
    summarise(
      n_affected_years = n(),
      first_affected = min(coverage_year),
      last_affected = max(coverage_year),
      .groups = "drop"
    )
  print(head(missing_summary, 10))
  message("\nNote: This typically affects early years (2000-2003) where we lack")
  message("births from 1997-1999. Can use later coverage years for analysis.")
} else {
  message("\n✓ All coverage years have corresponding birth cohort data")
}

# ------------------------------------------------------------------------------
# Visualization
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("CREATING VISUALIZATIONS")
message(rep("-", 60))

# Plot 1: Birth trends over time (selected countries)
high_birth_countries <- birth_stats %>%
  arrange(desc(mean_births)) %>%
  head(8) %>%
  pull(country)

p1_birth_trends <- birth_cohorts %>%
  filter(country %in% high_birth_countries) %>%
  ggplot(aes(x = year, y = births / 1e6, color = country_name)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Birth Trends in High-Population Countries",
    x = "Year",
    y = "Annual Births (millions)",
    color = "Country"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("output/diagnostics/birth_cohort_trends.png",
       p1_birth_trends, width = 12, height = 8)

# Plot 2: Distribution of birth cohort sizes
p2_birth_dist <- birth_cohorts %>%
  filter(year >= 2000) %>%  # Focus on analysis period
  ggplot(aes(x = births)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "Distribution of Annual Birth Cohort Sizes (2000-2025)",
    x = "Annual Births (log scale)",
    y = "Number of Country-Years"
  ) +
  theme_minimal()

ggsave("output/diagnostics/birth_cohort_size_distribution.png",
       p2_birth_dist, width = 10, height = 6)

# ------------------------------------------------------------------------------
# Calculate Cohort Sizes by Age
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("CALCULATING AGE-STRATIFIED COHORTS")
message(rep("-", 60))

# For each observation year, calculate population by age (birth cohort)
# This will be used to weight immunity gaps by cohort size

# Get range of years we need
coverage_year_range <- coverage_years %>%
  summarise(min_year = min(coverage_year), max_year = max(coverage_year))

message("\nCreating age-stratified cohorts for observation years ",
        coverage_year_range$min_year, "-", coverage_year_range$max_year)

# For each observation year, we need cohorts aged 0-14 years
# (WHO recommends tracking immunity in children and adolescents)

age_stratified_cohorts <- crossing(
  observation_year = coverage_year_range$min_year:coverage_year_range$max_year,
  age = 0:14
) %>%
  mutate(birth_year = observation_year - age) %>%
  crossing(
    birth_cohorts %>%
      distinct(country, country_name) %>%
      filter(country %in% coverage_countries$country)
  ) %>%
  left_join(
    birth_cohorts %>% select(country, birth_year = year, cohort_size = births),
    by = c("country", "birth_year")
  ) %>%
  filter(!is.na(cohort_size))

message("\nAge-stratified cohort data created:")
message("  Total records: ", nrow(age_stratified_cohorts))
message("  Countries: ", n_distinct(age_stratified_cohorts$country))
message("  Observation years: ",
        min(age_stratified_cohorts$observation_year), "-",
        max(age_stratified_cohorts$observation_year))

# Calculate total 0-14 population by country-year
total_014_pop <- age_stratified_cohorts %>%
  group_by(country, country_name, observation_year) %>%
  summarise(
    total_pop_0_14 = sum(cohort_size, na.rm = TRUE),
    n_cohorts_available = n(),
    .groups = "drop"
  )

message("\nTotal 0-14 population calculated for ",
        nrow(total_014_pop), " country-years")

# ------------------------------------------------------------------------------
# Save Results
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("SAVING RESULTS")
message(rep("-", 60))

# Save annual birth cohorts
write_csv(birth_cohorts,
          "data/processed/cohort_immunity_profiles/birth_cohorts.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/birth_cohorts.csv")

# Save birth statistics
write_csv(birth_stats,
          "data/processed/cohort_immunity_profiles/birth_cohort_statistics.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/birth_cohort_statistics.csv")

# Save age-stratified cohorts
write_csv(age_stratified_cohorts,
          "data/processed/cohort_immunity_profiles/age_stratified_cohorts.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/age_stratified_cohorts.csv")

# Save total 0-14 population
write_csv(total_014_pop,
          "data/processed/cohort_immunity_profiles/total_0_14_population.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/total_0_14_population.csv")

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

message("\n", rep("=", 60))
message("SUMMARY")
message(rep("=", 60))

message("\nBirth cohort data prepared:")
message("  Countries: ", n_distinct(birth_cohorts$country))
message("  Year range: ", min(birth_cohorts$year), "-", max(birth_cohorts$year))
message("  Total country-years: ", nrow(birth_cohorts))

message("\nAge-stratified cohorts (0-14 years):")
message("  Observation years: ",
        min(age_stratified_cohorts$observation_year), "-",
        max(age_stratified_cohorts$observation_year))
message("  Total records: ", nrow(age_stratified_cohorts))

message("\nFor cohort reconstruction, we need:")
message("  1. ✓ MCV1 and MCV2 coverage (COMPLETE)")
message("  2. ✓ Immunization schedules (COMPLETE)")
message("  3. ✓ Birth cohort populations (COMPLETE)")
message("  4. ⏳ SIA data (supplementary immunization activities)")

message("\nReady to proceed to:")
message("  Option A: Download SIA data (if available)")
message("  Option B: Begin cohort reconstruction without SIAs (baseline)")
message("\nNext script: python3 scripts/01_cohort/reconstruct_cohorts.py")
message("\n", rep("=", 60))

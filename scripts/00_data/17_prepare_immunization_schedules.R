# ==============================================================================
# Prepare Immunization Schedules for Cohort Mapping
# ==============================================================================
# Purpose: Extract MCV1 and MCV2 target ages from vaccine schedule data
#          to map calendar-year coverage to birth cohorts
#
# Input: data/raw/vaccine_schedule_data.csv
# Output: data/processed/cohort_immunity_profiles/measles_schedule.csv
# ==============================================================================

library(tidyverse)

message("=============================================================")
message("PHASE 2B: IMMUNIZATION SCHEDULE DATA PREPARATION")
message("=============================================================\n")

# ------------------------------------------------------------------------------
# Load Vaccine Schedule Data
# ------------------------------------------------------------------------------

message("Loading vaccine schedule data...")
schedule_raw <- read_csv("data/raw/vaccine_schedule_data.csv",
                         show_col_types = FALSE)

message("Raw schedule data: ", nrow(schedule_raw), " rows")
message("Columns: ", paste(names(schedule_raw), collapse = ", "))

# Check unique vaccines
unique_vaccines <- schedule_raw %>%
  distinct(vaccine) %>%
  arrange(vaccine)

message("\nUnique vaccines in data: ", nrow(unique_vaccines))

# ------------------------------------------------------------------------------
# Extract Measles Schedule
# ------------------------------------------------------------------------------

message("\nExtracting measles-containing vaccine schedules...")

# Filter for measles vaccines
mcv_schedule <- schedule_raw %>%
  filter(str_detect(vaccine, "(?i)measles|MCV|MMR|MR")) %>%
  select(
    country,
    country_name,
    who_region,
    year,
    vaccine = vaccine,
    schedule_round,
    target_population,
    geographic_area,
    age_administered,
    source_comment
  )

message("Measles schedule records: ", nrow(mcv_schedule))
message("Countries: ", n_distinct(mcv_schedule$country))

# Check what age formats we have
message("\nAge administration formats found:")
age_formats <- mcv_schedule %>%
  count(age_administered, sort = TRUE) %>%
  head(20)
print(age_formats)

# ------------------------------------------------------------------------------
# Parse Age Administered
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("PARSING TARGET AGES")
message(rep("-", 60))

parse_age_to_months <- function(age_str) {
  # Parse age strings like "M9", "Y1", "M15-18" to months
  if (is.na(age_str)) return(NA_real_)

  age_str <- trimws(age_str)

  # Months format: M9, M12, M15-18
  if (str_detect(age_str, "^M\\d+")) {
    age_num <- as.numeric(str_extract(age_str, "\\d+"))
    return(age_num)
  }

  # Years format: Y1, Y2
  if (str_detect(age_str, "^Y\\d+")) {
    age_num <- as.numeric(str_extract(age_str, "\\d+"))
    return(age_num * 12)
  }

  # Birth: BIRTH, B
  if (str_detect(age_str, "(?i)^B(IRTH)?$")) {
    return(0)
  }

  # Week format: W6
  if (str_detect(age_str, "^W\\d+")) {
    age_num <- as.numeric(str_extract(age_str, "\\d+"))
    return(age_num / 4.33)  # Approximate weeks to months
  }

  # Default: return NA
  return(NA_real_)
}

mcv_schedule_parsed <- mcv_schedule %>%
  mutate(
    target_age_months = map_dbl(age_administered, parse_age_to_months),
    target_age_years = target_age_months / 12
  )

# Check parsing success
parsing_summary <- mcv_schedule_parsed %>%
  summarise(
    total = n(),
    parsed = sum(!is.na(target_age_months)),
    unparsed = sum(is.na(target_age_months)),
    pct_parsed = 100 * parsed / total
  )

message("\nAge parsing results:")
print(parsing_summary)

# Show unparsed ages
unparsed_ages <- mcv_schedule_parsed %>%
  filter(is.na(target_age_months)) %>%
  count(age_administered, sort = TRUE)

if (nrow(unparsed_ages) > 0) {
  message("\nUnparsed age formats:")
  print(unparsed_ages)
}

# ------------------------------------------------------------------------------
# Classify by Dose
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("CLASSIFYING DOSES")
message(rep("-", 60))

mcv_schedule_classified <- mcv_schedule_parsed %>%
  mutate(
    dose = case_when(
      schedule_round == "1" ~ "MCV1",
      schedule_round == "2" ~ "MCV2",
      schedule_round == "3" ~ "MCV3",
      str_detect(vaccine, "(?i)^MCV1") ~ "MCV1",
      str_detect(vaccine, "(?i)^MCV2") ~ "MCV2",
      str_detect(vaccine, "(?i)1st") ~ "MCV1",
      str_detect(vaccine, "(?i)2nd") ~ "MCV2",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(dose %in% c("MCV1", "MCV2"))  # Keep only first two doses

dose_summary <- mcv_schedule_classified %>%
  group_by(dose) %>%
  summarise(
    n_records = n(),
    n_countries = n_distinct(country),
    mean_age_months = mean(target_age_months, na.rm = TRUE),
    median_age_months = median(target_age_months, na.rm = TRUE),
    min_age_months = min(target_age_months, na.rm = TRUE),
    max_age_months = max(target_age_months, na.rm = TRUE),
    .groups = "drop"
  )

message("\nDose classification summary:")
print(dose_summary)

# ------------------------------------------------------------------------------
# Create Country-Year Schedule Dataset
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("CREATING COUNTRY-YEAR SCHEDULE")
message(rep("-", 60))

# Get most recent schedule for each country (since schedules don't change often)
country_schedule <- mcv_schedule_classified %>%
  filter(!is.na(target_age_months)) %>%
  group_by(country, country_name, dose) %>%
  arrange(desc(year)) %>%
  slice(1) %>%  # Most recent schedule
  ungroup() %>%
  select(
    country,
    country_name,
    who_region,
    dose,
    target_age_months,
    target_age_years,
    vaccine,
    most_recent_year = year
  )

message("\nCountry schedules created:")
message("  Countries with MCV1 schedule: ",
        n_distinct(country_schedule$country[country_schedule$dose == "MCV1"]))
message("  Countries with MCV2 schedule: ",
        n_distinct(country_schedule$country[country_schedule$dose == "MCV2"]))

# Check for countries missing schedule data
coverage_countries <- read_csv("data/processed/cohort_immunity_profiles/measles_dose_coverage_wide.csv",
                               show_col_types = FALSE) %>%
  distinct(country, country_name)

missing_schedule <- coverage_countries %>%
  anti_join(
    country_schedule %>% filter(dose == "MCV1"),
    by = "country"
  )

if (nrow(missing_schedule) > 0) {
  message("\nWARNING: ", nrow(missing_schedule), " countries have coverage data but no schedule:")
  print(missing_schedule)
  message("\nWill use WHO default schedules:")
  message("  MCV1: 9 months (WHO recommendation for high-transmission settings)")
  message("  MCV2: 15-18 months (typical second dose timing)")
}

# ------------------------------------------------------------------------------
# Fill Missing Schedules with Defaults
# ------------------------------------------------------------------------------

# WHO default schedules
default_schedule <- tibble(
  dose = c("MCV1", "MCV2"),
  target_age_months = c(9, 15),  # Conservative defaults
  target_age_years = c(0.75, 1.25)
)

message("\nUsing default schedules for countries without data:")
print(default_schedule)

# Complete schedule for all countries
complete_schedule <- coverage_countries %>%
  crossing(default_schedule %>% select(dose)) %>%
  left_join(country_schedule, by = c("country", "country_name", "dose")) %>%
  mutate(
    target_age_months = coalesce(target_age_months,
                                 ifelse(dose == "MCV1", 9, 15)),
    target_age_years = target_age_months / 12,
    schedule_source = ifelse(is.na(vaccine), "WHO_default", "country_reported")
  ) %>%
  select(
    country,
    country_name,
    dose,
    target_age_months,
    target_age_years,
    schedule_source,
    vaccine,
    who_region,
    most_recent_year
  )

# ------------------------------------------------------------------------------
# Visualization
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("CREATING VISUALIZATIONS")
message(rep("-", 60))

# Plot 1: Distribution of target ages
p1_age_dist <- complete_schedule %>%
  ggplot(aes(x = target_age_months, fill = dose)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~dose, ncol = 1, scales = "free_y") +
  labs(
    title = "Distribution of MCV Target Ages Across Countries",
    x = "Target Age (months)",
    y = "Number of Countries",
    fill = "Dose"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("output/diagnostics/mcv_target_age_distribution.png",
       p1_age_dist, width = 10, height = 6)

# Plot 2: MCV1 vs MCV2 timing
schedule_wide <- complete_schedule %>%
  select(country, country_name, dose, target_age_months) %>%
  pivot_wider(names_from = dose, values_from = target_age_months)

p2_dose_timing <- schedule_wide %>%
  filter(!is.na(MCV1) & !is.na(MCV2)) %>%
  ggplot(aes(x = MCV1, y = MCV2)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "MCV1 vs MCV2 Target Ages Across Countries",
    subtitle = "Red line: MCV2 = MCV1 (impossible)",
    x = "MCV1 Target Age (months)",
    y = "MCV2 Target Age (months)"
  ) +
  theme_minimal()

ggsave("output/diagnostics/mcv_dose_timing_comparison.png",
       p2_dose_timing, width = 8, height = 6)

# ------------------------------------------------------------------------------
# Save Results
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("SAVING RESULTS")
message(rep("-", 60))

write_csv(complete_schedule,
          "data/processed/cohort_immunity_profiles/measles_schedule.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/measles_schedule.csv")

# Save schedule summary
schedule_summary <- complete_schedule %>%
  group_by(dose, schedule_source) %>%
  summarise(
    n_countries = n(),
    mean_age_months = mean(target_age_months),
    sd_age_months = sd(target_age_months),
    min_age_months = min(target_age_months),
    max_age_months = max(target_age_months),
    .groups = "drop"
  )

write_csv(schedule_summary,
          "data/processed/cohort_immunity_profiles/measles_schedule_summary.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/measles_schedule_summary.csv")

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

message("\n", rep("=", 60))
message("SUMMARY")
message(rep("=", 60))

message("\nSchedule data prepared for ", n_distinct(complete_schedule$country), " countries")
message("\nSchedule sources:")
source_summary <- complete_schedule %>%
  count(schedule_source, dose) %>%
  pivot_wider(names_from = dose, values_from = n)
print(source_summary)

message("\nTypical target ages:")
message("  MCV1: ", round(mean(complete_schedule$target_age_months[complete_schedule$dose == "MCV1"]), 1),
        " months (median: ",
        median(complete_schedule$target_age_months[complete_schedule$dose == "MCV1"]), ")")
message("  MCV2: ", round(mean(complete_schedule$target_age_months[complete_schedule$dose == "MCV2"]), 1),
        " months (median: ",
        median(complete_schedule$target_age_months[complete_schedule$dose == "MCV2"]), ")")

message("\nFor cohort reconstruction, we need:")
message("  1. ✓ MCV1 and MCV2 coverage (COMPLETE)")
message("  2. ✓ Immunization schedules (COMPLETE)")
message("  3. ⏳ Birth cohort populations")
message("  4. ⏳ SIA data")

message("\nNext script: 02c_prepare_birth_cohorts.R")
message("\n", rep("=", 60))

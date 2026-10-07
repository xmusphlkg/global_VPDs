# ==============================================================================
# Prepare Pertussis (DTP) Immunization Schedules
# ==============================================================================
# Purpose: Define DTP vaccination schedule for cohort mapping
#
# Output: data/processed/cohort_immunity_profiles/pertussis_schedule.csv
#
# Standard DTP Schedule:
#   - DTP1: 2 months
#   - DTP2: 4 months
#   - DTP3: 6 months (primary series completion)
#   - DTP4: 15-18 months (first booster)
# ==============================================================================

library(tidyverse)

message("=============================================================")
message("PERTUSSIS IMMUNIZATION SCHEDULE PREPARATION")
message("=============================================================\n")

# ------------------------------------------------------------------------------
# Load Coverage Countries
# ------------------------------------------------------------------------------

message("Loading coverage data to get country list...")
coverage_countries <- read_csv("data/processed/cohort_immunity_profiles/pertussis_dose_coverage_wide.csv",
                               show_col_types = FALSE) %>%
  distinct(country, country_name)

message("Countries: ", n_distinct(coverage_countries$country))

# ------------------------------------------------------------------------------
# WHO Standard DTP Schedule
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("WHO STANDARD DTP SCHEDULE")
message(rep("-", 60))

# WHO-recommended DTP schedule (most countries follow this)
default_schedule <- tibble(
  dose = c("DTP1", "DTP3"),
  target_age_months = c(2, 6),
  target_age_years = c(2/12, 6/12),
  schedule_note = c(
    "First dose of primary series",
    "Completion of primary series (3rd dose)"
  )
)

message("\nWHO standard DTP schedule:")
print(default_schedule)

message("\nNote: We track DTP1 and DTP3 for primary series.")
message("      DTP4 (booster) not available in WUENIC dataset.")
message("      Most countries administer DTP at 2, 4, 6 months.")

# ------------------------------------------------------------------------------
# Create Complete Schedule for All Countries
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("CREATING COUNTRY SCHEDULES")
message(rep("-", 60))

# Apply WHO standard schedule to all countries
# (In reality, schedules vary slightly, but WHO standard is widely adopted)
complete_schedule <- coverage_countries %>%
  crossing(default_schedule %>% select(dose, target_age_months, target_age_years)) %>%
  mutate(
    schedule_source = "WHO_standard"
  ) %>%
  select(
    country,
    country_name,
    dose,
    target_age_months,
    target_age_years,
    schedule_source
  )

message("\n✓ Schedules created for ", n_distinct(complete_schedule$country), " countries")
message("  DTP1: ", sum(complete_schedule$dose == "DTP1"), " countries")
message("  DTP3: ", sum(complete_schedule$dose == "DTP3"), " countries")

# ------------------------------------------------------------------------------
# Schedule Summary Statistics
# ------------------------------------------------------------------------------

schedule_summary <- complete_schedule %>%
  group_by(dose, schedule_source) %>%
  summarise(
    n_countries = n(),
    mean_age_months = mean(target_age_months),
    median_age_months = median(target_age_months),
    .groups = "drop"
  )

message("\nSchedule summary:")
print(schedule_summary)

# ------------------------------------------------------------------------------
# Save Results
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("SAVING RESULTS")
message(rep("-", 60))

write_csv(complete_schedule,
          "data/processed/cohort_immunity_profiles/pertussis_schedule.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/pertussis_schedule.csv")

write_csv(schedule_summary,
          "data/processed/cohort_immunity_profiles/pertussis_schedule_summary.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/pertussis_schedule_summary.csv")

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

message("\n", rep("=", 60))
message("SUMMARY")
message(rep("=", 60))

message("\nPertussis (DTP) schedule prepared for ", n_distinct(complete_schedule$country), " countries")
message("\nStandard timing:")
message("  DTP1: 2 months (first dose)")
message("  DTP3: 6 months (primary series completion)")

message("\nFor cohort reconstruction:")
message("  ✓ Coverage data prepared")
message("  ✓ Schedules prepared")
message("  ⏳ Next: Cohort reconstruction WITH waning immunity")

message("\n", rep("=", 60))

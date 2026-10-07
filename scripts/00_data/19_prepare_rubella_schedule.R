# ==============================================================================
# Prepare Rubella Immunization Schedules for Cohort Mapping
# ==============================================================================
# Purpose: Extract RCV1 and RCV2 target ages from vaccine schedule data
#          to map calendar-year coverage to birth cohorts
#
# Note: Rubella is often delivered via MMR vaccine, so we extract both
#       standalone rubella vaccines (MR) and MMR schedules
#
# Input: data/raw/vaccine_schedule_data.csv
# Output: data/processed/cohort_immunity_profiles/rubella_schedule.csv
# ==============================================================================

library(tidyverse)

message("=============================================================")
message("RUBELLA IMMUNIZATION SCHEDULE DATA PREPARATION")
message("=============================================================\n")

# ------------------------------------------------------------------------------
# Load Vaccine Schedule Data
# ------------------------------------------------------------------------------

message("Loading vaccine schedule data...")
schedule_raw <- read_csv("data/raw/vaccine_schedule_data.csv",
                         show_col_types = FALSE)

message("Raw schedule data: ", nrow(schedule_raw), " rows")

# ------------------------------------------------------------------------------
# Extract Rubella Schedule (including MMR and MR vaccines)
# ------------------------------------------------------------------------------

message("\nExtracting rubella-containing vaccine schedules...")

# Filter for rubella-containing vaccines (RCV, MMR, MR, MMRV)
rcv_schedule <- schedule_raw %>%
  filter(str_detect(vaccine, "(?i)rubella|RCV|MMR|^MR$|MMRV")) %>%
  # Focus on routine immunization
  filter(target_population == "General/routine") %>%
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

message("Rubella schedule records: ", nrow(rcv_schedule))
message("Countries: ", n_distinct(rcv_schedule$country))

# Check vaccine types
message("\nRubella-containing vaccine types found:")
vaccine_types <- rcv_schedule %>%
  count(vaccine, sort = TRUE) %>%
  head(10)
print(vaccine_types)

# ------------------------------------------------------------------------------
# Parse Age Administered
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("PARSING TARGET AGES")
message(rep("-", 60))

parse_age_to_months <- function(age_str) {
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

  return(NA_real_)
}

rcv_schedule_parsed <- rcv_schedule %>%
  mutate(
    target_age_months = map_dbl(age_administered, parse_age_to_months),
    target_age_years = target_age_months / 12
  )

parsing_summary <- rcv_schedule_parsed %>%
  summarise(
    total = n(),
    parsed = sum(!is.na(target_age_months)),
    unparsed = sum(is.na(target_age_months)),
    pct_parsed = 100 * parsed / total
  )

message("\nAge parsing results:")
print(parsing_summary)

# ------------------------------------------------------------------------------
# Classify by Dose
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("CLASSIFYING DOSES")
message(rep("-", 60))

rcv_schedule_classified <- rcv_schedule_parsed %>%
  mutate(
    dose = case_when(
      schedule_round == "1" ~ "RCV1",
      schedule_round == "2" ~ "RCV2",
      schedule_round == "3" ~ "RCV3",
      str_detect(vaccine, "(?i)^RCV1") ~ "RCV1",
      str_detect(vaccine, "(?i)^RCV2") ~ "RCV2",
      TRUE ~ paste0("RCV", schedule_round)
    )
  ) %>%
  filter(!is.na(target_age_months))

dose_distribution <- rcv_schedule_classified %>%
  count(dose, sort = TRUE)

message("\nDose distribution:")
print(dose_distribution)

# ------------------------------------------------------------------------------
# Take Most Recent Schedule Per Country
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("SELECTING CURRENT SCHEDULES")
message(rep("-", 60))

# For each country, take the most recent year's schedule
rcv_schedule_current <- rcv_schedule_classified %>%
  group_by(country, dose) %>%
  filter(year == max(year)) %>%
  # If multiple entries for same dose in same year, take median age
  summarise(
    country_name = first(country_name),
    who_region = first(who_region),
    year = first(year),
    target_age_months = median(target_age_months, na.rm = TRUE),
    vaccine = first(vaccine),
    .groups = "drop"
  ) %>%
  arrange(country, dose)

message("Current schedules: ", nrow(rcv_schedule_current), " country-dose combinations")

# Summary by dose
dose_age_summary <- rcv_schedule_current %>%
  group_by(dose) %>%
  summarise(
    n_countries = n(),
    mean_age_months = mean(target_age_months, na.rm = TRUE),
    median_age_months = median(target_age_months, na.rm = TRUE),
    min_age_months = min(target_age_months, na.rm = TRUE),
    max_age_months = max(target_age_months, na.rm = TRUE),
    .groups = "drop"
  )

message("\nAge distribution by dose:")
print(dose_age_summary)

# ------------------------------------------------------------------------------
# Fill Missing RCV2 with Default Values
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("HANDLING MISSING RCV2 SCHEDULES")
message(rep("-", 60))

# Countries with RCV1 but no RCV2
countries_with_rcv1 <- rcv_schedule_current %>%
  filter(dose == "RCV1") %>%
  pull(country)

countries_with_rcv2 <- rcv_schedule_current %>%
  filter(dose == "RCV2") %>%
  pull(country)

countries_missing_rcv2 <- setdiff(countries_with_rcv1, countries_with_rcv2)

message("Countries with RCV1: ", length(countries_with_rcv1))
message("Countries with RCV2: ", length(countries_with_rcv2))
message("Countries missing RCV2: ", length(countries_missing_rcv2))

# For countries missing RCV2, create a placeholder based on typical schedule
# (RCV2 typically given 3-6 months after RCV1, or at 15-18 months)
if (length(countries_missing_rcv2) > 0) {
  message("\nCreating default RCV2 schedules for countries without data...")
  message("Using default: RCV1 + 6 months or 18 months (whichever is later)")

  default_rcv2 <- rcv_schedule_current %>%
    filter(country %in% countries_missing_rcv2, dose == "RCV1") %>%
    mutate(
      dose = "RCV2",
      target_age_months = pmax(target_age_months + 6, 18),  # At least 18 months
      vaccine = paste0(vaccine, " (2nd dose, estimated)")
    )

  rcv_schedule_final <- bind_rows(
    rcv_schedule_current,
    default_rcv2
  ) %>%
    arrange(country, dose)
} else {
  rcv_schedule_final <- rcv_schedule_current
}

# ------------------------------------------------------------------------------
# Save Prepared Data
# ------------------------------------------------------------------------------

message("\n", rep("-", 60))
message("SAVING PREPARED DATA")
message(rep("-", 60))

write_csv(rcv_schedule_final, "data/processed/cohort_immunity_profiles/rubella_schedule.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/rubella_schedule.csv")

# Save summary
schedule_summary <- rcv_schedule_final %>%
  group_by(dose) %>%
  summarise(
    n_countries = n(),
    mean_age_months = mean(target_age_months),
    median_age_months = median(target_age_months),
    .groups = "drop"
  )

write_csv(schedule_summary, "data/processed/cohort_immunity_profiles/rubella_schedule_summary.csv")
message("✓ Saved: data/processed/cohort_immunity_profiles/rubella_schedule_summary.csv")

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

message("\n", rep("=", 60))
message("SUMMARY")
message(rep("=", 60))

message("\nRubella immunization schedules prepared:")
message("  • Countries with schedules: ", n_distinct(rcv_schedule_final$country))
message("  • RCV1 schedules: ", sum(rcv_schedule_final$dose == "RCV1"))
message("  • RCV2 schedules: ", sum(rcv_schedule_final$dose == "RCV2"))

message("\nTypical target ages:")
print(schedule_summary)

message("\nNext: python3 scripts/01_cohort/reconstruct_cohorts.py")
message(rep("=", 60), "\n")

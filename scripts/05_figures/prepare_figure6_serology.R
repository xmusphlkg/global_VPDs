# Plot-ready data for the external serology triangulation figure.
# NHANES and aggregate survey estimates remain separate estimands; this
# module only derives the declared within-source age contrasts. Kilifi DTP
# rows remain labelled external age profiles rather than pooled estimates.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

prepare_figure6_serology <- function(out = "output/revised") {
  read_result <- function(path, required) {
    x <- read_csv(file.path(out, path), show_col_types = FALSE, progress = FALSE)
    stopifnot(all(required %in% names(x)), nrow(x) > 0)
    x
  }

  endpoint_levels <- c("External age marker (%)", "Programme-derived gap (%)")
  disease_levels <- c("Measles", "Rubella", "Mumps", "Diphtheria", "Pertussis")

  sero <- read_result(
    "supplementary_tables/TableS15_nhanes_age_validation.csv",
    c("disease", "age", "seropositive", "seronegative",
      "seronegative_ci_low", "seronegative_ci_high", "model_gap",
      "profile_status", "n_valid_assay")
  ) %>%
    filter(disease %in% c("Measles", "Rubella", "Mumps"), age %in% 6:14) %>%
    mutate(
      disease = factor(disease, levels = disease_levels),
      profile_status = factor(profile_status, levels = c("known", "projected")),
      endpoint = "External age marker (%)"
    )
  stopifnot(
    nrow(sero) == 27,
    all(is.finite(as.matrix(sero[, c("age", "seronegative", "model_gap")])))
  )

  aligned <- bind_rows(
    sero %>% transmute(disease, age, profile_status, n_valid_assay,
      source = "NHANES 2017–March 2020",
      endpoint = factor("External age marker (%)", levels = endpoint_levels),
      estimate = seronegative, lower = seronegative_ci_low,
      upper = seronegative_ci_high),
    sero %>% transmute(disease, age, profile_status, n_valid_assay,
      source = "NHANES-aligned model profile, 2018",
      endpoint = factor("Programme-derived gap (%)", levels = endpoint_levels),
      estimate = model_gap, lower = NA_real_, upper = NA_real_)
  )

  kilifi_path <- file.path(out, "external_validation/kilifi_dtp_age_serology.csv")
  if (!file.exists(kilifi_path)) {
    kilifi_path <- file.path("data", "external", "kilifi_2009_2021",
      "kilifi_dtp_age_serology.csv")
  }
  kilifi <- read_csv(kilifi_path, show_col_types = FALSE, progress = FALSE)
  stopifnot(all(c("survey_year", "age_band", "age_mid", "n", "diph_minimal_pct",
    "diph_minimal_low", "diph_minimal_high", "pert_recent_12mo_pct",
    "pert_recent_6mo_pct") %in% names(kilifi)), nrow(kilifi) > 0)
  kilifi <- kilifi %>%
    filter(survey_year == 2021,
    age_band %in% c("1-4", "5-9", "10-14"))
  stopifnot(nrow(kilifi) == 3, all(is.finite(kilifi$age_mid)))

  dtp_top <- bind_rows(
    kilifi %>% transmute(
      disease = factor("Diphtheria", levels = disease_levels), age = age_mid,
      profile_status = factor("known", levels = c("known", "projected")),
      n_valid_assay = n, source = "Kilifi, Kenya, 2021",
      endpoint = factor("External age marker (%)", levels = endpoint_levels),
      estimate = diph_minimal_pct / 100, lower = diph_minimal_low / 100,
      upper = diph_minimal_high / 100),
    kilifi %>% transmute(
      disease = factor("Pertussis", levels = disease_levels), age = age_mid,
      profile_status = factor("known", levels = c("known", "projected")),
      n_valid_assay = n, source = "Kilifi, Kenya, 2021",
      endpoint = factor("External age marker (%)", levels = endpoint_levels),
      # The source has no validated protection cutoff. This sums its two
      # operational recent-exposure categories (62.5 EU/mL or higher).
      estimate = (pert_recent_12mo_pct + pert_recent_6mo_pct) / 100,
      lower = NA_real_, upper = NA_real_)
  )

  model_profiles <- read_result(
    "history_matching/eligible_history_panel.csv",
    c("country", "disease", "year", "scenario", "gap_waning_strict_1_4",
      "gap_waning_strict_5_9", "gap_waning_strict_10_14")
  ) %>% filter(country == "KEN", year == 2021,
    scenario == "component_strict", disease %in% c("Diphtheria", "Pertussis"))
  stopifnot(nrow(model_profiles) == 2)
  dtp_bottom <- model_profiles %>%
    pivot_longer(starts_with("gap_waning_strict_"),
      names_to = "age_band", values_to = "estimate") %>%
    mutate(
      age = recode(age_band, gap_waning_strict_1_4 = 2.5,
        gap_waning_strict_5_9 = 7, gap_waning_strict_10_14 = 12),
      disease = factor(disease, levels = disease_levels),
      profile_status = factor("known", levels = c("known", "projected")),
      n_valid_assay = NA_integer_, source = "Kenya model profile, 2021",
      endpoint = factor("Programme-derived gap (%)", levels = endpoint_levels),
      lower = NA_real_, upper = NA_real_) %>%
    select(disease, age, profile_status, n_valid_assay, source, endpoint,
      estimate, lower, upper)
  aligned <- bind_rows(aligned, dtp_top, dtp_bottom)

  known <- sero %>% filter(profile_status == "known") %>% arrange(disease, age)
  age_changes <- known %>%
    select(disease, age, profile_status, n_valid_assay, cutoff,
      model_gap, seronegative, seronegative_ci_low, seronegative_ci_high) %>%
    group_by(disease) %>% mutate(
    age_youngest = first(age), age_oldest = last(age),
    model_gap_baseline = first(model_gap),
    seronegative_baseline = first(seronegative)
  ) %>% ungroup() %>%
    pivot_longer(c(model_gap, seronegative),
      names_to = "endpoint", values_to = "estimate") %>%
    mutate(
      baseline = if_else(endpoint == "model_gap", model_gap_baseline,
        seronegative_baseline),
      change_pp = 100 * (estimate - baseline),
      endpoint = factor(recode(endpoint,
        model_gap = "Programme-derived gap", seronegative = "NHANES seronegativity"),
        levels = c("NHANES seronegativity", "Programme-derived gap"))
    )
  stopifnot(nrow(age_changes) == 42,
    all(is.finite(age_changes$change_pp)),
    all(filter(age_changes, age == age_youngest)$change_pp == 0))
  contrasts <- known %>% group_by(disease) %>% summarise(
    age_youngest = first(age), age_oldest = last(age),
    model_gap_change = 100 * (last(model_gap) - first(model_gap)),
    seronegativity_change = 100 * (last(seronegative) - first(seronegative)),
    spearman_rho = cor(model_gap, seronegative, method = "spearman"),
    n_known_ages = n(), .groups = "drop"
  ) %>%
    pivot_longer(c(model_gap_change, seronegativity_change),
      names_to = "endpoint", values_to = "change_pp") %>%
    mutate(endpoint = recode(endpoint,
      model_gap_change = "Programme-derived gap",
      seronegativity_change = "NHANES seronegativity"),
      endpoint = factor(endpoint,
        levels = c("Programme-derived gap", "NHANES seronegativity")),
      disease = factor(disease, levels = c("Measles", "Rubella", "Mumps")),
      rho_label = sprintf("Spearman ρ = %s%.2f; n = %d",
        ifelse(spearman_rho < 0, "−", "+"), abs(spearman_rho), n_known_ages))
  stopifnot(nrow(contrasts) == 6, all(is.finite(contrasts$change_pp)))

  read_external <- function(path, required) {
    x <- read_csv(file.path(out, "external_validation", path),
      show_col_types = FALSE, progress = FALSE)
    stopifnot(all(required %in% names(x)), nrow(x) > 0)
    x
  }
  esen2_measles <- read_external("esen2_measles_age_bands.csv",
    c("country", "age_lower", "age_upper", "seronegative_pct"))
  esen2_rubella <- read_external("esen2_rubella_age_bands.csv",
    c("country", "age_lower", "age_upper", "seronegative_pct"))
  esen2_mumps <- read_external("esen2_mumps_age_bands.csv",
    c("country", "age_lower", "age_upper", "seropositive_pct"))

  esen_contrast <- function(data, disease, young, old, value, positive = FALSE) {
    y <- data %>% filter(age_lower == young[1], age_upper == young[2]) %>%
      transmute(country, young = .data[[value]])
    o <- data %>% filter(age_lower == old[1], age_upper == old[2]) %>%
      transmute(country, old = .data[[value]])
    z <- inner_join(y, o, by = "country")
    if (!positive) z <- mutate(z, young = 100 - young, old = 100 - old)
    z %>% mutate(
      disease = factor(disease, levels = c("Measles", "Rubella", "Mumps")),
      source = "ESEN2", contrast_pp = old - young
    ) %>% filter(is.finite(contrast_pp))
  }
  esen <- bind_rows(
    esen_contrast(esen2_measles, "Measles", c(2, 4), c(10, 19), "seronegative_pct"),
    esen_contrast(esen2_rubella, "Rubella", c(2, 4), c(10, 14), "seronegative_pct"),
    esen_contrast(esen2_mumps, "Mumps", c(1, 4), c(10, 14), "seropositive_pct", TRUE)
  )
  stopifnot(nrow(filter(esen, disease == "Measles")) == 18,
    nrow(filter(esen, disease == "Rubella")) == 15,
    nrow(filter(esen, disease == "Mumps")) == 16)
  esen_summary <- esen %>% group_by(disease) %>% summarise(
    median_pp = median(contrast_pp), q25_pp = quantile(contrast_pp, .25),
    q75_pp = quantile(contrast_pp, .75), n = n(),
    n_positive = sum(contrast_pp > 0), .groups = "drop"
  ) %>% mutate(summary_label = sprintf("%d/%d positive", n_positive, n))

  aggregate_sero <- read_result(
    "supplementary_tables/TableS31_aggregate_serology_summary.csv",
    c("source", "disease", "contrast", "estimate_pp", "p25_pp", "p75_pp")
  )
  forest <- aggregate_sero %>% filter(
    (grepl("^ESEN2", source) & !grepl("standard age", contrast)) |
    (grepl("^KiGGS", source) & grepl("positive percentage", contrast)) |
    (grepl("^PIENTER2", source) & grepl("^mean ages", contrast))
  ) %>% mutate(
    disease = factor(disease, levels = c("Measles", "Rubella", "Mumps")),
    source_short = case_when(grepl("^ESEN2", source) ~ "ESEN2",
      grepl("^KiGGS", source) ~ "KiGGS", grepl("^PIENTER2", source) ~ "PIENTER2"),
    source_short = factor(source_short, levels = c("ESEN2", "KiGGS", "PIENTER2")),
    y = case_when(
      disease == "Measles" & source_short == "ESEN2" ~ 6,
      disease == "Measles" & source_short == "KiGGS" ~ 5,
      disease == "Rubella" & source_short == "ESEN2" ~ 3.5,
      disease == "Rubella" & source_short == "KiGGS" ~ 2.5,
      disease == "Mumps" & source_short == "ESEN2" ~ 1,
      disease == "Mumps" & source_short == "KiGGS" ~ 0,
      disease == "Mumps" & source_short == "PIENTER2" ~ -1
    ),
    row_label = case_when(
      source_short == "ESEN2" ~ paste0(as.character(disease), "\n  ", as.character(source_short)),
      TRUE ~ paste0("  ", as.character(source_short))
    )
  )
  stopifnot(nrow(forest) == 7, all(is.finite(forest$estimate_pp)))

  esen_rows <- forest %>% filter(source_short == "ESEN2") %>% select(disease, y)
  esen <- left_join(esen, esen_rows, by = "disease") %>%
    arrange(disease, contrast_pp, country) %>% group_by(disease) %>%
    mutate(y_jitter = y + rep(c(-.14, 0, .14), length.out = n())) %>% ungroup()
  esen_summary <- left_join(esen_summary, esen_rows, by = "disease")
  esen_reference <- filter(forest, source_short == "ESEN2")
  stopifnot(all(abs(esen_summary$median_pp - esen_reference$estimate_pp) < 1e-8),
    all(abs(esen_summary$q25_pp - esen_reference$p25_pp) < 1e-8),
    all(abs(esen_summary$q75_pp - esen_reference$p75_pp) < 1e-8))

  dtp <- aligned %>% filter(disease %in% c("Diphtheria", "Pertussis")) %>%
    mutate(
      age_band = factor(recode(as.character(age), `2.5` = "1-4",
        `7` = "5-9", `12` = "10-14"), levels = c("1-4", "5-9", "10-14")),
      estimate_pct = 100 * estimate, lower_pct = 100 * lower,
      upper_pct = 100 * upper,
      marker_construct = case_when(
        endpoint == "Programme-derived gap (%)" ~ "Programme-derived gap",
        disease == "Diphtheria" ~ "Minimum antitoxin protection",
        TRUE ~ "Recent anti-PT exposure (>=62.5 EU/mL)"
      )
    )
  stopifnot(nrow(dtp) == 12, !anyNA(dtp$age_band),
    all(is.finite(dtp$estimate_pct)),
    all(dtp$lower_pct <= dtp$estimate_pct, na.rm = TRUE),
    all(dtp$upper_pct >= dtp$estimate_pct, na.rm = TRUE))

  list(aligned = aligned, age_changes = age_changes, contrasts = contrasts,
    esen = esen, esen_summary = esen_summary, forest = forest, dtp = dtp)
}

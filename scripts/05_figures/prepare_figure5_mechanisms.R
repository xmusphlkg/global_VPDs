# Data preparation for the fixed-parameter age-structured mechanism figure.
# The outputs are structural summaries of the canonical age mechanism; they are
# not posterior draws or observed case predictions.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
})

prepare_figure5_mechanisms <- function(out = "output/revised") {
  age_path <- file.path(out, "age_seirv", "annual_age_seirv_results.csv")
  sens_path <- file.path(out, "age_seirv", "annual_age_seirv_sensitivity.csv")
  age <- read_csv(age_path, show_col_types = FALSE, progress = FALSE)
  sens <- read_csv(sens_path, show_col_types = FALSE, progress = FALSE)
  manifest <- jsonlite::fromJSON(file.path(out, "age_seirv", "manifest.json"))
  required_age <- c("country", "year", "disease", "mode", "gap_0_4", "gap_5_9",
                    "gap_10_14", "adult_gap", "age_gradient_10_14_minus_0_4",
                    "effective_R_initial")
  required_sens <- c("country", "year", "disease", "mode",
                     "adult_susceptibility_scenario", "transmission_blocking_scenario",
                     "effective_R_initial")
  stopifnot(all(required_age %in% names(age)), all(required_sens %in% names(sens)))
  disease_levels <- c("Measles", "Rubella", "Pertussis", "Diphtheria", "Mumps")
  mechanism_diseases <- c("Pertussis", "Diphtheria", "Mumps")
  stopifnot(
    nrow(age) > 0, nrow(sens) > 0,
    setequal(unique(age$disease), disease_levels),
    setequal(unique(age$mode), c("unwaning", "waning")),
    setequal(unique(age$year), 2014:2023),
    all(is.finite(as.matrix(age[, c("gap_0_4", "gap_5_9", "gap_10_14",
      "adult_gap", "age_gradient_10_14_minus_0_4", "effective_R_initial")]))),
    all(is.finite(as.matrix(sens[, c("adult_susceptibility_scenario",
      "transmission_blocking_scenario", "effective_R_initial")]))))
  # Use the assumptions that generated the completed mechanism outputs.
  blocking <- vapply(manifest$parameters[disease_levels],
    function(parameters) parameters$transmission_blocking, numeric(1))
  source_reduction <- vapply(manifest$parameters[disease_levels],
    function(parameters) if (is.null(parameters$source_reduction)) 0 else parameters$source_reduction,
    numeric(1))
  stopifnot(all(is.finite(blocking)), all(blocking >= 0 & blocking <= 1))
  age <- age %>% mutate(
    disease = factor(disease, levels = disease_levels),
    mode = factor(mode, levels = c("unwaning", "waning")))
  sens <- sens %>% mutate(disease = factor(disease, levels = disease_levels),
                          mode = factor(mode, levels = c("unwaning", "waning")))
  countries <- unique(age$country)

  gradient <- age %>%
    group_by(disease, year, mode) %>%
    summarise(gradient_pp = 100 * mean(age_gradient_10_14_minus_0_4), .groups = "drop")
  # Paired country means retain absolute age-band gaps and country heterogeneity.
  age_plane <- age %>% group_by(country, country_name, disease, mode) %>%
    summarise(younger_gap_pct = 100 * mean(gap_0_4),
      older_gap_pct = 100 * mean(gap_10_14),
      n_years = n_distinct(year), .groups = "drop")
  stopifnot(nrow(age_plane) == length(countries) * 5L * 2L,
    all(age_plane$n_years == 10L),
    all(is.finite(age_plane$younger_gap_pct)),
    all(is.finite(age_plane$older_gap_pct)))
  endpoint <- gradient %>% group_by(disease, mode) %>%
    summarise(gradient_pp = mean(gradient_pp), .groups = "drop")

  country <- age %>%
    group_by(country, disease, mode) %>%
    summarise(gradient_pp = 100 * mean(age_gradient_10_14_minus_0_4),
              Re = mean(effective_R_initial), .groups = "drop") %>%
    pivot_wider(names_from = mode, values_from = c(gradient_pp, Re)) %>%
    mutate(delta_gradient_pp = gradient_pp_waning - gradient_pp_unwaning,
           delta_Re = Re_waning - Re_unwaning) %>%
    left_join(age %>% distinct(country, country_name), by = "country")

  sensitivity <- sens %>%
    filter(disease %in% mechanism_diseases) %>%
    group_by(disease, adult_susceptibility_scenario,
             transmission_blocking_scenario, mode) %>%
    summarise(Re = mean(effective_R_initial), .groups = "drop") %>%
    pivot_wider(names_from = mode, values_from = Re) %>%
    mutate(delta_Re = waning - unwaning,
           relative_delta_Re_pct = 100 * delta_Re / unwaning)

  # Every central grid cell must reproduce the corresponding annual output.
  central <- sensitivity %>%
    filter(abs(adult_susceptibility_scenario - .30) < 1e-10,
      abs(transmission_blocking_scenario - blocking[as.character(disease)]) < 1e-10) %>%
    left_join(country %>% group_by(disease) %>%
      summarise(annual_unwaning = mean(Re_unwaning),
        annual_waning = mean(Re_waning), .groups = "drop"), by = "disease")
  stopifnot(nrow(central) == length(mechanism_diseases),
    all(abs(central$unwaning - central$annual_unwaning) < 1e-10),
    all(abs(central$waning - central$annual_waning) < 1e-10))

  # Single-age cohort inputs are shown before the model's population-weighted
  # aggregation to 0–4, 5–9 and 10–14; the infant scenario is not in this panel.
  read_cohort <- function(filename, diseases) {
    read_csv(file.path(Sys.getenv("VPDS_COHORT_DIR", "output/five_disease_revision"), filename),
      col_select = c(country, observation_year, age, disease,
        gap_no_waning_projected, gap_waning_projected),
      show_col_types = FALSE, progress = FALSE) %>%
      filter(country %in% countries, observation_year %in% 2014:2023,
        disease %in% diseases, age %in% 1:14)
  }
  # The shared table has strict-component mumps histories. Figure 4 and the
  # age mechanism use the separately labelled all-MCV proxy instead.
  cohort <- bind_rows(
    read_cohort("cohort_level_all_ages.csv", c("Pertussis", "Diphtheria")),
    read_cohort("mumps_all_mcv_proxy_cohort_level.csv", "Mumps"))
  cohort <- cohort %>%
    mutate(profile_fallback = is.na(gap_no_waning_projected) | is.na(gap_waning_projected))
  stopifnot(all(cohort$disease[cohort$profile_fallback] == "Mumps"))
  cohort <- cohort %>% mutate(
    gap_no_waning_projected = if_else(profile_fallback,
      0.30, gap_no_waning_projected),
    gap_waning_projected = if_else(profile_fallback,
      0.30, gap_waning_projected))
  stopifnot(nrow(cohort) == length(mechanism_diseases) * length(countries) * 10L * 14L,
    !anyDuplicated(cohort[c("country", "observation_year", "age", "disease")]),
    all(is.finite(cohort$gap_no_waning_projected)),
    all(is.finite(cohort$gap_waning_projected)))
  single_age <- cohort %>%
    pivot_longer(c(gap_no_waning_projected, gap_waning_projected),
      names_to = "mode", values_to = "gap") %>%
    mutate(mode = recode(mode, gap_no_waning_projected = "unwaning",
      gap_waning_projected = "waning"),
      susceptibility = 1 - blocking[as.character(disease)] * (1 - gap)) %>%
    group_by(disease, observation_year, age, mode) %>%
    summarise(susceptibility = mean(susceptibility),
      n_profile_fallback = sum(profile_fallback), .groups = "drop") %>%
    rename(year = observation_year) %>%
    mutate(disease = factor(disease, levels = mechanism_diseases))
  single_age_delta <- single_age %>%
    pivot_wider(names_from = mode, values_from = susceptibility) %>%
    mutate(delta = waning - unwaning)

  # Vaccine-action routing: map the band gaps used by the contact kernel below
  # to recipient acquisition s and source infectiousness t with the same
  # formulae, so each channel is the operator input that changes under waning.
  # Adult gaps are identical across modes and are omitted.
  channel_bands <- c(gap_0_4 = "0–4", gap_5_9 = "5–9", gap_10_14 = "10–14")
  channels <- age %>% filter(disease %in% mechanism_diseases) %>%
    select(country, year, disease, mode, all_of(names(channel_bands))) %>%
    pivot_longer(all_of(names(channel_bands)), names_to = "age_band", values_to = "gap") %>%
    mutate(b = blocking[as.character(disease)],
      r = source_reduction[as.character(disease)],
      s = 1 - b * (1 - gap),
      t = 1 - r * (1 - gap) * (1 - b) / pmax(s, 1e-12)) %>%
    group_by(country, disease, age_band, mode) %>%
    summarise(s = mean(s), t = mean(t), .groups = "drop") %>%
    pivot_wider(names_from = mode, values_from = c(s, t)) %>%
    transmute(country, disease = factor(as.character(disease), levels = mechanism_diseases),
      age_band = factor(channel_bands[age_band], levels = channel_bands),
      delta_acquisition_pp = 100 * (s_waning - s_unwaning),
      delta_infectiousness_pp = 100 * (t_waning - t_unwaning))
  stopifnot(nrow(channels) == length(mechanism_diseases) * length(countries) * 3L,
    all(is.finite(channels$delta_acquisition_pp)),
    all(is.finite(channels$delta_infectiousness_pp)))

  # The saved reciprocal matrices use 2019 demography as a contact reference.
  # Combine them with annual susceptibility inputs and retain the native 16
  # contact bands; this is a structural illustration, not a fitted matrix.
  age_groups <- c("0_4", "5_9", "10_14", "15_19", "20_24", "25_29", "30_34",
                  "35_39", "40_44", "45_49", "50_54", "55_59", "60_64",
                  "65_69", "70_74", "75_plus")
  age_labels <- c("0–4", "5–9", "10–14", "15–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75+")
  matrix_dir <- file.path(out, "age_seirv", "contact_matrices_reciprocal")
  matrices <- setNames(lapply(countries, function(country) {
    matrix_data <- read_csv(file.path(matrix_dir, paste0(country, ".csv")),
      show_col_types = FALSE, progress = FALSE, name_repair = "minimal")
    stopifnot(identical(names(matrix_data)[-1], age_groups),
      identical(matrix_data[[1]], age_groups))
    matrix <- as.matrix(matrix_data[, -1])
    stopifnot(identical(dim(matrix), c(16L, 16L)), all(is.finite(matrix)))
    matrix
  }), countries)
  kernel_rows <- list()
  for (i in seq_len(nrow(age))) {
    r <- age[i, ]
    if (!as.character(r$disease) %in% mechanism_diseases) next
    profile <- c(r$gap_0_4, r$gap_5_9, r$gap_10_14, rep(r$adult_gap, 13))
    b <- blocking[[as.character(r$disease)]]
    s <- 1 - b * (1 - profile)
    source <- 1 - source_reduction[[as.character(r$disease)]] *
      (1 - profile) * (1 - b) / pmax(s, 1e-12)
    kernel_rows[[length(kernel_rows) + 1L]] <- list(
      disease = as.character(r$disease), mode = as.character(r$mode),
      country = as.character(r$country),
      K = diag(s) %*% matrices[[as.character(r$country)]] %*% diag(source))
  }
  kernel_mean <- list()
  # Country-mean operators are the resampling units for the Fig. 5e bootstrap.
  # Every country contributes the same ten years, so the equal country-year
  # mean equals the mean of the country means.
  kernel_country_mean <- list()
  for (dis in mechanism_diseases) for (mode in c("unwaning", "waning")) {
    mats <- lapply(kernel_rows, function(x) if (x$disease == dis && x$mode == mode) x$K else NULL)
    mats <- mats[!vapply(mats, is.null, logical(1))]
    key <- paste(dis, mode, sep = "_")
    kernel_mean[[key]] <- Reduce("+", mats) / length(mats)
    kernel_country_mean[[key]] <- setNames(lapply(countries, function(cty) {
      own <- lapply(kernel_rows, function(x)
        if (x$disease == dis && x$mode == mode && x$country == cty) x$K else NULL)
      own <- own[!vapply(own, is.null, logical(1))]
      stopifnot(length(own) == 10L)
      Reduce("+", own) / length(own)
    }), countries)
    stopifnot(max(abs(Reduce("+", kernel_country_mean[[key]]) / length(countries) -
      kernel_mean[[key]])) < 1e-10)
  }
  k_scale <- max(vapply(kernel_mean, max, numeric(1)))
  k_delta <- map_dfr(mechanism_diseases, function(dis) {
    z <- (kernel_mean[[paste(dis, "waning", sep = "_")]] -
      kernel_mean[[paste(dis, "unwaning", sep = "_")]]) / k_scale
    expand_grid(recipient_age = age_labels, source_age = age_labels) %>%
      mutate(disease = dis, value = as.vector(t(z)))
  })
  k_recipient <- k_delta %>% group_by(disease, recipient_age) %>%
    summarise(value = sum(value), .groups = "drop")
  k_source <- k_delta %>% group_by(disease, source_age) %>%
    summarise(value = sum(value), .groups = "drop")

  dominant_shares <- function(K) {
    er <- eigen(K)
    el <- eigen(t(K))
    u <- abs(Re(er$vectors[, which.max(Re(er$values))]))
    v <- abs(Re(el$vectors[, which.max(Re(el$values))]))
    c_i <- v * (K %*% u)[, 1]
    100 * c_i / sum(c_i)
  }
  contribution <- map_dfr(mechanism_diseases, function(dis) map_dfr(c("unwaning", "waning"), function(mode) {
    # Evaluate outside tibble(): its data mask would recycle `mode` to 16 rows.
    shares <- dominant_shares(kernel_mean[[paste(dis, mode, sep = "_")]])
    tibble(disease = dis, mode = mode, age_band = age_labels, contribution = shares)
  })) %>% mutate(mode = factor(mode, levels = c("unwaning", "waning")))

  contribution_delta <- contribution %>%
    pivot_wider(names_from = mode, values_from = contribution) %>%
    mutate(delta_contribution_pp = waning - unwaning)

  # Fixed-parameter percentile bootstrap intervals for Fig. 5e and 5f. The
  # replicate count and seed match the Table 1 country bootstrap. Fig. 5e
  # resamples countries (paired across diseases and modes) and recomputes the
  # mean operator and its dominant mode; Fig. 5f resamples years within each
  # country, with the same year draws shared across countries and diseases.
  boot_B <- 2000L
  boot_seed <- 20261006L
  n_boot_countries <- length(countries)
  boot_years <- sort(unique(age$year))
  n_boot_years <- length(boot_years)
  set.seed(boot_seed)
  country_draws <- matrix(sample.int(n_boot_countries, boot_B * n_boot_countries,
    replace = TRUE), nrow = boot_B)
  year_draws <- matrix(sample.int(n_boot_years, boot_B * n_boot_years,
    replace = TRUE), nrow = boot_B)
  percentile <- function(x) unname(quantile(x, c(.025, .975), type = 7))

  contribution_boot <- map_dfr(mechanism_diseases, function(dis) {
    deltas <- t(vapply(seq_len(boot_B), function(b) {
      idx <- country_draws[b, ]
      shares <- lapply(c("unwaning", "waning"), function(mode) {
        own <- kernel_country_mean[[paste(dis, mode, sep = "_")]][idx]
        dominant_shares(Reduce("+", own) / n_boot_countries)
      })
      shares[[2]] - shares[[1]]
    }, numeric(length(age_labels))))
    stopifnot(identical(dim(deltas), c(boot_B, length(age_labels))),
      all(is.finite(deltas)))
    limits <- apply(deltas, 2, percentile)
    tibble(disease = dis, age_band = age_labels,
      delta_contribution_lower_95 = limits[1, ],
      delta_contribution_upper_95 = limits[2, ])
  })
  contribution_delta <- contribution_delta %>%
    left_join(contribution_boot, by = c("disease", "age_band")) %>%
    mutate(resampling_unit = "country", n_countries = n_boot_countries,
      bootstrap_replicates = boot_B)
  stopifnot(all(is.finite(contribution_delta$delta_contribution_lower_95)),
    all(contribution_delta$delta_contribution_lower_95 <=
      contribution_delta$delta_contribution_pp),
    all(contribution_delta$delta_contribution_pp <=
      contribution_delta$delta_contribution_upper_95))

  annual_delta <- age %>%
    filter(disease %in% mechanism_diseases) %>%
    select(country, year, disease, mode, effective_R_initial) %>%
    pivot_wider(names_from = mode, values_from = effective_R_initial) %>%
    mutate(delta_Re = waning - unwaning) %>%
    arrange(country, disease, year)
  country_boot <- annual_delta %>%
    group_by(country, disease) %>%
    summarise(limits = {
      stopifnot(n() == n_boot_years, identical(as.numeric(year), as.numeric(boot_years)))
      list(percentile(rowMeans(matrix(delta_Re[year_draws], nrow = boot_B))))
    }, .groups = "drop") %>%
    mutate(delta_Re_lower_95 = vapply(limits, `[`, numeric(1), 1),
      delta_Re_upper_95 = vapply(limits, `[`, numeric(1), 2)) %>%
    select(-limits)
  stopifnot(nrow(country_boot) == n_boot_countries * length(mechanism_diseases),
    all(is.finite(country_boot$delta_Re_lower_95)),
    all(country_boot$delta_Re_lower_95 <= country_boot$delta_Re_upper_95))
  country <- country %>%
    left_join(country_boot, by = c("country", "disease")) %>%
    mutate(resampling_unit = if_else(disease %in% mechanism_diseases, "year", NA_character_),
      n_years = n_boot_years,
      bootstrap_replicates = if_else(disease %in% mechanism_diseases, boot_B, NA_integer_))
  stopifnot(all(country$delta_Re_lower_95[country$disease %in% mechanism_diseases] <=
      country$delta_Re[country$disease %in% mechanism_diseases]),
    all(country$delta_Re[country$disease %in% mechanism_diseases] <=
      country$delta_Re_upper_95[country$disease %in% mechanism_diseases]))
  bootstrap <- list(
    analysis = "Figure 5 fixed-parameter percentile bootstrap intervals (panels e and f)",
    replicates = boot_B, seed = boot_seed,
    generator = "R set.seed() + sample.int() (Mersenne-Twister); country draws then year draws from one stream",
    interval = "2.5th and 97.5th percentiles; R type-7 (linear) quantiles",
    panel_b = list(statistic = "delta_contribution_pp", resampling_unit = "country",
      n_units = n_boot_countries,
      scheme = "Resample the ten countries with replacement, average their country-mean operators, recompute the dominant right/left eigenvectors and normalized shares for both modes, take waning minus unwaning; draws shared across diseases and modes"),
    panel_c = list(statistic = "delta_Re", resampling_unit = "year within country",
      n_units = n_boot_years,
      scheme = "Resample the ten paired annual waning-minus-unwaning contrasts within each country with replacement and average; draws shared across countries and diseases"),
    countries_in_resampling_order = countries, years = boot_years,
    interval_scope = "Variation from resampling the selected countries (b) or years (c) conditional on fixed R0, acquisition blocking, adult and infant protection states, cohort reconstruction, demography and contact matrices; not posterior uncertainty")

  list(blocking = blocking, mechanism_diseases = mechanism_diseases,
       surface = single_age, surface_delta = single_age_delta, channels = channels,
       gradient = gradient, endpoint = endpoint, age_plane = age_plane, country = country,
       kernel_delta = k_delta, kernel_recipient = k_recipient,
       kernel_source = k_source, contribution = contribution,
       contribution_delta = contribution_delta, sensitivity = sensitivity,
       bootstrap = bootstrap)
}

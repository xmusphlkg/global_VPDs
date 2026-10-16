# Plot-ready data for the cohort-memory atlas. No case or transmission refits.
prepare_figure1_atlas <- function(out = "output/revised") {
  read_result <- function(path) readr::read_csv(file.path(out, path), show_col_types = FALSE, progress = FALSE)
  exact <- read_result("history_matching/exact_matched_country_years.csv")
  history <- read_result("history_matching/eligible_history_panel.csv")
  matching <- read_result("history_matching/matching_summary.csv")
  cohorts <- readr::read_csv("output/five_disease_revision/cohort_level_all_ages.csv",
    col_select = c(country, observation_year, age, birth_year, cohort_size, disease,
                   gap_waning_strict, p0_strict, p1_strict, p2_strict,
                   effective_dose1_year, effective_dose2_year),
    show_col_types = FALSE, progress = FALSE)
  countries <- c("TZA", "UKR")
  levels <- pub_disease_levels
  primary <- matching %>% filter(configuration == "exact_equal_country_year")
  stopifnot(nrow(primary) == 5, !anyDuplicated(primary$disease))

  # A: observed Ukrainian programme history makes the low-coverage cohort band visible.
  a_heat <- cohorts %>% filter(country == "UKR", disease == "Pertussis",
    observation_year %in% 2014:2023, age %in% 1:14) %>%
    transmute(layer = "landscape", country, year = observation_year, age,
              birth_cohort = birth_year, gap_pp = gap_waning_strict * 100)
  a_coverage <- history %>% filter(country == "UKR", disease == "Pertussis", year %in% 2014:2023) %>%
    transmute(layer = "coverage", country, year, coverage)
  a_profile <- a_heat %>% filter(year == 2023) %>% mutate(layer = "current_age_profile")
  a_cohorts <- tibble(birth_cohort = seq(2002, 2020, 3)) %>%
    mutate(layer = "cohort_diagonal", x = pmax(2014, birth_cohort + 1),
           xend = pmin(2023, birth_cohort + 14), y = x - birth_cohort,
           yend = xend - birth_cohort) %>% filter(x < xend)
  a_bands <- a_coverage %>% filter(coverage < 30) %>%
    transmute(layer = "low_coverage_year", year, xmin = year - .5, xmax = year + .5)
  stopifnot(nrow(a_heat) == 140, nrow(a_profile) == 14, nrow(a_coverage) == 10,
            all(is.finite(a_heat$gap_pp)), setequal(a_bands$year, 2014:2016))

  # B: a reproducible ±2-percentage-point matched endpoint set, independently
  # of case-model eligibility. Figure 1b does not plot the separate 5-point
  # matching sensitivity used in panel D.
  endpoint_band_pp <- 2
  b_endpoints <- history %>% filter(disease == "Pertussis", year == 2023,
    abs(coverage - 83) <= endpoint_band_pp) %>% distinct(country, endpoint_coverage = coverage)
  b_coverage <- history %>% filter(disease == "Pertussis", year %in% 2014:2023) %>%
    inner_join(b_endpoints, by = "country") %>%
    transmute(layer = "coverage_history", country, year, coverage, endpoint_coverage,
              focal = country %in% countries)
  b_deficit <- b_coverage %>% filter(focal) %>%
    mutate(layer = "annual_deficit", deficit_pp = 100 - coverage)
  stopifnot(nrow(filter(b_coverage, focal)) == 20,
    all(filter(b_coverage, focal, year == 2023)$coverage == 83),
    all(abs(b_coverage$endpoint_coverage - 83) <= endpoint_band_pp))

  # C: repeat the established DTP reconstruction formula on its nine sensitivity
  # scenarios. The range is parameter sensitivity, never a sampling interval.
  c_input <- cohorts %>% filter(country %in% countries, disease == "Pertussis",
    observation_year == 2023, age %in% 1:14) %>%
    mutate(t1 = pmax(0, observation_year - coalesce(effective_dose1_year, as.numeric(observation_year))),
           t2 = pmax(0, observation_year - coalesce(effective_dose2_year, as.numeric(observation_year))))
  parameter_grid <- read_result("dtp_parameter_sensitivity.csv") %>% filter(disease == "Pertussis") %>%
    distinct(tau, partial_series_coefficient, full_series_coefficient)
  c_grid <- tidyr::crossing(c_input, parameter_grid) %>%
    mutate(layer = "parameter_scenario", gap_pp = 100 * (p0_strict +
      p1_strict * (1 - partial_series_coefficient * exp(-t1 / tau)) +
      p2_strict * (1 - full_series_coefficient * exp(-t2 / tau))))
  # Compute the envelope before naming the baseline to avoid sequential summarise masking.
  c_profile <- c_grid %>% group_by(country, age, birth_year, cohort_size) %>%
    summarise(lower_pp = min(gap_pp), upper_pp = max(gap_pp),
              gap_pp = first(gap_waning_strict) * 100, .groups = "drop") %>% mutate(layer = "current_profile")
  c_groups <- c_input %>% mutate(age_group = case_when(age <= 4 ~ "1–4", age <= 9 ~ "5–9", TRUE ~ "10–14")) %>%
    group_by(country, age_group) %>% summarise(age = mean(age),
      gap_pp = weighted.mean(gap_waning_strict * 100, cohort_size), .groups = "drop") %>% mutate(layer = "age_group_mean")
  c_total <- c_input %>% group_by(country) %>%
    summarise(gap_pp = weighted.mean(gap_waning_strict * 100, cohort_size), .groups = "drop") %>%
    mutate(layer = "population_weighted_gap")
  reference <- read_result("current_coverage_matched_examples.csv") %>% filter(country %in% countries, disease == "Pertussis", year == 2023)
  c_base <- filter(c_grid, tau == 5, partial_series_coefficient == .5)
  stopifnot(nrow(c_grid) == 252, nrow(c_profile) == 28,
    max(abs(c_base$gap_pp - c_base$gap_waning_strict * 100)) < 1e-8,
    all(c_profile$lower_pp <= c_profile$gap_pp), all(c_profile$upper_pp >= c_profile$gap_pp),
    max(abs(left_join(c_total, reference, by = "country") %>% transmute(e = gap_pp / 100 - model_gap_waning) %>% pull(e))) < 1e-8)

  # D: separate residual distributions, exact SD intervals and sensitivity glyphs.
  row_y <- function(d) length(levels) + 1 - match(d, levels)
  d_points <- exact %>% transmute(layer = "residual", disease, country, year,
    residual_pp = within_stratum_gap_residual_pp, y = row_y(disease))
  set.seed(20261006)
  d_points$point_y <- d_points$y - .19 + runif(nrow(d_points), -.065, .065)
  d_density <- bind_rows(lapply(levels, function(dis) {
    v <- filter(d_points, disease == dis)$residual_pp
    z <- density(v, from = min(v), to = max(v), n = 160)
    tibble(layer = "half_density", disease = dis, residual_pp = z$x,
           y = row_y(dis) + .04, upper_y = row_y(dis) + .04 + .30 * z$y / max(z$y))
  }))
  d_quartiles <- d_points %>% group_by(disease, y) %>% summarise(
    median_pp = median(residual_pp), q25_pp = quantile(residual_pp, .25),
    q75_pp = quantile(residual_pp, .75), n = n(), .groups = "drop") %>% mutate(layer = "residual_quartiles")
  configurations <- c("exact_equal_country_year", "bin_2pp_equal_country_year",
                       "bin_5pp_equal_country_year", "exact_birth_cohort_weighted")
  d_sd <- matching %>% filter(configuration %in% configurations) %>%
    mutate(layer = "sd_summary", y = row_y(disease),
           specification = factor(configuration, levels = configurations,
                                   labels = c("Exact", "2-percentage-point band", "5-percentage-point band", "Weighted")),
           glyph_y = y + c(.21, .07, -.07, -.21)[match(configuration, configurations)])
  stopifnot(nrow(d_sd) == 20, nrow(d_points) == nrow(exact), all(is.finite(d_sd$gap_sd_pp)))

  # E: country-cluster bootstrap, pooling resampled country-year rows exactly
  # as in the point estimand; resample all three age bands together.
  set.seed(20261006)
  e_summary <- bind_rows(lapply(levels, function(dis) {
    x <- filter(exact, disease == dis)
    values <- as.matrix(select(x, gap_waning_strict_1_4, gap_waning_strict_5_9, gap_waning_strict_10_14)) * 100
    stopifnot(all(is.finite(values)))
    sums <- rowsum(values, x$country)
    counts <- as.numeric(table(factor(x$country, levels = rownames(sums))))
    ncountry <- nrow(sums)
    multiplicity <- rmultinom(2000, ncountry, rep(1 / ncountry, ncountry))
    boot <- t(multiplicity) %*% sums / as.numeric(t(multiplicity) %*% counts)
    point <- colMeans(values)
    ci <- apply(boot, 2, quantile, c(.025, .975))
    gradient <- boot[, 3] - boot[, 1]
    bind_rows(tibble(layer = "age_cell", disease = dis, age_group = c("1–4", "5–9", "10–14"),
      x = 1:3, estimate_pp = point, lower_pp = ci[1, ], upper_pp = ci[2, ]),
      tibble(layer = "age_gradient", disease = dis, estimate_pp = point[3] - point[1],
             lower_pp = quantile(gradient, .025), upper_pp = quantile(gradient, .975))) %>%
      mutate(y = row_y(dis), n = nrow(x), n_countries = ncountry, bootstrap_replicates = 2000)
  }))

  # F: both coordinates are centred inside the identical disease/year/current-
  # coverage strata. This is the estimand behind the published within slopes.
  f_points <- exact %>% filter(is.finite(past_10y_deficit), is.finite(gap_waning_strict)) %>%
    group_by(disease, stratum) %>% filter(n_distinct(country) >= 2) %>%
    mutate(deficit_residual_pp = past_10y_deficit - mean(past_10y_deficit),
           gap_residual_pp = 100 * (gap_waning_strict - mean(gap_waning_strict))) %>% ungroup() %>%
    transmute(layer = "within_stratum_observation", disease, country, year, stratum,
              deficit_pp = past_10y_deficit, gap_pp = gap_waning_strict * 100,
              deficit_residual_pp, gap_residual_pp)
  f_slopes <- primary %>% transmute(layer = "within_stratum_slope", disease,
    estimate_pp = gap_per_10pp_past_deficit_pp, lower_pp = gap_per_10pp_past_deficit_pp_low,
    upper_pp = gap_per_10pp_past_deficit_pp_high, bootstrap_replicates)
  f_recalculated <- f_points %>% group_by(disease) %>% summarise(
    beta10 = 10 * sum(deficit_residual_pp * gap_residual_pp) / sum(deficit_residual_pp^2), .groups = "drop")
  stopifnot(max(abs(left_join(f_slopes, f_recalculated, by = "disease") %>% transmute(e = estimate_pp - beta10) %>% pull(e))) < 1e-8)
  f_lines <- bind_rows(lapply(levels, function(dis) {
    d <- filter(f_points, disease == dis)
    s <- filter(f_slopes, disease == dis)
    xx <- seq(min(d$deficit_residual_pp), max(d$deficit_residual_pp), length.out = 120)
    tibble(layer = "slope_band", disease = dis, deficit_residual_pp = xx,
      gap_residual_pp = xx * s$estimate_pp / 10,
      lower_pp = pmin(xx * s$lower_pp / 10, xx * s$upper_pp / 10),
      upper_pp = pmax(xx * s$lower_pp / 10, xx * s$upper_pp / 10))
  }))
  f_bins <- f_points %>% mutate(x = floor(deficit_residual_pp / 2) * 2 + 1,
                                y = floor(gap_residual_pp / 2) * 2 + 1) %>%
    count(disease, x, y, name = "n") %>% mutate(layer = "density_bin")
  stopifnot(sum(f_bins$n) == nrow(f_points))

  list(a = bind_rows(a_heat, a_coverage, a_cohorts, a_bands),
       b = bind_rows(b_coverage, b_deficit),
       c = bind_rows(c_profile, c_groups, c_total,
         c_grid %>% select(layer, country, age, birth_year, cohort_size, tau,
                            partial_series_coefficient, full_series_coefficient, gap_pp)),
       d = bind_rows(d_points, d_density, d_quartiles, d_sd), e = e_summary,
       f = bind_rows(f_points, f_slopes, f_lines, f_bins))
}

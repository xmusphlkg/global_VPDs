# Paired model-history contrasts. All aggregation precedes posterior summaries.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})

prepare_figure4_scenarios <- function(out = "output/revised") {
  read <- function(path, required) {
    types <- if (endsWith(path, "scenario_contrasts.csv")) cols(recovery_year = col_double()) else NULL
    z <- read_csv(file.path(out, path), col_types = types, show_col_types = FALSE, progress = FALSE)
    stopifnot(nrow(problems(z)) == 0)
    stopifnot(nrow(z) > 0, all(required %in% names(z)))
    z
  }
  qs <- function(z, groups, value) z %>% group_by(across(all_of(groups))) %>%
    summarise(median = median(.data[[value]]), lo = quantile(.data[[value]], .025),
      hi = quantile(.data[[value]], .975), lo50 = quantile(.data[[value]], .25),
      hi50 = quantile(.data[[value]], .75), n_draws = n(), .groups = "drop")
  annual_qs <- function(z, groups, value) qs(z, groups, value) %>%
    left_join(z %>% group_by(across(all_of(groups))) %>%
      summarise(mean = mean(.data[[value]]), .groups = "drop"), by = groups)
  sc <- read("bayesian_transmission/scenario_contrasts.csv",
    c("draw", "disease", "country", "history", "duration", "total_cases"))
  pp <- read("bayesian_transmission/posterior_predictive.csv",
    c("draw", "disease", "country", "year", "history", "duration", "expected_reports")) %>%
    rename(predicted_cases = expected_reports)
  holdout <- read("bayesian_transmission/time_holdout_posterior_predictive.csv",
    c("country", "disease", "year", "observed_cases", "nb_predictive_median",
      "nb_predictive_q025", "nb_predictive_q975", "nb_predictive_covered")) %>%
    transmute(country, disease, year, cases = observed_cases,
      pred_median = nb_predictive_median, pred_lo = nb_predictive_q025,
      pred_hi = nb_predictive_q975, covered_95 = as.logical(nb_predictive_covered))
  selection <- read("transmission_country_selection/selection_audit.csv", c("country", "selection_order")) %>%
    filter(!is.na(selection_order)) %>% arrange(selection_order)
  nc <- nrow(selection); nd <- n_distinct(sc$draw)
  key <- c("draw", "disease", "country", "history", "duration")
  stopifnot(nc == 10, setequal(sc$country, selection$country),
    setequal(sc$disease, c("Measles", "Rubella", "Pertussis", "Diphtheria", "Mumps")),
    nrow(sc) == nc * nd * 5 * 5, !anyDuplicated(sc[key]),
    !anyDuplicated(pp[c(key, "year")]), setequal(pp$year, 2014:2023),
    nrow(pp) == 10 * nrow(sc), all(is.finite(sc$total_cases)), all(sc$total_cases > 0),
    all(is.finite(pp$predicted_cases)), all(pp$predicted_cases > 0))
  # Saved scenario totals must equal the ten annual expected counts.
  totals <- pp %>% group_by(across(all_of(key))) %>% summarise(annual_sum = sum(predicted_cases), .groups = "drop") %>%
    inner_join(sc, by = key)
  stopifnot(nrow(totals) == nrow(sc), max(abs(totals$annual_sum / totals$total_cases - 1)) < 1e-10)
  duration <- sc %>% filter(history == "observed") %>%
    select(all_of(key), total_cases) %>% pivot_wider(names_from = duration, values_from = total_cases) %>%
    mutate(delta_cases = waning - unwaning, percent = 100 * (waning / unwaning - 1))
  stopifnot(nrow(duration) == nc * nd * 5, all(is.finite(duration$percent)))
  pooled <- sc %>% group_by(draw, disease, history, duration) %>%
    summarise(total = sum(total_cases), n_countries = n(), .groups = "drop")
  stopifnot(all(pooled$n_countries == nc))
  # a: annual expected reports on a common equal-country scale. Keep the
  # original ten-year summaries and interaction layers below because panel B
  # and the supplementary interaction table use those paired totals.
  ad <- pooled %>% filter(history != "recovery_5y") %>%
    mutate(cases_per_country = total / nc)
  annual_ad <- pp %>% filter(history != "recovery_5y") %>%
    group_by(draw, disease, history, duration, year) %>%
    summarise(cases_per_country = mean(predicted_cases), .groups = "drop")
  annual_a <- bind_rows(mutate(annual_ad, layer = "annual_draw"),
    annual_qs(annual_ad, c("disease", "history", "duration", "year"), "cases_per_country") %>%
      mutate(layer = "annual_summary"))
  # Country detail is summarized from the same draws without pooling countries.
  # The main panel displays posterior means; central intervals remain available
  # in the country supplements and the source data.
  country_annual <- annual_qs(filter(pp, history != "recovery_5y"),
    c("country", "disease", "history", "duration", "year"), "predicted_cases") %>%
    left_join(select(selection, country, country_name, selection_order), by = "country") %>%
    arrange(selection_order, disease, history, duration, year)
  a <- bind_rows(mutate(ad, layer = "draw"),
    qs(ad, c("disease", "history", "duration"), "cases_per_country") %>% mutate(layer = "summary"))
  interaction <- ad %>% select(draw, disease, history, duration, cases_per_country) %>%
    pivot_wider(names_from = c(history, duration), values_from = cases_per_country) %>%
    mutate(interaction_cases = (observed_waning - observed_unwaning) -
      (sustained_high_waning - sustained_high_unwaning))
  a <- bind_rows(a, mutate(interaction, layer = "interaction_draw"),
    qs(interaction, "disease", "interaction_cases") %>% mutate(layer = "interaction_summary"),
    annual_a)
  annual_summary <- filter(annual_a, layer == "annual_summary")
  stopifnot(nrow(annual_ad) == nd * 5 * 4 * 10,
    all(is.finite(annual_ad$cases_per_country)),
    all(is.finite(annual_summary$mean)), all(annual_summary$mean > 0),
    all(is.finite(annual_summary$median)),
    all(annual_summary$lo <= annual_summary$median),
    all(annual_summary$median <= annual_summary$hi))
  stopifnot(nrow(country_annual) == nc * 5 * 4 * 10,
    !anyDuplicated(country_annual[c("country", "disease", "history", "duration", "year")]),
    all(country_annual$n_draws == nd), all(is.finite(country_annual$mean)),
    all(country_annual$mean > 0), all(country_annual$lo > 0),
    all(country_annual$lo <= country_annual$median & country_annual$median <= country_annual$hi),
    !anyNA(country_annual$country_name), !anyNA(country_annual$selection_order))
  # b: paired duration contrast under observed history. Keep the original
  # report-count estimand; the renderer applies the signed logarithmic display.
  bd <- duration %>% group_by(draw, disease) %>%
    summarise(delta_cases = mean(delta_cases), .groups = "drop")
  b <- bind_rows(mutate(bd, layer = "draw"),
    qs(bd, "disease", "delta_cases") %>% mutate(layer = "summary"))
  stopifnot(nrow(bd) == nd * 5, all(is.finite(bd$delta_cases)))
  # c: ratios of country sums, never averages of country-specific ratios.
  baseline <- pooled %>% filter(history == "observed", duration == "waning") %>%
    select(draw, disease, baseline = total)
  cd <- pooled %>% filter(duration == "waning") %>% inner_join(baseline, by = c("draw", "disease")) %>%
    mutate(relative_percent = 100 * total / baseline,
      reduction_percent = 100 - relative_percent)
  c <- bind_rows(mutate(cd, layer = "draw"),
    qs(cd, c("disease", "history", "duration"), "relative_percent") %>%
      mutate(reduction_percent = 100 - median,
        reduction_lo = 100 - hi, reduction_hi = 100 - lo,
        reduction_lo50 = 100 - hi50, reduction_hi50 = 100 - lo50,
        layer = "summary"))
  # d: every country/disease estimate retains its paired 95% interval.
  d <- qs(duration, c("country", "disease"), "percent") %>%
    left_join(select(selection, country, country_name, selection_order), by = "country") %>%
    mutate(excludes_zero = lo > 0 | hi < 0)
  # e: cumulative ratios computed within paired draws at every year.
  ed <- pp %>% filter(duration == "waning", disease %in% c("Pertussis", "Diphtheria", "Mumps")) %>%
    group_by(draw, disease, history, year) %>% summarise(annual_cases = sum(predicted_cases), .groups = "drop") %>%
    arrange(draw, disease, history, year) %>% group_by(draw, disease, history) %>%
    mutate(cumulative_cases = cumsum(annual_cases)) %>% ungroup()
  eb <- ed %>% filter(history == "observed") %>% select(draw, disease, year, baseline = cumulative_cases)
  ed <- ed %>% inner_join(eb, by = c("draw", "disease", "year")) %>%
    mutate(relative_percent = 100 * cumulative_cases / baseline)
  e <- bind_rows(mutate(ed, layer = "draw"),
    qs(ed, c("disease", "history", "year"), "relative_percent") %>% mutate(layer = "summary"))
  endpoint <- filter(e, layer == "summary", year == 2023) %>%
    inner_join(filter(c, layer == "summary"), by = c("disease", "history"), suffix = c("_e", "_c"))
  stopifnot(nrow(endpoint) == 9, max(abs(endpoint$median_e - endpoint$median_c)) < 1e-10,
    max(abs(endpoint$lo_e - endpoint$lo_c)) < 1e-10, max(abs(endpoint$hi_e - endpoint$hi_c)) < 1e-10)
  # The matched-present high-history scenario must reduce the ten-year
  # expected total in every disease at the posterior median. Annual paths can
  # cross because lower early incidence leaves more susceptibles for later
  # transmission, so this check is deliberately applied to the cumulative
  # endpoint rather than every calendar year.
  high_endpoint <- filter(c, layer == "summary", history == "sustained_high")
  stopifnot(nrow(high_endpoint) == 5, all(high_endpoint$median < 100))
  # f/g: predictive (negative-binomial observation) intervals, not mean intervals.
  stopifnot(nrow(holdout) == 179L, setequal(holdout$year, 2020:2023),
    !anyDuplicated(holdout[c("country", "disease", "year")]),
    all(is.finite(holdout$cases)), all(holdout$cases >= 0),
    all(holdout$pred_lo <= holdout$pred_median & holdout$pred_median <= holdout$pred_hi),
    all(holdout$covered_95 == (holdout$cases >= holdout$pred_lo & holdout$cases <= holdout$pred_hi)))
  coverage <- holdout %>% group_by(disease) %>% summarise(n_observed = n(),
    n_covered = sum(covered_95), coverage_percent = 100 * mean(covered_95), .groups = "drop")
  f <- mutate(holdout, layer = "observation")
  g <- bind_rows(mutate(coverage, layer = "coverage"),
    tibble(layer = "overall", n_observed = nrow(holdout), n_covered = sum(holdout$covered_95),
      coverage_percent = 100 * mean(holdout$covered_95)))
  for (z in list(a, c, d, e)) {
    ss <- if ("layer" %in% names(z)) filter(z, layer == "summary") else z
    stopifnot(all(is.finite(ss$median)), all(ss$lo <= ss$median & ss$median <= ss$hi), all(ss$n_draws == nd))
  }
  bs <- filter(b, layer == "summary")
  stopifnot(nrow(bs) == 5, all(is.finite(bs$median)),
    all(bs$lo <= bs$median & bs$median <= bs$hi), all(bs$n_draws == nd))
  list(a = a, b = b, c = c, d = d, e = e, f = f, g = g,
    country_annual = country_annual)
}

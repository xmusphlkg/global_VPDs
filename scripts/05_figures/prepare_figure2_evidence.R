# Prepare empirical evidence from authoritative fits; never refit case models.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(glmmTMB)})
prepare_figure2_evidence <- function(out = "output/revised") {
  read_result <- function(path, required) {
    x <- read_csv(file.path(out, path), show_col_types = FALSE, progress = FALSE)
    stopifnot(nrow(x) > 0, all(required %in% names(x)))
    x
  }
  coeff <- read_result("disease_model_coefficients.csv", c("disease", "model", "rate_ratio", "rr_low", "rr_high", "n"))
  cmp <- read_result("disease_model_comparisons.csv", c("comparison", "delta_aic", "aic_base", "aic_extended", "n"))
  stopifnot(all(abs(cmp$delta_aic - (cmp$aic_extended - cmp$aic_base)) < 1e-8))
  a <- coeff %>% filter(model == "M3_Cohort_time") %>%
    left_join(cmp %>% filter(comparison == "M3 vs M1") %>%
      transmute(disease, comparison_n = n, delta_aic_m3_m1 = delta_aic,
        aic_gain = -delta_aic), by = "disease")
  stopifnot(nrow(a) == 5, !anyDuplicated(a$disease), all(a$n == a$comparison_n),
    all(is.finite(a$aic_gain)), all(a$rr_low > 0), all(a$rr_low <= a$rate_ratio), all(a$rr_high >= a$rate_ratio))
  sample <- read_result("analysis_data.csv", c("disease", "country")) %>%
    group_by(disease) %>% summarise(n_country = n_distinct(country), sample_n = n(), .groups = "drop")
  a <- left_join(a, sample, by = "disease")
  stopifnot(all(a$n == a$sample_n), all(a$n_country > 0),
    all(is.finite(a$p_value)), sum(!is.na(a$q_value)) == 4)
  b <- cmp %>% filter(comparison %in% c("M1 vs M0", "M2 vs M1", "M3 vs M1", "M3 vs M2")) %>%
    mutate(layer = "disease_comparison", aic_gain = -delta_aic)
  stopifnot(nrow(b) == 20)
  pooled <- read_result("pooled_model_comparisons.csv",
    c("analysis", "comparison", "delta_aic", "aic_base", "aic_extended", "lrt_p_value", "df", "n")) %>%
    filter(comparison == "M4 interaction vs M3") %>%
    mutate(layer = "pooled_interaction", aic_gain = -delta_aic)
  stopifnot(nrow(pooled) == 2, !anyDuplicated(pooled$analysis),
    setequal(pooled$analysis, c("four_disease_primary", "five_disease_proxy")),
    all(is.finite(pooled$aic_gain)), all(pooled$lrt_p_value > 0 & pooled$lrt_p_value <= 1),
    all(abs(pooled$aic_gain - (pooled$aic_base - pooled$aic_extended)) < 1e-8))
  b <- bind_rows(b, pooled)

  # Standardize over each disease's observed country-year covariate distribution.
  # The country random intercept is fixed at zero; this is not integration over
  # the random-effect distribution and not a causal intervention estimate.
  # Use the full fixed-effect covariance for pointwise log-delta-method intervals.
  curves <- list()
  for (dis in c("Measles", "Rubella", "Pertussis", "Diphtheria", "Mumps")) {
    models <- readRDS(file.path(out, "models", paste0(tolower(dis), ".rds")))
    fit <- models$M3
    beta <- fixef(fit)$cond
    gap_term <- if ("gap_10" %in% names(beta)) "gap_10" else "gap0_10"
    V <- as.matrix(vcov(fit)$cond)
    X <- model.matrix(fit, component = "cond")
    ref <- filter(a, disease == dis)
    stopifnot(fit$fit$convergence == 0, isTRUE(fit$sdr$pdHess),
      nobs(fit) == ref$n, abs(exp(beta[[gap_term]]) - ref$rate_ratio) < 1e-6,
      abs(AIC(models$M1) - AIC(fit) - ref$aic_gain) < 1e-6,
      identical(colnames(X), names(beta)), identical(colnames(V), names(beta)))
    observed_gap <- X[, gap_term] * 10
    support <- quantile(observed_gap, c(.025, .975))
    central_support <- quantile(observed_gap, c(.1, .9))
    xs <- seq(support[1], support[2], length.out = 101)
    rows <- lapply(xs, function(gap) {
      design <- X; design[, gap_term] <- gap / 10
      rates <- exp(drop(design %*% beta)) * 1e5
      rate <- mean(rates)
      gradient <- colSums(design * rates) / sum(rates)
      se <- sqrt(drop(t(gradient) %*% V %*% gradient))
      tibble(disease = dis, layer = "standardized_curve", gap_pp = gap,
        rate_per_100k = rate, rate_low = exp(log(rate) - qnorm(.975) * se),
        rate_high = exp(log(rate) + qnorm(.975) * se), n = nrow(X),
        random_intercept = 0, support_low_pp = support[1], support_high_pp = support[2],
        central_support_low_pp = central_support[1], central_support_high_pp = central_support[2])
    })
    curve <- bind_rows(rows)
    stopifnot(all(is.finite(curve$rate_high)), all(curve$rate_low > 0),
      abs(diff(log(curve$rate_per_100k[c(1, nrow(curve))])) / diff(range(xs)) * 10 - beta[[gap_term]]) < 1e-8)
    curves[[dis]] <- bind_rows(curve, tibble(disease = dis, layer = "observed_gap",
      country = as.character(fit$frame$country), year = as.integer(as.character(fit$frame[["factor(year)"]])),
      gap_pp = observed_gap))
  }
  c <- bind_rows(curves)
  diag <- read_result("five_disease_diagnostics/coefficients.csv", c("specification", "term", "rr", "ok")) %>%
    filter(specification %in% c("L5_primary_add_prior", "W_within_between")) %>%
    mutate(layer = "decomposition", component = case_when(term == "gap_10" ~ "Overall",
      term == "gap_10_within" ~ "Within", term == "gap_10_between" ~ "Between"))
  loco <- read_result("five_disease_diagnostics/leave_one_country_out.csv", c("disease", "estimate", "omitted_country", "ok")) %>%
    mutate(layer = "leave_one_country_out")
  stopifnot(nrow(diag) == 15, all(diag$ok), nrow(loco) == sum(a$n_country), all(loco$ok),
    sum(loco$disease == "Measles") == 68, sum(loco$disease == "Pertussis") == 176)
  d <- bind_rows(diag, loco)
  e <- read_result("dtp_parameter_sensitivity.csv", c("tau", "partial_series_coefficient", "rr", "delta_aic_vs_historical", "converged")) %>%
    mutate(layer = "parameter_grid", delta_aic_m3_m1 = delta_aic_vs_historical,
      aic_gain = -delta_aic_vs_historical,
      primary = (disease == "Pertussis" & tau == 5 & partial_series_coefficient == .5) |
        (disease == "Diphtheria" & tau == 30 & partial_series_coefficient == .71))
  stopifnot(nrow(e) == 18, all(e$converged), sum(e$primary) == 2, all(e$rr > 0),
    all(is.finite(e$aic_gain)))
  primary <- filter(e, primary) %>% left_join(select(a, disease, rate_ratio, primary_gain = aic_gain), by = "disease")
  stopifnot(all(abs(log(primary$rr / primary$rate_ratio)) < 1e-4), all(abs(primary$aic_gain - primary$primary_gain) < 1e-6))
  hold <- read_result("time_holdout_by_disease.csv", c("scope", "cutoff", "disease", "model", "RMSE", "n_test"))
  pooled_hold <- read_result("time_holdout_metrics.csv", c("scope", "cutoff", "model", "RMSE", "n_test"))
  hold <- bind_rows(mutate(hold, layer = "disease_holdout"), mutate(pooled_hold, disease = "All", layer = "pooled_holdout")) %>%
    filter(scope == "five_disease_proxy", model %in% c("M1", "M3")) %>%
    select(layer, scope, disease, cutoff, model, RMSE, n_test, n_train) %>%
    pivot_wider(names_from = model, values_from = c(RMSE, n_test, n_train))
  stopifnot(nrow(hold) == 12, all(hold$n_test_M1 == hold$n_test_M3), all(hold$n_train_M1 == hold$n_train_M3),
    all(hold$RMSE_M1 > 0))
  f <- hold %>% mutate(rmse_change_pct = 100 * (RMSE_M3 / RMSE_M1 - 1),
      test_start = cutoff + 1L, test_end = if_else(cutoff == 2016, 2019L, 2023L),
      fit_scope = "disease-specific; year fixed effects", prediction_scope = "pooled common gap; linear time") %>%
    left_join(select(a, disease, delta_aic_m3_m1), by = "disease")
  stopifnot(all(is.finite(f$rmse_change_pct)), all(is.finite(filter(f, layer == "disease_holdout")$delta_aic_m3_m1)))
  list(a = a, b = b, c = c, d = d, e = e, f = f)
}

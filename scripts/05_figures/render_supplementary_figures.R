#!/usr/bin/env Rscript

# Revised supplementary figures S1-S8.
# S1-S3 are six-panel narratives; S4-S8 give country-level scenario trajectories.
# Figure legends are maintained in manuscripts/SUPPLEMENTARY_MATERIAL.md.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

source("scripts/05_figures/publication_figure_style.R")
source("scripts/05_figures/helpers/figure1_layout_helpers.R")

project_root <- normalizePath(".", mustWork = TRUE)
result_root <- file.path(project_root, "output", "revised")
figure_dir <- file.path(result_root, "figures", "supplementary")
source_dir <- file.path(result_root, "source_data", "figure_redesign")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

disease_levels <- pub_disease_levels
disease_cols <- pub_disease_colours

theme_panel <- function(grid = "y") {
  theme_publication(
    base_size = 5.7, grid = grid, show_titles = TRUE,
    plot_margin = margin(2.5, 3.5, 2.5, 3.5)
  ) + theme_panel_letter() +
    theme(
      axis.text.x = element_text(size = 4.8),
      axis.text.y = element_text(size = 4.8),
      axis.title = element_text(size = 5.7, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 4.6),
      legend.title = element_text(size = 4.9, face = "bold"),
      legend.key.height = unit(2.1, "mm"),
      legend.key.width = unit(3.5, "mm"),
      panel.spacing.x = unit(2.5, "mm")
    )
}

read_required <- function(path, columns = character()) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  dat <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  missing <- setdiff(columns, names(dat))
  if (length(missing)) {
    stop("Input ", path, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  dat
}

num <- function(x) suppressWarnings(as.numeric(x))
fmt_disease <- function(x) factor(as.character(x), levels = disease_levels)

write_source_long <- function(stem, panels) {
  stopifnot(length(panels) == 6, all(vapply(panels, nrow, integer(1)) > 0))
  # Preserve heterogeneous panel tables in one auditable long-form source file.
  out <- purrr::imap_dfr(panels, function(dat, panel_name) {
    if (is.null(dat) || !nrow(dat)) return(tibble())
    dat <- as_tibble(dat)
    dat$row_id <- seq_len(nrow(dat))
    dat %>%
      pivot_longer(-row_id, names_to = "variable", values_to = "value",
                   values_transform = as.character) %>%
      mutate(panel = panel_name, .before = 1)
  })
  readr::write_csv(out, file.path(source_dir, paste0(stem, "_source_data.csv")))
  invisible(out)
}

s2_only <- "--s2-only" %in% commandArgs(trailingOnly = TRUE)
s3_only <- "--s3-only" %in% commandArgs(trailingOnly = TRUE)
s4_s8_only <- "--s4-s8-only" %in% commandArgs(trailingOnly = TRUE)
if (s4_s8_only) {
  source("scripts/05_figures/render_figure4_country_supplement.R")
  render_figure4_country_supplement(result_root)
  quit(save = "no", status = 0)
}
save_figure <- function(plot, stem, height_mm = 142, panels) {
  for (i in seq_along(panels)) {
    write_csv(panels[[i]]$data, file.path(source_dir, paste0(sub("_.*", "", stem), "_", letters[i], ".csv")))
  }
  save_publication_figure(
    plot, stem = stem, directory = figure_dir,
    width_mm = pub_full_width_mm, height_mm = height_mm, dpi = 600
  )
}

message("Loading supplementary figure inputs...")

# S1 uses paired childhood scenarios and their fixed-input sensitivities.
s1_sensitivity <- read_required(
  file.path(result_root, "bayesian_transmission", "scenario_sensitivity.csv"),
  c("high_history_reduction_pct", "recovery_reduction_pct", "duration_contrast_pct")
)
s1_child <- s1_sensitivity %>% filter(sensitivity == "high_history_cap") %>%
  group_by(disease, level) %>% summarise(
    high = median(high_history_reduction_pct), recovery = median(recovery_reduction_pct), .groups = "drop") %>%
  pivot_longer(c(high, recovery), names_to = "scenario", values_to = "reduction") %>%
  mutate(disease = fmt_disease(disease))
s1_paired <- s1_sensitivity %>%
  filter(sensitivity %in% c("adult_gap", "natural_immunity_duration_years", "covid_contact_scale_2020_2023")) %>%
  mutate(setting = paste0(recode(sensitivity, adult_gap = "Adult gap", natural_immunity_duration_years = "Natural immunity", covid_contact_scale_2020_2023 = "Contact scale"), ": ", level)) %>%
  group_by(disease, setting) %>% summarise(
    high = median(high_history_reduction_pct), recovery = median(recovery_reduction_pct), .groups = "drop") %>%
  pivot_longer(c(high, recovery), names_to = "scenario", values_to = "reduction") %>%
  mutate(disease = fmt_disease(disease))
s1_draws <- read_csv(file.path(result_root, "bayesian_transmission", "posterior_predictive.csv"), show_col_types = FALSE)
s1_totals <- read_csv(file.path(result_root, "bayesian_transmission", "scenario_contrasts.csv"), show_col_types = FALSE) %>%
  group_by(draw, disease, history, duration) %>% summarise(total = sum(total_cases), .groups="drop") %>%
  unite(scenario, history, duration, sep="__") %>% pivot_wider(names_from=scenario, values_from=total)
s1_pooled <- s1_totals %>% mutate(waning_vs_unwaning_pct=100*(observed__waning/observed__unwaning-1),
  waning_vs_flat_pct=100*(1-sustained_high__waning/observed__waning),
  unwaning_vs_flat_pct=100*(1-recovery_5y__waning/observed__waning)) %>%
  group_by(disease) %>% summarise(across(ends_with("_pct"), median), .groups="drop") %>%
  mutate(disease=fmt_disease(disease), period="2014_2023")
s1_pred <- read_csv(file.path(result_root, "bayesian_transmission", "time_holdout_metrics.csv"), show_col_types=FALSE) %>%
  select(disease, expected_interval_coverage_pct, nb_predictive_coverage_pct) %>%
  pivot_longer(-disease, names_to="model", values_to="coverage_95_pct") %>%
  mutate(disease=fmt_disease(disease), split="temporal_holdout")
s1_adult <- read_csv(file.path(result_root, "bayesian_transmission", "scenario_sensitivity.csv"), show_col_types=FALSE) %>%
  filter(sensitivity=="adult_gap") %>% group_by(draw,disease,level) %>% summarise(total=sum(total_expected_reports),.groups="drop") %>%
  group_by(draw,disease) %>% mutate(relative_impact_pct=100*(total/total[level==0.3]-1)) %>% ungroup() %>%
  group_by(disease,level) %>% summarise(relative_impact_pct=median(relative_impact_pct),.groups="drop") %>%
  mutate(disease=fmt_disease(disease),adult_baseline_susceptibility=level,period="2014_2023")
s1_heat <- s1_draws %>% filter(history=="observed") %>%
  group_by(draw,disease,year,duration) %>% summarise(total=sum(expected_reports),.groups="drop") %>%
  pivot_wider(names_from=duration,values_from=total) %>% mutate(value=100*(waning/unwaning-1)) %>%
  group_by(disease,year) %>% summarise(value=median(value),.groups="drop") %>% mutate(disease=fmt_disease(disease))
stopifnot(nrow(s1_heat)==50, all(is.finite(s1_heat$value)))

# S2 inputs: reconstruction, matching, recovery and observation audit.
s2_panel <- read_required(
  file.path(project_root, "output/five_disease_revision/country_year_14cohort.csv"),
  c("disease", "observation_year", "scenario", "eligible_14_strict", "eligible_80pct_strict")
) %>% filter(scenario == "component_strict") %>%
  mutate(disease = fmt_disease(disease), year = num(observation_year))
s2_status <- read_required(
  file.path(result_root, "case_data_quality_summary.csv"),
  c("disease", "case_status", "size", "observed_rows_all_statuses", "source_rows_absent")
) %>% mutate(disease = fmt_disease(disease), size = num(size))
s2_strata <- read_required(
  file.path(result_root, "history_matching", "matching_strata.csv"),
  c("disease", "configuration", "gap_sd_pp", "gap_range_pp")
) %>% mutate(disease = fmt_disease(disease), gap_sd_pp = num(gap_sd_pp), gap_range_pp = num(gap_range_pp))
s2_match <- read_required(
  file.path(result_root, "history_matching", "matching_summary.csv"),
  c("disease", "configuration", "gap_sd_pp", "constant_protection_gap_sd_pp", "age_standardized_gap_sd_pp")
) %>% mutate(disease = fmt_disease(disease))
s2_recovery <- read_required(
  file.path(result_root, "history_matching", "recovery_summary.csv"),
  c("disease", "definition", "n_events")
) %>% select(disease, definition, n_events) %>% mutate(disease = fmt_disease(disease))
s2_audit <- read_required(
  file.path(result_root, "case_quality_observation_audit.csv"),
  c("disease", "n_joined", "n_current_missing", "n_lag_missing", "n_covariate_missing", "n_model_complete")
) %>% mutate(disease = fmt_disease(disease))

# S3 inputs: Bayesian diagnostics, predictive checks, and scenario posterior.
s3_diag <- read_required(
  file.path(result_root, "bayesian_transmission", "diagnostics.csv"),
  c("disease", "parameter", "rhat_rank_normalized", "ess_bulk", "ess_tail")
) %>% mutate(disease = fmt_disease(disease))
s3_ppc <- read_required(
  file.path(result_root, "bayesian_transmission", "posterior_predictive_checks.csv"),
  c("disease", "covered_95")
) %>% mutate(disease = fmt_disease(disease), covered_95 = as.logical(covered_95))
s3_holdout <- read_required(
  file.path(result_root, "bayesian_transmission", "time_holdout_posterior_predictive.csv"),
  c("disease", "nb_predictive_covered")
) %>% mutate(disease = fmt_disease(disease), nb_predictive_covered = as.logical(nb_predictive_covered))
s3_audit <- read_required(
  file.path(result_root, "bayesian_transmission", "input_audit.csv"),
  c("country", "disease", "n_rows", "n_observed", "n_positive", "mumps_proxy")
) %>% mutate(disease = fmt_disease(disease))
s3_contact <- read_required(
  file.path(result_root, "bayesian_transmission", "contact_matrix_audit.csv"),
  c("country", "corrected_spectral_radius", "corrected_max_reciprocity_error")
) %>% mutate(corrected_spectral_radius = num(corrected_spectral_radius), corrected_max_reciprocity_error = num(corrected_max_reciprocity_error))

if (!s2_only && !s3_only) {
message("Building Supplementary Figure S1...")

s1_a <- ggplot(s1_child, aes(level, reduction, colour = scenario, group = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = pub_colours$reference, linewidth = .25) +
  geom_line(linewidth = .55) + geom_point(size = .8) +
  scale_y_continuous(limits = function(x) range(c(0, x)), breaks = breaks_pretty(n = 4),
    expand = expansion(mult = c(0, .05))) +
  facet_wrap(~disease, nrow = 1, labeller = labeller(disease = pub_disease_short_labels)) +
  scale_colour_manual(values = c(high = pub_colours$red, recovery = pub_colours$blue),
                      labels = c(high = "Child high protection", recovery = "5-year recovery")) +
  scale_x_continuous(breaks = c(0, .1, .2), labels = percent_format(accuracy = 1)) +
  labs(title = "A", x = "10–14-year gap cap (%)", y = "Reduction in expected reports (%)", colour = NULL) + theme_panel("y")

s1_b <- ggplot(s1_paired, aes(reduction, setting, colour = scenario, shape = scenario)) +
  geom_point(position = position_dodge(width = .5), size = .8) +
  facet_wrap(~disease, nrow = 1, labeller = labeller(disease = pub_disease_short_labels)) +
  scale_colour_manual(values = c(high = pub_colours$red, recovery = pub_colours$blue),
                      labels = c(high = "Child high protection", recovery = "5-year recovery")) +
  scale_shape_manual(values = c(high = 16, recovery = 17),
                     labels = c(high = "Child high protection", recovery = "5-year recovery")) +
  labs(title = "B", x = "Reduction in expected reports (%)", y = NULL, colour = NULL, shape = NULL) + theme_panel("x")

s1_c <- ggplot(s1_heat, aes(factor(year), disease, fill = value)) + geom_tile(colour = pub_colours$light_grey, linewidth = 0.15) +
  geom_text(aes(label = sprintf("%.1f", value)), size = 1.8) +
  scale_fill_gradient2(low = pub_colours$teal, mid = "white", high = pub_colours$red, midpoint = 0, na.value = pub_colours$light_grey) +
  guides(fill = figure1_colourbar("Change under waning versus no waning (%)", width = 24)) +
  labs(title = "C", x = "Calendar year", y = NULL, fill = "Change under waning versus no waning (%)") +
  theme_panel("none") + figure1_bottom_legend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

s1_d <- ggplot(s1_pooled %>% filter(period == "2014_2023") %>%
    pivot_longer(c(waning_vs_unwaning_pct, waning_vs_flat_pct, unwaning_vs_flat_pct), names_to = "comparison", values_to = "percent") %>%
    mutate(comparison = recode(comparison, waning_vs_unwaning_pct = "Waning vs unwaning", waning_vs_flat_pct = "10–14-year reference reduction", unwaning_vs_flat_pct = "Recovery reduction")),
  aes(disease, percent, fill = comparison)
) + geom_hline(yintercept = 0, linetype = "dashed", colour = pub_colours$reference, linewidth = 0.25) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  scale_fill_manual(values = c("Waning vs unwaning" = pub_colours$orange, "10–14-year reference reduction" = pub_colours$red, "Recovery reduction" = pub_colours$blue)) +
  labs(title = "D", x = NULL, y = "Contrast in expected reports (%)", fill = NULL) +
  theme_panel("y") + theme(axis.text.x = element_text(angle = 35, hjust = 1))

s1_e <- ggplot(s1_pred %>% filter(split == "temporal_holdout"), aes(disease, coverage_95_pct, colour = model, shape = model)) +
  geom_hline(yintercept = 95, linetype = "dashed", linewidth = 0.25, colour = pub_colours$reference) +
  geom_point(position = position_dodge(width = 0.5), size = 1.4) +
  scale_colour_manual(values=c(expected_interval_coverage_pct=pub_colours$blue, nb_predictive_coverage_pct=pub_colours$orange),
    labels=c(expected_interval_coverage_pct="Interval for expected reports", nb_predictive_coverage_pct="Negative-binomial observation-predictive interval")) +
  scale_shape_discrete(labels=c(expected_interval_coverage_pct="Interval for expected reports", nb_predictive_coverage_pct="Negative-binomial observation-predictive interval")) +
  scale_y_continuous(labels = label_number(suffix = "%"), limits = c(0, 100),
    breaks = seq(0, 100, 25), expand = expansion(mult = c(0, .05))) +
  labs(title = "E", x = NULL, y = "Coverage of 95% predictive intervals (%)", colour = NULL, shape = NULL) +
  theme_panel("y") + theme(axis.text.x = element_text(angle = 35, hjust = 1))

s1_f <- ggplot(s1_adult %>% filter(period == "2014_2023"), aes(adult_baseline_susceptibility, relative_impact_pct, colour = disease, group = disease)) +
  geom_hline(yintercept = 0, colour = pub_colours$reference, linewidth = 0.2) + geom_line(linewidth = 0.55) + geom_point(size = 1.0) +
  scale_colour_manual(values = disease_cols, drop = FALSE) + scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "F", x = "Adult protection gap (%)", y = "Change in expected reports\nrelative to 30% gap (%)", colour = "Disease") +
  theme_panel("y")

fig_s1 <- (s1_a + s1_b) / (s1_c + s1_d) / (s1_e + s1_f)
save_figure(fig_s1, "FigureS1_revised_paired_sensitivity", 165, list(s1_a,s1_b,s1_c,s1_d,s1_e,s1_f))
write_source_long("FigureS1_revised", list(child_cap_sensitivity = s1_child, paired_sensitivity = s1_paired, annual_heatmap = s1_heat,
  pooled_contrasts = s1_pooled, predictive_diagnostics = s1_pred, adult_susceptibility = s1_adult))

}

if (!s3_only) {
message("Building Supplementary Figure S2...")

s2_a <- s2_panel %>% group_by(disease, year) %>% summarise(
  complete = mean(as.logical(eligible_14_strict), na.rm = TRUE) * 100,
  coverage = mean(as.logical(eligible_80pct_strict), na.rm = TRUE) * 100, .groups = "drop"
) %>% pivot_longer(c(complete, coverage), names_to = "metric", values_to = "percent") %>%
  ggplot(aes(year, percent, colour = metric, group = metric)) + geom_line(linewidth = 0.55) + geom_point(size = 0.6) +
  facet_wrap(~disease, nrow = 1, labeller = labeller(disease = pub_disease_short_labels)) + scale_colour_manual(values = c(complete = pub_colours$green, coverage = pub_colours$blue), labels = c(complete = "All 14 ages", coverage = "At least 80%")) +
  scale_x_continuous(breaks = c(2014, 2023), labels = c("14", "23")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
    labels = label_number(suffix = "%"), expand = expansion(mult = c(0, .05))) + labs(title = "A", x = "Calendar year (2014–2023)", y = "Eligible rows (%)", colour = NULL) + theme_panel("y")

s2_b <- ggplot(s2_status %>% group_by(disease, case_status) %>% summarise(size = sum(size, na.rm = TRUE), .groups = "drop"), aes(disease, size, fill = case_status)) +
  geom_col(position = "fill", width = 0.72) + scale_y_continuous(breaks = seq(0, 1, .25),
    labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = c(reported_positive = pub_colours$red, reported_zero = pub_colours$blue, unreported = pub_colours$reference, not_applicable = pub_colours$light_grey), labels = c(reported_positive = "Positive", reported_zero = "Zero", unreported = "Unreported", not_applicable = "Not applicable"), drop = FALSE) +
  labs(title = "B", x = NULL, y = "Share of status rows (%)", fill = "Case status") + theme_panel("none") + theme(axis.text.x = element_text(angle = 35, hjust = 1))

s2_c <- ggplot(s2_strata %>% filter(configuration == "exact_equal_country_year", is.finite(gap_sd_pp)), aes(gap_sd_pp)) +
  geom_histogram(bins = 24, fill = pub_colours$blue, colour = "white", linewidth = 0.15) + facet_wrap(~disease, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), breaks = breaks_pretty(n = 4),
    expand = expansion(mult = c(0, .05))) +
  labs(title = "C", x = "Gap SD (percentage points)", y = "Number of strata") + theme_panel("y")

s2_d <- s2_match %>% filter(configuration == "exact_equal_country_year") %>%
  pivot_longer(c(gap_sd_pp, constant_protection_gap_sd_pp, age_standardized_gap_sd_pp), names_to = "estimand", values_to = "sd_pp") %>%
  mutate(estimand = recode(estimand, gap_sd_pp = "Waning gap", constant_protection_gap_sd_pp = "Constant protection", age_standardized_gap_sd_pp = "Age standardized")) %>%
  ggplot(aes(disease, sd_pp, fill = estimand)) + geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  scale_fill_manual(values = c("Waning gap" = pub_colours$orange, "Constant protection" = pub_colours$reference, "Age standardized" = pub_colours$teal)) +
  scale_y_continuous(limits = c(0, NA), breaks = breaks_pretty(n = 4),
    expand = expansion(mult = c(0, .05))) +
  labs(title = "D", x = NULL, y = "Gap SD (percentage points)", fill = NULL) + theme_panel("y") + theme(axis.text.x = element_text(angle = 35, hjust = 1))

# Former Figure 1g: retain both recovery definitions and event counts.
s2_recovery <- s2_recovery %>% mutate(
  definition_label = factor(definition,
    levels = c("prior_3y_mean_le85", "each_prior_3y_le85"), labels = c("Mean", "Strict"))
  )
stopifnot(nrow(s2_recovery) == 6, all(is.finite(s2_recovery$n_events)))
s2_e <- ggplot(s2_recovery, aes(disease, n_events, fill = definition_label)) +
  geom_col(position = position_dodge(width = .72), width = .62) +
  scale_fill_manual(values = c(Mean = pub_colours$orange, Strict = pub_colours$blue),
    name = "Recovery definition") +
  scale_y_continuous(breaks = pretty(s2_recovery$n_events, n = 4),
    expand = expansion(mult = c(0, .12))) +
  labs(title = "E", x = NULL, y = "Eligible recovery events") + theme_panel("y") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

s2_f <- s2_audit %>% select(disease, n_joined, n_current_missing, n_lag_missing, n_covariate_missing, n_model_complete) %>%
  pivot_longer(-disease, names_to = "audit", values_to = "n") %>% mutate(audit = recode(audit, n_joined = "Joined", n_current_missing = "Current missing", n_lag_missing = "Lag missing", n_covariate_missing = "Covariate missing", n_model_complete = "Model complete")) %>%
  ggplot(aes(disease, n, fill = audit)) + geom_col(position = "dodge", width = 0.72) +
  scale_fill_manual(values = c(Joined = pub_colours$blue, `Current missing` = pub_colours$orange, `Lag missing` = pub_colours$red, `Covariate missing` = pub_colours$gold, `Model complete` = pub_colours$green)) +
  scale_y_continuous(limits = c(0, NA), breaks = breaks_pretty(n = 4),
    expand = expansion(mult = c(0, .05))) +
  labs(title = "F", x = NULL, y = "Country–year rows", fill = NULL) + theme_panel("y") + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

fig_s2 <- (s2_a + s2_b) / (s2_c + s2_d) / (s2_e + s2_f)
save_figure(fig_s2, "FigureS2_revised_reconstruction_robustness", 165, list(s2_a,s2_b,s2_c,s2_d,s2_e,s2_f))
write_source_long("FigureS2_revised", list(cohort_completeness = s2_panel, case_status = s2_status, matched_strata = s2_strata,
  sensitivity = s2_match, recovery_events = s2_recovery, observation_audit = s2_audit))

}

if (!s2_only) {
message("Building Supplementary Figure S3...")

s3_a <- ggplot(s3_diag, aes(disease, rhat_rank_normalized, colour = parameter)) +
  geom_point(position = position_dodge(width = 0.5), size = 1.3) + scale_colour_manual(values = c(log_R0 = pub_colours$blue, log_nb_size = pub_colours$orange)) +
  labs(title = "A", x = NULL, y = "R-hat", colour = "Parameter") + theme_panel("y") + theme(axis.text.x = element_text(angle = 35, hjust = 1))

s3_b <- s3_diag %>% pivot_longer(c(ess_bulk, ess_tail), names_to = "ess_type", values_to = "ess") %>%
  ggplot(aes(disease, ess, colour = parameter, shape = ess_type)) +
  geom_point(position = position_dodge(width = 0.5), size = 1.3) + scale_colour_manual(values = c(log_R0 = pub_colours$blue, log_nb_size = pub_colours$orange)) +
  labs(title = "B", x = NULL, y = "Effective sample size", colour = "Parameter", shape = "ESS type") + theme_panel("y") + theme(axis.text.x = element_text(angle = 35, hjust = 1))

s3_c <- s3_ppc %>% group_by(disease) %>% summarise(coverage = mean(covered_95, na.rm = TRUE) * 100, n = n(), .groups = "drop") %>%
  ggplot(aes(disease, coverage, fill = disease)) + geom_col(width = 0.68) + scale_fill_manual(values = pub_evidence_colours, guide = "none") + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
    labels = label_number(suffix = "%"), expand = expansion(mult = c(0, .05))) +
  labs(title = "C", x = NULL, y = "Posterior predictive coverage (%)") + theme_panel("y") + theme(axis.text.x = element_text(angle = 35, hjust = 1))

s3_d <- s3_holdout %>% group_by(disease) %>% summarise(coverage = mean(nb_predictive_covered, na.rm = TRUE) * 100, n = n(), .groups = "drop") %>%
  ggplot(aes(disease, coverage, fill = disease)) + geom_col(width = 0.68) + scale_fill_manual(values = pub_evidence_colours, guide = "none") + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
    labels = label_number(suffix = "%"), expand = expansion(mult = c(0, .05))) +
  labs(title = "D", x = NULL, y = "Observed within 95%\nnegative-binomial predictive interval (%)") + theme_panel("y") + theme(axis.text.x = element_text(angle = 35, hjust = 1))

s3_e <- s3_audit %>% group_by(disease) %>% summarise(positive = sum(n_positive), zero = sum(n_observed-n_positive), missing = sum(n_rows-n_observed), .groups = "drop") %>%
  pivot_longer(c(positive, zero, missing), names_to = "count_type", values_to = "n") %>%
  ggplot(aes(disease, n, fill = count_type)) + geom_col(width = 0.65) +
  scale_fill_manual(values = c(positive = pub_colours$red, zero = pub_colours$blue, missing = pub_colours$reference), labels = c(missing = "Unobserved", positive = "Positive", zero = "Zero")) +
  scale_x_discrete(labels = identity) +
  scale_y_continuous(limits = c(0, NA), breaks = breaks_pretty(n = 4),
    expand = expansion(mult = c(0, .05))) +
  labs(title = "E", x = NULL, y = "Country–year rows", fill = NULL) + theme_panel("y") + theme(axis.text.x = element_text(angle = 35, hjust = 1))

s3_f <- ggplot(s3_contact, aes(corrected_spectral_radius, corrected_max_reciprocity_error)) +
  geom_point(size = 1.4, colour = pub_colours$blue) +
  geom_text_repel(aes(label = country), seed = 42, size = 2.0,
    colour = pub_colours$muted_text, box.padding = .25, point.padding = .15,
    min.segment.length = 0, max.overlaps = Inf, segment.size = .2) +
  scale_y_log10(labels = label_scientific(), limits = c(3e-17, 3e-16)) +
  labs(title = "F", x = "Corrected contact spectral radius", y = "Maximum reciprocity error") + theme_panel("y")

fig_s3 <- (s3_a + s3_b) / (s3_c + s3_d) / (s3_e + s3_f)
save_figure(fig_s3, "FigureS3_revised_bayesian_input_audit", 165, list(s3_a,s3_b,s3_c,s3_d,s3_e,s3_f))
write_source_long("FigureS3_revised", list(rhat = s3_diag, ess = s3_diag, posterior_predictive = s3_ppc, time_holdout = s3_holdout,
  input_audit = s3_audit, contact_matrix = s3_contact))

}

if (!s2_only && !s3_only) {
  source("scripts/05_figures/render_figure4_country_supplement.R")
  render_figure4_country_supplement(result_root)
}

message("Supplementary figures written to: ", figure_dir)
message("Source tables written to: ", source_dir)

if (!s2_only && !s3_only) {
  source("scripts/05_figures/render_figure5_mechanisms.R")
  render_age_mechanism_figures(result_root, figures = "S9")
}

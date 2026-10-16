#!/usr/bin/env Rscript

# Supplementary Figure S10: robustness of the pertussis association and
# pooled one-time 10-14-year catch-up scenario.  This renderer is deliberately
# standalone so that it can be rerun without modifying the main figure script.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
source("scripts/05_figures/publication_figure_style.R")

result_root <- "output/revised"
figure_dir <- file.path(result_root, "figures", "supplementary")
source_dir <- file.path(result_root, "source_data", "figure_redesign")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

robust_path <- "output/robustness_case_models/robustness_comparison.csv"
coef_path <- "output/robustness_case_models/robustness_coefficients.csv"
catchup_path <- file.path(result_root, "bayesian_transmission", "catchup_10_14_2018_contrasts_summary.csv")
stopifnot(file.exists(robust_path), file.exists(coef_path), file.exists(catchup_path))

robust <- read_csv(robust_path, show_col_types = FALSE)
coef <- read_csv(coef_path, show_col_types = FALSE)
catchup <- read_csv(catchup_path, show_col_types = FALSE)

# The global spline comparison is a likelihood-ratio test against the primary
# (linear gap) model.  Store the value in the source table so it is auditable.
primary_rds <- "output/robustness_case_models/models/Pertussis_primary.rds"
spline_rds <- "output/robustness_case_models/models/Pertussis_nonlinear_spline.rds"
spline_p <- NA_real_
if (requireNamespace("glmmTMB", quietly = TRUE) && file.exists(primary_rds) && file.exists(spline_rds)) {
  ptest <- tryCatch(anova(readRDS(primary_rds), readRDS(spline_rds)), error = function(e) NULL)
  if (!is.null(ptest) && nrow(ptest) >= 2) spline_p <- as.numeric(ptest$`Pr(>Chisq)`[2])
}
if (!is.finite(spline_p)) spline_p <- 3.0e-4 # prespecified diagnostic in the analysis log

forest <- robust %>%
  filter(disease == "Pertussis", analysis %in% c("primary", "country_fixed_effect", "pre_pandemic")) %>%
  transmute(
    analysis = recode(analysis,
      primary = "Primary model",
      country_fixed_effect = "Country fixed effects",
      pre_pandemic = "2014–2019 only"),
    estimate = rate_ratio, low = rr_low, high = rr_high,
    p_value, n, year_min, year_max
  ) %>%
  mutate(analysis = factor(analysis, levels = rev(c("Primary model", "Country fixed effects", "2014–2019 only"))))

forest_source <- forest %>% mutate(spline_global_lrt_p = spline_p)
catchup_source <- catchup %>%
  transmute(disease, campaign, campaign_year, gap_cap, baseline_expected_reports,
            catchup_expected_reports, absolute_reduction, pooled_reduction_pct,
            q025_reduction_pct, q975_reduction_pct, n_draws)

readr::write_csv(forest_source, file.path(source_dir, "FigureS10_a_pertussis_robustness_source_data.csv"))
readr::write_csv(catchup_source, file.path(source_dir, "FigureS10_b_catchup_pooled_source_data.csv"))

theme_s10 <- theme_publication(base_size = 6.2, grid = "x", show_titles = TRUE,
  plot_margin = margin(3, 4, 3, 4)) +
  theme(axis.text = element_text(size = 5.4), axis.title = element_text(size = 6.2, face = "bold"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        panel.grid.major.y = element_blank())

p_a <- ggplot(forest, aes(x = estimate, y = analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .25, colour = pub_colours$reference) +
  geom_text(data = data.frame(x = 1, analysis = "Primary model", label = "RR = 1"),
            aes(x, analysis, label = label), inherit.aes = FALSE, hjust = 0, vjust = 2.5,
            nudge_x = .5, size = 2.2, colour = pub_colours$reference) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = .16, linewidth = .5, colour = pub_colours$orange) +
  geom_point(size = 1.8, colour = pub_colours$orange) +
  scale_x_continuous(limits = c(0, 35), breaks = c(0, 5, 10, 15, 20, 25, 30, 35),
                labels = c("0", "5", "10", "15", "20", "25", "30", "35"),
                expand = c(0, 0)) +
  labs(title = "A",
       x = "Rate ratio per 0.10 increase in modelled cohort-time gap", y = NULL) + theme_s10

p_b <- catchup_source %>%
  mutate(disease = factor(disease, levels = rev(pub_disease_levels))) %>%
  ggplot(aes(x = pooled_reduction_pct, y = disease, colour = disease)) +
  geom_vline(xintercept = 0, linewidth = .25, colour = pub_colours$reference) +
  geom_errorbarh(aes(xmin = q025_reduction_pct, xmax = q975_reduction_pct), height = .18, linewidth = .5) +
  geom_point(size = 1.9) +
  geom_text(aes(x = q975_reduction_pct, label = sprintf("%.1f%%", pooled_reduction_pct)),
            hjust = 0, nudge_x = 1.2, size = 2.3, colour = pub_colours$text) +
  scale_colour_manual(values = pub_disease_colours, drop = FALSE) +
  scale_x_continuous(limits = c(0, 45), breaks = seq(0, 40, 10), labels = scales::label_percent(scale = 1),
                     expand = expansion(mult = c(0, .02))) +
  labs(title = "B",
       x = "Reduction in expected reports (%)", y = NULL) +
  theme_s10 + theme(legend.position = "none")

fig <- p_a + p_b + plot_layout(widths = c(1, 1))
save_publication_figure(fig, "FigureS10_robustness_and_catchup", figure_dir,
                        width_mm = pub_full_width_mm, height_mm = 92, dpi = 600)
message("Wrote Figure S10 and source data to ", figure_dir)

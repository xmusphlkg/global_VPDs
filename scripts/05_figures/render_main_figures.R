#!/usr/bin/env Rscript
# Render the approved result narratives from existing analysis outputs.
# This script never fits a model or pools posterior draws across countries.
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr)
  library(tidyr); library(readr); library(scales)
})
source("scripts/05_figures/publication_figure_style.R")
out <- Sys.getenv("VPDS_FIGURE_OUT", "output/revised")
fig_dir <- file.path(out, "figures/main_nc")
src_dir <- file.path(out, "source_data/figure_redesign")
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
read_input <- function(path, required) {
  d <- read_csv(file.path(out, path), show_col_types = FALSE, progress = FALSE)
  stopifnot(nrow(d) > 0, all(required %in% names(d)))
  if ("disease" %in% names(d)) d$disease <- factor(d$disease, levels = pub_disease_levels)
  d
}
cc <- function(guide = "none") scale_colour_manual(values = pub_disease_colours, guide = guide)
fc <- function(guide = "none") scale_fill_manual(values = pub_disease_colours, guide = guide)
panel <- function(p, x = NULL, y = NULL, grid = "y") {
  p + theme_publication(grid = grid, show_titles = TRUE) + theme_panel_letter() +
    labs(subtitle = NULL, x = x, y = y) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.margin = margin(3, 4, 3, 4))
}
zero <- function() geom_hline(yintercept = 0, colour = pub_colours$reference, linewidth = .25, linetype = 2)
qsummary <- function(d, groups, value) {
  d %>% group_by(across(all_of(groups))) %>% summarise(
    median = median(.data[[value]]), lo = quantile(.data[[value]], .025),
    hi = quantile(.data[[value]], .975), n_draws = n(), .groups = "drop")
}
# Figure legends are maintained in manuscripts/manuscript_nature_communications.md.
source("scripts/05_figures/render_figure1_atlas.R")
source("scripts/05_figures/render_figure2_evidence.R")
source("scripts/05_figures/render_figure3_measles_resurgence.R")
source("scripts/05_figures/render_figure4_scenarios.R")
source("scripts/05_figures/render_figure5_mechanisms.R")
source("scripts/05_figures/render_figure6_serology.R")

# Each --figureN-only entry point leaves the other main figures untouched.
renderers <- list(
  "1" = list(render_figure1_atlas, "Figure 1 cohort-memory atlas: six composite panels"),
  "2" = list(render_figure2_evidence, "Figure 2 empirical evidence: six composite panels"),
  "3" = list(render_figure3_measles_resurgence, "Figure 3 EU/EEA measles resurgence: seven panels, with Supplementary S12"),
  "4" = list(render_figure4_scenarios, "Figure 4 paired history and waning scenarios: seven composite panels"),
  "5" = list(render_figure5_mechanisms, "Figure 5 age and transmission mechanisms with Supplementary S9: six panels"),
  "6" = list(render_figure6_serology, "Figure 6 serology triangulation: three panels"))
only <- sub("^--figure([1-6])-only$", "\\1",
  grep("^--figure[1-6]-only$", commandArgs(trailingOnly = TRUE), value = TRUE))
for (number in if (length(only)) only else names(renderers)) {
  renderers[[number]][[1]](out)
  message("Rendered ", renderers[[number]][[2]], "; 600-dpi PNG and vector PDF.")
}

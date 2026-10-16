# Independent antibody age structure and its boundary with reconstructed gaps.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(scales)
})
source("scripts/05_figures/publication_figure_style.R")
source("scripts/05_figures/helpers/figure1_layout_helpers.R")
source("scripts/05_figures/prepare_figure6_serology.R")

render_figure6_serology <- function(out = "output/revised") {
  data <- prepare_figure6_serology(out)
  disease_cols <- pub_evidence_colours
  charcoal <- "#4B4D4F"
  signed <- function(x) sprintf("%s%g", ifelse(x > 0, "+", ""), x)
  mark <- function(p, letter) p + labs(title = letter)
  t <- function(grid = "none") theme_publication(grid = grid, show_titles = TRUE) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(1, 2, 1, 2),
      axis.title = figure1_axis_title_style(),
      axis.title.x = figure1_axis_title_style(margin(t = 2)),
      axis.title.y = figure1_axis_title_style(margin(r = 2)),
      axis.text = element_text(size = 5.4),
      axis.ticks.length = unit(.7, "mm"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0,
        margin = margin(b = 3)),
      plot.title.position = "plot",
      strip.text = figure1_facet_title_style(),
      legend.text = element_text(size = 5.1),
      legend.title = figure1_axis_title_style(),
      legend.key.width = unit(4, "mm"),
      legend.key.height = unit(2.6, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "mm"),
      panel.spacing = unit(3, "mm")
    )

  # Marginal proportion intervals cannot be reused as intervals for baseline
  # differences without survey covariance estimates; these changes are descriptive.
  a <- data$age_changes %>% mutate(
    series_colour = if_else(endpoint == "NHANES seronegativity",
      unname(disease_cols[as.character(disease)]), charcoal)
  )
  endpoint_labels <- a %>% filter(age == age_oldest) %>% mutate(
    label = sprintf("%s %+.1f pp",
      if_else(endpoint == "NHANES seronegativity", "NHANES", "Programme-derived gap"), change_pp),
    label_y = if_else(endpoint == "NHANES seronegativity", 19, 16.5)
  )
  pa <- ggplot(a, aes(age, change_pp, group = endpoint)) +
    geom_hline(yintercept = 0, colour = pub_colours$reference,
      linewidth = .25, linetype = 2) +
    geom_line(aes(colour = series_colour, linetype = endpoint), linewidth = .6) +
    geom_point(aes(colour = series_colour, shape = endpoint), size = 1.4) +
    geom_text(data = endpoint_labels,
      aes(x = 6.25, y = label_y, label = label, colour = series_colour),
      inherit.aes = FALSE, hjust = 0, size = 2.05, show.legend = FALSE) +
    facet_wrap(vars(disease), nrow = 1, drop = TRUE) +
    scale_colour_identity() +
    scale_linetype_manual(values = c("NHANES seronegativity" = "solid",
      "Programme-derived gap" = "22"),
      labels = c("Solid coloured line: NHANES seronegativity",
        "Grey dashed line: programme-derived gap"), name = NULL) +
    scale_shape_manual(values = c("NHANES seronegativity" = 16,
      "Programme-derived gap" = 1), name = NULL) +
    scale_x_continuous(limits = c(6, 14.3), breaks = c(6, 8, 10, 12, 14),
      expand = c(0, 0)) +
    scale_y_continuous(limits = c(-12, 21), breaks = c(-10, 0, 10, 20),
      labels = signed, expand = c(0, 0)) +
    labs(x = "Age (years)", y = "Change from youngest\nincluded age") +
    t("y") + figure1_bottom_legend() +
    guides(shape = "none", linetype = guide_legend(nrow = 1,
      override.aes = list(colour = c(disease_cols[["Mumps"]], charcoal))))

  # Country profiles and their summaries occupy a single ESEN2 source row;
  # source-level KiGGS/PIENTER2 estimates remain separate, with no pooling.
  source_points <- filter(data$forest, source_short != "ESEN2")
  pb <- ggplot() +
    geom_vline(xintercept = 0, colour = pub_colours$reference,
      linewidth = .3, linetype = 2) +
    geom_hline(yintercept = c(4.25, 2.25),
      colour = pub_colours$light_grey, linewidth = .2) +
    geom_point(data = data$esen, aes(contrast_pp, y_jitter, colour = disease),
      alpha = .4, size = 1.25) +
    geom_segment(data = data$esen_summary,
      aes(x = q25_pp, xend = q75_pp, y = y, yend = y, colour = disease),
      linewidth = .95) +
    geom_point(data = data$esen_summary, aes(median_pp, y, colour = disease),
      shape = 21, fill = "white", size = 2.4, stroke = .8) +
    geom_text(data = data$esen_summary,
      aes(x = 48, y = y + .37, label = sprintf("%d/%d positive", n_positive, n),
        fontface = if_else(disease == "Mumps", "bold", "plain")),
      hjust = 1, size = 2.1, colour = pub_colours$text) +
    geom_point(data = source_points,
      aes(estimate_pp, y, colour = disease, shape = source_short), size = 2) +
    scale_colour_manual(values = disease_cols, guide = "none") +
    scale_shape_manual(values = c(KiGGS = 15, PIENTER2 = 17), guide = "none") +
    scale_x_continuous(limits = c(-25, 50), breaks = c(-20, 0, 20, 40),
      labels = signed, expand = c(0, 0)) +
    scale_y_continuous(breaks = data$forest$y, labels = data$forest$row_label,
      limits = c(-1.6, 6.6), expand = c(0, 0)) +
    labs(x = "Older − younger seropositivity", y = NULL) +
    t("x") + theme(axis.text.y = element_text(size = 5.4),
      axis.ticks.y = element_blank(), axis.line.y = element_blank())

  c <- data$dtp %>% mutate(
    series_colour = if_else(endpoint == "External age marker (%)",
      unname(disease_cols[as.character(disease)]), charcoal),
    disease = factor(disease, levels = c("Diphtheria", "Pertussis")),
    construct = factor(if_else(endpoint == "External age marker (%)",
      "External marker", "Programme-derived gap"),
      levels = c("External marker", "Programme-derived gap"))
  )
  marker_labels <- c %>% filter(endpoint == "External age marker (%)") %>%
    distinct(disease, construct) %>% mutate(label = if_else(disease == "Diphtheria",
      "Antitoxin\nprotection", "Recent anti-PT\nexposure"))
  pc <- ggplot(c, aes(as.numeric(age_band), estimate_pct, colour = series_colour)) +
    geom_errorbar(data = filter(c, is.finite(lower_pct) & is.finite(upper_pct)),
      aes(ymin = lower_pct, ymax = upper_pct), width = .13, linewidth = .3) +
    geom_line(aes(group = endpoint, linetype = endpoint), linewidth = .5) +
    geom_point(aes(shape = endpoint), size = 1.45) +
    geom_text(aes(y = if_else(is.finite(upper_pct), upper_pct + 5, estimate_pct + 7),
      label = if_else(endpoint == "External age marker (%)",
        sprintf("%.0f", estimate_pct), sprintf("%.1f", estimate_pct))),
      size = 1.9, show.legend = FALSE) +
    geom_text(data = marker_labels, aes(x = 2, y = 104, label = label),
      inherit.aes = FALSE, hjust = .5, vjust = 1, lineheight = .9, size = 1.7,
      colour = pub_colours$text) +
    facet_grid(disease ~ construct, switch = "y", axes = "all_x",
      labeller = labeller(construct = c("External marker" = "External\nmarker",
        "Programme-derived gap" = "Programme-derived\ngap"))) +
    scale_colour_identity() +
    scale_shape_manual(values = c("External age marker (%)" = 16,
      "Programme-derived gap (%)" = 1), guide = "none") +
    scale_linetype_manual(values = c("External age marker (%)" = "solid",
      "Programme-derived gap (%)" = "22"), guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = levels(c$age_band),
      limits = c(.65, 3.35), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 108), breaks = c(0, 50, 100), expand = c(0, 0)) +
    labs(x = "Age band (years)", y = "Seropositivity or protection-gap estimate (%)") +
    t("y") + theme(strip.text = figure1_facet_title_style(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 5.8),
      panel.spacing.x = unit(3, "mm"), panel.spacing.y = unit(4, "mm"))

  fig <- wrap_plots(A = patchwork::free(mark(pa, "a"), type = "panel"),
    B = patchwork::free(mark(pb, "b"), type = "panel"),
    C = patchwork::free(mark(pc, "c"), type = "panel"), design = "AA\nBC",
    heights = c(1, 1.4), widths = c(1.8, 1), guides = "keep")
  fig_dir <- file.path(out, "figures/main_nc")
  src_dir <- file.path(out, "source_data/figure_redesign")
  dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
  save_publication_figure(fig, "Figure6_NC_redesigned", fig_dir, height_mm = 132, dpi = 600)
  displayed <- list(a = a, b = bind_rows(
    data$esen %>% mutate(layer = "country_profile"),
    data$esen_summary %>% mutate(layer = "ESEN2_summary", contrast_pp = median_pp),
    source_points %>% mutate(layer = "source_estimate", contrast_pp = estimate_pp)), c = c)
  for (i in names(displayed)) write_csv(displayed[[i]],
    file.path(src_dir, paste0("Figure6_", i, ".csv")))
  invisible(fig)
}

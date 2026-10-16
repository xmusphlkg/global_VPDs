# Figure 2: association, historical information, robustness and prediction limits.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(scales)})
source("scripts/05_figures/publication_figure_style.R")
source("scripts/05_figures/helpers/figure1_layout_helpers.R")
source("scripts/05_figures/prepare_figure2_evidence.R")
source("scripts/05_figures/helpers/figure2_forest_panel.R")

render_figure2_evidence <- function(out = "output/revised") {
  data <- prepare_figure2_evidence(out)
  a <- data$a; b <- data$b; c <- data$c; d <- data$d; e <- data$e; f <- data$f
  colours <- pub_evidence_colours
  fills <- colours
  cc <- function() scale_colour_manual(values = colours, guide = "none")
  fc <- function() scale_fill_manual(values = fills, guide = "none")
  t <- function(grid = "none") theme_publication(grid = grid, show_titles = TRUE) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(1, 2, 1, 2), axis.title = figure1_axis_title_style(),
      axis.title.x = figure1_axis_title_style(margin(t = 2)),
      axis.title.y = figure1_axis_title_style(margin(r = 2)),
      axis.text = element_text(size = 5.4), axis.ticks.length = unit(.7, "mm"),
      plot.title = element_text(size = 6.2, face = "bold"),
      plot.subtitle = figure1_facet_title_style(),
      strip.text = figure1_facet_title_style(), legend.text = element_text(size = 5.2),
      legend.title = figure1_axis_title_style(),
      panel.spacing = unit(3, "mm"), panel.spacing.x = unit(3.5, "mm"),
      legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
      legend.box.just = "center", legend.box.spacing = unit(1, "mm"),
      legend.spacing.y = unit(.4, "mm"),
      legend.key.height = unit(2.4, "mm"), legend.key.width = unit(3.5, "mm"))
  mark <- function(p, letter) p + labs(title = letter) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 3)), plot.title.position = "plot")
  # All d facets retain axes/ticks; labels appear only on the outer margins.
  null_x <- function(x = 0) geom_vline(xintercept = x, colour = "#9EA3A6", linetype = 2, linewidth = .25)
  pa <- mark(figure2_forest_panel(a, colours, fills), "a")

  # The three historical representations each add one alternative term to the
  # current-coverage baseline; they are not nested.
  bc <- filter(b, layer == "disease_comparison") %>% mutate(y = match(disease, rev(pub_disease_levels)),
    x = match(comparison, c("M1 vs M0", "M2 vs M1", "M3 vs M1", "M3 vs M2")))
  b_matrix <- ggplot(bc, aes(x, y, fill = delta_aic)) + geom_tile(width = .94, height = .87) +
    geom_text(aes(label = sprintf("%+.1f", ifelse(abs(delta_aic) < .05, 0, delta_aic))), size = 2) +
    scale_fill_gradientn(colours = pub_gradients$diverging_full,
      limits = c(-11, 11), guide = "none") +
    scale_x_continuous(breaks = 1:4, labels = c("Coverage\nhistory\nvs current\ncoverage", "Cohort gap\nvs coverage\nhistory",
        "Cohort-time\ngap vs coverage\nhistory", "Cohort-time\ngap vs\ncohort gap"),
      limits = c(.5, 4.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:5, labels = rev(pub_disease_levels),
      limits = c(.5, 5.5), expand = c(0, 0)) +
    labs(x = "ΔAIC (alternative − reference)", y = NULL) + t() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      axis.text.x = element_text(size = 5.0, lineheight = .95))
  pb <- mark(b_matrix, "b")
  top_row <- figure2_forest_panel(a, colours, fills, comparison = b_matrix)

  c_curve <- filter(c, layer == "standardized_curve") %>% mutate(disease = factor(disease, levels = pub_disease_levels))
  c_rug <- filter(c, layer == "observed_gap") %>% mutate(disease = factor(disease, levels = pub_disease_levels))
  c_support <- c_curve %>% group_by(disease, central_support_low_pp, central_support_high_pp) %>%
    summarise(band_low = min(rate_low), band_high = max(rate_high), .groups = "drop")
  pc <- ggplot(c_curve, aes(gap_pp, rate_per_100k, colour = disease)) +
    geom_rect(data = c_support, aes(xmin = central_support_low_pp, xmax = central_support_high_pp,
      ymin = band_low, ymax = band_high), inherit.aes = FALSE, fill = "#F0F1F2", colour = NA) +
    geom_ribbon(aes(ymin = rate_low, ymax = rate_high, fill = disease), alpha = .14, colour = NA) +
    geom_line(linewidth = .65) +
    geom_rug(data = c_rug, aes(x = gap_pp, y = NULL), sides = "b", length = unit(1, "mm"), alpha = .16, linewidth = .22) +
    facet_wrap(~disease, nrow = 1, scales = "free") + cc() + fc() +
    scale_y_log10(breaks = c(.001, .01, .1, 1, 10, 100, 1000), labels = label_number()) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3), expand = expansion(mult = .05)) +
    labs(x = "Cohort-time gap (%)", y = "Cases\nper 100,000 population") + t("y")
  pc <- mark(pc, "c") + theme(strip.text = element_blank(), plot.margin = margin(1, 5, 1, 2))

  dd <- filter(d, layer == "decomposition") %>% mutate(disease = factor(disease, levels = pub_disease_levels),
    y = match(component, c("Between", "Within", "Overall")))
  # Keep the log-rate-ratio axis sparse enough to remain legible in the
  # narrow stacked facets while covering the full observed interval.
  d_forest_x_limits <- 10^range(floor(log10(dd$rr_low)), ceiling(log10(dd$rr_high)))
  d_forest_x_breaks <- sort(unique(c(d_forest_x_limits, 1)))
  d_y_breaks <- pretty(range(dd$y), n = length(unique(dd$y)))
  d_y_breaks <- sort(unique(d_y_breaks[d_y_breaks %in% dd$y]))
  d_y_limits <- range(d_y_breaks)
  d_forest <- ggplot(dd, aes(rr, y, colour = disease)) + null_x(1) +
    geom_segment(aes(x = rr_low, xend = rr_high, yend = y), linewidth = .5) +
    geom_point(size = 1.3) + facet_wrap(~disease, nrow = 1) + cc() +
    scale_x_log10(breaks = d_forest_x_breaks, limits = d_forest_x_limits,
      labels = function(x) vapply(x, function(v) format(v, scientific = FALSE, trim = TRUE), character(1)), expand = expansion(mult = .08)) +
    scale_y_continuous(breaks = d_y_breaks, labels = c("Between", "Within", "Overall")[d_y_breaks],
      limits = d_y_limits, expand = expansion(mult = .12)) +
    labs(x = "Rate ratio (95% CI)", y = NULL) + t() + theme(strip.text = element_blank())
  dl <- filter(d, layer == "leave_one_country_out") %>% mutate(disease = factor(disease, levels = pub_disease_levels))
  d_loco_x_breaks <- pretty(range(c(dl$estimate, 0)), n = 3)
  d_loco_x_limits <- range(d_loco_x_breaks)
  d_loco_y_breaks <- pretty(c(0, 1), n = 2)
  d_loco_y_limits <- range(d_loco_y_breaks)
  d_loco <- ggplot(dl, aes(estimate, after_stat(scaled), colour = disease, fill = disease)) + null_x() +
    geom_density(alpha = .18, linewidth = .45) +
    geom_rug(aes(y = NULL), sides = "b", alpha = .18, linewidth = .2, length = unit(1, "mm")) +
    facet_wrap(~disease, nrow = 1) + cc() + fc() +
    scale_x_continuous(limits = d_loco_x_limits, breaks = d_loco_x_breaks, expand = c(0, 0)) +
    scale_y_continuous(limits = d_loco_y_limits, breaks = d_loco_y_breaks,
      expand = expansion(mult = .08)) +
    labs(x = "Leave-one-out log(RR)", y = NULL) + t() +
    theme(axis.title.y = element_blank(), strip.text = element_blank())
  pd <- wrap_plots(mark(d_forest, "d"), d_loco, ncol = 1, heights = c(.68, 1))

  # Discrete 3x3 sensitivity grids. Label the M3-minus-M1 AIC changes directly instead of implying
  # that a smooth, densely sampled response surface has been fitted.
  sensitivity_plot <- function(dis, tag = NULL) {
    z <- filter(e, disease == dis) %>% mutate(x = match(tau, sort(unique(tau))),
      y = match(partial_series_coefficient, sort(unique(partial_series_coefficient))))
    p <- ggplot(z, aes(x, y, fill = delta_aic_m3_m1)) + geom_tile(width = .96, height = .92) +
      geom_text(aes(label = sprintf("%+.1f", ifelse(abs(delta_aic_m3_m1) < .05, 0, delta_aic_m3_m1))), size = 2) +
      geom_text(data = filter(z, primary), aes(x = x + .29, y = y + .29), label = "*", size = 3.3) +
    scale_fill_gradientn(colours = pub_gradients$diverging_full,
      limits = c(-11, 11), breaks = c(-10, 0, 10), name = "ΔAIC (cohort-time gap − coverage history)",
        guide = figure1_colourbar("ΔAIC (cohort-time gap − coverage history)", width = 24)) +
      scale_x_continuous(breaks = 1:3, labels = sort(unique(z$tau)), limits = c(.5, 3.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = 1:3, labels = sort(unique(z$partial_series_coefficient)), limits = c(.5, 3.5), expand = c(0, 0)) +
      labs(x = "Duration (years)", y = if (dis == "Pertussis") "Partial-series\nprotection coefficient" else NULL, subtitle = dis) + t()
    if (!is.null(tag)) p <- mark(p, tag)
    p
  }
  pe <- wrap_plots(sensitivity_plot("Pertussis", "e"), sensitivity_plot("Diphtheria"),
    nrow = 1, guides = "collect") & figure1_bottom_legend()

  ff <- filter(f, layer == "disease_holdout") %>%
    mutate(period = paste0("Train through ", cutoff, "\nTest ", test_start, "–", test_end))
  pf <- ggplot(ff, aes(delta_aic_m3_m1, rmse_change_pct)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "#FBF2E9") +
    geom_hline(yintercept = 0, colour = "#9EA3A6", linewidth = .25, linetype = 2) + null_x() +
    geom_point(aes(colour = disease, fill = disease), shape = 21, size = 1.9, stroke = .55) +
    ggrepel::geom_text_repel(aes(label = disease, colour = disease),
      size = 1.95, seed = 20261006, min.segment.length = 0, segment.size = .2,
      box.padding = .35, point.padding = .22, max.overlaps = Inf, show.legend = FALSE) +
    facet_wrap(~period, nrow = 1) + cc() + fc() +
    scale_x_continuous(limits = c(-12, 14), breaks = c(-10, 0, 10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-20, 22), breaks = c(-20, 0, 20), expand = c(0, 0)) +
    labs(x = "ΔAIC (cohort-time gap − coverage history)", y = "Held-out RMSE change (%)") + t()
  pf <- mark(pf, "f")
  blocks <- lapply(list(pa, pb, pc, pd, pe, pf), function(p) wrap_elements(full = p))
  fig <- wrap_plots(
    wrap_elements(full = top_row),
    blocks[[3]], blocks[[4]],
    wrap_plots(blocks[[5]], blocks[[6]], nrow = 1, widths = c(1, 1.05)),
    ncol = 1, heights = c(.9, .82, 1.3, 1), guides = "keep")
  save_publication_figure(fig, "Figure2_NC_redesigned", file.path(out, "figures/main_nc"), height_mm = 180, dpi = 600)
  src <- file.path(out, "source_data/figure_redesign")
  dir.create(src, recursive = TRUE, showWarnings = FALSE)
  publication_comparison_names <- c(
    "M1 vs M0" = "Coverage history vs current coverage",
    "M2 vs M1" = "Constant-protection cohort gap vs coverage history",
    "M3 vs M1" = "Cohort-time gap vs coverage history",
    "M3 vs M2" = "Cohort-time gap vs constant-protection cohort gap",
    "M4 interaction vs M3" = "Disease-by-gap interaction vs cohort-time gap"
  )
  for (i in names(data)) {
    exported <- data[[i]]
    if (i == "b") exported$comparison <- unname(publication_comparison_names[exported$comparison])
    write_csv(exported, file.path(src, paste0("Figure2_", i, ".csv")))
  }
  invisible(fig)
}

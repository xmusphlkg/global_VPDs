# Figure 4: history and duration as separate axes of paired scenarios.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(scales)})
source("scripts/05_figures/publication_figure_style.R")
source("scripts/05_figures/helpers/figure1_layout_helpers.R")
source("scripts/05_figures/prepare_figure4_scenarios.R")
source("scripts/05_figures/helpers/figure4_country_forest_panel.R")
source("scripts/05_figures/helpers/figure4_scenario_style.R")

render_figure4_scenarios <- function(out = "output/revised") {
  data <- prepare_figure4_scenarios(out)
  colours <- pub_evidence_colours
  fills <- colours
  cc <- function() scale_colour_manual(values = colours, guide = "none")
  fc <- function() scale_fill_manual(values = fills, guide = "none")
  local_inside <- function(x = .5, y = .9, direction = "vertical") theme(
    legend.position = "inside", legend.position.inside = c(x, y),
    legend.justification.inside = c("center", "top"), legend.direction = direction,
    legend.background = element_rect(fill = scales::alpha("white", .9), colour = NA),
    legend.key.height = unit(2.3, "mm"), legend.spacing.y = unit(.3, "mm"))
  t <- function(grid = "none") theme_publication(grid = grid, show_titles = TRUE) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(1, 2, 1, 2), axis.title = figure1_axis_title_style(),
      axis.title.x = figure1_axis_title_style(margin(t = 2)),
      axis.title.y = figure1_axis_title_style(margin(r = 2)),
      axis.text = element_text(size = 5.4), axis.ticks.length = unit(.7, "mm"),
      plot.subtitle = figure1_facet_title_style(), strip.text = figure1_facet_title_style(),
      panel.spacing.x = unit(3.5, "mm"), legend.text = element_text(size = 5.2),
      legend.title = figure1_axis_title_style(), legend.key.height = unit(2.4, "mm"),
      legend.key.width = unit(3.5, "mm"), legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"))
  mark <- function(p, letter) p + labs(title = letter) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 3)),
      plot.title.position = "plot")
  rows <- function() scale_y_continuous(breaks = 1:5, labels = rev(pub_disease_levels),
    limits = c(.5, 6.25), expand = c(0, 0))
  null <- function(x) geom_vline(xintercept = x, linetype = 2, colour = "#9EA3A6", linewidth = .25)
  cs <- filter(data$c, layer == "summary") %>% mutate(y = match(disease, rev(pub_disease_levels)))

  # a: posterior means of annual expected reports, equally averaged over
  # countries within each draw. Country trajectories and intervals are in S4–S8.
  scenario_cols <- figure4_scenario_colours
  scenario_labels <- figure4_scenario_labels
  aa <- filter(data$a, layer == "annual_summary") %>% mutate(
    disease = factor(disease, levels = pub_disease_levels),
    disease_label = factor(as.character(disease), levels = pub_disease_levels),
    scenario = paste(history, duration, sep = "_"),
    scenario = factor(scenario, levels = names(scenario_cols)))
  # Per-facet y breaks: free_y lets each disease facet pick pretty breaks
  # from its own data range, so small-count diseases aren't squeezed by
  # a shared scale designed for a large-count disease.  limits are set to
  # range(breaks) per facet so the axis sits exactly on the tick endpoints.
  pa <- ggplot(aa, aes(year, mean, group = scenario, colour = scenario)) +
    geom_line(aes(linetype = duration), linewidth = .55) +
    facet_wrap(~disease_label, nrow = 1, scales = "free_y", axes = "all", axis.labels = "margins") +
    scale_colour_manual(values = scenario_cols, breaks = names(scenario_cols),
      labels = unname(scenario_labels), name = NULL) +
    scale_linetype_manual(values = figure4_scenario_linetypes, guide = "none") +
    scale_x_continuous(limits = c(2013.7, 2023.3), breaks = c(2014, 2016, 2018, 2020, 2022),
      expand = c(0, 0)) +
    scale_y_continuous(limits = function(values) range(pretty(values)),
      labels = label_number(scale_cut = cut_short_scale()), expand = c(0, 0)) +
    labs(x = "Calendar year", y = "Annual expected reports\nper country") + t("y") +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = .8))) +
    figure1_bottom_legend()
  pa <- mark(pa, "a")

  # b: posterior duration contrasts on a declared signed log axis with a linear core.
  # Exact-zero durable references remain explicit point summaries without a
  # fabricated density.
  trans <- list(transform = figure4_signed_log)
  densities <- filter(data$b, layer == "draw") %>% group_by(disease) %>% group_modify(~{
    x <- trans$transform(.x$delta_cases / 1000)
    if (diff(range(x)) < 1e-12) return(tibble(x = numeric(), height = numeric()))
    den <- density(x, from = min(x), to = max(x), n = 256)
    tibble(x = den$x, height = .48 * den$y / max(den$y))
  }) %>% ungroup() %>% mutate(y = match(disease, rev(pub_disease_levels)))
  bs <- filter(data$b, layer == "summary") %>%
    mutate(y = match(disease, rev(pub_disease_levels)),
      x = trans$transform(median / 1000), xl = trans$transform(lo / 1000),
      xh = trans$transform(hi / 1000),
      label = if_else(median > 0, paste0("+", comma(round(median))), comma(round(median))),
      # Positive contrasts are labelled left of their interval at row height so
      # the text clears both the posterior density and the interval whiskers.
      label_x = if_else(median > 0, xl, xh),
      label_y = if_else(median > 0, y, y - .32),
      label_hjust = if_else(median > 0, 1.12, -.08))
  # Leave room to the right of the largest interval for its text annotation;
  # the previous 5% pad placed the pertussis label beyond the panel boundary.
  b_upper <- max(1100, max(bs$hi, na.rm = TRUE) / 1000 * 1.3)
  b_breaks <- c(-1, 0, 1, 10, 100, 1000)
  pb <- ggplot(densities, aes(x, colour = disease, fill = disease)) + null(0) +
    geom_ribbon(aes(ymin = y, ymax = y + height), alpha = .28, linewidth = .3) +
    geom_segment(data = bs, aes(x = xl, xend = xh, y = y, yend = y), linewidth = .55) +
    geom_point(data = bs, aes(x = x, y = y), shape = 21, size = 1.7, stroke = .5) +
    geom_text(data = bs, aes(x = label_x, y = label_y,
      label = label, hjust = label_hjust),
      inherit.aes = FALSE, size = 2.05, colour = pub_colours$text, show.legend = FALSE) +
    cc() + fc() + rows() +
    scale_x_continuous(breaks = trans$transform(b_breaks),
      labels = c("−1", "0", "1", "10", "100", "1000"),
      limits = trans$transform(c(-1, b_upper)),
      expand = expansion(mult = c(0, .08))) +
    labs(x = "Waning − no waning\n(10³ reports)", y = NULL) + t() +
    theme(plot.margin = margin(1, 4, 1, 2),
      axis.text.y = element_text(size = 5.4))
  pb <- mark(pb, "b")

  # c: history effects under waning, expressed as reductions from the observed
  # history baseline so zero is the common left-hand reference.
  ckey <- c(recovery_5y = "5-year recovery", sustained_high = "10–14-year reference")
  cbase <- filter(cs, history == "observed") %>% select(disease, y) %>%
    left_join(filter(cs, history == "sustained_high") %>% select(disease, low = reduction_percent), by = "disease")
  cplot <- filter(cs, history != "observed") %>% mutate(
    estimate = reduction_percent,
    lo_reduction = reduction_lo, hi_reduction = reduction_hi)
  pc <- ggplot(cplot, aes(estimate, y, colour = disease, fill = disease)) + null(0) +
    geom_segment(data = cbase, aes(x = 0, xend = low, yend = y), colour = "#DADCDD", linewidth = .45) +
    geom_segment(aes(x = lo_reduction, xend = hi_reduction, yend = y), linewidth = .6) +
    geom_point(aes(shape = history), size = 1.7, stroke = .5) +
    cc() + fc() + scale_shape_manual(values = c(recovery_5y = 24, sustained_high = 22),
      breaks = names(ckey), labels = ckey, name = NULL) +
    guides(shape = guide_legend(ncol = 1, override.aes = list(colour = "#444444", fill = "#444444"))) +
    rows() + scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
      expand = expansion(mult = c(0, .08))) +
    labs(x = "Reduction relative to\nobserved history (%)", y = NULL,
      subtitle = NULL) + t() + local_inside(.5, 1, direction = "vertical") +
      theme(axis.text.y = element_blank(), plot.margin = margin(1, 4, 1, 2))
  pc <- mark(pc, "c")

  # d: country heterogeneity, focused on the three duration-sensitive diseases.
  pd <- mark(figure4_country_forest_panel(data$d, colours,
    diseases = c("Pertussis", "Diphtheria", "Mumps")), "d")

  # d: cumulative trajectories retain the recovery information removed from
  # the former one-dimensional history panel.
  es <- filter(data$e, layer == "summary") %>%
    mutate(disease = factor(disease, levels = c("Pertussis", "Diphtheria", "Mumps")))
  scenario <- c(observed = "Observed history", recovery_5y = "5-year recovery", sustained_high = "10–14-year reference")
  endpoint_labels <- es %>% filter(year == 2023, history %in% c("recovery_5y", "sustained_high")) %>%
    select(disease, history, median) %>% pivot_wider(names_from = history, values_from = median) %>%
    mutate(x = 2022.85, y = 12,
      label = paste0("5-year −", round(100 - recovery_5y, 1), "%\n10–14-year −", round(100 - sustained_high, 1), "%"))
  pe <- ggplot(es, aes(year, median, group = history)) +
    annotate("rect", xmin = 2018, xmax = 2022, ymin = 0, ymax = 112, fill = "#F4F5F5") +
    geom_vline(xintercept = 2018, colour = "#9EA3A6", linewidth = .25, linetype = 2) +
    geom_ribbon(data = filter(es, history != "observed"), aes(ymin = lo, ymax = hi, fill = disease),
      alpha = .12, colour = NA) +
    geom_line(data = filter(es, history == "observed"), aes(linetype = history), colour = "#52575A", linewidth = .45) +
    geom_line(data = filter(es, history != "observed"), aes(colour = disease, linetype = history, alpha = history), linewidth = .65) +
    geom_text(data = endpoint_labels, aes(x = x, y = y, label = label, colour = disease),
      inherit.aes = FALSE, hjust = 1, vjust = 0, size = 2.05, lineheight = .9, show.legend = FALSE) +
    facet_wrap(~disease, nrow = 1) + cc() + fc() +
    scale_alpha_manual(values = c(recovery_5y = 1, sustained_high = .6), guide = "none") +
    scale_linetype_manual(values = c(observed = "solid", recovery_5y = "longdash", sustained_high = "dotted"),
      breaks = names(scenario), labels = scenario, name = NULL) +
    guides(linetype = guide_legend(nrow = 1, override.aes = list(colour = "#52575A", alpha = 1))) +
    scale_x_continuous(limits = c(2013.65, 2023.35), breaks = c(2014, 2018, 2023), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 112), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
    labs(x = "Calendar year", y = "Cumulative expected reports\n(% of observed-history total)") + t("y") +
    theme(panel.spacing.x = unit(5, "mm")) +
    figure1_bottom_legend()
  pe <- mark(pe, "e")

  # e/f: separate the observed-versus-predictive adequacy plot from its
  # disease-level coverage summary so each panel has enough horizontal room.
  ff <- data$f
  cov <- filter(data$g, layer == "coverage") %>% mutate(y = match(disease, rev(pub_disease_levels)))
  f_breaks <- c(0, 100, 10000, 1000000)
  f_limit <- 1500000
  scatter <- ggplot(ff, aes(cases, pred_median, colour = disease, fill = disease)) +
    geom_abline(slope = 1, intercept = 0, linewidth = .3, colour = "#9EA3A6", linetype = 2) +
    geom_point(shape = 21, size = .9, stroke = .2, alpha = .5) + cc() + fc() +
    scale_x_continuous(transform = "log1p", breaks = f_breaks,
      labels = c("0", "100", "10k", "1m"), limits = c(0, f_limit), expand = c(0, 0)) +
    scale_y_continuous(transform = "log1p", breaks = f_breaks,
      labels = c("0", "100", "10k", "1m"), limits = c(0, f_limit), expand = c(0, 0)) +
    labs(x = "Observed reports",
      y = "Held-out negative-binomial\npredictive median") + t() +
    # Identical transformed limits keep the equality reference at 45 degrees.
    coord_equal()
  coverage <- ggplot(cov, aes(coverage_percent, y, colour = disease, fill = disease)) + null(95) +
    geom_segment(aes(x = 80, xend = coverage_percent, yend = y), alpha = .45, linewidth = .55) +
    geom_point(shape = 21, size = 1.35, stroke = .4) +
    cc() + fc() +
    scale_y_continuous(breaks = cov$y, labels = cov$disease,
      limits = c(.5, 5.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(80, 101), breaks = c(80, 90, 100), expand = c(0, 0)) +
    labs(x = "95% predictive-interval\ncoverage (%)", y = NULL) + t() +
    theme(axis.text.y = element_text(size = 5))
  pf <- mark(scatter, "f")
  pg <- mark(coverage, "g")

  row2 <- wrap_plots(patchwork::free(pb),
    patchwork::free(pc), pd, nrow = 1,
    widths = c(1.25, 1, 3.5), guides = "keep")
  row3 <- wrap_plots(patchwork::free(pe, type = "panel"), pf, pg, nrow = 1,
    widths = c(3, 1, 1), guides = "keep")
  blocks <- list(
    patchwork::free(pa, type = "panel"),
    row2, row3
  )
  fig <- wrap_plots(blocks, ncol = 1, heights = c(1.2, 1.1, 0.85), guides = "keep")
  save_publication_figure(fig, "Figure4_NC_redesigned", file.path(out, "figures/main_nc"), height_mm = 145, dpi = 600)
  src <- file.path(out, "source_data/figure_redesign")
  dir.create(src, recursive = TRUE, showWarnings = FALSE)
  displayed_data <- list(a = data$a, b = data$b, c = data$c,
    d = filter(data$d, disease %in% c("Pertussis", "Diphtheria", "Mumps")),
    e = data$e, f = data$f, g = data$g)
  for (i in names(displayed_data)) write_csv(displayed_data[[i]], file.path(src, paste0("Figure4_", i, ".csv")))
  interaction <- filter(data$a, layer == "interaction_summary") %>%
    select(disease, median, lo50, hi50, lo, hi, n_draws) %>%
    mutate(estimand = "(Observed waning - Observed no waning) - (High waning - High no waning)",
      unit = "Expected reports per country, 2014-2023", n_countries = 10)
  table_dir <- file.path(out, "supplementary_tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(interaction, file.path(out, "bayesian_transmission/history_waning_interaction.csv"))
  write_csv(interaction, file.path(table_dir, "TableS43_history_waning_interaction.csv"))
  invisible(fig)
}

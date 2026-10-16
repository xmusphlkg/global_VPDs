# Six composite panels; recovery-event counts are shown in Supplementary Figure S2E.
source("scripts/05_figures/prepare_figure1_atlas.R")
source("scripts/05_figures/helpers/figure1_layout_helpers.R")

render_figure1_atlas <- function(out = "output/revised") {
  data <- prepare_figure1_atlas(out)
  a <- data$a; b <- data$b; c <- data$c; d <- data$d; e <- data$e
  f <- data$f
  ax <- figure1_axis_factory()
  country_colours <- c(TZA = pub_colours$blue, UKR = pub_colours$orange)
  country_names <- c(TZA = "Tanzania", UKR = "Ukraine")
  # Display labels are local to Figure 1; reconstruction scenarios stay intact.
  disease_labels <- setNames(pub_disease_levels, pub_disease_levels)
  t <- function(grid = "none") theme_publication(grid = grid, show_titles = TRUE) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(1, 2, 1, 2), axis.title = figure1_axis_title_style(),
      axis.title.x = figure1_axis_title_style(margin(t = 2)),
      axis.title.y = figure1_axis_title_style(margin(r = 2)),
      axis.text = element_text(size = 5.4), axis.ticks.length = unit(.7, "mm"),
      plot.subtitle = element_text(size = 5.5, margin = margin(b = 2)),
      legend.title = figure1_axis_title_style(), legend.text = element_text(size = 5.2),
      legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
      legend.box.just = "center",
      legend.box.spacing = unit(1, "mm"), legend.spacing.y = unit(.4, "mm"),
      legend.key.height = unit(2.4, "mm"), legend.key.width = unit(3.5, "mm"))
  local_inside <- function(x, y, justification = c(0, 0)) theme(
    legend.position = "inside", legend.position.inside = c(x, y),
    legend.justification.inside = justification,
    legend.background = element_rect(fill = scales::alpha("white", .88), colour = NA))
  gap_fill <- function(guide = "none", breaks = pretty(c(0, 100), n = 2), title = "Programme-derived gap (%)") scale_fill_gradientn(
    colours = pub_gradients$sequential_red,
    limits = range(breaks), breaks = breaks,
    name = title, guide = guide)
  # Integer calendar labels stay sparse; every underlying break comes from pretty.
  year_labels <- function(x) ifelse(x %% 2 == 0 & x <= 2023,
    format(x, scientific = FALSE, trim = TRUE), "")
  year_x <- function(id) ax$x(id, c(2014, 2023), n = 10, labels = year_labels,
    limits = c(2013.5, 2023.5))
  age_y <- function(id, labels = waiver()) ax$y(id, c(.5, 14.5), n = 3, labels = labels)
  hide_y <- theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
  row_y <- function(id, labels = TRUE) scale_y_continuous(
    breaks = 1:5, limits = c(.5, 5.5), expand = c(0, 0),
    labels = if (labels) rev(pub_disease_levels) else NULL)

  # A: retain alignment of year and age scales, while freeing the margins that
  # previously made the small coverage plot reserve space for the large heatmap.
  band <- filter(a, layer == "low_coverage_year")
  a_cov <- ggplot(filter(a, layer == "coverage"), aes(year, coverage)) +
    geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = pub_colours$light_grey, alpha = .45) +
    geom_line(colour = country_colours[["UKR"]], linewidth = .6) +
    geom_point(size = .65, colour = country_colours[["UKR"]]) +
    annotate("text", x = 2015, y = 90, label = "<30%", size = 1.95,
      colour = pub_colours$muted_text) +
    annotate("text", x = 2022.65, y = 58, label = "83.0%", size = 1.95,
      colour = country_colours[["UKR"]]) +
    year_x("A_coverage_x") + ax$y("A_coverage_y", c(0, 100), n = 2) +
    labs(x = NULL, y = "DTP3\ncoverage (%)") + t() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          plot.margin = margin(0, 2, 8, 2))

  a_heat <- ggplot(filter(a, layer == "landscape"), aes(year, age, fill = gap_pp)) +
    geom_tile(width = 1, height = 1) +
    geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = .5, ymax = 14.5),
      inherit.aes = FALSE, fill = pub_colours$reference, alpha = .12) +
    geom_segment(data = filter(a, layer == "cohort_diagonal"), aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, colour = "white", linewidth = .35, alpha = .85) +
    annotate("segment", x = 2019.5, xend = 2020.35, y = 12.1, yend = 12.35,
      colour = "white", linewidth = .3, arrow = arrow(length = unit(1, "mm"))) +
    annotate("text", x = 2019.4, y = 12.1, label = "Birth cohorts", hjust = 1,
      colour = "white", size = 2.1) +
    gap_fill(figure1_colourbar("Programme-derived gap (%)", width = 30)) +
    year_x("A_heatmap_x") + age_y("A_heatmap_y") +
    labs(x = "Calendar year", y = "Age (years)") + t() +
    figure1_bottom_legend()

  a_inner <- wrap_plots(
    a_cov,
    a_heat,
    design = "A\nB",heights = c(.85, 3.8), guides = "keep")

  # B: keep the country key next to its trajectories and the deficit key under
  # its two strips. The single endpoint label avoids two labels colliding at 83%.
  b_cov <- filter(b, layer == "coverage_history")
  b_end <- filter(b_cov, focal, year == 2023)
  b_start <- filter(b_cov, focal, year == min(year)) %>% mutate(
    coverage_label = sprintf("%.1f%%", coverage),
    label_y = coverage + if_else(coverage > 80, -8, 8))
  b_plot <- ggplot(b_cov, aes(year, coverage, group = country)) +
    geom_line(data = filter(b_cov, !focal), colour = "#CDD1D3", linewidth = .26, alpha = .65) +
    geom_line(data = filter(b_cov, focal), aes(colour = country), linewidth = .65) +
    geom_point(data = b_start, aes(colour = country), size = 1.4, show.legend = FALSE) +
    geom_text(data = b_start, aes(x = year + .16, y = label_y, label = coverage_label, colour = country),
      hjust = 0, size = 1.9, show.legend = FALSE) +
    geom_point(data = b_end, aes(colour = country), size = 1.4, show.legend = FALSE) +
    annotate("text", x = 2023.42, y = 63, label = "83.0%", hjust = 1, size = 2.0) +
    scale_colour_manual(values = country_colours, labels = country_names, name = NULL) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = .65, alpha = 1))) +
    ax$y("B_coverage_y", c(0, 100), n = 2) + year_x("B_coverage_x") +
    labs(x = NULL, y = "DTP3 coverage (%)") + t("y") +
    local_inside(.63, .01) + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
  b_strip_data <- filter(b, layer == "annual_deficit") %>% mutate(row = if_else(country == "TZA", 2, 1))
  b_strip <- ggplot(b_strip_data, aes(year, row, fill = deficit_pp)) +
    geom_tile(width = .94, height = .88) +
    scale_fill_gradient(low = "#FAEFE4", high = "#A8431E", limits = range(pretty(c(0, 100), n = 2)),
      breaks = pretty(c(0, 100), n = 2), guide = figure1_colourbar("Annual coverage deficit (%)", width = 23)) +
    year_x("B_strip_x") + ax$y("B_strip_y", c(.5, 2.5), n = 3,
      labels = figure1_carrier_labels(c(1, 2), c("UKR", "TZA"))) +
    labs(x = "Calendar year", y = NULL) + t() + hide_y +
    figure1_bottom_legend() + theme(plot.margin = margin(0, 2, 0, 2))
  b_inner <- b_plot / b_strip + plot_layout(heights = c(2.5, .8), guides = "keep")

  # C: distinguish country totals from the single-year, age-group and sensitivity
  # encodings in a local key below the trajectories.
  c_profile <- filter(c, layer == "current_profile")
  c_groups <- filter(c, layer == "age_group_mean")
  c_totals <- filter(c, layer == "population_weighted_gap")
  c_key <- setNames(sprintf("%s  %.1f%%", country_names[c_totals$country], c_totals$gap_pp), c_totals$country)
  c_marks <- c("Single-year", "Age-group mean", "Sensitivity range")
  c_plot <- ggplot(c_profile, aes(age, gap_pp, colour = country, fill = country)) +
    annotate("rect", xmin = .5, xmax = 4.5, ymin = 0, ymax = 100, fill = "#EDF3F8", alpha = .65) +
    annotate("rect", xmin = 9.5, xmax = 14.5, ymin = 0, ymax = 100, fill = "#FBF3EB", alpha = .8) +
    geom_ribbon(aes(ymin = lower_pp, ymax = upper_pp, alpha = "Sensitivity range"),
      colour = NA, show.legend = c(alpha = NA, fill = FALSE)) +
    geom_line(aes(linetype = "Single-year"), linewidth = .65) +
    geom_point(data = c_groups, aes(shape = "Age-group mean"),
      fill = "white", size = 1.6, stroke = .5,
      show.legend = c(shape = NA, colour = FALSE)) +
    geom_rug(sides = "b", length = unit(.75, "mm"), alpha = .65, show.legend = FALSE) +
    scale_colour_manual(values = country_colours, labels = c_key,
      name = "Country-weighted gap (%)") +
    scale_fill_manual(values = country_colours, guide = "none") +
    # Identical titles and breaks merge these aesthetics into one mark guide.
    # Each layer contributes only its own key: line, open point or ribbon.
    scale_linetype_manual(values = setNames(rep("solid", 3), c_marks),
      limits = c_marks, name = "Mark type") +
    scale_shape_manual(values = setNames(rep(21, 3), c_marks),
      limits = c_marks, name = "Mark type") +
    scale_alpha_manual(values = setNames(rep(.13, 3), c_marks),
      limits = c_marks, name = "Mark type") +
    guides(
      colour = guide_legend(ncol = 1, order = 1),
      linetype = guide_legend(ncol = 1, order = 2,
        override.aes = list(colour = c(pub_colours$text, pub_colours$text, NA),
          fill = c(NA, "white", pub_colours$reference), alpha = c(1, 1, .13))),
      shape = guide_legend(ncol = 1, order = 2),
      alpha = guide_legend(ncol = 1, order = 2)) +
    ax$x("C_profile_x", c(0, 15), n = 3) + ax$y("C_profile_y", c(0, 100), n = 4) +
    labs(x = "Age (years)", y = "Programme-derived\nprotection gap (%)") + t("y") +
    local_inside(.97, .09, c(1, 0)) +
    theme(legend.box = "horizontal", legend.box.just = "top",
      legend.spacing.x = unit(3, "mm"), legend.title.position = "top")

  # D: shared row positions, but independent horizontal margins for the residual,
  # SD and matching-rule axes. Connected symbols group each disease comparison.
  d_rain <- ggplot() +
    geom_vline(xintercept = 0, colour = "#999999", linewidth = .23, linetype = 2) +
    geom_ribbon(data = filter(d, layer == "half_density"), aes(residual_pp, ymin = y, ymax = upper_y, group = disease),
      fill = pub_colours$blue, alpha = .23) +
    geom_point(data = filter(d, layer == "residual"), aes(residual_pp, point_y), size = .28, alpha = .22, colour = "#55595D") +
    geom_segment(data = filter(d, layer == "residual_quartiles"), aes(x = q25_pp, xend = q75_pp, y = y, yend = y), linewidth = .7) +
    geom_point(data = filter(d, layer == "residual_quartiles"), aes(median_pp, y), size = 1.1) +
    ax$x("D_residual_x", d$residual_pp, n = 7) + row_y("D_residual_y") +
    labs(x = "Within-stratum gap residual\n(percentage points)", y = "Disease") + t() +
    theme(plot.margin = margin(1, 5, 1, 2))
  d_sd <- filter(d, layer == "sd_summary", configuration == "exact_equal_country_year")
  d_exact <- ggplot(d_sd, aes(gap_sd_pp, y)) +
    geom_segment(aes(x = gap_sd_pp_low, xend = gap_sd_pp_high, yend = y), linewidth = .5, colour = pub_colours$blue) +
    geom_point(size = 1.25, colour = pub_colours$blue) +
    geom_text(aes(x = 3.9, y = y + .23, label = sprintf("%.2f", gap_sd_pp)), hjust = 1, size = 1.85) +
    ax$x("D_sd_x", c(0, d_sd$gap_sd_pp_high, 4), n = 2) + row_y("D_sd_y", FALSE) +
    labs(x = "Residual SD\n(percentage points)", y = NULL) + t() + hide_y +
    theme(plot.margin = margin(1, 5, 1, 4))
  d_sensitivity <- ggplot(filter(d, layer == "sd_summary") %>% arrange(disease, specification),
    aes(gap_sd_pp, glyph_y, shape = specification)) +
    geom_path(aes(group = disease), colour = scales::alpha(pub_colours$orange, .55), linewidth = .35) +
    geom_point(size = 1.25, colour = pub_colours$orange) +
    scale_shape_manual(values = c("Exact" = 16, "2-percentage-point band" = 17, "5-percentage-point band" = 15, "Weighted" = 18),
      breaks = c("Exact", "2-percentage-point band", "5-percentage-point band", "Weighted"),
      labels = c("Exact", "2-percentage-point band", "5-percentage-point band", "Weighted"), name = NULL) +
    guides(shape = guide_legend(nrow = 1)) +
    ax$x("D_sensitivity_x", c(0, filter(d, layer == "sd_summary")$gap_sd_pp), n = 2) + row_y("D_sensitivity_y", FALSE) +
    labs(x = "Residual SD\n(percentage points)", y = NULL) + t() + hide_y +
    # One-row key right-aligned to the panel edge; it may extend left under
    # the neighbouring SD axis, where the shared legend row is otherwise empty.
    figure1_bottom_legend() + theme(legend.text = element_text(size = 4.8),
      legend.key.width = unit(2.2, "mm"), legend.spacing.x = unit(.6, "mm"),
      legend.justification = "right")
  d_inner <- wrap_plots(d_rain, patchwork::free(d_exact, side = "lr"),
    patchwork::free(d_sensitivity, side = "lr"), nrow = 1,
    widths = c(2.75, 1.1, 1.65), guides = "keep")

  # E: show the age-band means as heatmap cell labels, with a local colour key.
  # The older-minus-younger gradient rows stay in the source data but are not drawn.
  e_cells <- filter(e, layer == "age_cell") %>% mutate(
    text_colour = if_else(estimate_pp > 78, "white", pub_colours$text))
  # The title sits above the centred colourbar so the key fits the narrower panel e.
  e_matrix <- ggplot(e_cells, aes(x, y, fill = estimate_pp)) +
    geom_tile(width = .94, height = .86) +
    gap_fill(figure1_colourbar("Programme-derived gap (%)", width = 22, title_position = "top")) +
    geom_text(aes(label = sprintf("%.1f", estimate_pp), colour = text_colour), size = 2.0) +
    scale_colour_identity() +
    scale_x_continuous(breaks = 1:3, limits = c(.5, 3.5),
      labels = c("1–4", "5–9", "10–14"), expand = c(0, 0), minor_breaks = NULL) +
    row_y("E_matrix_y") + labs(x = "Age group (years)", y = "Disease") + t() +
    theme(axis.text.y = element_text(size = 5.1)) +
    figure1_bottom_legend()

  # F: historical deficit residuals and disease-specific within-stratum slopes.
  f_levels <- pub_disease_levels
  f <- f %>% mutate(disease = factor(disease, levels = f_levels))
  f_bins <- filter(f, layer == "density_bin")
  f_slopes <- filter(f, layer == "within_stratum_slope") %>%
    mutate(annotation = sprintf("%.3f (%.3f–%.3f)",
      estimate_pp / 100, lower_pp / 100, upper_pp / 100))
  f_xvalues <- c(f_bins$x - 1, f_bins$x + 1, f$deficit_residual_pp)
  f_yvalues <- c(f_bins$y - 1, f_bins$y + 1, f$lower_pp[f$layer == "slope_band"], f$upper_pp[f$layer == "slope_band"])
  f_xbreaks <- pretty(range(f_xvalues, na.rm = TRUE), n = 5)
  f_ybreaks <- pretty(range(f_yvalues, na.rm = TRUE), n = 4)
  f_slopes <- f_slopes %>% mutate(label_x = min(f_xbreaks) + .06 * diff(range(f_xbreaks)),
    label_y = max(f_ybreaks) - .035 * diff(range(f_ybreaks)))
  f_plot <- ggplot() +
    geom_hline(yintercept = 0, colour = "#D2D5D6", linewidth = .22) +
    geom_vline(xintercept = 0, colour = "#D2D5D6", linewidth = .22) +
    geom_tile(data = f_bins, aes(x, y, fill = n), width = 2, height = 2) +
    geom_ribbon(data = filter(f, layer == "slope_band"), aes(deficit_residual_pp, ymin = lower_pp, ymax = upper_pp),
      fill = pub_colours$orange, alpha = .18) +
    geom_line(data = filter(f, layer == "slope_band"), aes(deficit_residual_pp, gap_residual_pp),
      colour = "#B35426", linewidth = .55) +
    scale_fill_gradient(low = "#DFE7EB", high = "#3F708C", trans = "log10",
      limits = c(1, 1000), breaks = c(1, 10, 100, 1000),
      guide = figure1_colourbar("Count", width = 20, direction = "vertical", title_position = "top")) +
    geom_text(data = f_slopes, aes(label_x, label_y, label = annotation), hjust = 0, vjust = 1,
      colour = "#A84821", size = 1.65) +
    facet_wrap(~disease, nrow = 1, labeller = as_labeller(disease_labels)) +
    ax$x("F_density_x", f_xvalues, n = 5, labels = function(x) ifelse(x %% 20 == 0, x, "")) +
    ax$y("F_density_y", f_yvalues, n = 4) +
    labs(x = "Preceding-decade coverage-deficit residual\n(percentage points)",
      y = "Protection-gap residual\n(percentage points)") + t() +
    figure1_right_legend() +
    theme(strip.text = element_text(size = 5.8), panel.spacing.x = unit(3.5, "mm"),
      plot.caption = element_blank(), plot.margin = margin(1, 2, 5, 2))
  # Each named panel is a native patchwork, freed from unrelated outer alignment.
  # Explicit keep prevents a parent from hoisting any panel's local guides.
  mark <- function(p, letter) p + labs(title = letter) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 3)),
          plot.title.position = "plot")
  a_inner[[1]] <- mark(a_inner[[1]], "a")
  b_inner[[1]] <- mark(b_inner[[1]], "b")
  c_plot <- mark(c_plot, "c")
  d_inner[[1]] <- mark(d_inner[[1]], "d")
  e_matrix <- mark(e_matrix, "e")
  f_plot <- mark(f_plot, "f")
  blocks <- lapply(list(a_inner, b_inner, c_plot, d_inner, e_matrix, f_plot),
    function(p) patchwork::free(p, type = "panel"))
  # Columns 3-6 still sum to 2.25, so rows a-c keep their proportions; the
  # heatmap-only panel e takes 1.65 and panel d gains the released width.
  fig <- wrap_plots(blocks, design = "AABBBB\nAACCCC\nDDDEEE\nFFFFFF",
    widths = c(1.65, 1.65, .60, .55, .55, .55),
    heights = c(1, 1, 1.25, .85), guides = "keep")
  fig_dir <- file.path(out, "figures/main_nc")
  src_dir <- file.path(out, "source_data/figure_redesign")
  dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
  save_publication_figure(fig, "Figure1_NC_redesigned", fig_dir, height_mm = 185, dpi = 600)
  output_data <- list(a = a, b = b, c = c, d = d, e = e, f = f)
  for (i in seq_along(output_data)) write_csv(output_data[[i]], file.path(src_dir, paste0("Figure1_", letters[i], ".csv")))
  invisible(NULL)
}

# Figure 5: age redistribution and contact-weighted transmission potential.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(cowplot)
})
source("scripts/05_figures/publication_figure_style.R")
source("scripts/05_figures/helpers/figure1_layout_helpers.R")
source("scripts/05_figures/prepare_figure5_mechanisms.R")

render_age_mechanism_figures <- function(out = "output/revised", figures = c("5", "S9"),
                                         delta_transform = "identity", fig5_suffix = "") {
  data <- prepare_figure5_mechanisms(out)
  mechanism_diseases <- data$mechanism_diseases
  disease_cols <- pub_evidence_colours
  t <- function(grid = "none") theme_publication(grid = grid, show_titles = TRUE) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(1, 2, 1, 2), axis.title = figure1_axis_title_style(),
      axis.title.x = figure1_axis_title_style(margin(t = 2, b = 8)),
      axis.title.y = figure1_axis_title_style(margin(r = 2)),
      axis.text = element_text(size = 5.4), axis.ticks.length = unit(.7, "mm"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0,
        margin = margin(b = 3)), plot.title.position = "plot",
      strip.text = figure1_facet_title_style(), legend.text = element_text(size = 5.1),
      legend.title = figure1_axis_title_style(), legend.title.position = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
      panel.spacing = unit(3, "mm"))
  mark <- function(p, letter) p + labs(title = letter)
  # Match Figure 2e's horizontal, left-title guide while leaving room for
  # the zero tick, which extends left of each colourbar's first rectangle.
  colourbar <- function(title, width = 24) {
    guide_colourbar(title = title, display = "rectangles", nbin = 64,
      direction = "horizontal", ticks.colour = "white", ticks.linewidth = .3,
      theme = theme(
        legend.key.width = unit(width, "mm"), legend.key.height = unit(1.7, "mm"),
        legend.title.position = "left",
        legend.title = figure1_legend_title_style(margin(r = 10)),
        legend.margin = margin(0, 0, 0, 0)))
  }
  colourbar_vertical <- function(title, height = 32) {
    guide_colourbar(title = title, display = "rectangles", nbin = 64,
      direction = "vertical", ticks.colour = "white", ticks.linewidth = .3,
      theme = theme(
        legend.key.width = unit(2.5, "mm"), legend.key.height = unit(height, "mm"),
        legend.title.position = "left",
        legend.title = element_text(size = 6, face = "bold", colour = pub_colours$text,
          angle = 90, margin = margin(r = 5))))
  }
  year_scale <- scale_x_continuous(breaks = c(2014, 2018, 2023),
    limits = c(2013.5, 2023.5), expand = c(0, 0))
  zero <- geom_hline(yintercept = 0, colour = pub_colours$reference,
    linewidth = .25, linetype = 2)

  # a: paired age-gap phase planes expose absolute gaps, direction and countries.
  a_dat <- data$age_plane %>% mutate(
    disease = factor(as.character(disease), levels = pub_disease_levels),
    mode = factor(mode, levels = c("unwaning", "waning")),
    year_start = 2014L, year_end = 2023L)
  a_pair <- a_dat %>% select(country, disease, mode, younger_gap_pct, older_gap_pct) %>%
    pivot_wider(names_from = mode, values_from = c(younger_gap_pct, older_gap_pct)) %>%
    mutate(identical_modes = near(younger_gap_pct_waning, younger_gap_pct_unwaning) &
      near(older_gap_pct_waning, older_gap_pct_unwaning))
  a_mean <- a_dat %>% group_by(disease, mode) %>%
    summarise(younger_gap_pct = mean(younger_gap_pct),
      older_gap_pct = mean(older_gap_pct), .groups = "drop")
  a_mean_pair <- a_mean %>%
    pivot_wider(names_from = mode, values_from = c(younger_gap_pct, older_gap_pct)) %>%
    mutate(identical_modes = near(younger_gap_pct_waning, younger_gap_pct_unwaning) &
      near(older_gap_pct_waning, older_gap_pct_unwaning))
  a_counts <- a_pair %>% group_by(disease) %>% summarise(
    n_countries = n(),
    n_reversed = sum(older_gap_pct_unwaning < younger_gap_pct_unwaning &
      older_gap_pct_waning > younger_gap_pct_waning), .groups = "drop")
  stopifnot(nrow(a_mean) == 10L, nrow(a_mean_pair) == 5L,
    !anyNA(a_pair), all(a_dat$younger_gap_pct >= 0 & a_dat$younger_gap_pct <= 100),
    all(a_dat$older_gap_pct >= 0 & a_dat$older_gap_pct <= 100))
  a_compare <- ggplot() +
    annotate("polygon", x = c(0, 0, 100), y = c(0, 100, 100),
      fill = "#F3F5F6", colour = NA) +
    geom_abline(slope = 1, intercept = 0, colour = pub_colours$reference,
      linewidth = .3, linetype = "dashed") +
    geom_segment(data = filter(a_pair, !identical_modes),
      aes(x = younger_gap_pct_unwaning, y = older_gap_pct_unwaning,
        xend = younger_gap_pct_waning, yend = older_gap_pct_waning, colour = disease),
      linewidth = .25, alpha = .3) +
    geom_point(data = filter(a_dat, mode == "unwaning"),
      aes(younger_gap_pct, older_gap_pct, colour = disease),
      shape = 21, fill = "white", size = .9, stroke = .3, alpha = .65) +
    geom_point(data = filter(a_dat, mode == "waning"),
      aes(younger_gap_pct, older_gap_pct, colour = disease), size = .65, alpha = .7) +
    geom_segment(data = filter(a_mean_pair, !identical_modes),
      aes(x = younger_gap_pct_unwaning, y = older_gap_pct_unwaning,
        xend = younger_gap_pct_waning, yend = older_gap_pct_waning, colour = disease),
      linewidth = .5) +
    geom_point(data = filter(a_mean, mode == "unwaning"),
      aes(younger_gap_pct, older_gap_pct, colour = disease, shape = mode),
      size = 2.1, stroke = .7, fill = "white") +
    geom_point(data = filter(a_mean, mode == "waning"),
      aes(younger_gap_pct, older_gap_pct, colour = disease, shape = mode), size = 1.45) +
    facet_wrap(~disease, nrow = 1, axes = "all", axis.labels = "margins") +
    scale_colour_manual(values = disease_cols, guide = "none") +
    scale_shape_manual(values = c(unwaning = 21, waning = 16),
      labels = c(unwaning = "No waning", waning = "Waning"), name = NULL) +
    scale_x_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100),
      oob = scales::oob_keep, expand = c(0, 0)) +
    coord_fixed(ylim = c(0, 100), clip = "off") +
    labs(x = "Programme-derived protection gap at ages 0–4 (%)",
      y = "Programme-derived\nprotection gap at ages 10–14 (%)") + t() +
      guides(shape = guide_legend(ncol = 1,
      override.aes = list(colour = pub_colours$text, fill = "white", size = 1.8))) +
    figure1_bottom_legend() +
    theme(legend.position = c(.01, 1), legend.justification = c(0, 1),
      legend.background = element_rect(fill = scales::alpha("white", .82), colour = NA),
      legend.box = "horizontal", legend.direction = "horizontal",
      legend.key.size = unit(2.5, "mm"), legend.spacing.x = unit(2, "mm"),
      # The wider right margin keeps the last facet's centred "100" x tick
      # label inside the canvas; expand = 0 places it on the panel edge.
      legend.margin = margin(1, 1, 1, 1), plot.margin = margin(2, 8, 2, 2),
      strip.text = element_text(margin = margin(b = 10), size = 5.8, face = "bold"))
  pa <- mark(a_compare, "a")

  # b: true single-age cohort inputs before aggregation to the contact bands.
  b_abs <- data$surface %>% mutate(view = factor(if_else(mode == "waning",
    "Waning", "No waning"), levels = c("No waning", "Waning")))
  b_delta <- data$surface_delta %>% mutate(view = "Waning − no waning",
    delta_pp = 100 * delta)
  b_dat <- bind_rows(b_abs %>% transmute(disease, year, age,
      view = as.character(view), susceptibility, delta_pp = NA_real_, n_profile_fallback),
    b_delta %>% transmute(disease, year, age, view,
      susceptibility = NA_real_, delta_pp, n_profile_fallback)) %>%
    mutate(transmission_blocking = unname(data$blocking[as.character(disease)]),
      cohort_input = if_else(disease == "Mumps",
        "mumps_all_mcv_proxy_cohort_level.csv", "cohort_level_all_ages.csv"))

  # b: waning − no waning difference retains the year × age surface.
  b_delta_limit <- max(pretty(c(0, b_delta$delta_pp), n = 3))
  pb_delta <- ggplot(b_delta, aes(year, age, fill = delta_pp)) +
    geom_tile(width = .98, height = .97) + facet_wrap(~disease, nrow = 1) +
    scale_fill_gradientn(colours = pub_gradients$sequential_blue_deep,
      limits = c(0, b_delta_limit), breaks = seq(0, b_delta_limit, by = 10),
      transform = delta_transform, oob = scales::squish,
      name = "Waning-induced\nsusceptibility increase",
      guide = colourbar_vertical("Waning-induced\nsusceptibility increase", height = 30)) +
    year_scale + scale_y_reverse(breaks = c(1, 5, 10, 14), limits = c(14.5, .5),
      expand = c(0, 0)) +
    labs(x = "Calendar year", y = "Age (years)") + t() +
    theme(legend.position = "right", legend.box = "vertical", legend.direction = "vertical",
      strip.text = element_text(face = "bold"))
  pb <- mark(pb_delta, "b")
  # pb <- pb + theme(plot.margin = margin(0, 0, 0, 0))

  # c: vaccine-action routing. The same waning-induced gap increase (a) enters
  # the contact operator on orthogonal channels: recipient acquisition for
  # pertussis and mumps, source infectiousness for diphtheria.
  c_country <- data$channels
  c_mean <- c_country %>% group_by(disease, age_band) %>%
    summarise(across(c(delta_acquisition_pp, delta_infectiousness_pp), mean),
      .groups = "drop")
  channel_labels <- c(delta_acquisition_pp = "Recipient acquisition\nsusceptibility",
    delta_infectiousness_pp = "Source onward\ninfectiousness")
  c_long <- c_country %>%
    pivot_longer(names(channel_labels), names_to = "channel", values_to = "delta_pp") %>%
    mutate(channel = factor(channel_labels[channel], levels = channel_labels))
  c_bar <- c_long %>% group_by(disease, channel, age_band) %>%
    summarise(mean = mean(delta_pp), lo = min(delta_pp), hi = max(delta_pp),
      .groups = "drop")
  c_zero <- c_bar %>% group_by(disease, channel) %>%
    summarise(zero = all(abs(hi) < 1e-8), .groups = "drop") %>% filter(zero) %>%
    mutate(label = "No change")
  c_off_axis <- max(pmin(c_country$delta_acquisition_pp, c_country$delta_infectiousness_pp))
  stopifnot(nrow(c_bar) == 18L, c_off_axis < 1e-8, nrow(c_zero) == 3L)
  c_limit <- max(pretty(c(0, c_bar$hi), n = 3))
  pc_route <- ggplot(c_bar, aes(mean, age_band)) +
    geom_col(aes(fill = disease, alpha = age_band), width = .72, orientation = "y") +
    geom_vline(xintercept = 0, colour = pub_colours$light_grey, linewidth = .3) +
    geom_errorbar(data = filter(c_bar, hi > 1e-8),
      aes(xmin = lo, xmax = hi), width = .25, linewidth = .3,
      colour = pub_colours$text, orientation = "y") +
    geom_text(data = c_zero, aes(x = c_limit / 2, y = "5–9", label = label),
      size = 1.75, colour = pub_colours$muted_text, fontface = "italic") +
    facet_grid(disease ~ channel, switch = "y") +
    scale_fill_manual(values = disease_cols, guide = "none") +
    scale_alpha_manual(values = c(`0–4` = .45, `5–9` = .7, `10–14` = 1), guide = "none") +
    scale_x_continuous(limits = c(0, c_limit), breaks = seq(0, c_limit, by = 10),
      expand = expansion(mult = c(0, .04))) +
    scale_y_discrete(limits = rev(levels(c_bar$age_band))) +
    labs(x = "Waning-induced change", y = NULL) + t() +
    theme(strip.placement = "outside",
      strip.text.x = element_text(face = "bold", lineheight = .9),
      strip.text.y.left = element_text(face = "bold", angle = 90),
      panel.spacing.x = unit(2.5, "mm"),
      panel.spacing.y = unit(1.5, "mm"))
  pc4 <- mark(pc_route, "c")
  # pc4 <- pc4 + theme(plot.margin = margin(0, 0, 0, 0))

  # Figure 5d: show active recipient rows and all sixteen source-age bands.
  cdat <- data$kernel_delta %>% mutate(
    disease = factor(disease, levels = mechanism_diseases),
    source_age = factor(source_age, levels = unique(source_age)))
  # Source protection (diphtheria) can change adult recipient rows as well.
  child_ages <- if (any(vapply(data$blocking, function(x) x == 0, logical(1)))) {
    unique(cdat$recipient_age)
  } else c("0–4", "5–9", "10–14")
  c_active <- cdat %>% filter(recipient_age %in% child_ages) %>%
    mutate(recipient_age = factor(recipient_age, levels = rev(child_ages)))
  c_recipient_dat <- data$kernel_recipient %>% mutate(
    disease = factor(disease, levels = mechanism_diseases),
    recipient_age = factor(recipient_age, levels = child_ages)) %>%
    filter(recipient_age %in% child_ages)
  c_source_dat <- data$kernel_source %>% mutate(
    disease = factor(disease, levels = mechanism_diseases),
    source_age = factor(source_age, levels = levels(cdat$source_age)))
  # Legend: pretty breaks are computed first and also fix the colour limits,
  # so the bar ends on a labelled tick. A square-root transform keeps the many
  # small entries (median non-zero cell ~0.006 against a ~0.19 maximum)
  # distinguishable instead of collapsing them onto the first colour. The
  # caption states the square-root scale; the uneven tick spacing shows it.
  c_breaks <- pretty(c(0, max(abs(cdat$value))), n = 3)
  c_limits <- range(c_breaks)
  stopifnot(min(c_limits) == 0, all(cdat$value >= 0), all(cdat$value <= max(c_limits)))
  c_title <- "Waning-induced\ntransmission-weight increase"
  # Show only the five sentinel age bands on the heatmap axes (and the
  # matching S9 source-age marginal) to keep tick labels readable at 1:1.
  source_ticks <- levels(c_active$source_age)[c(1, 5, 9, 13, 16)]
  recipient_ticks <- if (length(child_ages) >= 16) child_ages[c(1, 5, 9, 13, 16)] else child_ages
  recipient_limit <- max(pretty(c(0, c_recipient_dat$value), n = 3))
  source_limit <- max(pretty(c(0, c_source_dat$value), n = 3))
  c_heat <- ggplot(c_active, aes(source_age, recipient_age, fill = value)) +
    geom_tile(colour = "white", linewidth = .12) + facet_wrap(~disease, nrow = 1) +
    scale_fill_gradientn(colours = pub_gradients$sequential_red,
      limits = c_limits, breaks = c_breaks, transform = "sqrt",
      labels = label_number(accuracy = .01),
      name = c_title, guide = colourbar_vertical(c_title, height = 30)) +
    scale_x_discrete(breaks = source_ticks) +
    # Sparse y-axis labels keep the age labels readable.
    scale_y_discrete(breaks = recipient_ticks) +
    labs(x = "Source age (years)", y = "Recipient age (years)") + t() +
    figure1_bottom_legend() + theme(legend.position = "right",
      legend.box = "vertical", legend.direction = "vertical") +
    theme(strip.text = element_text(face = "bold"),
      plot.margin = margin(0, 0, 0, 0))
  c_recipient <- ggplot(c_recipient_dat, aes(recipient_age, value)) +
    geom_col(fill = pub_colours$orange, width = .65) +
    facet_wrap(~disease, nrow = 1) +
    scale_y_continuous(limits = c(0, recipient_limit),
      breaks = pretty(c(0, recipient_limit), n = 3), expand = c(0, 0)) +
    scale_x_discrete(breaks = if (length(child_ages) > 3) child_ages[c(1, 6, 11, 16)] else child_ages) +
    labs(x = "Recipient age group", y = NULL) + t("y") +
    theme(strip.text = element_blank(), axis.text.x = element_text(size = 4.6,
      angle = if (length(child_ages) > 3) 35 else 0,
      hjust = if (length(child_ages) > 3) 1 else .5))
  c_source <- ggplot(c_source_dat, aes(source_age, value)) +
    geom_col(fill = pub_colours$orange, width = .72) +
    facet_wrap(~disease, nrow = 1) +
    scale_x_discrete(breaks = source_ticks) +
    scale_y_continuous(limits = c(0, source_limit),
      breaks = pretty(c(0, source_limit), n = 3), expand = c(0, 0)) +
    labs(x = "Source age group", y = NULL) + t("y") +
    theme(strip.text = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1, size = 4.6))
  pc <- wrap_plots(mark(c_heat, "c"), c_recipient, c_source,
    ncol = 1, heights = c(2.1, 1, 1), guides = "keep")

  # Figure 5e: signed redistribution of normalized dominant-mode shares as
  # Fig. 2d-style point estimates with 95% paired country-bootstrap intervals.
  d_dat <- data$contribution_delta %>%
    mutate(disease = factor(disease, levels = mechanism_diseases),
      age_band = factor(age_band, levels = rev(unique(age_band))))
  d_range <- range(c(0, d_dat$delta_contribution_lower_95,
    d_dat$delta_contribution_upper_95))
  d_breaks <- pretty(d_range, n = 3)
  d_breaks <- d_breaks[d_breaks >= d_range[1] & d_breaks <= d_range[2]]
  # Create dummy data for legend
  legend_dummy <- data.frame(x = 0, y = 0)
  pd <- ggplot(d_dat, aes(delta_contribution_pp, age_band, colour = disease)) +
    geom_vline(xintercept = 0, colour = pub_colours$reference,
      linewidth = .25, linetype = 2) +
    geom_segment(aes(x = delta_contribution_lower_95, xend = delta_contribution_upper_95,
      yend = age_band, linetype = "95% CI"), linewidth = .5) +
    geom_point(aes(shape = "Estimate"), size = 1.5) + facet_wrap(~disease, nrow = 1) +
    scale_colour_manual(values = disease_cols, guide = "none") +
    scale_shape_manual(values = c("Estimate" = 16), name = NULL) +
    scale_linetype_manual(values = c("95% CI" = 1), name = NULL) +
    scale_x_continuous(limits = d_range, breaks = d_breaks,
      labels = function(x) paste0(ifelse(x > 0, "+", ""), x),
      expand = expansion(mult = .08)) +
    scale_y_discrete() +
    labs(x = "Change in normalised dominant-mode\nage contribution",
      y = "Age (years)") + t("x") +
    guides(shape = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(colour = pub_colours$text, fill = "white")),
      linetype = guide_legend(nrow = 1, byrow = TRUE)) +
    figure1_bottom_legend() +
    theme(legend.position = "bottom", legend.box = "horizontal", legend.direction = "horizontal",
      legend.key.size = unit(2, "mm"))
  pd <- mark(pd, "d")

  # Figure 5f: country responses ranked by the pertussis contrast, as Fig. 2d-style
  # point estimates with 95% within-country year-bootstrap intervals.
  # Use same country name mapping as Figure 4d for consistency
  country_name_map_fig4 <- c(
    IND = "India", CHN = "China", USA = "USA", IDN = "Indonesia",
    PAK = "Pakistan", BRA = "Brazil", RUS = "Russia", DEU = "Germany",
    THA = "Thailand", GBR = "UK"
  )
  country_order <- data$country %>% filter(disease == "Pertussis") %>%
    arrange(delta_Re) %>% pull(country)
  e_dat <- data$country %>% filter(disease %in% mechanism_diseases) %>%
    mutate(country = factor(country, levels = country_order),
           country_display = country_name_map_fig4[as.character(country)],
           country_display = factor(country_display, 
                                    levels = country_name_map_fig4[country_order]),
           disease = factor(disease, levels = mechanism_diseases))
  pe <- ggplot(e_dat, aes(delta_Re, country_display, colour = disease)) +
    geom_vline(xintercept = 0, colour = pub_colours$reference,
      linewidth = .25, linetype = 2) +
    geom_segment(aes(x = delta_Re_lower_95, xend = delta_Re_upper_95, yend = country_display,
      linetype = "95% CI"), linewidth = .5) +
    geom_point(aes(shape = "Estimate"), size = 1.5) + facet_wrap(~disease, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = disease_cols, guide = "none") +
    scale_shape_manual(values = c("Estimate" = 16), name = NULL) +
    scale_linetype_manual(values = c("95% CI" = 1), name = NULL) +
    scale_x_continuous(breaks = function(limits) pretty(limits, n = 3),
      labels = label_number(),
      expand = expansion(mult = .08)) +
    labs(x = "Change in initial Rₑ", y = "Country") + t("x") +
    guides(shape = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(colour = pub_colours$text, fill = "white")),
      linetype = guide_legend(nrow = 1, byrow = TRUE)) +
    figure1_bottom_legend() +
    theme(legend.position = "bottom", legend.box = "horizontal", legend.direction = "horizontal",
      legend.key.size = unit(2, "mm"),
      panel.spacing.x = unit(3, "mm"))
  pe <- mark(pe, "e")

  fig_dir <- file.path(out, "figures/main_nc")
  src_dir <- file.path(out, "source_data/figure_redesign")
  dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
  free_panel <- function(p) patchwork::free(p, type = "panel")
  # Figure 5: panels a–c show age redistribution and vaccine-action
  # routing; panels d–f show contact-weighted transmission. The former
  # sensitivity grid is supplied as a supplementary table and omitted here.
  # Compose rows 2 and 3 as a standalone patchwork: b/d use the wider
  # left column, and c spans the narrower right column.
  row_bc_dc <- wrap_plots(
    B = mark(pb, "b"), C = free_panel(mark(pc4, "c")),
    D = mark(c_heat, "d"), design = "BC\nDC",
    widths = c(1.7, 1.1), guides = "keep")
  # Cowplot combines the three rows; e/f share the final row equally.
  row_a <- free_panel(mark(pa, "a"))
  row_e_f <- plot_grid(mark(pd, "e"), mark(pe, "f"), ncol = 2,
    rel_widths = c(1, 1), align = "none")
  fig_combined <- plot_grid(row_a, row_bc_dc, row_e_f,
    ncol = 1, rel_heights = c(1.1, 1.65, 1.4), align = "none")
  fig_s9 <- wrap_plots(mark(c_recipient, "A"), mark(c_source, "B"),
    ncol = 1, heights = c(1, 1.1))
  if ("5" %in% figures)
    save_publication_figure(fig_combined, paste0("Figure5_NC_redesigned", fig5_suffix), fig_dir,
      height_mm = 190, dpi = 600)
  if ("S9" %in% figures)
    save_publication_figure(fig_s9, "FigureS9_revised_contact_marginals",
      file.path(out, "figures/supplementary"), height_mm = 95, dpi = 600)
  a_source <- bind_rows(
    a_dat %>% mutate(layer = "country_period_mean"),
    a_mean %>% mutate(layer = "equal_country_mean", year_start = 2014L, year_end = 2023L)) %>%
    mutate(gradient_pp = older_gap_pct - younger_gap_pct) %>%
    left_join(select(a_counts, disease, n_countries, n_reversed), by = "disease")
  if ("5" %in% figures) {
    write_csv(data$gradient, file.path(src_dir, "Figure5_a_annual_gradients.csv"))
    b_source <- b_dat %>% filter(view == "Waning − no waning") %>% select(-susceptibility)
    c_source <- c_country %>% mutate(layer = "country_period_mean") %>%
      bind_rows(c_mean %>% mutate(layer = "equal_country_mean")) %>%
      mutate(transmission_blocking = unname(data$blocking[as.character(disease)]))
    displayed <- list(a = a_source, b = b_source, c = c_source,
      d = data$kernel_delta %>% mutate(layer = "matrix"), e = d_dat, f = e_dat)
    for (i in names(displayed))
      write_csv(displayed[[i]], file.path(src_dir, paste0("Figure5_", i, ".csv")))
    # The former sensitivity grid is retained as a supplementary table.
    write_csv(data$sensitivity, file.path(src_dir, "Figure5_sensitivity_grid.csv"))
    jsonlite::write_json(data$bootstrap,
      file.path(out, "age_seirv", "figure5_bootstrap_manifest.json"),
      auto_unbox = TRUE, pretty = TRUE)
  }
  if ("S9" %in% figures) {
    write_csv(data$kernel_recipient, file.path(src_dir, "FigureS9_a.csv"))
    write_csv(data$kernel_source, file.path(src_dir, "FigureS9_b.csv"))
  }
  invisible(list(figure5 = fig_combined, supplementary9 = fig_s9))
}

render_figure5_mechanisms <- function(out = "output/revised") {
  render_age_mechanism_figures(out)
}

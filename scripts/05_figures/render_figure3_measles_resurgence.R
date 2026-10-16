# Figure 3: EU/EEA measles resurgence in the ten highest-burden countries.
# Monthly panels use the WHO/UNICEF archive; ECDC supplies the case-based age
# and vaccination layers. The logic runs from burden (a–d),
# through where the cohort gaps sit (e) and whether they rank incidence beyond
# current coverage (f), to whether they locate the ages of cases within
# countries (g). Supplementary Fig. S12 gives the added-variable view of the
# 1–4-year gap (A) and where national cohorts stop describing the outbreak (B).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
})

source("scripts/05_figures/publication_figure_style.R")
source("scripts/05_figures/helpers/figure1_layout_helpers.R")

render_figure3_measles_resurgence <- function(out = "output/revised") {
  data_dir <- file.path(out, "measles_resurgence_full")
  figure_dir <- file.path(out, "figures/main_nc")
  src_dir <- file.path(out, "source_data/figure_redesign")
  dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
  read_fig <- function(name, required) {
    d <- read_csv(file.path(data_dir, name), show_col_types = FALSE)
    stopifnot(nrow(d) > 0, all(required %in% names(d)))
    d
  }
  hot <- read_fig("measles_hotspot_country_summary.csv", c("country", "hotspot", "rank_region", "cases_2024_26", "wuenic_mcv1_2025", "wuenic_mcv2_2025"))
  monthly <- read_fig("measles_hotspot_monthly.csv", c("country", "year", "month", "cases", "incidence_per_100k"))
  monthly_all <- read_fig("measles_country_monthly_2016_2026.csv",
    c("country", "year", "month", "cases", "incidence_per_100k"))
  region_monthly <- read_fig("measles_region_monthly.csv", c("year", "month", "region_cases", "region_population"))
  age_gaps <- read_fig("measles_hotspot_age_gaps.csv", c("country", "band", "gap", "hotspot", "rank_region"))
  adj <- read_fig("measles_coverage_adjusted_gap.csv", c("term", "rho", "lower", "upper"))
  adj_pts <- read_fig("measles_coverage_adjusted_points.csv", c("country", "hotspot", "resid_gap_1_4", "resid_incidence", "annual_incidence"))
  fe <- read_fig("measles_age_alignment_fe.csv", c("model", "irr", "lower", "upper", "n_countries", "n_country_years", "n_cases"))
  fe_pts <- read_fig("measles_age_alignment_fe_points.csv", c("country", "top10", "gap_resid", "log2_oe", "expected"))
  screen <- read_fig("measles_screening_coverage.csv", c("country", "outcome_year", "hotspot", "national_mcv1", "screening_equivalent_coverage", "screening_equivalent_lower", "screening_equivalent_upper", "known_status"))
  screen_sum <- read_fig("measles_screening_coverage_summary.csv", c("group", "weighted_deficit", "lower", "upper"))

  top_order <- hot %>% filter(hotspot) %>% arrange(rank_region) %>% pull(country)
  top_cols <- c(`TRUE` = pub_colours$purple, `FALSE` = pub_colours$teal)
  top_labels <- c(`TRUE` = "Top 10 EU/EEA", `FALSE` = "Other EU/EEA")
  
  # Extract single-word country names from full names
  # Strategy: use ISO code to look up standardized single-word names
  country_name_map <- c(
    AUT = "Austria",
    BEL = "Belgium",
    BGR = "Bulgaria",
    DEU = "Germany",
    ESP = "Spain",
    FRA = "France",
    ITA = "Italy",
    NLD = "Netherlands",
    POL = "Poland",
    ROU = "Romania"
  )
  
  # Alternative: extract first word from country_name column if available
  extract_single_word_name <- function(iso_code, full_name) {
    # Predefined mapping for consistency
    if (iso_code %in% names(country_name_map)) {
      return(country_name_map[iso_code])
    }
    # Fallback: extract first meaningful word
    # Remove common prefixes/suffixes and parenthetical notes
    clean_name <- gsub("\\s*\\(.*?\\)\\s*", "", full_name)
    # Take first word
    first_word <- strsplit(clean_name, "\\s+")[[1]][1]
    return(first_word)
  }
  
  # Apply to top order countries
  top_order_display <- sapply(top_order, function(code) {
    full_name <- hot %>% filter(country == code) %>% pull(country_name) %>% first()
    extract_single_word_name(code, full_name)
  })
  names(top_order_display) <- top_order
  last_complete <- function(d) filter(d, year < 2026 | (year == 2026 & month <= 6))

  base_theme <- theme_publication(grid = "y", show_titles = TRUE) + theme(
    plot.margin = margin(1, 2, 1, 2),
    axis.title = figure1_axis_title_style(),
    axis.title.x = figure1_axis_title_style(margin(t = 2)),
    axis.title.y = figure1_axis_title_style(margin(r = 2)),
    axis.ticks.length = unit(.7, "mm"),
    legend.position = "bottom", legend.title.position = "left",
    legend.title = figure1_axis_title_style(),
    legend.text = element_text(size = pub_text$legend_text),
    legend.key.width = unit(4, "mm"), legend.key.height = unit(2.6, "mm"),
    legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    panel.grid.major.y = element_line(colour = pub_colours$grid, linewidth = .18),
    panel.grid.major.x = element_blank()
  )
  tag <- function(p, letter) p + labs(title = letter) + theme_panel_letter()

  # a. EU/EEA monthly incidence.  The country archive is collapsed to the
  # three displayed population-weighted series; individual country traces are
  # reserved for panels b–d and are not overlaid here.
  region_trend <- region_monthly %>% last_complete() %>%
    mutate(time = year + (month - 1) / 12, incidence_per_100k = 100000 * region_cases / region_population,
      series = "EU/EEA total")
  country_trend <- monthly_all %>%
    select(country, year, month, cases, population_2025) %>%
    left_join(select(hot, country, hotspot), by = "country") %>%
    mutate(series = if_else(coalesce(hotspot, FALSE), "Top 10 EU/EEA", "Other EU/EEA")) %>%
    group_by(series, year, month) %>%
    summarise(cases = sum(cases, na.rm = TRUE), population = sum(population_2025, na.rm = TRUE), .groups = "drop") %>%
    mutate(time = year + (month - 1) / 12, incidence_per_100k = 100000 * cases / population) %>%
    last_complete()
  a_breaks <- seq(2016, 2026, 2)
  a_dat <- bind_rows(region_trend, country_trend)
  a_y <- pretty(c(0, a_dat$incidence_per_100k), n = 5)
  panel_a <- ggplot(a_dat, aes(time, incidence_per_100k, colour = series)) +
    geom_line(aes(group = series), linewidth = .55) +
    scale_x_continuous(breaks = a_breaks, limits = c(2016, 2026.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = a_y, limits = range(a_y), expand = c(0, 0)) +
    labs(x = "Calendar month", y = "EU/EEA measles incidence\nper 100,000 population") +
    scale_colour_manual(values = c("EU/EEA total" = pub_colours$reference,
      "Top 10 EU/EEA" = top_cols[["TRUE"]], "Other EU/EEA" = top_cols[["FALSE"]]),
      name = "Countries") +
    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = .8, alpha = 1))) +
    base_theme

  # b–c. Burden ranking with current national coverage: coverage of the top ten is
  # mostly high, so annual coverage alone does not identify them.
  b_dat <- hot %>% filter(hotspot) %>% 
    mutate(country = factor(country, levels = rev(top_order)),
           country_display = country_name_map[as.character(country)],
           country_display = factor(country_display, levels = rev(top_order_display)))
  b_cases <- ggplot(b_dat, aes(cases_2024_26, country_display)) +
    # Bars start at the axis floor because a log axis has no zero.
    geom_segment(aes(x = 100, xend = cases_2024_26, yend = country_display), colour = pub_colours$purple, linewidth = 2.4) +
    geom_text(aes(label = label_number(big.mark = ",")(cases_2024_26)), hjust = -.1, size = 2.0, colour = pub_colours$text) +
    scale_x_log10(limits = c(100, 3e5), breaks = c(1e2, 1e3, 1e4, 1e5), labels = label_number(big.mark = ","), expand = c(0, 0)) +
    labs(x = "Reported measles cases,\nJan 2024–Jun 2026", y = NULL) +
    base_theme + theme(panel.grid.major.y = element_blank(),
      plot.margin = margin(1, 2, 1, 2))
  b_cov_dat <- b_dat %>% select(country_display, MCV1 = wuenic_mcv1_2025, MCV2 = wuenic_mcv2_2025) %>%
    pivot_longer(c(MCV1, MCV2), names_to = "dose", values_to = "coverage")
  b_cov <- ggplot(b_cov_dat, aes(coverage * 100, country_display)) +
    geom_vline(xintercept = 95, linetype = 2, colour = pub_colours$reference, linewidth = .25) +
    geom_line(aes(group = country_display), colour = pub_colours$light_grey, linewidth = .5) +
    geom_point(aes(colour = dose), size = 1.25) +
    scale_colour_manual(values = c(MCV1 = pub_colours$blue, MCV2 = pub_colours$orange), name = "2025 WUENIC vaccine dose") +
    scale_x_continuous(limits = c(60, 100), breaks = c(60, 80, 95), expand = c(0, 0)) +
    labs(x = "WUENIC coverage (%)", y = NULL) +
    base_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(), plot.margin = margin(1, 2, 1, 2))
  # d. Monthly incidence of the same ten countries, standardised within country.
  # Missing country-months stay NA and are shown in grey.
  c_dat <- tidyr::crossing(country = top_order, year = 2016:2026, month = 1:12) %>% last_complete() %>%
    left_join(monthly %>% select(country, year, month, cases, incidence_per_100k), by = c("country", "year", "month")) %>%
    mutate(time = year + (month - 1) / 12, country = factor(country, levels = rev(top_order))) %>%
    group_by(country) %>%
    mutate(z = (incidence_per_100k - mean(incidence_per_100k, na.rm = TRUE)) / sd(incidence_per_100k, na.rm = TRUE)) %>%
    ungroup()
  panel_c <- ggplot(c_dat, aes(time, country, fill = z)) +
    geom_tile(width = 1 / 12, height = .88) +
    scale_fill_gradientn(colours = rev(pub_gradients$diverging_full), limits = c(-6, 6), breaks = c(-6, -3, 0, 3, 6),
      labels = c("-6σ", "-3σ", "mean", "+3σ", "+6σ"), name = "Z-score", oob = squish, na.value = "#BDBDBD") +
    scale_x_continuous(breaks = a_breaks, limits = c(2016 - 1 / 24, 2026.5 - 1 / 24), expand = c(0, 0)) +
    scale_y_discrete(labels = NULL) +
    labs(x = "Calendar month", y = NULL) +
    guides(fill = figure1_colourbar("Z-score", width = 30, height = 1.6)) +
    base_theme + theme(axis.line.y = element_line(linewidth = .26, colour = pub_colours$text),
      axis.ticks.y = element_line(linewidth = .22, colour = pub_colours$text), panel.grid.major.y = element_blank(),
      plot.margin = margin(1, 2, 1, 2))

  # e. Where the cohort gap sits: 2023 reconstruction (exposure fixed before the
  # outcome window) by age band, with the EU/EEA median for reference.
  band_levels <- c("1_4", "5_9", "10_14", "15_19", "20_29", "30_39")
  band_labels <- c("1–4", "5–9", "10–14", "15–19", "20–29", "30–39")
  d_ref <- age_gaps %>% group_by(band) %>% summarise(gap = median(gap, na.rm = TRUE), .groups = "drop") %>% mutate(country = "EU/EEA\nmedian")
  d_dat <- bind_rows(age_gaps %>% filter(hotspot) %>% select(country, band, gap), d_ref) %>%
    mutate(country_display = ifelse(country == "EU/EEA\nmedian", "EU/EEA\nmedian", country_name_map[as.character(country)]),
           country = factor(country, levels = rev(c(top_order, "EU/EEA\nmedian"))),
           country_display = factor(country_display, levels = rev(c(top_order_display, "EU/EEA\nmedian"))),
           band = factor(band, levels = band_levels, labels = band_labels))
  panel_d <- ggplot(d_dat, aes(band, country_display, fill = gap * 100)) +
    geom_tile(colour = "white", linewidth = .35) +
    geom_text(aes(label = ifelse(is.na(gap), "–", sprintf("%.0f", gap * 100)), colour = !is.na(gap) & gap > .22), size = 2.0) +
    geom_hline(yintercept = 1.5, colour = pub_colours$text, linewidth = .3) +
    scale_fill_gradientn(colours = pub_gradients$sequential_red, limits = c(0, 40), oob = squish, na.value = "#EDEDED",
      name = "Reconstructed cohort gap (%)") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = pub_colours$text), guide = "none") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Age band in 2023 (years)", y = NULL) +
    guides(fill = figure1_colourbar("Reconstructed\ncohort gap (%)", width = 22, height = 1.6)) +
    base_theme + theme(legend.location = "plot", plot.margin = margin(1, 5, 1, 2))

  # f. Across all EU/EEA countries: crude correlation of 2024–2026 incidence with
  # 2023 coverage (shaded), and rank partial correlations with each 2023 cohort
  # gap conditional on 2023 MCV1 and MCV2.
  term_labels <- c(current_mcv1 = "MCV1 2023 (crude)", current_mcv2 = "MCV2 2023 (crude)",
    gap_1_4 = "Gap 1–4 y", gap_5_9 = "Gap 5–9 y", gap_10_14 = "Gap 10–14 y",
    gap_15_19 = "Gap 15–19 y", gap_20_29 = "Gap 20–29 y", gap_30_39 = "Gap 30–39 y")
  e_dat <- adj %>% filter(term %in% names(term_labels)) %>%
    mutate(crude = term %in% c("current_mcv1", "current_mcv2"),
      term = factor(term, levels = rev(names(term_labels)), labels = rev(term_labels)))
  n_e <- max(adj$n_countries)
  panel_e <- ggplot(e_dat, aes(rho, term, colour = crude)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = length(term_labels) - 1.5, ymax = length(term_labels) + .5, fill = "#F2F2F2") +
    geom_vline(xintercept = 0, linetype = 2, colour = pub_colours$reference, linewidth = .3) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .25, linewidth = .42) +
    geom_point(size = 1.35) +
    scale_colour_manual(values = c(`TRUE` = pub_colours$reference, `FALSE` = pub_colours$blue), guide = "none") +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -.5, 0, .5, 1), expand = c(0, 0)) +
    labs(x = sprintf("Correlation with 2024–2026 incidence\n(ρ; n = %d)", n_e), y = NULL) +
    base_theme + theme(panel.grid.major.y = element_blank(), plot.margin = margin(1, 5, 1, 2))

  # S12A. Added-variable view of the 1–4-year gap in f, with the top ten
  # labelled. ρ and its interval are the country-bootstrap values printed in f;
  # the top ten are chosen on the outcome, so no within-group fits are drawn.
  s12a_dat <- adj_pts %>% filter(is.finite(resid_gap_1_4)) %>% mutate(hotspot = factor(hotspot, levels = c(FALSE, TRUE)))
  s12a_rho <- adj %>% filter(term == "gap_1_4")
  stopifnot(nrow(s12a_rho) == 1,
    abs(cor(s12a_dat$resid_gap_1_4, s12a_dat$resid_incidence) - s12a_rho$rho) < 1e-6)
  s12a_x_breaks <- pretty(s12a_dat$resid_gap_1_4, n = 4)
  s12a_y_breaks <- pretty(s12a_dat$resid_incidence, n = 4)
  s12a <- ggplot(s12a_dat, aes(resid_gap_1_4, resid_incidence)) +
    geom_hline(yintercept = 0, colour = pub_colours$light_grey, linewidth = .25) +
    geom_vline(xintercept = 0, colour = pub_colours$light_grey, linewidth = .25) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = pub_colours$reference,
      fill = pub_colours$reference, linewidth = .45, alpha = .12) +
    geom_point(aes(colour = hotspot, size = annual_incidence), alpha = .9) +
    geom_text_repel(data = filter(s12a_dat, hotspot == "TRUE"), aes(label = country), size = 2.0, seed = 7,
      box.padding = .25, min.segment.length = 0, segment.size = .2, colour = top_cols[["TRUE"]]) +
    annotate("text", x = min(s12a_x_breaks), y = max(s12a_y_breaks), hjust = -.05, vjust = 1.2, size = 2.1,
      colour = pub_colours$text, label = sprintf("All %d countries:\nρ = %.2f (%.2f–%.2f)",
        s12a_rho$n_countries, s12a_rho$rho, s12a_rho$lower, s12a_rho$upper)) +
    scale_colour_manual(values = top_cols, labels = top_labels, name = "Country group") +
    scale_size_continuous(range = c(.6, 3.2), breaks = c(.1, 1, 10), labels = c("0.1", "1", "10"),
      name = "Annual incidence per 100,000 population") +
    scale_x_continuous(breaks = s12a_x_breaks, limits = range(s12a_x_breaks), expand = c(0, 0)) +
    scale_y_continuous(breaks = s12a_y_breaks, limits = range(s12a_y_breaks), expand = c(0, 0)) +
    labs(x = "Age 1–4 gap adjusted for MCV1 and MCV2\n(rank residual)",
      y = "2024–2026 incidence\n(rank residual)") +
    guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(size = 1.5)),
      size = guide_legend(order = 2, nrow = 1, override.aes = list(colour = pub_colours$text))) +
    base_theme + theme(legend.position = c(.99, .01), legend.justification = c(1, 0),
      legend.direction = "vertical", legend.box = "vertical",
      legend.title.position = "top",
      legend.key.height = unit(3.5, "mm"), legend.spacing.y = unit(1.2, "mm"),
      legend.background = element_rect(fill = scales::alpha("white", .82), colour = NA),
      legend.margin = margin(1.5, 1.5, 1.5, 1.5), plot.margin = margin(1, 5, 1, 2))

  # g. Within-country age alignment: case ages (ECDC 2023–25) against the
  # preceding-year gap for the same birth cohorts (age a in year t maps to
  # age a-1 in t-1); the primary fit uses closed 1–29-year bands.
  fe_main <- fe %>% filter(model == "eu_eea_band_fe")
  fe_top <- fe %>% filter(model == "top10_band_fe")
  g_dat <- fe_pts %>% mutate(top10 = factor(top10, levels = c(FALSE, TRUE)))
  g_x_breaks <- pretty(g_dat$gap_resid * 100, n = 4)
  g_y_breaks <- pretty(g_dat$log2_oe, n = 4)
  g_x <- seq(min(g_dat$gap_resid * 100), max(g_dat$gap_resid * 100), length.out = 40)
  g_model_labels <- c(
    "All EU/EEA" = sprintf("All EU/EEA (n = %d):\nIRR %.2f (%.2f–%.2f)", fe_main$n_countries, fe_main$irr, fe_main$lower, fe_main$upper),
    "Top 10" = sprintf("Top 10 (n = %d):\nIRR %.2f (%.2f–%.2f)", fe_top$n_countries, fe_top$irr, fe_top$lower, fe_top$upper))
  g_line <- bind_rows(
    tibble(model = "All EU/EEA", x = g_x, y = log2(fe_main$irr) / 10 * g_x,
      ymin = log2(fe_main$lower) / 10 * g_x, ymax = log2(fe_main$upper) / 10 * g_x),
    tibble(model = "Top 10", x = g_x, y = log2(fe_top$irr) / 10 * g_x,
    ymin = log2(fe_top$lower) / 10 * g_x, ymax = log2(fe_top$upper) / 10 * g_x)
  ) %>% mutate(model = factor(model, levels = names(g_model_labels)),
    lo = pmin(ymin, ymax), hi = pmax(ymin, ymax))
  g_big <- g_dat %>% filter(top10 == "TRUE") %>% group_by(country) %>% slice_max(expected, n = 1) %>% ungroup()
  panel_g <- ggplot(g_dat, aes(gap_resid * 100, log2_oe)) +
    geom_hline(yintercept = 0, colour = pub_colours$light_grey, linewidth = .25) +
    geom_ribbon(data = g_line, aes(x = x, ymin = lo, ymax = hi, fill = model), inherit.aes = FALSE, alpha = .12) +
    geom_line(data = g_line, aes(x, y, colour = model), linewidth = .5) +
    geom_point(aes(size = expected, colour = top10), alpha = .7, stroke = .3) +
  scale_colour_manual(values = c(top_cols, "All EU/EEA" = pub_colours$reference,
    "Top 10" = top_cols[["TRUE"]]),
    breaks = names(g_model_labels), labels = unname(g_model_labels), name = "IRR") +
  scale_fill_manual(values = c("All EU/EEA" = pub_colours$reference,
    "Top 10" = top_cols[["TRUE"]]), guide = "none") +
    scale_size_continuous(range = c(.4, 3), breaks = c(10, 100, 1000), name = "Expected case count") +
    scale_x_continuous(breaks = g_x_breaks, limits = range(g_x_breaks), expand = c(0, 0)) +
    scale_y_continuous(breaks = g_y_breaks, limits = range(g_y_breaks), expand = c(0, 0)) +
    labs(x = "Age-aligned gap residual\n(percentage points)", y = "log₂(observed/expected cases)") +
  guides(colour = guide_legend(order = 1, nrow = 2, byrow = TRUE),
    size = guide_legend(order = 2, nrow = 1,
      override.aes = list(shape = 16, colour = pub_colours$text, fill = NA,
        linetype = 0, linewidth = 0))) +
  base_theme + theme(legend.position = c(.99, .01), legend.justification = c(1, 0),
    legend.direction = "vertical", legend.box = "vertical",
    legend.title.position = "top",
    legend.key.height = unit(5.5, "mm"), legend.spacing.y = unit(1.2, "mm"),
    legend.box.spacing = unit(1, "mm"),
      legend.background = element_rect(fill = scales::alpha("white", .82), colour = NA),
      legend.margin = margin(1.5, 1.5, 1.5, 1.5), legend.spacing.x = unit(1.5, "mm"),
      plot.margin = margin(1, 5, 1, 2))

  # S12B. Screening-equivalent coverage under assumed vaccine effectiveness
  # (VE 0.93–0.97), implied by the vaccinated share of cases, against national
  # MCV1 coverage of the same ages. Points below the diagonal indicate that the
  # case composition is compatible with lower local coverage under the model.
  s12b_dat <- screen %>% mutate(hotspot = factor(hotspot, levels = c(FALSE, TRUE)))
  s12b_sum <- screen_sum %>% filter(group %in% c("top10", "other")) %>%
    mutate(name = c(top10 = "Top 10", other = "Other")[group]) %>% arrange(desc(group)) %>%
    mutate(txt = sprintf("%s: %.0f pp (%.0f–%.0f)", name, 100 * weighted_deficit, 100 * lower, 100 * upper)) %>%
    pull(txt) %>% paste(collapse = "\n")
  s12b <- ggplot(s12b_dat, aes(national_mcv1 * 100, screening_equivalent_coverage * 100, colour = hotspot)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = pub_colours$reference, linewidth = .3) +
    geom_linerange(aes(ymin = screening_equivalent_lower * 100, ymax = screening_equivalent_upper * 100), linewidth = .3, alpha = .6) +
    geom_point(aes(size = known_status), alpha = .85) +
    annotate("text", x = 61, y = 99, hjust = 0, vjust = 1, size = 2.1, colour = pub_colours$text,
      label = paste0("Case-weighted deficit\n", s12b_sum)) +
    scale_colour_manual(values = top_cols, labels = top_labels, name = "Country group") +
    scale_size_continuous(range = c(.6, 3.4), breaks = c(100, 1000, 10000), labels = label_number(big.mark = ","),
      name = "Cases with known status") +
    scale_x_continuous(limits = c(60, 100), breaks = c(60, 70, 80, 90, 100), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
    labs(x = "National MCV1 coverage at case ages (%)",
      y = "Screening-equivalent coverage\nunder assumed vaccine effectiveness (%)") +
    guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(size = 1.5)),
      size = guide_legend(order = 2, nrow = 1)) +
    base_theme + theme(legend.position = c(0.01, .001), legend.justification = c(0, 0),
      legend.direction = "vertical", legend.box = "vertical",
      legend.title.position = "top",
      legend.key.height = unit(5.5, "mm"), legend.spacing.y = unit(1.2, "mm"),
      legend.box.spacing = unit(1, "mm"),
      legend.background = element_rect(fill = scales::alpha("white", .82), colour = NA),
      legend.margin = margin(1.5, 1.5, 1.5, 1.5), plot.margin = margin(1, 5, 1, 2))

  # Source data follow the printed panel letters.
  displayed <- list(a = a_dat, b = b_dat,
    c = b_cov_dat, d = c_dat, e = d_dat, f = e_dat, g = g_dat)
  for (i in names(displayed)) write_csv(displayed[[i]], file.path(src_dir, paste0("Figure3_", i, ".csv")))
  write_csv(s12a_dat, file.path(src_dir, "FigureS12_a.csv"))
  write_csv(s12b_dat, file.path(src_dir, "FigureS12_b.csv"))

  # cowplot measures grobs on its null device; the default PDF device
  # under-measures bold legend titles and clips them. ggsave writes PNG via
  # ragg, so measure with the matching AGG device.
  old_null <- cowplot::set_null_device("agg")
  on.exit(cowplot::set_null_device(old_null), add = TRUE)
  row1 <- wrap_plots(tag(panel_a, "a"), tag(b_cases, "b"),
    tag(b_cov, "c"), tag(panel_c, "d"), widths = c(1, .95, .8, 1.1), guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
      legend.title.position = "left")
  # Keep f and g as regular grobs because freeing a panel that carries a
  # multi-row legend triggers a patchwork free-table simplification error.
  row2 <- cowplot::plot_grid(patchwork::free(tag(panel_d, "e"), type = "panel"),
    tag(panel_e, "f"), tag(panel_g, "g"), ncol = 3,
    rel_widths = c(1, 1, 1.2), align = "none")
  figure <- cowplot::plot_grid(row1, row2, ncol = 1,
    rel_heights = c(.9, .9), align = "none")
  save_publication_figure(figure, "Figure3_NC_redesigned", figure_dir, width_mm = 183, height_mm = 140)

  figure_s12 <- cowplot::plot_grid(tag(s12a, "A"), tag(s12b, "B"), ncol = 2,
    rel_widths = c(1, 1), align = "none")
  save_publication_figure(figure_s12, "FigureS12_measles_ranking_and_screening",
    file.path(out, "figures/supplementary"), width_mm = 183, height_mm = 95)
  invisible(figure)
}

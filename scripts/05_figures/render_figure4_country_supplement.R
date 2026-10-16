# Figures S4–S8: annual history–duration trajectories for every model country.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})
source("scripts/05_figures/publication_figure_style.R")
source("scripts/05_figures/prepare_figure4_scenarios.R")
source("scripts/05_figures/helpers/figure4_scenario_style.R")

render_figure4_country_supplement <- function(out = "output/revised", data = NULL) {
  if (is.null(data)) data <- prepare_figure4_scenarios(out)$country_annual
  required <- c("country", "country_name", "selection_order", "disease", "history",
    "duration", "year", "mean", "median", "lo", "hi", "lo50", "hi50", "n_draws")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Country trajectory data missing: ", paste(missing, collapse = ", "), call. = FALSE)

  key <- c("country", "disease", "history", "duration", "year")
  summaries <- c("mean", "median", "lo", "hi", "lo50", "hi50")
  countries <- data %>% distinct(country, country_name, selection_order) %>% arrange(selection_order)
  stopifnot(
    nrow(countries) == 10L, !anyDuplicated(countries$country),
    !anyDuplicated(countries$country_name), setequal(countries$selection_order, 1:10),
    !anyNA(data[required]), !anyDuplicated(data[key]), nrow(data) == 2000L,
    setequal(data$disease, pub_disease_levels), setequal(data$year, 2014:2023),
    setequal(data$history, c("observed", "sustained_high")),
    setequal(data$duration, c("waning", "unwaning")),
    n_distinct(data$n_draws) == 1L, all(data$n_draws >= 2L),
    all(data$n_draws == as.integer(data$n_draws)),
    all(vapply(data[summaries], function(x) all(is.finite(x) & x > 0), logical(1))),
    all(data$lo <= data$lo50 & data$lo50 <= data$median &
      data$median <= data$hi50 & data$hi50 <= data$hi)
  )
  # Unique keys, the allowed values above and this complete row count ensure
  # that every disease contains all 10 countries × 4 scenarios × 10 years.
  counts <- data %>% count(disease)
  stopifnot(all(counts$n == 400L))

  figure_dir <- file.path(out, "figures", "supplementary")
  source_dir <- file.path(out, "source_data", "figure_redesign")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  scenario_levels <- names(figure4_scenario_colours)
  country_labels <- stats::setNames(
    paste0(letters[seq_len(nrow(countries))], "  ", countries$country_name), countries$country)

  plots <- vector("list", length(pub_disease_levels))
  names(plots) <- pub_disease_levels
  for (i in seq_along(pub_disease_levels)) {
    disease <- pub_disease_levels[[i]]
    figure_id <- paste0("FigureS", i + 3L)
    stem <- paste0(figure_id, "_revised_", tolower(disease), "_country_trajectories")
    table <- data %>% filter(.data$disease == .env$disease) %>%
      mutate(
        country = factor(country, levels = countries$country),
        scenario = factor(paste(history, duration, sep = "_"), levels = scenario_levels)
      ) %>% arrange(selection_order, scenario, year)
    stopifnot(!anyNA(table$scenario))

    p <- ggplot(table, aes(year, mean, colour = scenario, fill = scenario, group = scenario)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .12, colour = NA, show.legend = FALSE) +
      geom_line(aes(linetype = duration), linewidth = .55) +
      facet_wrap(~country, ncol = 2, scales = "free_y", axes = "all",
        axis.labels = "all", labeller = labeller(country = country_labels)) +
      scale_colour_manual(values = figure4_scenario_colours, breaks = scenario_levels,
        labels = unname(figure4_scenario_labels[scenario_levels]), name = NULL) +
      scale_fill_manual(values = figure4_scenario_colours, guide = "none") +
      scale_linetype_manual(values = figure4_scenario_linetypes, guide = "none") +
      scale_x_continuous(limits = c(2013.7, 2023.3), breaks = c(2014, 2017, 2020, 2023),
        expand = c(0, 0)) +
      scale_y_log10(labels = label_number(scale_cut = cut_short_scale()),
        breaks = figure4_log_breaks,
        limits = function(values) range(figure4_log_breaks(values)),
        expand = c(0, 0)) +
      labs(x = "Calendar year", y = "Annual expected reports per country\n(log₁₀ scale)") +
      guides(colour = guide_legend(nrow = 1, byrow = TRUE,
        override.aes = list(linewidth = .8,
          linetype = unname(figure4_scenario_linetypes[c("waning", "unwaning", "waning", "unwaning")])))) +
      theme_publication(grid = "y", plot_margin = margin(3, 4, 3, 4)) +
      theme(
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6.8, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 3)),
        axis.title.y = element_text(margin = margin(r = 4)),
        axis.ticks.length = unit(.8, "mm"),
        strip.text = element_text(size = 7, face = "bold", hjust = 0,
          margin = margin(t = 1, b = 3)),
        panel.spacing.x = unit(6, "mm"), panel.spacing.y = unit(4, "mm"),
        legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 5.8), legend.key.width = unit(5, "mm"),
        legend.key.height = unit(2.8, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(3, "mm"), legend.spacing.x = unit(2, "mm")
      )

    write_csv(select(table, all_of(required), scenario),
      file.path(source_dir, paste0(figure_id, "_a.csv")))
    save_publication_figure(p, stem, figure_dir, width_mm = pub_full_width_mm,
      height_mm = 210, dpi = 600)
    plots[[disease]] <- p
    message("Rendered ", figure_id, ": ", disease,
      "; 10 countries, four posterior mean trajectories and central 95% intervals.")
  }
  invisible(plots)
}

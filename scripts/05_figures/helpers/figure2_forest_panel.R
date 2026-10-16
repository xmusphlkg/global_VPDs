# A compact forest/table panel. ggdist draws the supplied Wald intervals;
# cowplot aligns headers, shared disease rows and the two quantitative axes.
# Package references: https://cran.r-project.org/package=ggdist
#                     https://cran.r-project.org/package=cowplot
figure2_forest_panel <- function(a, colours, fills, comparison = NULL) {
  for (pkg in c("ggdist", "cowplot")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Figure 2 requires ", pkg, "; install.packages(\"", pkg, "\")", call. = FALSE)
    }
  }
  a <- dplyr::mutate(a, y = match(disease, rev(pub_disease_levels)))
  shared_y <- function() scale_y_continuous(limits = c(.5, 5.5), breaks = NULL, expand = c(0, 0))
  band <- function(log_axis = FALSE) geom_rect(data = data.frame(y = c(2, 4),
      left = if (log_axis) .06 else -Inf, right = if (log_axis) 40 else Inf),
    aes(xmin = left, xmax = right, ymin = y - .5, ymax = y + .5),
    inherit.aes = FALSE, fill = "#F4F5F5", colour = NA)
  base <- function(numeric_axis = FALSE) {
    theme_publication(show_titles = TRUE, plot_margin = margin(0, 0, 0, 0)) +
      theme(plot.subtitle = figure1_facet_title_style(),
        plot.title = element_blank(), plot.background = element_rect(fill = "white", colour = NA),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.x = element_text(size = 5.4),
        axis.title.x = figure1_axis_title_style(margin(t = 2)),
        axis.line.x = if (numeric_axis) element_line(linewidth = .26) else element_blank(),
        axis.ticks.x = if (numeric_axis) element_line(linewidth = .22) else element_blank(),
        axis.ticks.length = unit(.7, "mm"), legend.position = "none")
  }
  text_axis <- function() scale_x_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0))
  disease <- ggplot(a, aes(.08, y)) + band() +
    geom_text(aes(label = disease),
      hjust = 0, size = 2.05) +
    shared_y() + text_axis() + labs(x = NULL, y = NULL, subtitle = "Disease") + base()
  counts <- ggplot(a, aes(.5, y)) + band() +
    geom_text(aes(label = scales::comma(n, accuracy = 1)), size = 2.0) +
    shared_y() + text_axis() + labs(x = NULL, y = NULL, subtitle = "n") + base()
  forest <- ggplot(a, aes(rate_ratio, y, colour = disease, fill = disease)) + band(TRUE) +
    geom_vline(xintercept = 1, colour = "#91989B", linetype = 2, linewidth = .24) +
    ggdist::geom_pointinterval(aes(xmin = rr_low, xmax = rr_high), orientation = "horizontal",
      shape = 21, point_size = 1.8, stroke = .45, linewidth = .45,
      interval_size_domain = c(0, 1), interval_size_range = c(0, 1), show.legend = FALSE) +
    scale_colour_manual(values = colours) + scale_fill_manual(values = fills) +
    scale_x_log10(limits = c(.06, 40), breaks = c(.1, 1, 10),
      labels = c("0.1", "1", "10"), expand = c(0, 0)) +
    shared_y() + labs(x = "Rate ratio per\n0.10 increase in cohort-time gap", y = NULL, subtitle = "Association") +
    base(TRUE)
  estimates <- ggplot(a, aes(.5, y)) + band() +
    geom_text(aes(label = sprintf("%.2f (%.2f–%.2f)", rate_ratio, rr_low, rr_high)), size = 1.95) +
    shared_y() + text_axis() + labs(x = NULL, y = NULL, subtitle = "RR (95% CI)") + base()
  stat_column <- function(values, title) {
    ggplot(a, aes(.5, y)) + band() +
      geom_text(aes(label = values), size = 2.0) +
      shared_y() + text_axis() + labs(x = NULL, y = NULL, subtitle = title) + base()
  }
  format_p <- function(x) ifelse(is.na(x), "—", ifelse(x < .001, sprintf("%.4f", x), sprintf("%.3f", x)))
  countries <- stat_column(as.character(a$n_country), "Countries")
  pvalues <- stat_column(format_p(a$p_value), "P")
  qvalues <- stat_column(format_p(a$q_value), "q")
  plots <- list(disease, countries, counts, forest, estimates, pvalues, qvalues)
  widths <- c(1.2, .9, .65, 1.9, 2.05, .7, .7)
  if (!is.null(comparison)) {
    plots <- c(plots, list(comparison + labs(subtitle = " ")))
    widths <- c(widths, sum(widths) / 1.65)
  }
  cowplot::plot_grid(plotlist = plots, nrow = 1, align = "h", axis = "tb",
    rel_widths = widths, labels = if (is.null(comparison)) NULL else c("a", rep("", 6), "b"),
    label_size = 8, label_x = 0, hjust = 0)

}

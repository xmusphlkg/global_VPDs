# Same ggdist point-interval and cowplot forest/table grammar as Figure 2a.
figure4_country_forest_panel <- function(d, colours,
                                         diseases = c("Pertussis", "Diphtheria", "Mumps")) {
  for (pkg in c("ggdist", "cowplot")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Figure 4d requires ", pkg, call. = FALSE)
  }
  d <- dplyr::mutate(d, y = 11 - selection_order,
    country_label = dplyr::recode(country, IND = "India", CHN = "China", USA = "USA", IDN = "Indonesia",
      PAK = "Pakistan", BRA = "Brazil", RUS = "Russia", DEU = "Germany", THA = "Thailand", GBR = "UK"))
  shared_y <- function() scale_y_continuous(limits = c(.5, 10.5), breaks = NULL, expand = c(0, 0))
  band <- function() geom_rect(data = data.frame(y = seq(1, 9, 2)),
    aes(xmin = -Inf, xmax = Inf, ymin = y - .5, ymax = y + .5),
    inherit.aes = FALSE, fill = "#F3F5F5", colour = NA)
  base <- function(numeric_axis = FALSE) theme_publication(plot_margin = margin(0, 0, 0, 0)) +
    theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.y = element_blank(), axis.text.x = element_text(size = 5.1),
      axis.line.x = if (numeric_axis) element_line(linewidth = .26) else element_blank(),
      axis.ticks.x = if (numeric_axis) element_line(linewidth = .22) else element_blank(),
      axis.ticks.length = unit(.7, "mm"), legend.position = "none")
  text_x <- function() scale_x_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0))
  head <- function(label) cowplot::ggdraw() + cowplot::draw_label(label, size = 5.8,
    fontfamily = "sans", fontface = "bold", colour = pub_colours$text)
  country <- ggplot(dplyr::distinct(d, y, country_label), aes(.04, y)) + band() +
    geom_text(aes(label = country_label), hjust = 0, size = 2.0) + shared_y() + text_x() + base()
  body <- list(country); widths <- 1.15
  # Separate axes expose small reference-disease changes; units remain percent.
  d <- dplyr::filter(d, disease %in% diseases)
  diseases <- diseases[diseases %in% unique(as.character(d$disease))]
  for (dis in diseases) {
    z <- dplyr::filter(d, disease == dis)
    # Manual breaks per disease
    ticks <- list(Measles = c(-1, 0, 1), Rubella = c(-1, 0, 1), Pertussis = c(0, 100, 200),
      Diphtheria = c(0, 20, 40, 60, 80), Mumps = c(0, 500, 1000, 1500))
    breaks <- ticks[[dis]]
    limits <- range(breaks)
    plot <- ggplot(z, aes(median, y)) + band() +
      geom_vline(xintercept = 0, colour = "#9EA3A6", linetype = 2, linewidth = .24) +
      ggdist::geom_pointinterval(aes(xmin = lo, xmax = hi), orientation = "horizontal",
        shape = 21, point_size = 1.4, stroke = .4, linewidth = .4,
        colour = colours[[dis]], fill = colours[[dis]],
        interval_size_domain = c(0, 1), interval_size_range = c(0, 1)) +
      scale_x_continuous(limits = limits, breaks = breaks, expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      shared_y() + base(TRUE)
    # Numeric columns are data, not in-panel explanatory annotations.
    labels <- ggplot(z, aes(.5, y)) + band() +
      geom_text(aes(label = sprintf("%+.2f", median)),
        size = 1.65) + shared_y() + text_x() + base()
    body <- c(body, list(plot, labels)); widths <- c(widths, 1.5, 1.03)
  }
  # One alignment pass keeps every country row and interval exactly level.
  aligned <- cowplot::plot_grid(plotlist = body, nrow = 1, align = "h", axis = "tb", rel_widths = widths)
  headers <- cowplot::plot_grid(plotlist = c(list(head("Country")),
    lapply(diseases, function(x) head(sub(" \\(proxy\\)$", "", unname(pub_disease_labels[x]))))),
    nrow = 1, rel_widths = c(1.15, rep(1.5 + 1.03, length(diseases))))
  footer <- cowplot::ggdraw() + cowplot::draw_label("Waning-associated increase vs no waning (%)", size = 6,
    fontfamily = "sans", fontface = "bold", colour = pub_colours$text)
  cowplot::plot_grid(headers, aligned, footer, ncol = 1, rel_heights = c(.08, 1, .08))
}

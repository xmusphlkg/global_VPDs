# Figure 1 quantitative axes share one contract: pretty() breaks, no expansion,
# limits default to the extrema of those breaks; calendar axes use tile edges.
# Categorical row/column
# labels use pretty-generated carrier axes with unused carrier ticks hidden.
figure1_axis_factory <- function() {
  axis <- function(aesthetic, id, values, n = 4, labels = waiver(), limits = NULL) {
    values <- values[is.finite(values)]
    stopifnot(length(values) > 0)
    breaks <- pretty(range(values), n = n)
    scale_limits <- if (is.null(limits)) range(breaks) else limits
    scale <- if (aesthetic == "x") scale_x_continuous else scale_y_continuous
    scale(breaks = breaks, limits = scale_limits, labels = labels, expand = c(0, 0),
          minor_breaks = NULL)
  }
  list(x = function(id, values, n = 4, labels = waiver(), limits = NULL)
    axis("x", id, values, n, labels, limits),
       y = function(id, values, n = 4, labels = waiver(), limits = NULL)
         axis("y", id, values, n, labels, limits))
}

figure1_axis_title_style <- function(margins = margin()) {
  element_text(size = 6, face = "bold", colour = pub_colours$text, margin = margins)
}

# Colourbar titles sit in a separate column to the left of a horizontal bar.
# Aligning the title's right and top edges to that column keeps long, wrapped
# titles clear of the bar labels when legends are placed below a panel.
figure1_legend_title_style <- function(margins = margin(r = 2)) {
  element_text(size = 6, face = "bold", colour = pub_colours$text,
    hjust = 1, vjust = 1, margin = margins)
}

figure1_bottom_legend <- function() {
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title.position = "left",
    legend.justification = "center",
    legend.title = figure1_axis_title_style(margin(r = 2)),
    legend.text = element_text(size = 5.2, colour = pub_colours$text),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    legend.spacing.x = unit(1.2, "mm")
  )
}

# A vertical key to the right of the panels, with its title above the bar.
figure1_right_legend <- function() {
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.title.position = "top",
    legend.justification = "center",
    legend.title = figure1_axis_title_style(margin(b = 1.5)),
    legend.text = element_text(size = 5.2, colour = pub_colours$text),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1.5, "mm")
  )
}

# `width` is the bar's long axis and `height` its thickness, whichever way the
# bar runs. A top title is left-aligned over the bar; a left title keeps the
# right/top alignment described above.
figure1_colourbar <- function(title, width = 22, height = 1.7, direction = "horizontal",
                              title_position = "left") {
  long <- unit(width, "mm"); thick <- unit(height, "mm")
  title_style <- if (title_position == "top") {
    element_text(size = 6, face = "bold", colour = pub_colours$text,
      hjust = 0, margin = margin(b = 1.5))
  } else figure1_legend_title_style()
  guide_colourbar(title = title, display = "rectangles", nbin = 64,
    direction = direction, theme = theme(
      legend.key.width = if (direction == "horizontal") long else thick,
      legend.key.height = if (direction == "horizontal") thick else long,
      legend.title.position = title_position,
      legend.title = title_style,
      legend.margin = margin(0, 0, 0, 0)))
}

figure1_carrier_labels <- function(positions, labels) {
  force(positions); force(labels)
  function(x) { ans <- unname(labels[match(x, positions)]); ans[is.na(ans)] <- ""; ans }
}

# Resolve the Figure 1f facet style once for matching headers in other figures.
# This includes the inherited sans family, dark text, centred bold 5.8-pt type,
# line height and strip padding rather than approximating it with subtitles.
figure1_facet_title_style <- function() {
  ggplot2::calc_element("strip.text.x.top",
    theme_publication() + theme(strip.text = element_text(size = 5.8)))
}

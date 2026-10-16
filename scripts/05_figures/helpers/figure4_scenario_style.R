# Shared four-scenario identities for Figure 4a and its country supplements.
figure4_scenario_colours <- c(
  observed_waning = "#0072B2", observed_unwaning = "#56B4E9",
  sustained_high_waning = "#D55E00", sustained_high_unwaning = "#E69F00")
figure4_scenario_labels <- c(
  observed_waning = "Observed history · waning",
  observed_unwaning = "Observed history · no waning",
  sustained_high_waning = "10–14-year reference · waning",
  sustained_high_unwaning = "10–14-year reference · no waning")
figure4_scenario_linetypes <- c(waning = "solid", unwaning = "dashed")

# Constant ratios give equally spaced ticks on a log axis. Use powers of two
# for narrow ranges and powers of ten for ranges spanning several decades;
# never append linear pretty endpoints to logarithmic ticks.
figure4_log_breaks <- function(values) {
  values <- values[is.finite(values) & values > 0]
  stopifnot(length(values) > 0)
  bounds <- range(values)
  ratio <- if (bounds[2] / bounds[1] <= 16) 2 else 10
  anchor <- 10^floor(log10(bounds[1]))
  exponents <- log(bounds / anchor, base = ratio)
  anchor * ratio^seq(floor(exponents[1] + 1e-10), ceiling(exponents[2] - 1e-10))
}

# Symmetric logarithmic display with a linear core of +/- one unit.
# With units of thousand reports, 0, 1, 10, 100 and 1000 have equal spacing.
figure4_signed_log <- function(x) {
  sign(x) * ifelse(abs(x) <= 1, abs(x), 1 + log10(pmax(abs(x), 1)))
}

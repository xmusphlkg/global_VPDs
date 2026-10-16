#!/usr/bin/env Rscript

# Supplementary Figure S11: nonlinear pertussis association and matched
# pre-pandemic model comparisons.
suppressPackageStartupMessages({library(tidyverse); library(glmmTMB); library(patchwork)})
source("scripts/05_figures/publication_figure_style.R")

out_dir <- "output/revised/figures/supplementary"
src_dir <- "output/revised/source_data/figure_redesign"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(src_dir, recursive=TRUE, showWarnings=FALSE)

d <- read_csv("output/robustness_case_models/analysis_data.csv", show_col_types=FALSE) %>%
  filter(disease == "Pertussis")
fit <- readRDS("output/robustness_case_models/models/Pertussis_nonlinear_spline.rds")
lo <- quantile(d$gap_10, .025, na.rm=TRUE); hi <- quantile(d$gap_10, .975, na.rm=TRUE)
grid <- tibble(
  population=100000, country=first(d$country), coverage_10=median(d$coverage_10),
  prior=median(d$prior), gdp=median(d$gdp), uhc=median(d$uhc),
  urban=median(d$urban), health=median(d$health), year=2018,
  gap_10=seq(as.numeric(lo), as.numeric(hi), length.out=121))
pr <- predict(fit, newdata=grid, type="link", re.form=NA, se.fit=TRUE)
ref_gap <- median(d$gap_10, na.rm=TRUE)
ref <- tibble(gap_10=ref_gap, population=100000, country=first(d$country), coverage_10=median(d$coverage_10),
              prior=median(d$prior), gdp=median(d$gdp), uhc=median(d$uhc), urban=median(d$urban),
              health=median(d$health), year=2018)
eta0 <- as.numeric(predict(fit, newdata=ref, type="link", re.form=NA))
# Form the fixed-effect contrast against the median reference.  Using the
# covariance of the spline coefficients keeps uncertainty in the reference
# value; predict(..., se.fit=TRUE) alone gives a pointwise interval around each
# fitted curve value and omits that shared-reference covariance.
tt <- delete.response(terms(fit))
fixed_design <- function(x, year_value = 2018) {
  years <- sort(unique(d$year))
  xx <- x[rep(seq_len(nrow(x)), each = length(years)), , drop = FALSE]
  xx$year <- rep(years, times = nrow(x))
  mm <- model.matrix(tt, xx)
  mm[xx$year == year_value, , drop = FALSE]
}
grid_mm <- fixed_design(grid)
ref_mm <- fixed_design(ref)
delta_mm <- sweep(grid_mm, 2, ref_mm[1, ], FUN = "-")
beta_vcov <- as.matrix(vcov(fit)$cond)
contrast_se <- sqrt(rowSums((delta_mm %*% beta_vcov) * delta_mm))
curve <- grid %>% transmute(gap_10, rr=exp(as.numeric(pr$fit)-eta0),
                            low=exp(as.numeric(pr$fit)-eta0-1.96*contrast_se),
                            high=exp(as.numeric(pr$fit)-eta0+1.96*contrast_se))
readr::write_csv(curve, file.path(src_dir, "FigureS11_a_pertussis_spline_curve_source_data.csv"))

pre <- read_csv("output/robustness_case_models/pre_pandemic_model_comparison.csv", show_col_types=FALSE)
readr::write_csv(pre, file.path(src_dir, "FigureS11_b_pre_pandemic_model_comparison_source_data.csv"))

th <- theme_publication(base_size=6.2, grid="x", show_titles=TRUE, plot_margin=margin(3,4,3,4)) +
  theme(axis.text=element_text(size=5.4), axis.title=element_text(size=6.2, face="bold"),
        plot.title=element_text(size=8, face="bold", hjust=0),
        panel.grid.major.y=element_blank())

p1 <- ggplot(curve, aes(gap_10, rr)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=.25, colour=pub_colours$reference) +
  geom_ribbon(aes(ymin=low, ymax=high), fill=pub_colours$orange, alpha=.18) +
  geom_line(linewidth=.65, colour=pub_colours$orange) +
  geom_rug(data=d, aes(x=gap_10), sides="b", inherit.aes=FALSE, alpha=.18, length=unit(1.5,"mm")) +
  scale_y_log10() +
  labs(title="A",
       x="Scaled cohort-time gap (10× absolute gap)", y="Rate ratio") + th

p2 <- ggplot(pre, aes(x=delta_AIC_M3_minus_M1, y=reorder(disease, delta_AIC_M3_minus_M1))) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=.25, colour=pub_colours$reference) +
  geom_point(size=1.8, colour=pub_colours$blue) +
  geom_text(aes(label=sprintf("%.2f", delta_AIC_M3_minus_M1)), hjust=-.25, size=2.2) +
  scale_x_continuous(expand=expansion(mult=c(.05,.18))) +
  labs(title="B",
       x="ΔAIC (cohort-time gap − coverage history)", y=NULL) + th

fig <- p1 + p2 + plot_layout(widths=c(1.25,.75))
save_publication_figure(fig, "FigureS11_nonlinear_and_pre_pandemic", out_dir,
                        width_mm=pub_full_width_mm, height_mm=88, dpi=600)
message("Wrote Figure S11 and source data")

# Redesigned figure source data

The active set contains six main figures and twelve supplementary figures.
Figures 1 and 2 contain six panels (a–f), Figure 3 seven (a–g), Figure 4 seven
(a–g), Figure 5 six (a–f) and Figure 6 three (a–c). Supplementary Figures
S4–S8 each contain ten country facets.
The numbered `Figure*_*.csv` files contain the data used by the
corresponding panels.
Panel tables retain input units; transformations such as proportions to
percentage points are explicit in the R renderer and figure legend.

Figure legends are maintained only in `manuscripts/manuscript_nature_communications.md`
and `manuscripts/SUPPLEMENTARY_MATERIAL.md`. Renderers export PNG/PDF figures
and source-data CSVs, without separate captions, audit reports or session logs.
`FigureS*_revised_source_data.csv` retains additional input provenance.

Figure 6 includes NHANES MMRV age profiles (a) and Kilifi, Kenya 2021 DTP
age profiles (c). The diphtheria series is minimum anti-toxin protection; the
pertussis series is the sum of source anti-PT categories at or above 62.5
EU/mL, used only as an operational recent-exposure marker. The Kilifi rows are
single-population external evidence and are aligned with Kenya's 2021 model
profile. Source files and threshold limitations are documented in
`output/revised/external_validation/KILIFI_DTP_README.md` and
`output/revised/external_validation/PIENTER_DTP_README.md`.

Reproduce from the repository root:

```sh
Rscript scripts/05_figures/render_main_figures.R
Rscript scripts/05_figures/render_supplementary_figures.R
```

Both renderers read existing analytical results and use the shared style in
`scripts/05_figures/publication_figure_style.R`; they do not refit models. Bayesian scenario
contrasts preserve posterior draw IDs and compare matched scenarios before
calculating intervals. Figure 4 reports expected reported-case contrasts,
whereas Figure 5 reports deterministic susceptibility, contact contributions
and initial effective reproduction numbers. The late-period
Bayesian calibration check is explicitly distinguished from a refitted forecast.

Earlier main and supplementary graphics remain recoverable from Git history.

### Figure 1 atlas

Regenerate only this figure with:

```sh
Rscript scripts/05_figures/render_main_figures.R --figure1-only
```

`prepare_figure1_atlas.R` separates data preparation from rendering. Figure 1
source tables use `layer` to identify plotted observations, summaries, ribbons
and sensitivity scenarios. Panel A uses Ukraine and the supported ages 1–14.
Panel C contains the actual nine-scenario parameter envelope, not a time-series
quantile band or a confidence interval. Panel E uses a reproducible country-cluster
bootstrap of age means and paired age contrasts. Panel F centres both variables
within the same matching strata; its slope bands come from the existing 2,000
country-bootstrap results. Recovery contrasts formerly in panel G now appear
in Supplementary Figure S2E (`FigureS2_e.csv`), with both recovery definitions,
conditional intervals and event counts. Regenerate S2 with
`Rscript scripts/05_figures/render_supplementary_figures.R --s2-only`. No case or transmission model is refitted.

Figure 1 layout uses native `patchwork::free()` compositions and local legends
(`guides = "keep"`). Display labels omit the Mumps suffix; the
underlying source rows are unchanged.

### Figure 2 empirical evidence

```sh
Rscript scripts/05_figures/render_main_figures.R --figure2-only
```

Panel a uses `ggdist::geom_pointinterval` (currently `ggdist` 3.3.3) and
`cowplot` to align the forest and table columns. No package version is pinned.
Install these dependencies if missing:

```sh
Rscript -e 'install.packages(c("ggdist", "cowplot"), repos="https://cloud.r-project.org")'
```

`prepare_figure2_evidence.R` reads the authoritative coefficient, AIC,
within/between, LOCO, biological-sensitivity and temporal-validation tables.
It uses the saved DTP Cohort-time gap fits to calculate covariate-standardized rates and
pointwise log-delta-method intervals from the full fixed-effect covariance;
country random intercepts are fixed at zero, not integrated out. No models
are refitted. The curves stop at the observed 2.5th/97.5th gap percentiles;
all observations remain in the rug source layer. Pale support bands use the
10th/90th percentiles calculated from those same observed gaps, stored as
`central_support_low_pp` and `central_support_high_pp`.

Displayed model names are Current coverage, Coverage history, Constant-protection
cohort gap and Cohort-time gap. Each alternative adds one historical representation
to the Current coverage baseline; the constant-protection cohort gap assumes
constant protection and the Cohort-time gap uses disease-specific protection
durations. The Figure 2b source table reports these formal names; internal source
identifiers remain in the analytical outputs. AIC comparisons retain ΔAIC = alternative AIC minus reference AIC; negative values favour the alternative.

`Figure2_a.csv` stores the primary rate ratios, confidence intervals,
country–year sample sizes and Cohort-time gap ΔAIC values over Coverage history.
`Figure2_b.csv` contains 20 `disease_comparison` rows for the matrix and two
`pooled_interaction` rows for its inset. The latter preserve the four-primary
and five-disease scopes, AIC, likelihood-ratio statistic, degrees of
freedom, P value and sample size from `pooled_model_comparisons.csv`
(Supplementary Table S5). All gains use reference minus alternative AIC.
`Figure2_c.csv` distinguishes standardized curves and observed-gap rugs. `Figure2_c_p10_p90.csv` provides the derived disease-specific P10–P90 gap contrasts and RRs reported in the manuscript (Supplementary Table S43).
`Figure2_d.csv` retains all decomposition and LOCO results.
`Figure2_e.csv` retains the 18 original parameter fits, including confidence
intervals and primary-scenario flags; tiles use equally spaced parameter
positions, with no interpolation. Both tile colours and signed labels are
ΔAIC (cohort-time gap − coverage history) over Coverage history, with the same zero-centred diverging scale
as panel b. Rate ratios and intervals are retained as source data only.
`Figure2_f.csv` retains both baseline and alternative RMSE, split dates and
sample sizes. Its x axis is disease-specific contemporaneous ΔAIC (cohort-time gap − coverage history); its
y axis comes from pooled temporal models with a common slope, observed lagged
cases and linear time. Both axes compare Cohort-time gap with Coverage history.
Facet labels show both training cutoffs and held-out test years. Overall RMSE
changes of +12.6% and +7.0% are reported in the caption. Prediction with
disease-specific gap slopes is not evaluated by these pooled models. These changes use overall
RMSE, not a mean of disease-specific percentage changes. Sign counts for the
LOCO fits, sensitivity-cell meanings and primary-scenario markers are also
explained in the caption. The displayed disease name is Mumps; its all-MCV mapping assumption is
documented in Methods (Component mapping) and retained in source provenance.

Figure 2 uses `AA / BC / DE / DF`: panel d spans the two lower rows, with e
above f on the right. Column headers and panel facet titles use the resolved
Figure 1f facet style (5.8-pt dark bold sans, centred, with matching padding).

Axes and legends follow Figure 1: 0.7-mm ticks, 0.26-mm axis lines,
0.22-mm tick lines, 6-pt axis titles with 2-pt padding, 5.4-pt tick labels,
and 3.5-mm horizontal facet gaps. Panel d draws both axes and ticks in every
facet. Only the left column has y-axis labels, and only the bottom row has
x-axis labels and titles; internal axes retain their lines and ticks.
Panels b and e retain their axes, tick marks and labels. The different data
scales (including log rate ratios and discrete parameter grids) are intentional.

### Figure 3 EU/EEA measles resurgence

```sh
Rscript scripts/05_figures/render_main_figures.R --figure3-only
```

`render_figure3_measles_resurgence.R` reads the completed extension in
`output/revised/measles_resurgence_full/` and writes `Figure3_a.csv` through
`Figure3_g.csv`, one per printed panel: WHO/UNICEF EU/EEA monthly incidence
(a), cases (b) and 2025 WUENIC coverage (c) of the ten leading countries,
within-country monthly Z-scores (d), 2023 age-band gaps (e), coverage-adjusted
rank correlations (f) and primary closed-band same-birth-cohort age-alignment
points (g). Supplementary Fig. S12 source data contain the added-variable
ranking view (a) and screening-equivalent coverage under assumed vaccine
effectiveness (b).

### Figure 4 paired history and waning scenarios

```sh
Rscript scripts/05_figures/render_main_figures.R --figure4-only
```

`prepare_figure4_scenarios.R` reads the corrected ten-country Bayesian
scenario, annual prediction and predictive-check tables. It verifies unique
keys, a complete 2014–2023 grid, agreement between annual sums and scenario
totals, and agreement between cumulative endpoints and ten-year ratios. The renderer uses the completed corrected-input refit and does not fit models itself. `render_figure4_scenarios.R` inherits the Figures 1–2 font,
axis and facet styles, and exports 183 × 165 mm, 600-dpi PNG and vector PDF.

`Figure4_a.csv` contains annual equal-country expected reports for all five
diseases × four scenario states, with annual draw and summary layers. Panel a
shows the posterior mean after equal-country averaging within each draw;
central 95% intervals remain in the source table but are not drawn. The renderer
uses four scenario colours and free log-scaled y axes by disease facet.
The original ten-year draw/summary layers
remain in the same source table for the paired downstream contrasts. The
interaction_draw and interaction_summary layers retain the paired original-
scale difference-in-differences. Its summary is also written to
`bayesian_transmission/history_waning_interaction.csv` and exported as
`TableS43_history_waning_interaction.csv`. Log-axis visual slopes do not
estimate this additive interaction.

`Figure4_b.csv` retains equal-country duration contrasts in cases.
`Figure4_c.csv` retains country-summed report ratios for observed, recovery
and high history under waning. `Figure4_d.csv` contains all 50 country/disease
percentage contrasts and intervals. `figure4_country_forest_panel.R` uses
`ggdist::geom_pointinterval` and `cowplot::plot_grid`, as in Figure 2a, to align
country rows, forests and two-decimal median columns and graphical 95% intervals. Disease
axes have different limits, with the same percentage units.
`Figure4_e.csv` retains cumulative expected-report ratios for pertussis and
diphtheria, including all three scenarios and all paired draws. The recovery
updates occur in 2018–2022; annual standardized trajectories can cross the
observed-history path even when the cumulative endpoint remains lower.
`Figure4_f.csv` retains the 179 observed rows used in the held-out negative-
binomial predictive-median scatter. `Figure4_g.csv` retains disease-level and
overall negative-binomial predictive coverage, distinct from mean-report
uncertainty.

Figures 2–3 share `pub_evidence_colours`: purple Measles, teal Rubella, orange
Pertussis, blue Diphtheria and magenta Mumps. Every disease has the same symbol
rules; open marks only encode no waning in a. In c the observed-history reference is a dashed line at 100%, without points.
Explanatory annotations and repeated result labels are moved into the
manuscript caption. Axes, disease headings, numeric interval columns and
necessary scenario keys remain visible.

Figure 4d has no significance stars; full interval numbers remain in Source
Data. Figure 4b uses a signed pseudo-log transform with sigma 0.25 thousand
reports and explicit negative/zero/positive ticks. Point-mass zero contrasts
have no invented density. Expected-report axes explicitly distinguish modelled
observed-history reference scenarios from surveillance observations. Recovery
starts at a vertical dashed line in 2018; explanations remain in the caption.

The age-mapping and missing-value audit is resolved by a complete refit. Model
bands now include ages 1–4, 5–9 and 10–14 correctly; the WPP infant share combines
an age-0 gap of 1.0 with ages 1–4. Ages 15+ retain the declared 0.30 assumption,
and missing childhood gaps use the same 0.30 in both duration scenarios. Fixed
nuisance inputs remain exactly identical across chains. Seven optimization
starts per disease precede four chains of 8,000 iterations (2,000 burn-in).
Scenario summaries use 300 paired draws per disease from the 12,000-iteration
mode-guarded chains. Pertussis tail ESS is below the publication target and is
labelled exploratory in the manuscript. Figure 4 f/g use the 179 observations
held out from calibration (2020–2023), not the fitted-period predictive table.
`model_manifest.json`, `input_audit.csv` and `map_multistart_audit.json`
document the active inputs and refit. Superseded outputs
remain recoverable from Git history.
Measles/rubella expected duration contrasts are exactly zero. Intervals for
expected reports condition on fixed structural inputs; observation noise enters
only the posterior predictive checks, and event replicates are separate.


### Supplementary Figures S4–S8 country trajectories

`FigureS4_a.csv` through `FigureS8_a.csv` contain measles, rubella, pertussis,
diphtheria and mumps trajectories, respectively. Each table has 400 rows:
ten selected countries × four history-by-duration scenarios × 2014–2023.
The fields `country`, `name` and `selection_order` preserve country identity
and selection order; `disease`, `history`, `duration` and `year` identify the
annual scenario. `mean` is the posterior mean of expected annual reports,
`median` is the posterior median, `lo`/`hi` are central 95% bounds,
`lo50`/`hi50` are central 50% bounds, and `n_draws` records the draw count.
The country summaries are computed directly from annual expected-report draws,
without country pooling or a model refit.

Lines show `mean` and ribbons show `lo`–`hi` using the Figure 4a scenario
colours. Each figure has ten independently log-scaled country facets.
Intervals exclude observation noise and uncertainty in fixed nuisance inputs.
Mumps follows the main mapping; Indonesia and Pakistan have no observed mumps
series and use the shared disease posterior without country-specific calibration.
Regenerate all five with
`Rscript scripts/05_figures/render_supplementary_figures.R --s4-s8-only`.

### Figure 5 age and transmission mechanisms

```sh
Rscript scripts/05_figures/render_main_figures.R --figure5-only
```

`Figure5_a.csv`–`Figure5_c.csv` hold age-gap gradients, the single-age
waning-minus-no-waning susceptibility grid and vaccine-action routing;
`Figure5_a_annual_gradients.csv` keeps the annual gradients behind panel a.
`Figure5_d.csv`–`Figure5_f.csv` hold the contact-weighted operator change,
dominant-mode contributions and country ΔRₑ. `Figure5_e_summary.csv`
summarises panel e by age band, and `Figure5_sensitivity_grid.csv` retains
the former sensitivity grid. Contact marginals are in `FigureS9_a.csv` and
`FigureS9_b.csv`.

### Figure 6 external serology

`Figure6_a.csv`–`Figure6_c.csv` hold the NHANES age changes, cross-survey MMR
contrasts and Kilifi 2021 DTP comparisons.

Files superseded by the 14 October 2026 renumbering (measles moved from
Figure 6 to Figure 3) are kept in `output/revised/archive/figure_renumber_20261014/`.

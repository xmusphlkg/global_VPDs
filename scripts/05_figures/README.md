# Nature Communications figure system

Run from the repository root:

```bash
Rscript scripts/05_figures/render_main_figures.R
Rscript scripts/05_figures/render_supplementary_figures.R
```

The active renderers read completed model outputs and write six main PNG/PDF
pairs to `output/revised/figures/main_nc/` plus nine supplementary PNG/PDF pairs to
`output/revised/figures/supplementary/`. PNG files use 600 dpi RGB output and
PDFs use Cairo vector graphics. Exact plotted values are written to the
matching source-data directories.

| Figure | Renderer | Scientific role |
|---|---|---|
| 1 | `render_figure1_atlas.R` | Annual coverage, historical cohort states and the three linked analyses |
| 2 | `render_figure2_evidence.R` | Disease-specific case associations, robustness and prediction boundaries |
| 3 | `render_figure3_measles_resurgence.R` | EU/EEA measles resurgence located by country and age (panels a–g; S12 contains the moved ranking/screening panels) |
| 4 | `render_figure4_scenarios.R` | Paired history and waning scenarios for expected reports |
| 5 | `render_figure5_mechanisms.R` | Age redistribution and contact-weighted transmission potential (panels a–f) |
| 6 | `render_figure6_serology.R` | External serology triangulation |

Each figure can be regenerated alone without refitting any model, for example:

```bash
Rscript scripts/05_figures/render_main_figures.R --figure3-only
```

Figure 3 reads the completed measles extension in
`output/revised/measles_resurgence_full/` (see `scripts/README.md` for the
four upstream steps). Its source tables `Figure3_a.csv`–`Figure3_g.csv`
follow the printed panel letters. The same call renders Supplementary Figure
S12 (added-variable view of the 1–4-year gap and screening-equivalent
coverage) with source tables `FigureS12_a.csv` and `FigureS12_b.csv`.

Figure 4a shows posterior mean annual expected reports after equal-country
averaging within each draw, without interval ribbons. Supplementary Figures
S4–S8 show the corresponding four scenarios separately for all ten countries,
with central 95% posterior intervals, one disease per figure. Their source
tables are `FigureS4_a.csv` through `FigureS8_a.csv`. To update only these
country trajectories from completed model outputs:

```bash
Rscript scripts/05_figures/render_supplementary_figures.R --s4-s8-only
```

Figure 5 shows age-gap gradients (a), single-age susceptibility (b) and
vaccine-action routing (c), then follows pertussis, diphtheria and mumps
through the contact operator (d), dominant-mode contributions (e) and country
responses (f). The former sensitivity grid is written to
`Figure5_sensitivity_grid.csv`; S9 contains contact marginals. It reads the
current blocking assumptions from the mechanism manifest. To refresh Figure 5
and the country-bootstrap summary used in Table 1:

```bash
python3 scripts/03_mechanism/summarize_figure5_country_bootstrap.py --replicates 2000 --seed 20261006
Rscript scripts/05_figures/render_main_figures.R --figure5-only
```

Figure 6 combines NHANES MMR markers, cross-survey MMR contrasts and Kilifi
2021 diphtheria minimum-protection and pertussis recent-exposure markers,
aligned to the corresponding Kenya 2021 programme profiles.

Audit the active figure set with:

```bash
python3 scripts/06_package/audit_figure_set.py \
  --main-markdown manuscripts/manuscript_nature_communications.md \
  --main-dir output/revised/figures/main_nc --main-formats png,pdf \
  --supplement-markdown manuscripts/SUPPLEMENTARY_MATERIAL.md \
  --supplement-dir output/revised/figures/supplementary --supplement-formats png,pdf
```

Superseded plotting scripts and exports remain recoverable from Git history.

# Main Figures

Six Nature Communications figures are rendered here as 600-dpi RGB PNG and
Cairo PDF pairs, `Figure1_NC_redesigned` through `Figure6_NC_redesigned`.
Their panel source tables are in `output/revised/source_data/figure_redesign/`.

Render from the repository root:

```bash
Rscript scripts/05_figures/render_main_figures.R
```

`--figureN-only` (N = 1–6) regenerates a single figure.

| Figure | Content | Panels |
|---|---|---|
| 1 | Cohort-memory atlas | a–f |
| 2 | Disease-specific case associations and prediction boundaries | a–f |
| 3 | EU/EEA measles resurgence located by country and age | a–i |
| 4 | Paired history and waning scenarios for expected reports | a–g |
| 5 | Age redistribution and contact-weighted transmission potential | a–f |
| 6 | External serology triangulation | a–c |

Figure 4a displays posterior mean annual expected reports after equal-country
averaging, without interval ribbons; country-specific curves and central 95%
posterior intervals are in Supplementary Figures S4–S8. Contact marginals for
Figure 5 are in Supplementary Figure S9. External markers and
programme-derived gaps in Figure 6 remain distinct estimands.

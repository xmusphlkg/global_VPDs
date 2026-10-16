# ESEN2 historical measles, mumps and rubella serology

These CSV files transcribe the age-band tables from the WHO Bulletin ESEN2 reports:

- Rubella: [Nardone et al., 2008, PMC2647372](https://pmc.ncbi.nlm.nih.gov/articles/PMC2647372/), Table 1 (sample frame) and Table 3 (age-band seronegativity).
- Measles: [Andrews et al., 2008, PMC2647410](https://pmc.ncbi.nlm.nih.gov/articles/PMC2647410/), Table 1 (sample frame) and Table 3 (age-band seronegativity).

The records are country-by-age-band aggregate estimates. `seronegative_pct` is the published percentage with a negative IgG result; conversion to seropositivity is `1 - seronegative_pct / 100`. `n_young_source` is the published count for the broad young-age source band (1–19 years for rubella and <20 years for measles). The reports omit a count for every age band, so `n_age_band_reported` remains empty. Blank seronegative values indicate an age band with no published estimate.

The serum banks were collected during 1996–2004. Eleven collections used residual clinical sera, five used population samples, and one combined both approaches. Assays differed among countries and were standardized by the ESEN protocol. These data support historical age-gradient calibration and meta-regression. They are separate from the NHANES estimand and should be analysed with country and collection-year strata; individual-level records and design weights are unavailable.

- Mumps: [Miller et al., 2012, PMC9152080](https://pmc.ncbi.nlm.nih.gov/articles/PMC9152080/), Table 4 (18 European countries, 1996–2004 serum banks). The table reports age-specific mumps IgG positivity and the tested denominator for most 1–4, 5–9, 10–14, and 15–19 year bands. Sweden used selected ages, Spain used ages 2–4 in the first band, and Luxembourg had an age-4-only first band; these rows carry an age annotation.

Mumps seropositivity is an antibody endpoint. The ESEN paper explicitly describes antibody levels as an imperfect proxy for protection during outbreaks; retain mumps as exploratory and pair it with vaccine schedule, coverage, strain, and incidence data.

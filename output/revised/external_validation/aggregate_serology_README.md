# Aggregate serology calibration summary

`aggregate_serology_summary.csv` reports within-survey age contrasts from three published aggregate sources. Each survey remains a separate assay and sampling frame. The calculations use the published age-band percentages; they create no synthetic individual records and perform no cross-era prevalence pooling.

## Sources and contrasts

- ESEN2 serum banks (1996–2004) use country-level contrasts between the oldest and youngest comparable child bands. Measles uses 10–19 minus 2–4 years, rubella uses 10–14 minus 2–4 years, and mumps uses 10–14 minus 1–4 years. Measles and rubella source values are seronegative percentages and are converted to seropositivity by `100 - seronegative_pct`. The mumps source reports seropositivity directly. The source contains 18 countries for measles and mumps and 17 for rubella; comparable contrasts are available for 18 measles, 16 mumps, and 15 rubella countries because of source-specific age-band omissions. Sweden's mumps values represent selected ages, and Spain/Luxembourg have nonstandard youngest bands; these records remain flagged in the input CSV and are retained in the country count only when the requested bands exist. A standard-band mumps sensitivity excluding annotated records is also reported. Sources: [measles ESEN2 Table 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC2647410/), [rubella ESEN2 Table 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC2647372/), and [mumps ESEN2 Table 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC9152080/).
- KiGGS Basiserhebung (2003–2006) reports national survey-weighted age-band percentages for 1–2 and 14–17 years. Positive, borderline and negative titre categories are retained separately; denominators remain as published. Source: [KiGGS Table 2](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0042867).
- PIENTER2 (2006–2007) reports national mumps seropositivity by age. The summary contrasts the unweighted mean across ages 10–14 with the mean across ages 5–8, bracketing the routine second dose at age 9, and reports the age-10 minus age-8 jump. These are descriptive age-band comparisons; source population weights remain unavailable for reconstruction. Source: [PIENTER2 Table 1](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0058234).

## Audit corrections

The ESEN2 rubella Romania values were corrected to the source values 82.2%, 62.6% and 36.7% seronegative for ages 2–4, 5–9 and 10–14 years. The source table labels were corrected to measles Table 3 and mumps Table 4. Sweden's three selected-age mumps records were added with an age annotation. The script checks these fields before writing the summary.

Run `python3 scripts/04_validation/summarize_external_serology.py` to reproduce the two CSV outputs in this directory.

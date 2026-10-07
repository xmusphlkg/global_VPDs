# Independent NHANES serological construct calibration

This analysis compares age-specific NHANES MMRV IgG serology with a programme-derived cohort gap. These endpoints remain distinct: the model gap is a vaccine-history protection scenario, while IgG reflects detected antibody after vaccination or infection.

The design frame retains all WTSSMMPP-positive serum records aged 6–59 years. Age domains are represented by indicators so strata and PSUs outside a domain contribute zero residuals, following CDC Taylor-linearization guidance. The primary interval is a stratified PSU Taylor ratio variance with logit-t confidence limits and domain degrees of freedom equal to represented PSUs minus represented strata. Kish intervals are retained as sensitivity columns.

Aged 6–14 years includes 2,541 age-eligible records and 1,592 records with both valid positive weight and valid assay interpretation. Missing assays are excluded rather than assigned a negative result. Single-year estimates with fewer than ten negatives are flagged. Broad age bands are in `nhanes_age_band_serology.csv`; model years 2017–2020 are in `nhanes_model_year_sensitivity.csv`; waning and no-waning scenarios are descriptive in `nhanes_waning_sensitivity.csv`.

USA 2018 projected cohort ages are flagged in `profile_status`; `known_only` metrics exclude them. Mumps uses the labelled all-MCV proxy and remains exploratory because the assay threshold is a seropositivity threshold without an established correlate of protection.

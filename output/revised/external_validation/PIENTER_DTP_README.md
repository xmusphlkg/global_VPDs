# Dutch PIENTER DTP serology inputs

These two files preserve public aggregate results that can be used to extend the external-serology context of Figure 5:

- `pienter2_diphtheria_age_serology.csv` contains age 0-14 rows from the individual-level supplementary workbook for the Dutch PIENTER2 2006/2007 diphtheria survey (PMC4749226, `pone.0148605.s002.xls`, `table4`). It reports raw national-sample counts and exact binomial confidence intervals for anti-diphtheria toxin IgG at the minimum-protection threshold (`>=0.01 IU/mL`) and full-protection threshold (`>=0.1 IU/mL`).
- `pienter2_diphtheria_age_bands.csv` contains the same raw records aggregated to ages 1-4, 5-8 and 10-14 for an older-minus-younger age contrast.
- `pienter3_pertussis_age_exposure.csv` contains the published national-sample weighted PIENTER3 2016/2017 estimates for anti-PT IgG `>=100 IU/mL` in age groups 7-11 through 80+ (PMC8589709, main text section 4.1/Table 3). The 7-11 estimate is 8.7% in the main-text paragraph and 8.6% in Table 3; the CSV keeps 8.7% and records this small published discrepancy.

Diphtheria assay values come from the published individual-level workbook and are **unweighted raw proportions**. The article's headline values (91% above 0.01 IU/mL in the national sample) are survey-weighted; these files do not recreate those weights. Use the raw age rows only as a transparent descriptive age profile or label them explicitly as unweighted. The source article defines `<0.01 IU/mL` as non-protective, `0.01-<0.1 IU/mL` as basic protection and `>=0.1 IU/mL` as full protection.

Pertussis anti-PT `>=100 IU/mL` is an operational indicator of recent infection/exposure when there has been no recent vaccination. It is not a protection threshold and there is no validated single serological correlate of protection for pertussis. Do not plot `1 - anti-PT prevalence` as susceptibility or combine these rows with measles/rubella/mumps seronegativity as one estimand.

Source URLs:

- Diphtheria: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4749226/>; supplementary archive: <https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4749226/supplementaryFiles>
- Pertussis: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8589709/>; supplementary archive: <https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8589709/supplementaryFiles>

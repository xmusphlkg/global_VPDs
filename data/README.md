# Data

The public data release has three layers.

- curated/ contains the active five-disease parameter specification and the documented mumps component evidence table.
- processed/ contains cohort, coverage, demographic and model-ready snapshots used by the active workflow.
- external/ contains public-use NHANES, ESEN2, KiGGS, Kilifi and PIENTER records with local provenance notes.
- raw/ contains selected public source snapshots, including the single-age WPP population table, Prem contact matrices, WUENIC coverage, WHO case/schedule extracts and World Bank indicators.

The source downloads are public data products. The preparation scripts in scripts/00_data/ can recreate omitted large downloads when the corresponding source files are available locally. The analysis uses country-level and de-identified public-use records; no private or identifying records are included.

Primary provenance links:

- WHO Immunization Data Portal and WHO/UNICEF WUENIC: https://immunizationdata.who.int/
- UN World Population Prospects: https://population.un.org/wpp/
- Prem, Cook and Jit contact matrices: https://doi.org/10.1371/journal.pcbi.1005697
- NHANES public-use files: https://www.cdc.gov/nchs/nhanes/
- ESEN2, KiGGS, Kilifi and PIENTER source notes: see the README files in data/external/.

# Public analysis scripts

Run commands from the repository root. The numbered directories follow the
public data-to-results workflow.

- `00_data/` downloads and harmonizes public source snapshots.
- `01_cohort/` reconstructs the five disease age-cohort protection tables.
- `02_models/` fits the count models, daily transmission scenarios and released robustness analyses.
- `03_mechanism/` runs the age-structured deterministic mechanism.
- `04_validation/` runs external validation, matched-history analyses and the EU/EEA measles-resurgence extension.

Core commands are:

```bash
python3 scripts/01_cohort/reconstruct_cohorts.py
Rscript scripts/02_models/run_five_disease_models.R output/five_disease_revision output/revised data/processed/case_data_quality_preserved.csv
.venv-transmission/bin/python scripts/02_models/run_bayesian_transmission.py --iterations 12000 --burn 2000 --chains 4 --draws 300 --workers 20 --seed 20261006
python3 scripts/03_mechanism/run_age_seirv_mechanism.py
python3 scripts/04_validation/validate_external_surveillance.py
python3 scripts/04_validation/run_measles_resurgence_full.py
Rscript scripts/04_validation/fit_measles_resurgence_models.R
python3 scripts/04_validation/build_measles_hotspot_analysis.py
Rscript scripts/04_validation/fit_measles_age_alignment_fe.R
Rscript scripts/02_models/robustness_case_models.R output/five_disease_revision output/robustness_case_models data/processed/case_data_quality_preserved.csv
python3 scripts/02_models/run_transmission_catchup_scenarios.py --draws 120 --campaign-year 2018 --cap 0.10
```

The public snapshot excludes the main manuscript, tests, package assembly and
internal audit or diagnostic utilities. Publication figure rendering, figure
source data, supplementary tables and final supplementary attachments are
included.

# Global VPDS public reproducibility release

This repository is the public code and data release for the global five-disease vaccine-protection analysis. It contains the analysis workflow, curated parameter inputs, processed country and cohort data, selected public source snapshots, and numerical model results.

The associated manuscript title is **Historical vaccination trajectories preserve age-specific differences in programme-derived protection after coverage convergence**.

The release is organized so that the main calculations can be rerun from the repository root:

    # Cohort reconstruction
    python3 scripts/01_cohort/reconstruct_cohorts.py

    # Deterministic age-structured mechanism
    python3 scripts/03_mechanism/run_age_seirv_mechanism.py

    # External serology validation
    python3 scripts/04_validation/validate_external_surveillance.py

The Bayesian transmission fit is computationally intensive. Install requirements-transmission.txt in an isolated environment before running scripts/02_models/run_bayesian_transmission.py. The cohort count models require an R installation with the packages used by the corresponding scripts. Figure rendering and manuscript assembly are maintained outside this public snapshot.

## Repository layout

- scripts/00_data/ downloads and prepares source inputs.
- scripts/01_cohort/ reconstructs age 1–14 birth-cohort protection states.
- scripts/02_models/ contains count-model and daily transmission-model code.
- scripts/03_mechanism/ contains the age-structured mechanism analysis.
- scripts/04_validation/ contains the public external-serology validation code.
- data/curated/ stores the active parameter and mumps component specifications.
- data/processed/ stores the reproducible cohort and model input snapshots.
- data/external/ stores public-use external validation data and derived tables.
- data/raw/ stores selected public source snapshots; larger downloads are documented by the data-preparation scripts.
- output/ stores reconstructed cohorts and numerical analysis results.

## Data provenance

The input layer uses WHO/UNICEF WUENIC coverage, WHO reported cases and schedules, UN World Population Prospects 2024, World Bank indicators, Prem et al. contact matrices, and public NHANES/ESEN2/KiGGS/Kilifi/PIENTER datasets. Source URLs and data-handling notes are recorded in data/README.md and the README files beside the external snapshots.

Results are scenario and model outputs conditioned on the documented parameter specifications. They are provided for reproducibility and should be interpreted with the assumptions in data/curated/five_disease_parameters.json and scripts/02_models/transmission_model_spec.json.

## Public-release boundary

This repository is a filtered public snapshot. Internal test suites, editable submission documents, package and figure-rendering work products, audit/diagnostic logs, parameter-review records, virtual environments, and large raw downloads are intentionally excluded. The exclusion rules are recorded in PUBLIC_RELEASE_SCOPE.md and .gitignore.

# Public analysis scripts

Run commands from the repository root. The numbered directories follow the data-to-results workflow.

1. 00_data/ downloads and harmonizes public source snapshots.
2. 01_cohort/ reconstructs the five disease age-cohort protection tables.
3. 02_models/ fits the count models and the daily transmission scenarios.
4. 03_mechanism/ runs the age-structured deterministic mechanism.
5. 04_validation/ runs the public external-serology validation.

The public release keeps code needed to regenerate the released data and results. Working-tree tests, submission assembly helpers, figure rendering, and audit-only utilities are not part of this repository.

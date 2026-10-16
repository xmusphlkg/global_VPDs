# Public release scope

This repository is a curated public snapshot of the VPDS analysis. It contains reproducible workflow code, curated and processed inputs, selected public source snapshots, numerical analysis results, publication figures, figure source data, supplementary tables and final supplementary attachments.

The source working tree was filtered before publication. The snapshot intentionally excludes:

- the main manuscript and editable manuscript working files;
- internal test suites;
- audit, diagnostic, review, parameter-audit and log artefacts;
- manuscript package/assembly helpers and private review exports;
- virtual environments and large raw downloads that exceed the public snapshot scope.

Final supplementary attachments, publication figures, figure source data and publication supplementary tables are included. The repository `.gitignore` repeats the exclusion rules so later local runs do not accidentally add internal work products. Raw inputs not included here can be regenerated with the download scripts in `scripts/00_data/`.

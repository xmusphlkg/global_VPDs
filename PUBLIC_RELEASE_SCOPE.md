# Public release scope

This repository is a curated public snapshot of the VPDS analysis. It contains reproducible workflow code, curated and processed inputs, selected public source snapshots (including the single-age population and contact matrices), and numerical analysis results.

The source working tree was filtered before publication. The snapshot intentionally excludes:

- editable submission documents and manuscript working files;
- internal test suites;
- audit, diagnostic, review and log artefacts;
- manuscript package/assembly helpers;
- publication figure exports and supplementary table exports;
- virtual environments and large raw downloads that exceed the public snapshot scope.

The repository .gitignore repeats these boundaries so later local runs do not accidentally add them. Raw inputs not included here can be regenerated with the download scripts in scripts/00_data/.

# Repaired Bayesian transmission layer

The restored model uses a daily solver, 2000-2013 infection-immunity burn-in, WPP annual age-population anchoring, and 2014-2019 calibration. Years 2020-2023 are reserved for time-out validation. The discrete transmission coefficient is normalized so the disease-free next-generation number equals the sampled R0.

Only R0 and negative-binomial size are sampled. Vaccine effects, clinical fractions, importation and reporting remain fixed at disease-specific prior centres, so report estimates are conditional scenarios, not fully identified burden posteriors. `scenario_sensitivity.csv` varies adult gap, high-history cap, 2020-2023 contact scale, natural-immunity duration, vaccine effects, clinical fraction, importation, and reporting on paired posterior draws. These grids quantify conditional variation; they do not resolve non-identifiability or replace external validation. Seasonality remains absent.

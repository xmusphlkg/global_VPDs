# Age-structured SEIRV mechanism layer

This is an independent mechanism analysis for 10 selected countries and five diseases.
It uses WPP 2024 single-age populations, Prem et al. (2017) all-location contact matrices,
and the reconstructed cohort protection index. Aggregate case counts are not used to fit beta
or reporting parameters; missing case rows therefore cannot become zero observations.

The cohort gap describes a protection-state scenario. Infection susceptibility is obtained
with `s = 1 - B * (1 - gap)`, where B is an explicit transmission-blocking scenario. This
separates clinical protection from acquisition or onward-transmission protection. Adult
susceptibility and B are sensitivity dimensions. Outputs are epidemic potentials and
structural comparisons, not causal effects or observed incidence forecasts.

Review mapping: the cohort protection mass already includes clinical VE.
Use an acquisition-to-clinical ratio for recipients and a separate source
factor for protected people among infections. The operator is `diag(s) C diag(t)`.
Diphtheria can have zero acquisition contrast and a positive source contrast.
Source protection is fixed at the central value in the acquisition grid;
the separate operator sensitivity table varies the source effect.

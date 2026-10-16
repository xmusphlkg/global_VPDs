# Systematic current-coverage matching and observed recovery histories

The analysis uses every eligible country–disease–year with current coverage, independently of reported cases. Its estimand is heterogeneity in an already reconstructed protection-gap index conditional on current coverage and calendar year. It does not estimate a causal recovery effect or independently validate the reconstruction.

## Design

- Primary strata are exact disease × year × current WUENIC percentage (MCV1 for measles and the mumps proxy, RCV1 for rubella, DTP3 for pertussis/diphtheria). Only strata containing at least two distinct countries contribute.
- The primary statistic is the pooled within-stratum standard deviation: sqrt(sum_i w_i (gap_i − weighted_mean_gap_stratum)^2 / sum_i w_i), reported in percentage points. Countries–years have equal weight; each eligible unit appears once. This is residual dispersion, not a mean pairwise contrast.
- Prespecified sensitivities use 2- and 5-percentage-point bins anchored at zero, exact complete current dose vectors (MCV1/MCV2 or DTP1/DTP3), and birth-cohort-size weights across countries. Bins are left-closed and right-open; 100% has its own terminal bin.
- All configurations also report equal-age-standardized gaps (1/14 per age), removing differences in the birth-cohort age composition, and constant-protection gaps.
- A within-stratum slope relates the gap to the mean coverage deficit over the 10 preceding calendar years, excluding the current year. It is scaled per 10 percentage points of past deficit, requires all 10 observations, and is a decomposition of the index, not an independent outcome validation.
- Percentile intervals use 2,000 bootstrap samples of whole countries, preserving their full time series. Repeated copies of a country are represented as multiplicity weights and never supply their own matched partner. Strata losing all but one distinct country do not contribute in that replicate. The intervals condition on observed coverage, eligibility, efficacy coefficients, duration assumptions and the selected index specification.

## Primary results

| Disease | Matched/eligible country–years | Countries | Exact strata | Gap SD, percentage points (95% country-bootstrap interval) | Equal-age gap SD, percentage points | Gap difference per 10-percentage-point past-deficit increase, percentage points (95% interval) |
|---|---:|---:|---:|---|---:|---|
| Measles | 271/368 | 70 | 96 | 2.25 (1.42–2.57) | 2.26 | 9.10 (8.31–10.25) |
| Rubella | 1078/1177 | 127 | 178 | 3.06 (2.31–3.36) | 3.04 | 8.80 (8.41–9.23) |
| Pertussis | 1782/1932 | 192 | 280 | 1.31 (1.11–1.34) | 0.83 | 1.55 (1.28–1.74) |
| Diphtheria | 1782/1932 | 192 | 280 | 2.41 (1.86–2.39) | 2.41 | 4.82 (4.41–5.10) |
| Mumps (all-MCV proxy) | 271/368 | 70 | 96 | 2.13 (1.40–2.43) | 1.99 | 8.23 (7.54–9.10) |

These contrasts establish systematically that the current marginal coverage percentage does not determine the reconstructed protection state. The equal-age and complete-dose-vector results distinguish historical information from current additional-dose coverage and population age composition. The degree of heterogeneity need not increase when waning is introduced: rapid decay can compress between-country differences near a high gap index.

## Observed recovery trajectories

Recovery is defined without looking at gaps or cases: the first year during 2014–2022 with current and following-year coverage ≥90%, previous-year coverage <90%, and mean coverage over the previous three years ≤85%. A stricter sensitivity requires every preceding year ≤85%. Each country contributes at most one event per disease and definition. Following-year coverage confirms persistence; this descriptive selection must not be used as a prospective predictor. All events, including those without eligible reconstructed gaps, are retained in the audit. Trajectories extend from −3 to +5 years where observed through 2023; unavailable gap years remain missing, with no extrapolation.

At the recovery year, eligible recovering countries are compared descriptively with every eligible country in the same year at exactly the same current percentage, maintaining ≥90% in each preceding year and the next year. The reported event contrast intervals are country-cluster percentile bootstrap intervals with the same conditional interpretation as the matching intervals; they are not causal intervention intervals, and no hypothesis test is assigned to these small event sets.

- Measles: 7 coverage-defined recoveries; 1 with an eligible gap; no eligible exact benchmark comparison.
- Rubella: 5 coverage-defined recoveries; 3 with an eligible gap; 3 exact benchmark comparisons, mean recovered-minus-stable gap 6.22 percentage points (-0.41–17.26, country bootstrap).
- Pertussis: 6 coverage-defined recoveries; 6 with an eligible gap; 6 exact benchmark comparisons, mean recovered-minus-stable gap 0.76 percentage points (-0.21–1.63, country bootstrap).
- Diphtheria: 6 coverage-defined recoveries; 6 with an eligible gap; 6 exact benchmark comparisons, mean recovered-minus-stable gap 2.98 percentage points (0.93–5.16, country bootstrap).
- Mumps: 7 coverage-defined recoveries; 1 with an eligible gap; no eligible exact benchmark comparison.

## Interpretation and limitations

All estimated gaps are vaccine-derived primary-series indices. Infection-acquired protection, boosters, catch-up campaigns, migration and individual joint-dose histories are not observed here. Global dose timing and protection duration remain scenario inputs. Birth-cohort weights are original cohort sizes, not surviving resident populations. Equal-age standardization addresses age composition only; it does not remove those other assumptions. WUENIC percentage rounding means 'exact' matches are exact matches of published estimates, not proof of identical true coverage.

The mumps analysis is an explicitly labelled all-MCV component proxy. The strict component subset is too small for comparable country matching and is not pooled into the proxy. Coverage recoveries for the mumps proxy refer to the MCV1 platform. The same DTP or MCV programme events recur across their associated diseases; they are not independent programme recoveries.

The matching is descriptive and conditions only on current coverage and year; health-system quality, past infection, product, dose programmes and reporting may differ. Positive within-stratum gap–past-deficit slopes partly follow from construction. Neither the slopes nor the residual dispersion prove that the index predicts case rates. Observed recovery trajectories and duration perturbations must be distinguished from a longitudinal, explicitly timed recovery intervention in a transmission model.

## Manuscript-ready methods

We evaluated the information lost by current coverage using all age-complete reconstructed country–years, independently of case-report availability. Within each disease and calendar year, countries were stratified by their exact current primary coverage percentage. We quantified the pooled within-stratum standard deviation of the protection-gap index and estimated 95% percentile intervals by resampling countries with their complete time series (2,000 replicates). Sensitivities used 2- and 5-percentage-point coverage bins, complete current dose vectors, birth-cohort weights and equal-age-standardized indices. We separately enumerated coverage-defined recovery episodes and displayed their observed histories without outcome-based selection. These analyses characterize conditional heterogeneity and programme trajectories; they are not causal intervention estimates.

#!/usr/bin/env python3
"""Country-centred EU/EEA measles-resurgence analyses for Figure 3 and Supplementary Fig. S12.

Inputs are the data layers written by ``run_measles_resurgence_full.py``.
Outputs, all in ``output/revised/measles_resurgence_full``:

* ``measles_hotspot_country_summary.csv`` / ``measles_hotspot_monthly.csv``:
  the ten highest-burden EU/EEA countries (cumulative reported cases,
  January 2024–June 2026) and their 2016–2026 monthly series (descriptive);
* ``measles_hotspot_age_gaps.csv``: 2023 cohort gap by age band for every
  EU/EEA country (the top ten are shown in the figure);
* ``measles_hotspot_gap_contrast.csv``: top-ten-minus-other gap contrast
  (descriptive bootstrap; sensitivity);
* ``measles_coverage_adjusted_gap.csv`` / ``_points.csv``: does the 2023
  cohort gap rank 2024–2026 incidence beyond 2023 WUENIC MCV1 and MCV2?
* ``measles_screening_coverage.csv``: screening-equivalent coverage under
  assumed vaccine effectiveness, implied by case vaccination status, against
  national coverage of the same ages.

All analyses are ecological and resample countries, not individuals.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "output" / "revised" / "measles_resurgence_full"
OUT = SRC
TOP_N = 10
AGE_BANDS = ["1_4", "5_9", "10_14", "15_19", "20_29", "30_39"]
AGE_LABELS = {"1_4": "1–4", "5_9": "5–9", "10_14": "10–14", "15_19": "15–19", "20_29": "20–29", "30_39": "30–39"}
ADJ_TERMS = [f"gap_{b}" for b in AGE_BANDS] + ["pandemic_shortfall"]
N_BOOT = 2000
# Outcome window: January 2024 to June 2026, the last complete month.
WINDOW_YEARS = 2.5
MIN_POPULATION = 1e5
# Screening method: VE of any dose against measles, uniform over this range.
VE_RANGE = (0.93, 0.97)
MIN_KNOWN_STATUS = 20


def wuenic_wide() -> pd.DataFrame:
    wuenic = pd.read_csv(SRC / "measles_mcv_wuenic_1980_2025.csv")
    wide = wuenic.pivot_table(index=["country", "year"], columns="antigen", values="coverage", aggfunc="first").reset_index()
    wide[["MCV1", "MCV2"]] = wide[["MCV1", "MCV2"]] / 100
    return wide


def hotspot_summary(panel: pd.DataFrame, profiles: pd.DataFrame, coverage: pd.DataFrame) -> pd.DataFrame:
    totals = (panel.groupby(["region", "country", "country_name"], as_index=False)
              .agg(cases_2024_26=("cases", "sum"), peak_monthly_cases=("cases", "max"),
                   active_months=("cases", lambda x: int((x > 0).sum())),
                   outbreak_months=("outbreak_month", "sum"), sustained_months=("sustained_2month", "sum")))
    totals["rank_region"] = totals["cases_2024_26"].rank(method="first", ascending=False).astype(int)
    totals["hotspot"] = totals["rank_region"].le(TOP_N)
    keep = ["country", "gap_1_14", "gap_1_29", "historical_gap_5y", "historical_gap_10y"] + [f"gap_{b}" for b in AGE_BANDS]
    current = coverage.loc[coverage.year.eq(2025), ["country", "MCV1", "MCV2"]].rename(
        columns={"MCV1": "wuenic_mcv1_2025", "MCV2": "wuenic_mcv2_2025"})
    summary = (totals.merge(profiles.loc[profiles.year.eq(2025), keep], on="country", how="left", validate="many_to_one")
               .merge(current, on="country", how="left", validate="many_to_one"))
    gaps = summary[[f"gap_{b}" for b in AGE_BANDS]]
    summary["dominant_gap_label"] = gaps.fillna(-1).idxmax(axis=1).str.replace("gap_", "", regex=False).map(AGE_LABELS).where(gaps.notna().any(axis=1))
    summary.to_csv(OUT / "measles_hotspot_country_summary.csv", index=False)
    return summary


def hotspot_gap_contrast(summary: pd.DataFrame, rng: np.random.Generator) -> None:
    rows = []
    for term in ["gap_1_14"] + [f"gap_{b}" for b in AGE_BANDS]:
        sub = summary.dropna(subset=[term])
        hot = sub.loc[sub.hotspot, term].to_numpy(float)
        other = sub.loc[~sub.hotspot, term].to_numpy(float)
        draws = (rng.choice(hot, (N_BOOT, len(hot))).mean(axis=1) - rng.choice(other, (N_BOOT, len(other))).mean(axis=1))
        rows.append({"term": term, "n_hotspot": len(hot), "n_other": len(other),
                     "hotspot_minus_other_gap": hot.mean() - other.mean(),
                     "lower": np.quantile(draws, .025), "upper": np.quantile(draws, .975)})
    pd.DataFrame(rows).to_csv(OUT / "measles_hotspot_gap_contrast.csv", index=False)


def age_gap_profiles(summary: pd.DataFrame, profiles: pd.DataFrame) -> None:
    """2023 cohort gap by age band, the exposure fixed before the outcome window."""
    gaps = profiles.loc[profiles.year.eq(2023), ["country"] + [f"gap_{b}" for b in AGE_BANDS]]
    long = gaps.melt(id_vars="country", var_name="band", value_name="gap")
    long["band"] = long["band"].str.replace("gap_", "", regex=False)
    long = long.merge(summary[["country", "country_name", "hotspot", "rank_region", "cases_2024_26"]], on="country")
    long.sort_values(["rank_region", "band"]).to_csv(OUT / "measles_hotspot_age_gaps.csv", index=False)


def _rank_residual(y: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    """Residual of rank(y) on an intercept and the ranks of the covariates."""
    X = np.column_stack([np.ones(len(y))] + [rankdata(c) for c in covariates.T])
    ry = rankdata(y)
    beta, *_ = np.linalg.lstsq(X, ry, rcond=None)
    return ry - X @ beta


def _partial_rho(x: np.ndarray, y: np.ndarray, covariates: np.ndarray) -> float:
    ok = np.isfinite(x)
    if ok.sum() < 8 or len(np.unique(x[ok])) < 3:
        return np.nan
    return float(np.corrcoef(_rank_residual(x[ok], covariates[ok]), _rank_residual(y[ok], covariates[ok]))[0, 1])


def coverage_adjusted_gap(summary, profiles, coverage, archive, rng) -> None:
    """Rank partial correlation of incidence with each 2023 cohort gap.

    Spearman ranks of the exposure and of log incidence are residualised on
    the ranks of 2023 WUENIC MCV1 and MCV2 across EU/EEA countries; intervals
    resample countries. Exposures are fixed before the outcome window.
    """
    window = archive.loc[(archive.year >= 2024) & ((archive.year < 2026) | (archive.month <= 6))]
    incidence = window.groupby("country").agg(cases=("cases", "sum"), population=("offset_population", "first"))
    incidence["annual_incidence"] = 1e5 * incidence["cases"] / incidence["population"] / WINDOW_YEARS
    cov23 = coverage.loc[coverage.year.eq(2023)].set_index("country")[["MCV1", "MCV2"]]
    baseline = coverage.loc[coverage.year.between(2015, 2019)].groupby("country")["MCV1"].mean()
    pandemic = coverage.loc[coverage.year.between(2020, 2022)].groupby("country")["MCV1"].mean()
    gaps23 = profiles.loc[profiles.year.eq(2023)].set_index("country")[[f"gap_{b}" for b in AGE_BANDS]]
    d = (summary.set_index("country")[["region", "country_name", "hotspot"]]
         .join(cov23).join((baseline - pandemic).rename("pandemic_shortfall")).join(gaps23).join(incidence))
    # Microstates are excluded: single imported clusters dominate their rates.
    d = d.loc[d.population.gt(MIN_POPULATION)].dropna(subset=["MCV1", "MCV2", "annual_incidence"]).copy()
    d["log_incidence"] = np.log10(d["annual_incidence"] + 0.01)

    rows = []
    y = d["log_incidence"].to_numpy(float)
    cov = d[["MCV1", "MCV2"]].to_numpy(float)
    stats = {"current_mcv1": lambda i: spearmanr(y[i], cov[i, 0]).statistic,
             "current_mcv2": lambda i: spearmanr(y[i], cov[i, 1]).statistic}
    stats.update({t: (lambda i, x=d[t].to_numpy(float): _partial_rho(x[i], y[i], cov[i])) for t in ADJ_TERMS})
    full = np.arange(len(d))
    idx = rng.integers(0, len(d), (N_BOOT, len(d)))
    for term, fn in stats.items():
        draws = np.array([fn(i) for i in idx])
        rows.append({"term": term, "adjusted": term.startswith("gap_") or term == "pandemic_shortfall", "rho": fn(full),
                     "lower": np.nanquantile(draws, .025), "upper": np.nanquantile(draws, .975),
                     "n_countries": int(d[term].notna().sum()) if term in d else len(d)})
    for term in ["gap_1_4"]:
        ok = d[term].notna().to_numpy()
        d.loc[d.index[ok], f"resid_{term}"] = _rank_residual(d[term].to_numpy(float)[ok], cov[ok])
        d.loc[d.index[ok], "resid_incidence"] = _rank_residual(y[ok], cov[ok])
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "measles_coverage_adjusted_gap.csv", index=False)
    d.reset_index().to_csv(OUT / "measles_coverage_adjusted_points.csv", index=False)
    print(result.round(2).to_string(index=False))


def screening_coverage(summary: pd.DataFrame, rng: np.random.Generator) -> None:
    """Screening-equivalent coverage under assumed vaccine effectiveness.

    With PCV the vaccinated share of cases (known status) and VE the vaccine
    effectiveness, cases drawn from a population with coverage p satisfy
    PCV = p(1 - VE) / (p(1 - VE) + 1 - p), so p = PCV / (PCV + (1 - PCV)(1 - VE)).
    The result is compatible with the observed case vaccination composition
    under the screening model; it is not a directly observed community
    coverage denominator. The comparison is with ``national_mcv1``, the
    reconstructed coverage of the same ages nationally. Intervals combine a
    Jeffreys beta for PCV with VE drawn uniformly from ``VE_RANGE``.
    """
    vacc = pd.read_csv(SRC / "measles_case_vaccination_panel.csv")
    vacc = vacc.loc[vacc.known_status.ge(MIN_KNOWN_STATUS) & vacc.national_mcv1.notna()].reset_index(drop=True)
    vacc["hotspot"] = vacc["country"].isin(summary.loc[summary.hotspot, "country"])
    k, n = vacc["vaccinated"].to_numpy(float), vacc["known_status"].to_numpy(float)
    pcv = rng.beta(k + .5, n - k + .5, (N_BOOT, len(vacc)))
    ve = rng.uniform(*VE_RANGE, (N_BOOT, 1))
    implied = pcv / (pcv + (1 - pcv) * (1 - ve))
    point_pcv = k / n
    vacc["pcv"] = point_pcv
    vacc["screening_equivalent_coverage"] = point_pcv / (point_pcv + (1 - point_pcv) * (1 - np.mean(VE_RANGE)))
    vacc["screening_equivalent_lower"] = np.quantile(implied, .025, axis=0)
    vacc["screening_equivalent_upper"] = np.quantile(implied, .975, axis=0)
    # Keep the old column names as aliases for downstream reproducibility.
    vacc["implied_coverage"] = vacc["screening_equivalent_coverage"]
    vacc["implied_lower"] = vacc["screening_equivalent_lower"]
    vacc["implied_upper"] = vacc["screening_equivalent_upper"]
    deficit = vacc["national_mcv1"].to_numpy() - implied
    vacc["coverage_deficit"] = vacc["national_mcv1"] - vacc["screening_equivalent_coverage"]
    vacc["deficit_lower"] = np.quantile(deficit, .025, axis=0)
    vacc["deficit_upper"] = np.quantile(deficit, .975, axis=0)
    vacc.to_csv(OUT / "measles_screening_coverage.csv", index=False)

    # Case-weighted summaries (each country-year weighted by known-status cases).
    rows = []
    for group, sub in [("top10", vacc.loc[vacc.hotspot]), ("other", vacc.loc[~vacc.hotspot]), ("all", vacc)]:
        w = sub["known_status"].to_numpy(float)
        cols = sub.index.to_numpy()
        draws = (deficit[:, cols] * w).sum(axis=1) / w.sum()
        rows.append({"group": group, "country_years": len(sub), "cases_known_status": int(w.sum()),
                     "weighted_deficit": float(np.average(sub["coverage_deficit"], weights=w)),
                     "lower": np.quantile(draws, .025), "upper": np.quantile(draws, .975),
                     "share_country_years_deficit_gt0": float((sub["deficit_lower"] > 0).mean())})
    pd.DataFrame(rows).to_csv(OUT / "measles_screening_coverage_summary.csv", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


def main() -> None:
    rng = np.random.default_rng(20261014)
    panel = pd.read_csv(SRC / "measles_monthly_model_panel.csv")
    archive = pd.read_csv(SRC / "measles_country_monthly_2016_2026.csv")
    profiles = pd.read_csv(SRC / "measles_age_profiles_1980_2025.csv")
    coverage = wuenic_wide()

    # Rank on the outcome window only; months after June 2026 are incomplete.
    panel = panel.loc[(panel.year < 2026) | (panel.month <= 6)]
    summary = hotspot_summary(panel, profiles, coverage)
    # Hotspots are ranked on 2024–2026 burden; their trajectories use 2016–2026.
    monthly = archive.loc[archive.country.isin(summary.loc[summary.hotspot, "country"]),
                          ["region", "country", "country_name", "year", "month", "cases", "incidence_per_100k"]].copy()
    monthly["time"] = monthly.year + (monthly.month - 1) / 12
    monthly.to_csv(OUT / "measles_hotspot_monthly.csv", index=False)

    hotspot_gap_contrast(summary, rng)
    age_gap_profiles(summary, profiles)
    coverage_adjusted_gap(summary, profiles, coverage, archive, rng)
    screening_coverage(summary, rng)
    print(summary.loc[summary.hotspot, ["rank_region", "country", "cases_2024_26", "wuenic_mcv1_2025", "wuenic_mcv2_2025", "dominant_gap_label"]]
          .sort_values("rank_region").to_string(index=False))


if __name__ == "__main__":
    main()

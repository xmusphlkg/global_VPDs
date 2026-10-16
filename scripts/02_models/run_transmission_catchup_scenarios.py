#!/usr/bin/env python3
"""Evaluate implementable catch-up vaccination counterfactuals.

This script reuses the fitted posterior in output/revised/bayesian_transmission/
without re-fitting MCMC.  It applies a one-time campaign to the 10--14 birth
cohort in ``campaign_year`` and carries that cohort forward through successive
five-year age bands.  The campaign reduces the reconstructed gap to a declared
cap (default 0.10) only for the affected cohort; it does not reset every later
10--14 cohort, as sustained_high does.

The resulting contrasts are conditional on the fitted transmission model and
are expected-report contrasts, not causal estimates of real-world cases.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

# Keep imports local to the existing model implementation and posterior kernel.
from run_bayesian_transmission import (
    AGE_LABELS, COUNTRIES, DISEASES, YEARS, SIM_YEARS, build_input,
    simulate_country, DISEASE_PRIORS, SPEC, _unpack,
)

ROOT = Path(__file__).resolve().parents[2]
POST = ROOT / "output" / "revised" / "bayesian_transmission" / "posterior_samples.csv"
OUT = ROOT / "output" / "revised" / "bayesian_transmission"


def apply_one_time_catchup(gaps: np.ndarray, campaign_year: int, cap: float = .10,
                           source_band: int = 2) -> np.ndarray:
    """Apply a one-time catch-up to the 10--14 cohort and age it annually.

    The transmission model uses five-year age bands, so the affected cohort is
    represented as a moving fraction across adjacent bands. At one year after
    the campaign, 80% remains in ages 10--14 and 20% has entered ages 15--19;
    after five years the full cohort occupies the next band. This avoids an
    artificial five-year delay at the band boundary.
    """
    out = np.array(gaps, dtype=float, copy=True)
    years = np.asarray(SIM_YEARS)
    for i, year in enumerate(years):
        if year < campaign_year:
            continue
        elapsed = int(year) - int(campaign_year)
        band = source_band + elapsed // 5
        remainder = elapsed % 5
        frac_current = (5 - remainder) / 5.0
        frac_next = remainder / 5.0
        if 0 <= band < out.shape[1]:
            out[i, band] = out[i, band] + frac_current * (float(cap) - out[i, band])
        if remainder and 0 <= band + 1 < out.shape[1]:
            out[i, band + 1] = out[i, band + 1] + frac_next * (float(cap) - out[i, band + 1])
    return np.clip(out, 0.0, 1.0)


def simulate_catchup(inp: dict, disease: str, country: str, theta: np.ndarray,
                     campaign_year: int, cap: float, duration: str = "waning",
                     adult_gap: float = .30) -> np.ndarray:
    """Run the shared kernel with a one-time 10--14 catch-up profile."""
    # Reproduce simulate_country's setup, then replace only the gap profile.
    data = inp["numerical"][country]
    gaps = (data["gap_w"][disease] if duration == "waning" else data["gap_0"][disease]).copy()
    gaps[:, 3:] = adult_gap
    gaps = apply_one_time_catchup(gaps, campaign_year, cap=cap, source_band=2)
    # Calling the existing function with observed history would rebuild gaps, so
    # use its shared daily kernel directly through a temporary input profile.
    # This keeps parameter handling exactly aligned with the fitted model.
    from transmission_kernel import simulate_daily as shared_simulate_daily
    par = _unpack(theta, disease)
    q = DISEASE_PRIORS[disease]
    contact_scale = np.where(SIM_YEARS >= 2020, 1.0, 1.0)
    annual, _, _, _ = shared_simulate_daily(
        data["population"], data["births"], data["contacts"], data["contact_rho"],
        gaps, par["R0"], q["latent"], q["infectious"], par["acquisition_ratio"],
        par["clinical"], par["importation"],
        par["report"][COUNTRIES.index(country)] * par.get("reporting_multiplier", 1.0),
        contact_scale, 0.0, SPEC["diseases"][disease].get("source_reduction", 0.0))
    return annual[SIM_YEARS >= YEARS[0]]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=120, help="posterior draws per disease")
    ap.add_argument("--campaign-year", type=int, default=2018)
    ap.add_argument("--cap", type=float, default=.10, help="post-catch-up gap cap")
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--output", type=Path, default=OUT / "catchup_10_14_2018_contrasts.csv")
    args = ap.parse_args()
    if not POST.exists():
        raise FileNotFoundError(f"posterior file not found: {POST}; run fitted model first")
    inp = build_input()
    posterior = pd.read_csv(POST)
    rng = np.random.default_rng(args.seed)
    rows = []
    for disease in DISEASES:
        z = posterior[posterior.disease.eq(disease)].copy()
        if z.empty:
            continue
        idx = rng.choice(len(z), size=min(args.draws, len(z)), replace=False)
        for draw, (_, r) in enumerate(z.iloc[idx].iterrows()):
            theta = np.array([float(r.log_R0), float(r.log_nb_size)])
            for country in COUNTRIES:
                baseline = simulate_country(inp, disease, country, theta,
                                            duration="waning", history="observed")
                catchup = simulate_catchup(inp, disease, country, theta,
                                           args.campaign_year, args.cap)
                for j, year in enumerate(YEARS):
                    rows.append({
                        "disease": disease, "country": country, "draw": draw,
                        "year": int(year), "campaign": "one_time_10_14_catchup",
                        "campaign_year": args.campaign_year, "gap_cap": args.cap,
                        "baseline_expected_reports": float(baseline[j]),
                        "catchup_expected_reports": float(catchup[j]),
                        "absolute_reduction": float(baseline[j] - catchup[j]),
                        "relative_reduction_pct": float(100 * (1 - catchup[j] / max(baseline[j], 1e-12))),
                    })
    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    # First pool country-years within each posterior draw, then summarise the
    # pooled contrast across draws. Summing all rows directly would multiply
    # expected reports by the number of posterior draws.
    by_draw = (out.groupby(["disease", "campaign", "campaign_year", "gap_cap", "draw"], as_index=False)
               .agg(baseline_expected_reports=("baseline_expected_reports", "sum"),
                    catchup_expected_reports=("catchup_expected_reports", "sum"),
                    absolute_reduction=("absolute_reduction", "sum")))
    by_draw["pooled_reduction_pct"] = 100.0 * by_draw["absolute_reduction"] / by_draw["baseline_expected_reports"].clip(lower=1e-12)
    summary = (by_draw.groupby(["disease", "campaign", "campaign_year", "gap_cap"], as_index=False)
               .agg(baseline_expected_reports=("baseline_expected_reports", "median"),
                    catchup_expected_reports=("catchup_expected_reports", "median"),
                    absolute_reduction=("absolute_reduction", "median"),
                    pooled_reduction_pct=("pooled_reduction_pct", "median"),
                    q025_reduction_pct=("pooled_reduction_pct", lambda x: x.quantile(.025)),
                    q975_reduction_pct=("pooled_reduction_pct", lambda x: x.quantile(.975)),
                    n_draws=("draw", "nunique")))
    summary_path = args.output.with_name(args.output.stem + "_summary.csv")
    summary.to_csv(summary_path, index=False)
    print(f"wrote {args.output}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()

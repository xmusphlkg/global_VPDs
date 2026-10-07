#!/usr/bin/env python3
"""Likelihood-conditioned, age-structured Bayesian transmission model.

This module is separate from the older fixed-parameter epidemic potential
script. It uses reconstructed protection as a time-varying vaccination input
and estimates disease-specific transmission (R0) and observation overdispersion
from annual case records. Acquisition blocking, clinical blocking, reporting
and importation are fixed at declared disease-specific centres. Four random-walk
Metropolis chains are run for every disease.

The model is a parsimonious state-space model rather than an incidence forecast:
the input data contain national totals, not age-specific cases.  Age structure,
births, aging and a country-specific Prem contact matrix are retained; natural
infection immunity is represented by R. Diphtheria toxoid is assigned near-zero
acquisition blocking and strong clinical blocking, while pertussis has separate,
weaker acquisition and clinical effects. Mumps uses the project
all-MCV proxy and is flagged as prior-dominated when observations are sparse.

Outputs are written under output/revised/bayesian_transmission/.
"""
from __future__ import annotations

import json
import math
import os
import re
import subprocess
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence
from xml.etree import ElementTree as ET

# Keep BLAS from creating one thread per process.  The posterior jobs are
# parallelized over disease × chain; oversubscribing all cores made MAP
# initialization slower and did not improve the numerical result.
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"
os.environ.setdefault("NUMBA_NUM_THREADS", "4")

import numpy as np
import pandas as pd
from scipy.special import gammaln, logsumexp
from scipy.optimize import minimize
from scipy.stats import nbinom

from select_transmission_countries import selected_countries
from model_spec import SPEC
from transmission_kernel import (simulate_daily as shared_simulate_daily,
                                  simulate_daily_panel as shared_simulate_daily_panel)

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
COHORT = Path(os.environ.get("VPDS_COHORT_DIR", ROOT / "output" / "five_disease_revision"))
OUT = Path(os.environ.get("VPDS_BAYES_OUT", ROOT / "output" / "revised" / "bayesian_transmission"))
OUT.mkdir(parents=True, exist_ok=True)

COUNTRIES, COUNTRY_SELECTION = selected_countries()
PREM_NAMES = COUNTRY_SELECTION.set_index("country").prem_sheet.to_dict()
DISEASES = ["Measles", "Rubella", "Pertussis", "Diphtheria", "Mumps"]
SIM_YEARS = np.arange(int(SPEC["demography"]["burn_in_start"]),
                      int(SPEC["time"]["analysis_years"][1]) + 1, dtype=int)
YEARS = np.arange(int(SPEC["time"]["analysis_years"][0]),
                  int(SPEC["time"]["analysis_years"][1]) + 1, dtype=int)
CALIBRATION_YEARS = np.arange(int(YEARS[0]), int(SPEC["time"]["calibration_end"]) + 1, dtype=int)
AGE_LABELS = ["0_4", "5_9", "10_14", "15_19", "20_24", "25_29", "30_34",
              "35_39", "40_44", "45_49", "50_54", "55_59", "60_64",
              "65_69", "70_74", "75_plus"]

# Acquisition reduction is applied per unit reconstructed clinical protection
# mass. Clinical risk is derived as gap/s; it has no second vaccine coefficient.
RECALIBRATED_MAPPING = SPEC.get("protection_mapping") == "clinical_protection_mass_v2"
if not RECALIBRATED_MAPPING:
    raise ValueError("Only the revised clinical-protection mapping is supported")
DISEASE_PRIORS = {
        d: dict(R0=p["R0"], latent=p["latent_days"], infectious=p["infectious_days"],
                acquisition_ratio=float(np.clip(p["acquisition_ratio"], 1e-8, 1-1e-8)),
                clinical=p["clinical_fraction"], report=p["reporting"])
        for d, p in SPEC["diseases"].items()
    }

# Annual national totals identify only transmission intensity and observation
# overdispersion in this conditional model. All biological and reporting inputs
# therefore remain explicit fixed disease-specific centres in ``DISEASE_PRIORS``.
PARAM_NAMES = ["log_R0", "log_nb_size"]
ACTIVE_INDEX = np.arange(len(PARAM_NAMES), dtype=int)
ACTIVE_PARAM_NAMES = PARAM_NAMES.copy()


def _column_number(ref: str) -> int:
    m = re.match(r"([A-Z]+)", ref)
    n = 0
    for ch in m.group(1):
        n = n * 26 + ord(ch) - 64
    return n - 1


def _read_xlsx(path: Path) -> Dict[str, np.ndarray]:
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    relns = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(ns + "si"):
                shared.append("".join(t.text or "" for t in si.iter(ns + "t")))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rel = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap = {r.attrib["Id"]: r.attrib["Target"] for r in rel}
        out = {}
        for sh in wb.find(ns + "sheets"):
            target = relmap[sh.attrib[relns + "id"]].lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
            root = ET.fromstring(z.read(target)); rows = {}; mr = mc = 0
            for row in root.iter(ns + "row"):
                ri = int(row.attrib["r"]) - 1
                for cell in row.findall(ns + "c"):
                    ci = _column_number(cell.attrib.get("r", "A1")); v = cell.find(ns + "v")
                    val = v.text if v is not None else None
                    if val is not None and cell.attrib.get("t") == "s":
                        val = shared[int(val)]
                    elif val is not None:
                        try: val = float(val)
                        except ValueError: pass
                    rows.setdefault(ri, {})[ci] = val; mr = max(mr, ri); mc = max(mc, ci)
            arr = np.full((mr + 1, mc + 1), np.nan, dtype=object)
            for i, cols in rows.items():
                for j, val in cols.items(): arr[i, j] = val
            out[sh.attrib["name"]] = arr
        return out


def load_contact_matrices() -> Dict[str, np.ndarray]:
    paths = [RAW / "contact_matrices_prem2017" / "MUestimates_all_locations_1.xlsx",
             RAW / "contact_matrices_prem2017" / "MUestimates_all_locations_2.xlsx"]
    sheets = {}
    for p in paths: sheets.update(_read_xlsx(p))
    out = {}
    for c in COUNTRIES:
        name = PREM_NAMES[c]
        raw = sheets[name]; start = 1 if raw.shape[0] >= 17 and isinstance(raw[0, 0], str) else 0
        out[c] = np.asarray(raw[start:start + 16, :16], dtype=float)
    return out


def load_demography() -> tuple[dict, dict, dict]:
    x = pd.read_csv(RAW / "population_by_single_age.csv")
    x = x[x.country.isin(COUNTRIES) & x.year.isin(SIM_YEARS)].copy()
    x.age = pd.to_numeric(x.age, errors="coerce"); x.population = pd.to_numeric(x.population, errors="coerce")
    def band(a):
        if a >= 75: return 15
        return int(a // 5)
    x["band"] = x.age.map(band)
    pop = {}; infant_share = {}
    for c in COUNTRIES:
        for y in SIM_YEARS:
            z = x[(x.country == c) & (x.year == y)]
            v = np.zeros(16)
            for b, g in z.groupby("band"): v[int(b)] = g.population.sum()
            pop[c, int(y)] = v
            infant_population = float(z.loc[z.age.eq(0), "population"].sum())
            infant_share[c, int(y)] = infant_population / max(v[0], 1.0)
    b = pd.read_csv(RAW / "births_data.csv")
    births = {(r.country, int(r.year)): float(r.births) for r in b.itertuples()
              if r.country in COUNTRIES and int(r.year) in SIM_YEARS}
    return pop, births, infant_share


def ages_in_band(index: int) -> range:
    if index == 0:
        return range(1, 5)
    if index < 15:
        return range(index * 5, index * 5 + 5)
    return range(75, 100)


def combine_infant_gap(infant_share: float, age_1_4_gap: float) -> float:
    return infant_share + (1.0 - infant_share) * age_1_4_gap


def _gap_profiles(infant_share: dict) -> tuple[dict, dict]:
    # Mumps does not have a complete strict antigen history in the primary
    # reconstruction.  It must use the labelled all-MCV proxy file; silently
    # filling its strict missing values would erase the intended proxy boundary.
    base = pd.read_csv(COHORT / "cohort_level_all_ages.csv", low_memory=False)
    proxy = pd.read_csv(COHORT / "mumps_all_mcv_proxy_cohort_level.csv", low_memory=False)
    base = base[base.disease.ne("Mumps")]
    proxy["disease"] = "Mumps"
    x = pd.concat([base, proxy], ignore_index=True)
    x = x[x.country.isin(COUNTRIES) & x.observation_year.isin(SIM_YEARS)].copy()
    x["gap_w"] = pd.to_numeric(x.gap_waning_projected, errors="coerce")
    x["gap_0"] = pd.to_numeric(x.gap_no_waning_projected, errors="coerce")
    # Missing reconstructed states receive the same neutral fallback in both
    # duration scenarios; filling no-waning with an already averaged waning
    # band induces artificial differences even for durable-protection diseases.
    outw, out0 = {}, {}
    for c in COUNTRIES:
        for y in SIM_YEARS:
            for d in DISEASES:
                z = x[(x.country == c) & (x.observation_year == y) & (x.disease == d)].sort_values("age")
                g = np.full(16, .30); g0vct = np.full(16, .30)
                for i in range(16):
                    q = z[z.age.isin(ages_in_band(i))]
                    if len(q):
                        g[i] = float(np.average(q.gap_w.fillna(.30), weights=q.cohort_size.fillna(1)))
                        g0vct[i] = float(np.average(q.gap_0.fillna(.30), weights=q.cohort_size.fillna(1)))
                infant = infant_share[c, int(y)]
                g[0] = combine_infant_gap(infant, g[0])
                g0vct[0] = combine_infant_gap(infant, g0vct[0])
                # A missing profile remains a missing reconstruction state.  A
                # neutral value is used only for the structural simulation and
                # is audited in the manifest; it is never treated as a case zero.
                outw[c, int(y), d] = np.clip(g, 0, 1); out0[c, int(y), d] = np.clip(g0vct, 0, 1)
    return outw, out0


def load_cases() -> pd.DataFrame:
    p = ROOT / "data" / "processed" / "case_data_quality_preserved.csv"
    x = pd.read_csv(p)
    x = x[x.country.isin(COUNTRIES) & x.year.isin(YEARS) & x.disease.isin(DISEASES)].copy()
    x["cases"] = pd.to_numeric(x.cases, errors="coerce").clip(lower=0)
    # Only source rows with a numeric count enter the likelihood.  Explicit
    # reported zeros remain observations; absent source rows are missing.
    x["observed"] = x.case_status.ne("absent_from_source") & x.cases.notna()
    return x[["country", "year", "disease", "cases", "observed", "case_status"]]


def build_input() -> dict:
    pop, births, infant_share = load_demography(); w, n = _gap_profiles(infant_share); contacts = load_contact_matrices(); cases = load_cases()
    annual_contacts = {}; contact_rho = {}; contact_audit = []
    for c, matrix in contacts.items():
        for year in SIM_YEARS:
            p = np.maximum(pop[c, int(year)], 1.0)
            corrected = .5 * (matrix + matrix.T * (p[None, :] / p[:, None]))
            annual_contacts[c, int(year)] = np.clip(corrected, 0, None)
            contact_rho[c, int(year)] = max(float(np.max(np.linalg.eigvals(annual_contacts[c, int(year)]).real)), 1e-9)
    for c in COUNTRIES:
        errors = []; radii = []
        for year in SIM_YEARS:
            p = np.maximum(pop[c, int(year)], 1.0)
            matrix = annual_contacts[c, int(year)]
            flows = matrix * p[:, None]
            errors.append(float(np.max(np.abs(flows - flows.T)) / max(float(np.max(flows)), 1.0)))
            radii.append(contact_rho[c, int(year)])
        contact_audit.append(dict(country=c, corrected_spectral_radius=float(np.mean(radii)),
            corrected_max_reciprocity_error=max(errors)))
    inp = dict(pop=pop, births=births, gap_w=w, gap_0=n, contacts=annual_contacts,
               contact_rho=contact_rho, contact_audit=contact_audit, cases=cases)
    inp["numerical"] = {}
    for c in COUNTRIES:
        inp["numerical"][c] = {
            "population": np.stack([pop[c, int(y)] for y in SIM_YEARS]),
            "births": np.array([births.get((c, int(y)), 0.0) for y in SIM_YEARS]),
            "contacts": np.stack([annual_contacts[c, int(y)] for y in SIM_YEARS]),
            "contact_rho": np.array([contact_rho[c, int(y)] for y in SIM_YEARS]),
            "gap_w": {d: np.stack([w[c, int(y), d] for y in SIM_YEARS]) for d in DISEASES},
            "gap_0": {d: np.stack([n[c, int(y), d] for y in SIM_YEARS]) for d in DISEASES},
        }
    inp["panel_numerical"] = {
        "population": np.stack([inp["numerical"][c]["population"] for c in COUNTRIES]),
        "births": np.stack([inp["numerical"][c]["births"] for c in COUNTRIES]),
        "contacts": np.stack([inp["numerical"][c]["contacts"] for c in COUNTRIES]),
        "contact_rho": np.stack([inp["numerical"][c]["contact_rho"] for c in COUNTRIES]),
        "gap_w": {d: np.stack([inp["numerical"][c]["gap_w"][d] for c in COUNTRIES]) for d in DISEASES},
        "gap_0": {d: np.stack([inp["numerical"][c]["gap_0"][d] for c in COUNTRIES]) for d in DISEASES},
    }
    return inp


def _unpack(theta: np.ndarray, d: str, overrides: Mapping[str, object] | None = None) -> dict:
    q = DISEASE_PRIORS[d]
    values = {"R0": float(np.exp(theta[0])), "nb_size": float(np.exp(theta[1])),
              "acquisition_ratio": q["acquisition_ratio"],
              "clinical": q["clinical"], "importation": 2e-4 / 365.0,
              "report": np.full(len(COUNTRIES), q["report"], dtype=float),
              "reporting_multiplier": 1.0}
    if overrides:
        values.update(overrides)
    return values


def _log_prior(theta: np.ndarray, d: str) -> float:
    q = DISEASE_PRIORS[d]
    def normal(x, m, s): return -.5 * ((x - m) / s) ** 2 - math.log(s * math.sqrt(2 * math.pi))
    lp = normal(theta[0], math.log(q["R0"]), .35)
    return float(lp + normal(theta[1], math.log(25), .8))


def _nb_loglik(y: float, mu: float, size: float) -> float:
    mu = max(float(mu), 1e-12); size = max(float(size), 1e-4); y = float(y)
    return float(gammaln(y + size) - gammaln(size) - gammaln(y + 1) + size * math.log(size / (size + mu)) + y * math.log(mu / (size + mu)))


def simulate_country(inp: dict, d: str, c: str, theta: np.ndarray, duration="waning", history="observed", stochastic=False, rng=None, recovery_year=None, adult_gap=.30, high_history_cap=.10, covid_contact_scale=1.0, natural_immunity_years=0.0, source_reduction=None, fixed_overrides: Mapping[str, object] | None = None) -> np.ndarray:
    par = _unpack(theta, d, fixed_overrides); q = DISEASE_PRIORS[d]
    data = inp["numerical"][c]
    gaps = (data["gap_w"][d] if duration == "waning" else data["gap_0"][d]).copy()
    # Reconstructed profiles stop at age 14; use the declared structural
    # adult gap for 15+ in the observed scenario. The high-history contrast
    # changes only the 10–14-year child band, preserving 0–9 years and adults.
    gaps[:, 3:] = adult_gap
    if history in ("sustained_high", "recovery_5y"):
        high = gaps.copy()
        analysis = SIM_YEARS >= YEARS[0]
        high[analysis, 2] = np.minimum(high[analysis, 2], high_history_cap)
        if history == "sustained_high":
            gaps = high
        elif recovery_year is not None:
            for i, year in enumerate(SIM_YEARS):
                if year >= recovery_year:
                    fraction = min(1.0, (year - recovery_year + 1) / 5.0)
                    gaps[i] = gaps[i] * (1.0 - fraction) + high[i] * fraction
    if stochastic:
        raise ValueError("Daily stochastic replicates are not generated by the expected-state MCMC kernel")
    contact_scale = np.where(SIM_YEARS >= 2020, covid_contact_scale, 1.0)
    annual, _, _, _ = shared_simulate_daily(
        data["population"], data["births"], data["contacts"], data["contact_rho"],
        gaps, par["R0"], q["latent"], q["infectious"], par["acquisition_ratio"],
        par["clinical"], par["importation"],
        par["report"][COUNTRIES.index(c)] * par.get("reporting_multiplier", 1.0),
        contact_scale, natural_immunity_years,
        SPEC["diseases"][d].get("source_reduction", 0.0) if source_reduction is None else source_reduction)
    return annual[SIM_YEARS >= YEARS[0]]


def simulate_panel(inp: dict, d: str, theta: np.ndarray, duration="waning",
                   history="observed", recovery_year=None, adult_gap=.30,
                   high_history_cap=.10, covid_contact_scale=1.0,
                   natural_immunity_years=0.0, source_reduction=None,
                   fixed_overrides: Mapping[str, object] | None = None) -> np.ndarray:
    """Evaluate all selected countries in one call to the shared kernel."""
    par = _unpack(theta, d, fixed_overrides)
    data = inp["panel_numerical"]
    gaps = (data["gap_w"][d] if duration == "waning" else data["gap_0"][d]).copy()
    gaps[:, :, 3:] = adult_gap
    if history in ("sustained_high", "recovery_5y"):
        high = gaps.copy()
        analysis = SIM_YEARS >= YEARS[0]
        high[:, analysis, 2] = np.minimum(high[:, analysis, 2], high_history_cap)
        if history == "sustained_high":
            gaps = high
        elif recovery_year is not None:
            for i, year in enumerate(SIM_YEARS):
                if year >= recovery_year:
                    fraction = min(1.0, (year - recovery_year + 1) / 5.0)
                    gaps[:, i, :] = gaps[:, i, :] * (1.0 - fraction) + high[:, i, :] * fraction
    contact_scale = np.where(SIM_YEARS >= 2020, covid_contact_scale, 1.0)
    report = np.asarray(par["report"], dtype=float) * float(par.get("reporting_multiplier", 1.0))
    annual = shared_simulate_daily_panel(
        data["population"], data["births"], data["contacts"], data["contact_rho"],
        gaps, par["R0"], DISEASE_PRIORS[d]["latent"], DISEASE_PRIORS[d]["infectious"],
        par["acquisition_ratio"], par["clinical"], par["importation"], report, contact_scale,
        natural_immunity_years,
        SPEC["diseases"][d].get("source_reduction", 0.0) if source_reduction is None else source_reduction)
    return annual[:, SIM_YEARS >= YEARS[0]]

def log_posterior(theta: np.ndarray, inp: dict, d: str, duration="waning", history="observed") -> float:
    if not np.all(np.isfinite(theta)): return -np.inf
    lp = _log_prior(theta, d)
    if not np.isfinite(lp): return -np.inf
    par = _unpack(theta, d); cases = inp["cases"]
    predictions = simulate_panel(inp, d, theta, duration=duration, history=history)
    for ci, c in enumerate(COUNTRIES):
        pred = predictions[ci]
        z = cases[(cases.country == c) & (cases.disease == d) & cases.observed &
                  cases.year.isin(CALIBRATION_YEARS)]
        for r in z.itertuples(): lp += _nb_loglik(r.cases, pred[int(r.year) - int(YEARS[0])], par["nb_size"])
    return float(lp)


def initial_theta(d: str) -> np.ndarray:
    q = DISEASE_PRIORS[d]
    return np.array([math.log(q["R0"]), math.log(25)], dtype=float)


def find_map(inp: dict, d: str) -> dict:
    """Profile several dispersed starts and preserve distinct optima."""
    base = initial_theta(d)
    def objective(active: np.ndarray) -> float:
        theta = base.copy(); theta[ACTIVE_INDEX] = active
        value = -log_posterior(theta, inp, d)
        return float(value) if np.isfinite(value) else 1e100
    starts = []
    for offset in (-1.0, -.5, 0.0, .5, 1.0):
        initial = base[ACTIVE_INDEX].copy(); initial[0] += offset
        result = minimize(objective, initial, method="Nelder-Mead",
            options={"maxiter": 180, "xatol": 0.01, "fatol": 0.05})
        value = float(objective(result.x)) if np.all(np.isfinite(result.x)) else 1e100
        starts.append({"start": initial.tolist(), "optimum": result.x.tolist(),
            "objective": value, "success": bool(result.success)})
    best = min(starts, key=lambda row: row["objective"])
    # The profile audit retains alternate modes, but the primary posterior is
    # a single highest-density basin.  Starting chains in different modes and
    # pooling them would produce a false R-hat failure rather than a valid
    # mixture posterior for this random-walk sampler.
    best_log_r0 = float(best["optimum"][0])
    alternates = [row for row in starts if abs(float(row["optimum"][0]) - best_log_r0) > .25]
    if alternates:
        alternate = min(alternates, key=lambda row: row["objective"])
        # Keep the primary chain in the high-density basin.  The midpoint can
        # still include the shoulder of a nearby alternative mode (as seen for
        # Pertussis), so place the guard one quarter of the way from the best
        # mode toward that alternate.
        boundary = .75 * best_log_r0 + .25 * float(alternate["optimum"][0])
        mode_bounds = (boundary, np.inf) if boundary < best_log_r0 else (-np.inf, boundary)
    else:
        mode_bounds = (-np.inf, np.inf)
    return {"best": np.asarray(best["optimum"]),
            "chain_starts": np.repeat(np.asarray(best["optimum"])[None, :], 4, axis=0),
            "mode_bounds": mode_bounds,
            "profile_audit": starts}


def run_chain(args):
    inp, d, chain, n_iter, burn, seed, duration, history, map_active, mode_bounds = args
    rng = np.random.default_rng(seed + chain * 1009)
    theta = initial_theta(d)
    theta[ACTIVE_INDEX] = map_active
    # All primary chains share the selected basin; independent jitter and
    # burn-in adaptation give each chain its own proposal distribution.
    theta[ACTIVE_INDEX] += rng.normal(0, .01, len(ACTIVE_INDEX))
    current = log_posterior(theta, inp, d, duration, history)
    scale = np.full(len(ACTIVE_INDEX), .06); draws = []; accepted = 0; acc_window = 0
    proposal_history = [theta[ACTIVE_INDEX].copy()]; chol = np.diag(scale)
    for it in range(n_iter):
        # Adapt the two-parameter block proposal only during burn-in.
        if it < burn and it >= 100 and (it + 1) % 50 == 0 and len(proposal_history) >= 50:
            cov = np.cov(np.asarray(proposal_history).T) + np.eye(len(ACTIVE_INDEX)) * 1e-4
            chol = np.linalg.cholesky(cov) * (2.38 / math.sqrt(len(ACTIVE_INDEX)))
        prop = theta.copy()
        prop[ACTIVE_INDEX] += chol @ rng.normal(size=len(ACTIVE_INDEX))
        val = (-np.inf if prop[0] < mode_bounds[0] or prop[0] > mode_bounds[1]
               else log_posterior(prop, inp, d, duration, history))
        if math.log(rng.random()) < val - current:
            theta, current = prop, val; accepted += 1; acc_window += 1
        proposal_history.append(theta[ACTIVE_INDEX].copy())
        if it < burn and (it + 1) % 50 == 0:
            ar = acc_window / 50; scale *= np.exp(np.clip(ar - .25, -.25, .25)); acc_window = 0
            if it < 100: chol = np.diag(scale)
            else: chol *= np.exp(np.clip(ar - .25, -.25, .25))
        if it >= burn: draws.append(np.r_[theta, current])
    return np.asarray(draws), accepted / n_iter


def find_map_task(args):
    inp, d = args
    return d, find_map(inp, d)


def posterior_scenarios(inp: dict, d: str, samples: np.ndarray, rng: np.random.Generator, n_draws: int = 120) -> tuple[pd.DataFrame, pd.DataFrame]:
    if len(samples) > n_draws: samples = samples[rng.choice(len(samples), n_draws, replace=False)]
    rows = []; pp = []
    scenarios = [("observed", "waning", None), ("observed", "unwaning", None),
                 ("sustained_high", "waning", None), ("sustained_high", "unwaning", None),
                 ("recovery_5y", "waning", 2018)]
    for k, th in enumerate(samples):
        for history, duration, ry in scenarios:
            for c in COUNTRIES:
                deterministic = simulate_country(inp, d, c, th, duration, history, False, rng, ry)
                # These are expected reports; only the observation layer adds
                # negative-binomial predictive variation.
                nb_size = _unpack(th, d)["nb_size"]
                nb_prob = nb_size / (nb_size + np.maximum(deterministic, 1e-12))
                nb_cases = rng.negative_binomial(nb_size, np.clip(nb_prob, 1e-12, 1 - 1e-12))
                for yi, year in enumerate(YEARS):
                    pp.append(dict(draw=k, disease=d, country=c, year=int(year), history=history, duration=duration,
                                   expected_reports=float(deterministic[yi]),
                                   nb_predictive_reports=float(nb_cases[yi]), nb_size=nb_size))
                rows.append(dict(draw=k, disease=d, country=c, history=history, duration=duration,
                                 recovery_year=ry or np.nan, total_cases=float(deterministic.sum()),
                                 expected_reports_per_capita=float(deterministic.sum() / max(inp["pop"][c, 2014].sum(), 1))))
    return pd.DataFrame(rows), pd.DataFrame(pp)


def posterior_scenarios_task(args):
    inp, d, samples, seed, n_draws = args
    return posterior_scenarios(inp, d, samples, np.random.default_rng(seed), n_draws)


def sensitivity_analysis(inp: dict, d: str, samples: np.ndarray, n_draws: int = 6) -> pd.DataFrame:
    """Evaluate structural settings as paired Figure 3 scenario contrasts."""
    indices = np.linspace(0, len(samples) - 1, min(n_draws, len(samples)), dtype=int)
    selected = samples[indices]
    settings = []
    adults = SPEC["scenarios"]["adult_gap"]
    settings.extend(("adult_gap", x, x, .10, 1.0, 0.0, {}) for x in adults)
    settings.extend(("high_history_cap", x, .30, x, 1.0, 0.0, {}) for x in (0.0, .10, .20))
    settings.extend(("covid_contact_scale_2020_2023", x, .30, .10, x, 0.0, {}) for x in (.50, .75, 1.0))
    settings.extend(("natural_immunity_duration_years", x, .30, .10, 1.0, x, {}) for x in (0.0, 15.0, 30.0))
    settings.extend(("acquisition_ratio_offset", x, .30, .10, 1.0, 0.0, {"acquisition_ratio_offset": x}) for x in (-.05, 0.0, .05))
    if d == "Diphtheria":
        settings.extend(("source_reduction_ratio", x, .30, .10, 1.0, 0.0, {"source_reduction": x})
                        for x in SPEC["diseases"][d]["source_reduction_grid"])
    settings.extend(("clinical_fraction_offset", x, .30, .10, 1.0, 0.0, {"clinical_fraction_offset": x}) for x in (-.05, 0.0, .05))
    settings.extend(("importation_multiplier", x, .30, .10, 1.0, 0.0, {"importation_multiplier": x}) for x in (.5, 1.0, 2.0))
    settings.extend(("reporting_multiplier", x, .30, .10, 1.0, 0.0, {"reporting_multiplier": x}) for x in (.5, 1.0, 2.0))
    scenario_specs = {
        "observed_waning": ("observed", "waning", None),
        "observed_unwaning": ("observed", "unwaning", None),
        "high_waning": ("sustained_high", "waning", None),
        "high_unwaning": ("sustained_high", "unwaning", None),
        "recovery_waning": ("recovery_5y", "waning", 2018),
    }
    rows = []
    for draw, theta in enumerate(selected):
        base = _unpack(theta, d)
        for name, level, adult_gap, high_cap, covid_scale, immune_years, change in settings:
            overrides = {}
            if "acquisition_ratio_offset" in change:
                overrides["acquisition_ratio"] = float(np.clip(base["acquisition_ratio"] + change["acquisition_ratio_offset"], 1e-5, 1 - 1e-5))
            if "clinical_fraction_offset" in change:
                overrides["clinical"] = float(np.clip(base["clinical"] + change["clinical_fraction_offset"], 1e-5, 1.0))
            if "importation_multiplier" in change:
                overrides["importation"] = base["importation"] * change["importation_multiplier"]
            if "reporting_multiplier" in change:
                overrides["reporting_multiplier"] = change["reporting_multiplier"]
            values = {}
            for scenario, (history, duration, recovery_year) in scenario_specs.items():
                totals = [simulate_country(inp, d, c, theta, duration, history,
                    recovery_year=recovery_year, adult_gap=adult_gap, high_history_cap=high_cap,
                    covid_contact_scale=covid_scale, natural_immunity_years=immune_years,
                    source_reduction=change.get("source_reduction"), fixed_overrides=overrides).sum()
                    for c in COUNTRIES]
                values[scenario] = float(np.mean(totals))
            observed = values["observed_waning"]
            row = dict(draw=draw, disease=d, sensitivity=name, level=level,
                       history="observed", duration="waning",
                       total_expected_reports=observed * len(COUNTRIES),
                       equal_country_mean_reports=observed, n_countries=len(COUNTRIES),
                       observed_waning_reports=observed,
                       observed_unwaning_reports=values["observed_unwaning"],
                       high_waning_reports=values["high_waning"],
                       high_unwaning_reports=values["high_unwaning"],
                       recovery_waning_reports=values["recovery_waning"])
            row["high_history_reduction_pct"] = 100.0 * (1.0 - values["high_waning"] / max(observed, 1e-12))
            row["recovery_reduction_pct"] = 100.0 * (1.0 - values["recovery_waning"] / max(observed, 1e-12))
            row["duration_contrast_pct"] = 100.0 * (values["observed_unwaning"] / max(observed, 1e-12) - 1.0)
            rows.append(row)
    return pd.DataFrame(rows)


def sensitivity_analysis_task(args):
    return sensitivity_analysis(*args)


def holdout_diagnostics(inp: dict, posterior: pd.DataFrame, n_draws: int = 160) -> pd.DataFrame:
    """Posterior predictive time holdout audit (2020–23).

    The posterior is fit only through 2019. These 2020–2023 predictions are
    generated without using their case observations in calibration.
    """
    rows = []
    for d in DISEASES:
        z = posterior[posterior.disease.eq(d)]
        draws = z.sample(min(n_draws, len(z)), random_state=20261006)
        for c in COUNTRIES:
            yobs = inp["cases"][(inp["cases"].country == c) & (inp["cases"].disease == d) & (inp["cases"].year >= 2020) & inp["cases"].observed]
            if yobs.empty: continue
            pred = []; pred_nb = []; nb_sizes = []
            for r in draws.itertuples(index=False):
                th = np.array([getattr(r, p) for p in PARAM_NAMES], dtype=float)
                sim = simulate_country(inp, d, c, th, "waning", "observed", False,
                                       np.random.default_rng(20261006))
                pred.append(sim[6:10])
                size = math.exp(float(r.log_nb_size)); nb_sizes.append(size)
                prob = size / (size + np.maximum(sim[6:10], 1e-12))
                pred_nb.append(np.random.default_rng(20261006 + int(r.draw)).negative_binomial(size, prob))
            pp = np.asarray(pred)
            pp_nb = np.asarray(pred_nb); sizes = np.asarray(nb_sizes)
            obs = np.array([float(yobs[yobs.year.eq(y)].cases.iloc[0]) if (yobs.year == y).any() else np.nan for y in range(2020, 2024)])
            for j, year in enumerate(range(2020, 2024)):
                if np.isfinite(obs[j]):
                    expected_lo, expected_hi = np.quantile(pp[:, j], [.025, .975])
                    nb_lo, nb_hi = np.quantile(pp_nb[:, j], [.025, .975])
                    log_mass = nbinom.logpmf(obs[j], sizes, sizes / (sizes + np.maximum(pp[:, j], 1e-12)))
                    rows.append(dict(disease=d, country=c, year=year, observed_cases=obs[j],
                                     predictive_median=float(np.nanmedian(pp[:, j])),
                                     predictive_q025=float(expected_lo), predictive_q975=float(expected_hi),
                                     expected_interval_covered=bool(expected_lo <= obs[j] <= expected_hi),
                                     nb_predictive_median=float(np.median(pp_nb[:, j])),
                                     nb_predictive_q025=float(nb_lo), nb_predictive_q975=float(nb_hi),
                                     nb_predictive_covered=bool(nb_lo <= obs[j] <= nb_hi),
                                     predictive_log_score=float(logsumexp(log_mass) - math.log(len(log_mass))),
                                     holdout_type="time_out_2020_2023"))
    return pd.DataFrame(rows)


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--iterations", type=int, default=12000); ap.add_argument("--burn", type=int, default=2000)
    ap.add_argument("--chains", type=int, default=4); ap.add_argument("--draws", type=int, default=300); ap.add_argument("--workers", type=int, default=min(20, os.cpu_count() or 1)); ap.add_argument("--seed", type=int, default=20261006)
    args = ap.parse_args(); inp = build_input(); all_post = []; acceptance = []; all_sc = []; all_pp = []
    # Run disease × chain jobs in one process pool.  The previous implementation
    # parallelized only the four chains within one disease and left four fifths
    # of the available cores idle while the other diseases waited.  Each job is
    # independent conditional on the shared input object, so this is an exact
    # parallelization of the posterior calculation rather than a change in the
    # statistical target.
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=min(args.workers, len(DISEASES))) as ex:
            maps = dict(ex.map(find_map_task, [(inp, d) for d in DISEASES]))
    else:
        maps = {d: find_map(inp, d) for d in DISEASES}
    tasks = [(inp, d, ch, args.iterations, args.burn, args.seed, "waning", "observed",
              maps[d]["chain_starts"][ch % len(maps[d]["chain_starts"])], maps[d]["mode_bounds"])
             for d in DISEASES for ch in range(args.chains)]
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=min(args.workers, len(tasks))) as ex:
            job_results = list(ex.map(run_chain, tasks))
    else:
        job_results = [run_chain(t) for t in tasks]
    by_disease = {d: job_results[i * args.chains:(i + 1) * args.chains] for i, d in enumerate(DISEASES)}
    scenario_jobs = []
    for d in DISEASES:
        result = by_disease[d]
        full_chains = np.stack([r[0][:, :len(PARAM_NAMES)] for r in result]); samples = full_chains.reshape(-1, len(PARAM_NAMES))
        acceptance.extend(dict(disease=d, chain=ch, acceptance=result[ch][1])
                          for ch in range(args.chains))
        for ch in range(args.chains):
            for it, row in enumerate(full_chains[ch]): all_post.append(dict(disease=d, chain=ch, draw=it, **{p: float(row[i]) for i, p in enumerate(PARAM_NAMES)}))
        scenario_jobs.append((inp, d, samples, args.seed + 991 + len(scenario_jobs), args.draws))
    # Scenario simulations are independent by disease and retain their own
    # random seeds.  Running them in parallel makes the stochastic interval
    # calculation use the available cores instead of serially delaying output.
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=min(args.workers, len(scenario_jobs))) as ex:
            scenario_results = list(ex.map(posterior_scenarios_task, scenario_jobs))
    else:
        scenario_results = [posterior_scenarios_task(x) for x in scenario_jobs]
    all_sc.extend(x[0] for x in scenario_results); all_pp.extend(x[1] for x in scenario_results)
    sensitivity_jobs = [(inp, d,
        np.stack([r[0][:, :len(PARAM_NAMES)] for r in by_disease[d]]).reshape(-1, len(PARAM_NAMES)))
        for d in DISEASES]
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=min(args.workers, len(sensitivity_jobs))) as ex:
            sensitivities = list(ex.map(sensitivity_analysis_task, sensitivity_jobs))
    else:
        sensitivities = [sensitivity_analysis_task(x) for x in sensitivity_jobs]
    post = pd.DataFrame(all_post); post.to_csv(OUT / "posterior_samples.csv", index=False)
    pd.DataFrame(acceptance).to_csv(OUT / "chain_acceptance.csv", index=False)
    subprocess.run(["Rscript", str(Path(__file__).with_name("summarize_transmission_diagnostics.R")),
                    str(OUT / "posterior_samples.csv"), str(OUT / "chain_acceptance.csv"),
                    str(OUT / "diagnostics.csv")], check=True)
    (OUT / "map_multistart_audit.json").write_text(json.dumps({d: {
        "best": maps[d]["best"].tolist(), "starts": maps[d]["profile_audit"],
        "primary_mode_bounds_log_R0": list(maps[d]["mode_bounds"])
    } for d in DISEASES}, indent=2) + "\n", encoding="utf-8")
    sc = pd.concat(all_sc, ignore_index=True); pp = pd.concat(all_pp, ignore_index=True); sc.to_csv(OUT / "scenario_contrasts.csv", index=False); pp.to_csv(OUT / "posterior_predictive.csv", index=False)
    # Quantiles are the interval-bearing scientific output used by the paper;
    # the draw-level table remains available for posterior contrasts.
    summary = (sc.groupby(["disease", "country", "history", "duration"], as_index=False)
               .agg(total_cases_median=("total_cases", "median"),
                    total_cases_q025=("total_cases", lambda z: z.quantile(.025)),
                    total_cases_q975=("total_cases", lambda z: z.quantile(.975)),
                    expected_reports_per_capita_median=("expected_reports_per_capita", "median"),
                    expected_reports_per_capita_q025=("expected_reports_per_capita", lambda z: z.quantile(.025)),
                    expected_reports_per_capita_q975=("expected_reports_per_capita", lambda z: z.quantile(.975)),
                    n_posterior_draws=("draw", "nunique")))
    summary.to_csv(OUT / "scenario_summary.csv", index=False)
    holdout = holdout_diagnostics(inp, post)
    holdout.to_csv(OUT / "time_holdout_posterior_predictive.csv", index=False)
    holdout_metrics = (holdout.groupby("disease", as_index=False).apply(lambda z: pd.Series({
        "n_observations": len(z), "MAE": np.mean(np.abs(z.observed_cases - z.predictive_median)),
        "RMSE": np.sqrt(np.mean((z.observed_cases - z.predictive_median) ** 2)),
        "mean_log_score": z.predictive_log_score.mean(),
        "expected_interval_coverage_pct": 100 * z.expected_interval_covered.mean(),
        "nb_predictive_coverage_pct": 100 * z.nb_predictive_covered.mean()}), include_groups=False)
        .reset_index())
    holdout_metrics.to_csv(OUT / "time_holdout_metrics.csv", index=False)
    pd.concat(sensitivities, ignore_index=True).to_csv(OUT / "scenario_sensitivity.csv", index=False)
    ppc = (pp.query("history == 'observed' and duration == 'waning'")
        .groupby(["country", "disease", "year"], as_index=False)
        .agg(pred_lo=("nb_predictive_reports", lambda z: z.quantile(.025)),
             pred_median=("nb_predictive_reports", "median"),
             pred_hi=("nb_predictive_reports", lambda z: z.quantile(.975)),
             n_draws=("draw", "nunique")))
    observed = inp["cases"].query("observed")[['country', 'disease', 'year', 'cases']]
    ppc = ppc.merge(observed, on=["country", "disease", "year"], how="inner", validate="one_to_one")
    ppc["covered_95"] = (ppc.cases >= ppc.pred_lo) & (ppc.cases <= ppc.pred_hi)
    ppc.to_csv(OUT / "posterior_predictive_checks.csv", index=False)
    audit = inp["cases"].query("country in @COUNTRIES and year in @YEARS").groupby(["country", "disease"], as_index=False).agg(n_rows=("cases", "size"), n_observed=("observed", "sum"), n_positive=("cases", lambda z: int((z > 0).sum())))
    audit["mumps_proxy"] = audit.disease.eq("Mumps"); audit.to_csv(OUT / "input_audit.csv", index=False)
    pd.DataFrame(inp["contact_audit"]).to_csv(OUT / "contact_matrix_audit.csv", index=False)
    manifest = dict(model="16-age daily expected-state SEIRV; annual NB observation likelihood", countries=COUNTRIES, country_selection="output/revised/transmission_country_selection/selection_manifest.json", diseases=DISEASES, years=YEARS.tolist(), burn_in_years=[int(SIM_YEARS[0]), int(YEARS[0] - 1)], calibration_years=CALIBRATION_YEARS.tolist(), validation_years=[2020, 2021, 2022, 2023], time_step_days=1, seed=args.seed, seed_policy="all stochastic streams derive deterministically from this root seed", beta_definition="R0 * daily recovery probability / contact spectral radius", demography="WPP annual age-band anchoring; proportional state scaling represents deaths and net migration in aggregate", chains=args.chains, iterations=args.iterations, burn=args.burn, posterior_draws=args.draws, max_workers=min(args.workers, len(tasks)), parameters=PARAM_NAMES, estimated_parameters=ACTIVE_PARAM_NAMES, fixed_nuisance_parameters=["acquisition_ratio", "clinical_fraction", "importation", "reporting"] , nuisance_policy="acquisition protection, clinical fraction, importation and country reporting are fixed at disease-specific centres because annual totals do not identify them separately from transmission", map_initialization={d: maps[d]["best"].tolist() for d in DISEASES}, chain_mode_policy="all primary chains start in the best profiled basin and proposals are restricted to the midpoint boundary separating it from the nearest alternate profile mode; alternate modes are retained in map_multistart_audit and not pooled", mode_bounds={d: list(maps[d]["mode_bounds"]) for d in DISEASES}, map_multistart_audit="map_multistart_audit.json", outputs=["posterior_samples.csv", "chain_acceptance.csv", "diagnostics.csv", "posterior_predictive.csv", "scenario_contrasts.csv", "scenario_summary.csv", "scenario_sensitivity.csv", "input_audit.csv", "contact_matrix_audit.csv", "time_holdout_posterior_predictive.csv", "time_holdout_metrics.csv", "posterior_predictive_checks.csv", "map_multistart_audit.json"], disease_priors=DISEASE_PRIORS, prior_widths={"log_R0": .35, "log_nb_size": .8}, natural_immunity="2000-2013 daily burn-in; age-specific historical infection records are unavailable", vaccination="cohort-derived gap state; acquisition and clinical effects are fixed at disease-specific centres", scenarios="observed versus matched-present high-history profile; 10-14-year child protection only; waning versus no-waning; five-year profile interpolation from 2018", likelihood="negative binomial on reported 2014-2019 totals only; 2020-2023 held out", limitations="posterior is conditional on fixed biological/reporting/importation parameters and scenario assumptions; validation must pass before figure release", contact_matrix="Prem et al. 2017 all-location matrices, reciprocity corrected for each WPP year", scenario_uncertainty="posterior expected-report contrasts and negative-binomial observation predictive replicates; no process-noise interval", sensitivity_output="scenario_sensitivity.csv varies structural assumptions and fixed centres on paired posterior draws")
    manifest.update({"model": "16-age daily expected-state SEIRV; annual NB observation likelihood",
        "protection_mapping": SPEC.get("protection_mapping", "legacy_index_mapping"),
        "parameter_specification": os.environ.get("VPDS_TRANSMISSION_SPEC", "scripts/02_models/transmission_model_spec.json"),
        "source_reduction": {d: p.get("source_reduction", 0.0) for d, p in SPEC["diseases"].items()},
        "burn_in_years": [int(SIM_YEARS[0]), int(YEARS[0] - 1)],
        "calibration_years": CALIBRATION_YEARS.tolist(), "validation_years": [2020, 2021, 2022, 2023],
        "time_step_days": 1,
        "seed": args.seed,
        "seed_policy": "all stochastic streams derive deterministically from this root seed",
        "beta_definition": "R0 * daily recovery probability / contact spectral radius",
        "demography": "WPP annual age-band anchoring; proportional state scaling represents deaths and net migration in aggregate",
        "map_initialization": {d: maps[d]["best"].tolist() for d in DISEASES},
        "map_multistart_audit": "map_multistart_audit.json",
        "natural_immunity": "2000-2013 daily burn-in; age-specific historical infection records are unavailable",
        "likelihood": "negative binomial on 2014-2019 reported annual totals only; 2020-2023 held out",
        "contact_matrix": "Prem et al. 2017 all-location matrices, reciprocity corrected by WPP age population each year",
        "scenario_uncertainty": "posterior expected-report contrasts and negative-binomial observation predictive replicates; no process-noise interval",
        "outputs": ["posterior_samples.csv", "chain_acceptance.csv", "diagnostics.csv", "posterior_predictive.csv", "scenario_contrasts.csv", "scenario_summary.csv", "scenario_sensitivity.csv", "input_audit.csv", "contact_matrix_audit.csv", "time_holdout_posterior_predictive.csv", "time_holdout_metrics.csv", "posterior_predictive_checks.csv", "map_multistart_audit.json"]})
    (OUT / "model_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    (OUT / "README.md").write_text(f"# Repaired Bayesian transmission layer\n\nThe restored model uses a daily solver, 2000-2013 infection-immunity burn-in, WPP annual age-population anchoring, and 2014-2019 calibration. Years 2020-2023 are reserved for time-out validation. The discrete transmission coefficient is normalized so the disease-free next-generation number equals the sampled R0.\n\nOnly R0 and negative-binomial size are sampled. Vaccine effects, clinical fractions, importation and reporting remain fixed at disease-specific prior centres, so report estimates are conditional scenarios, not fully identified burden posteriors. `scenario_sensitivity.csv` varies adult gap, high-history cap, 2020-2023 contact scale, natural-immunity duration, vaccine effects, clinical fraction, importation, and reporting on paired posterior draws. These grids quantify conditional variation; they do not resolve non-identifiability or replace external validation. Seasonality remains absent.\n")
    print("wrote", OUT)


if __name__ == "__main__": main()

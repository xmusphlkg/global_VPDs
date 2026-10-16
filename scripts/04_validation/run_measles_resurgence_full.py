#!/usr/bin/env python3
"""EU/EEA measles-resurgence extension.

Layers built here (model fitting is in the companion R scripts):
  1. a single-age durable-protection gap for every country-year (MCV2 set to
     zero before national introduction) and age-band profiles to 30–39 years;
  2. a 2016–2026 country-month panel with lagged cohort gaps, lagged annual
     WUENIC coverage and an EU/EEA external-pressure proxy;
  3. a case-age panel from ECDC TESSy case-based age distributions (EU/EEA,
     2023–2025);
  4. case vaccination status with the national MCV1 coverage of the same
     ages, for the screening-method comparison;
  5. conditional measles R_e contrasts for age-targeted catch-up histories.

The analysis is restricted to the 30 EU/EEA countries reporting to ECDC so
that burden ranking, cohort gaps, case ages and case vaccination status come
from one surveillance system. The outputs are ecological and exploratory.
External pressure is a lagged EU/EEA incidence proxy, not a direct measurement of importation. The 2026
monthly series is provisional and uses 2025 population as the denominator.
"""
from __future__ import annotations

import io
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "revised" / "measles_resurgence_full"
OUT.mkdir(parents=True, exist_ok=True)

EPI_FILE = ROOT / "data" / "raw" / "measles_rubella_epi_curve_data.csv"
POP_FILE = ROOT / "data" / "raw" / "population_by_single_age.csv"
CROSSWALK_FILE = ROOT / "data" / "processed" / "country_code_crosswalk.csv"
COVERAGE_XLSX = ROOT / "data" / "raw" / "who_immunization_portal_coverage.xlsx"
INTRO_FILE = ROOT / "data" / "raw" / "vaccine_introduction.csv"
CONTACT_DIR = ROOT / "output" / "revised" / "age_seirv" / "contact_matrices_reciprocal"

# EU/EEA (ECDC reporting countries). Liechtenstein is not a WHO member and has
# no WHO monthly series, so the monthly layers contain 29 countries.
EU_EEA = ["AUT", "BEL", "BGR", "CYP", "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA", "GRC", "HRV", "HUN", "IRL",
          "ISL", "ITA", "LIE", "LTU", "LUX", "LVA", "MLT", "NLD", "NOR", "POL", "PRT", "ROU", "SVK", "SVN", "SWE"]
FORMAL_YEARS = [2024, 2025, 2026]
ARCHIVE_YEARS = list(range(2016, 2027))
EXPOSURE_YEARS = list(range(2014, 2026))
VE1, VE2 = 0.93, 0.97
# The reconstruction reaches age 44 in 2025 (MCV1 from 1980); 30–39 is the
# oldest band that is complete in every analysis year.
AGE_BAND_RANGES = {"1_4": (1, 4), "5_9": (5, 9), "10_14": (10, 14), "15_19": (15, 19),
                   "20_29": (20, 29), "30_39": (30, 39)}
AGE_BANDS = ["0_0"] + list(AGE_BAND_RANGES)
AGE_MODEL_BANDS = list(AGE_BAND_RANGES)
MAX_AGE = 49
MIN_BAND_COMPLETENESS = 0.90


def require_columns(frame: pd.DataFrame, columns: set[str], name: str) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise ValueError(f"{name} missing columns: {sorted(missing)}")


def load_extended_coverage() -> pd.DataFrame:
    out = OUT / "measles_mcv_wuenic_1980_2025.csv"
    if out.exists():
        return pd.read_csv(out)
    from openpyxl import load_workbook

    wb = load_workbook(COVERAGE_XLSX, read_only=True, data_only=True)
    rows = wb["Data"].iter_rows(values_only=True)
    header = next(rows)
    idx = {name: i for i, name in enumerate(header)}
    records = []
    for row in rows:
        if len(row) <= idx["COVERAGE"]:
            continue
        if row[idx["COVERAGE_CATEGORY"]] != "WUENIC" or row[idx["ANTIGEN"]] not in {"MCV1", "MCV2"}:
            continue
        records.append({
            "country": row[idx["CODE"]], "country_name": row[idx["NAME"]],
            "year": int(row[idx["YEAR"]]), "antigen": row[idx["ANTIGEN"]],
            "coverage": pd.to_numeric(row[idx["COVERAGE"]], errors="coerce"),
        })
    coverage = pd.DataFrame(records).drop_duplicates(["country", "year", "antigen"])
    coverage.to_csv(out, index=False)
    return coverage


def load_population() -> pd.DataFrame:
    pop = pd.read_csv(POP_FILE)
    require_columns(pop, {"country", "year", "age", "population"}, "single-age population")
    pop = pop.loc[pop["year"].between(2014, 2025) & pop["age"].between(0, 100)].copy()
    pop["population"] = pd.to_numeric(pop["population"], errors="coerce")
    return pop


def load_mcv2_introduction() -> pd.Series:
    intro = pd.read_csv(INTRO_FILE)
    require_columns(intro, {"country", "vaccine_code", "intro_year"}, "vaccine introduction")
    intro = intro.loc[intro["vaccine_code"].eq("MCV2")].dropna(subset=["intro_year"])
    return intro.groupby("country")["intro_year"].min().astype(int)


def build_single_age_gaps() -> pd.DataFrame:
    """Durable-protection gap for every country, calendar year and single age.

    Dose coverage is taken from the WUENIC year in which the cohort reached the
    scheduled age (MCV1 at 1, MCV2 at 4). MCV2 is set to zero for dose years
    before national introduction; for dose years after introduction but before
    the WUENIC MCV2 series starts it is also set to zero and flagged. Because
    G = 1 - 0.93 C1 - 0.04 C2, that choice can overstate a cohort's gap by at
    most 4 percentage points. Cohorts without MCV1 history remain missing.
    """
    coverage = load_extended_coverage()
    population = load_population()
    population = population.loc[population["age"].le(MAX_AGE)].copy()
    wide = coverage.pivot_table(index=["country", "year"], columns="antigen", values="coverage", aggfunc="first")
    wide = (wide.reindex(columns=["MCV1", "MCV2"]) / 100.0).reset_index()
    intro = load_mcv2_introduction()

    base = population[["country", "year", "age", "population"]].copy()
    base["birth_year"] = base["year"] - base["age"]
    for dose, target_age in (("MCV1", 1), ("MCV2", 4)):
        base["dose_year"] = base["birth_year"] + target_age
        base = base.merge(wide[["country", "year", dose]].rename(columns={"year": "dose_year", dose: f"{dose}_coverage"}),
                          on=["country", "dose_year"], how="left", validate="many_to_one")
        base.loc[base["dose_year"].gt(base["year"]), f"{dose}_coverage"] = 0.0
    dose_year_2 = base["birth_year"] + 4
    intro_year = base["country"].map(intro)
    pre_intro = base["MCV2_coverage"].isna() & intro_year.notna() & dose_year_2.lt(intro_year)
    imputed = base["MCV2_coverage"].isna() & ~pre_intro & base["MCV1_coverage"].notna()
    base["mcv2_status"] = np.select([pre_intro, imputed, base["MCV2_coverage"].notna()],
                                    ["pre_introduction", "imputed_zero", "observed"], "missing")
    base.loc[pre_intro | imputed, "MCV2_coverage"] = 0.0
    c1 = base["MCV1_coverage"].where(base["MCV1_coverage"].between(0, 1))
    c2 = base["MCV2_coverage"].where(base["MCV2_coverage"].between(0, 1)).clip(upper=c1)
    base["gap"] = (1 - c1) + (c1 - c2) * (1 - VE1) + c2 * (1 - VE2)
    base["dropout"] = (c1 - c2).clip(lower=0)
    base = base.drop(columns=["dose_year"])
    base.to_csv(OUT / "measles_single_age_gaps.csv.gz", index=False)
    return base


def band_summary(single: pd.DataFrame, lo: int, hi: int, keys=("country", "year"), include_age0: bool = False) -> pd.DataFrame:
    """Population-weighted gap over ages lo..hi (vaccine-eligible ages >= 1).

    A band is reported only when at least ``MIN_BAND_COMPLETENESS`` of its
    population has a reconstructed gap, so that a band is never summarised
    from an unrepresentative subset of its cohorts.
    """
    sub = single.loc[single["age"].between(0 if include_age0 else max(lo, 1), hi)].copy()
    sub["valid_pop"] = sub["population"].where(sub["gap"].notna(), 0.0)
    sub["w_gap"] = (sub["gap"] * sub["population"]).fillna(0.0)
    sub["w_dropout"] = (sub["dropout"] * sub["population"]).fillna(0.0)
    sub["w_imputed"] = sub["population"].where(sub["mcv2_status"].eq("imputed_zero"), 0.0)
    out = sub.groupby(list(keys), as_index=False)[["population", "valid_pop", "w_gap", "w_dropout", "w_imputed"]].sum()
    out["completeness"] = out["valid_pop"] / out["population"]
    ok = out["completeness"].ge(MIN_BAND_COMPLETENESS)
    out["gap"] = (out["w_gap"] / out["valid_pop"]).where(ok)
    out["dropout"] = (out["w_dropout"] / out["valid_pop"]).where(ok)
    out["mcv2_imputed_share"] = (out["w_imputed"] / out["population"]).where(ok)
    return out[list(keys) + ["population", "gap", "dropout", "completeness", "mcv2_imputed_share"]]


def reconstruct_age_profiles(single: pd.DataFrame) -> pd.DataFrame:
    names = pd.read_csv(CROSSWALK_FILE)[["canonical_iso3", "country_name_wpp"]].rename(
        columns={"canonical_iso3": "country", "country_name_wpp": "country_name"}).drop_duplicates("country")
    profiles = single[["country", "year"]].drop_duplicates().merge(names, on="country", how="left")
    for band, (lo, hi) in {**AGE_BAND_RANGES, "1_14": (1, 14), "1_29": (1, 29)}.items():
        summary = band_summary(single, lo, hi).rename(columns={
            "population": f"population_{band}", "gap": f"gap_{band}", "dropout": f"dropout_{band}",
            "completeness": f"completeness_{band}", "mcv2_imputed_share": f"mcv2_imputed_{band}"})
        profiles = profiles.merge(summary, on=["country", "year"], how="left", validate="one_to_one")
    # Infants are not vaccine-eligible; keep their band as all-susceptible.
    profiles["gap_0_0"] = 1.0
    coverage = single.dropna(subset=["MCV1_coverage"])
    for dose in ("MCV1", "MCV2"):
        cov = (coverage.assign(w=coverage[f"{dose}_coverage"] * coverage["population"])
               .groupby(["country", "year"], as_index=False)[["w", "population"]].sum())
        cov[f"cohort_weighted_{dose.lower()}"] = cov["w"] / cov["population"]
        profiles = profiles.merge(cov[["country", "year", f"cohort_weighted_{dose.lower()}"]], on=["country", "year"], how="left")
    profiles["n_valid_age_bands"] = profiles[[f"gap_{b}" for b in AGE_BANDS]].notna().sum(axis=1)
    profiles = profiles.sort_values(["country", "year"])
    for window in (5, 10):
        profiles[f"historical_gap_{window}y"] = profiles.groupby("country")["gap_1_14"].transform(
            lambda s: s.shift(1).rolling(window, min_periods=max(3, window // 2)).mean())
    profiles.to_csv(OUT / "measles_age_profiles_1980_2025.csv", index=False)
    return profiles


def load_who_monthly() -> tuple[pd.DataFrame, pd.DataFrame]:
    epi = pd.read_csv(EPI_FILE, low_memory=False)
    require_columns(epi, {"region", "country", "country_name", "year", "month", "measles_total"}, "WHO monthly surveillance")
    # Retain the longer 2016–2026 archive for country-level plots. The formal
    # hurdle models are restricted to 2024–2026 in ``main``.
    epi = epi.loc[epi["country"].isin(EU_EEA) & epi["year"].between(min(ARCHIVE_YEARS), max(ARCHIVE_YEARS))].copy()
    epi["region"] = "EU/EEA"
    epi["cases"] = pd.to_numeric(epi["measles_total"], errors="coerce").fillna(0.0)
    epi["year"] = epi["year"].astype(int); epi["month"] = epi["month"].astype(int)
    pop = pd.read_csv(ROOT / "data" / "raw" / "population_data.csv")
    pop = pop.loc[pop["year"].eq(2025), ["country", "population"]].rename(columns={"population": "population_2025"})
    epi = epi.merge(pop, on="country", how="left", validate="many_to_one")
    epi["offset_population"] = epi["population_2025"].fillna(1.0)
    epi["incidence_per_100k"] = 100_000 * epi["cases"] / epi["offset_population"]
    # Lagged external pressure is a transparent importation-pressure proxy.
    regional = epi.groupby(["region", "year", "month"], as_index=False).agg(region_cases=("cases", "sum"), region_population=("offset_population", "sum"))
    epi = epi.merge(regional, on=["region", "year", "month"], how="left")
    epi["external_cases"] = (epi["region_cases"] - epi["cases"]).clip(lower=0)
    epi["external_pressure"] = np.log1p(100_000 * epi["external_cases"] / epi["region_population"].replace(0, np.nan)).fillna(0)
    epi = epi.sort_values(["country", "year", "month"])
    epi["outbreak_month"] = epi["cases"].ge(3).astype(int)
    epi["sustained_2month"] = ((epi["cases"].gt(0)) & (epi.groupby("country")["cases"].shift(-1).gt(0))).astype(int)
    epi["external_pressure_lag"] = epi.groupby("country")["external_pressure"].shift(1).fillna(epi["external_pressure"])
    regional.to_csv(OUT / "measles_region_monthly.csv", index=False)
    epi.to_csv(OUT / "measles_country_monthly_2016_2026.csv", index=False)
    return epi, regional


def annual_wuenic() -> pd.DataFrame:
    wide = load_extended_coverage().pivot_table(index=["country", "year"], columns="antigen", values="coverage", aggfunc="first")
    return (wide.reindex(columns=["MCV1", "MCV2"]) / 100.0).rename(columns={"MCV1": "mcv1_coverage", "MCV2": "mcv2_coverage"}).reset_index()


def add_lagged_exposures(epi: pd.DataFrame, profiles: pd.DataFrame) -> pd.DataFrame:
    """Attach preceding-year cohort gaps and preceding-year annual WUENIC coverage."""
    keep = ["country", "year", "gap_1_14", "gap_1_29", "historical_gap_5y", "historical_gap_10y"] + [f"gap_{b}" for b in AGE_MODEL_BANDS]
    expo = profiles[keep].merge(annual_wuenic(), on=["country", "year"], how="left")
    expo["year"] = expo["year"] + 1
    expo = expo.rename(columns={c: f"{c}_lag" for c in expo.columns if c not in {"country", "year"}})
    out = epi.merge(expo, on=["country", "year"], how="left", validate="many_to_one")
    out["month_sin"] = np.sin(2 * np.pi * out["month"] / 12)
    out["month_cos"] = np.cos(2 * np.pi * out["month"] / 12)
    out.to_csv(OUT / "measles_monthly_model_panel.csv", index=False)
    return out


# ECDC Surveillance Atlas, health topic 35 (measles). Measure 1250846 is the
# case distribution by age group (% of N) and 1250850 by vaccination status.
# Note: the 12074xx measures belong to health topic 37 (mumps) and must not
# be used here.
ECDC_MEASLES_MEASURES = {1250846: "age", 1250850: "vaccination"}
ECDC_YEARS = (2023, 2025)


def fetch_ecdc_case_distributions() -> pd.DataFrame:
    out = OUT / "ecdc_measles_case_distributions_2023_2025.csv"
    if out.exists():
        return pd.read_csv(out)
    url = "https://atlas.ecdc.europa.eu/public/AtlasService/rest/post/GetMeasureResultsForTimePeriodAndGeoLevel"
    frames = []
    for measure, kind in ECDC_MEASLES_MEASURES.items():
        for year in range(ECDC_YEARS[0], ECDC_YEARS[1] + 1):
            payload = {"measureIds": str(measure), "timeCodes": None, "startTimeCode": str(year),
                       "endTimeCodeExcl": str(year + 1), "geoLevel": "2"}
            response = requests.post(url, json=payload, timeout=120)
            response.raise_for_status()
            rows = pd.DataFrame(response.json()["GetMeasureResultsForTimePeriodAndGeoLevelResult"]["MeasureResults"])
            frames.append(rows.assign(kind=kind))
    data = pd.concat(frames, ignore_index=True)
    data["TimeCode"] = pd.to_numeric(data["TimeCode"], errors="coerce")
    data = data.loc[data["TimeCode"].between(*ECDC_YEARS) & data["GeoCountry"].notna()]
    data = data.drop_duplicates(["MeasureId", "GeoCountry", "TimeCode", "XLabel"])
    crosswalk = pd.read_csv(CROSSWALK_FILE)[["iso2", "canonical_iso3"]].dropna().drop_duplicates("iso2")
    data = data.merge(crosswalk, left_on="GeoCountry", right_on="iso2", how="left")
    # Greece is reported as EL in TESSy.
    data.loc[data["GeoCountry"].eq("EL"), "canonical_iso3"] = "GRC"
    keep = ["kind", "MeasureId", "GeoCountry", "canonical_iso3", "TimeCode", "XLabel", "YValue", "N"]
    data[keep].to_csv(out, index=False)
    return data[keep]


ECDC_BANDS = {"<1": (0, 0), "1-4": (1, 4), "5-9": (5, 9), "10-14": (10, 14), "15-19": (15, 19), "20-29": (20, 29), "30+": (30, None)}
# Open-ended adult bands are truncated here for the denominator and the gap;
# cases above the cap stay in the numerator (documented in Methods).
OPEN_BAND_CAP = 39


def _cohort_aligned_band_summary(single: pd.DataFrame, lo: int, hi: int) -> pd.DataFrame:
    """Gap for the same birth cohorts one year before the outcome.

    A case aged ``a`` in outcome year ``t`` was aged ``a-1`` in ``t-1``.
    The age-0 state is retained for the 1–4 band because those children had
    not yet reached the standardised MCV1 anchor in the preceding year.
    """
    out = band_summary(single, lo - 1, hi - 1, include_age0=(lo - 1 == 0)).rename(columns={"year": "exposure_year"})
    return out


def build_case_age_panel(single: pd.DataFrame) -> pd.DataFrame:
    """ECDC case-based age distributions joined to the cohort reconstruction.

    Each row is one country, outcome year and reported age band, with cases,
    the band population in the outcome year and the gap for the same birth
    cohorts one year earlier (age ``a`` in year ``t`` maps to age ``a-1`` in
    ``t-1``).
    """
    ecdc = fetch_ecdc_case_distributions()
    ecdc = ecdc.loc[ecdc["kind"].eq("age") & ecdc["canonical_iso3"].notna() & ecdc["N"].gt(0)].copy()
    ecdc[["age_lo", "age_hi"]] = pd.DataFrame(ecdc["XLabel"].map(ECDC_BANDS).tolist(), index=ecdc.index)
    ecdc = ecdc.assign(source="ECDC", region="EU/EEA", country=ecdc["canonical_iso3"], outcome_year=ecdc["TimeCode"].astype(int),
                       cases=(ecdc["YValue"].fillna(0) / 100 * ecdc["N"]).round())
    panel = ecdc[["source", "region", "country", "outcome_year", "age_lo", "age_hi", "cases"]].reset_index(drop=True)
    panel["open_ended"] = panel["age_hi"].isna()
    panel["age_hi_used"] = panel["age_hi"].fillna(OPEN_BAND_CAP).clip(upper=OPEN_BAND_CAP).astype(int)
    panel["age_lo"] = panel["age_lo"].astype(int)
    # Infant-only bands and bands starting above the cap carry no cohort contrast.
    panel = panel.loc[panel["age_hi_used"].ge(1) & panel["age_lo"].le(OPEN_BAND_CAP)].copy()
    max_pop_year = int(single["year"].max())
    rows = []
    for (lo, hi), idx in panel.groupby(["age_lo", "age_hi_used"]).groups.items():
        sub = panel.loc[idx]
        expo_cohort = _cohort_aligned_band_summary(single, lo, hi).rename(columns={"gap": "gap_cohort_aligned",
                                                                                     "completeness": "cohort_completeness",
                                                                                     "mcv2_imputed_share": "cohort_mcv2_imputed_share"})
        pop_all = single.loc[single["age"].between(lo, hi)].groupby(["country", "year"], as_index=False)["population"].sum()
        sub = sub.assign(exposure_year=sub["outcome_year"] - 1, population_year=sub["outcome_year"].clip(upper=max_pop_year))
        sub = sub.merge(expo_cohort[["country", "exposure_year", "gap_cohort_aligned", "cohort_completeness", "cohort_mcv2_imputed_share"]], on=["country", "exposure_year"], how="left")
        sub["gap"] = sub["gap_cohort_aligned"]
        sub["completeness"] = sub["cohort_completeness"]
        sub["mcv2_imputed_share"] = sub["cohort_mcv2_imputed_share"]
        sub = sub.merge(pop_all.rename(columns={"year": "population_year"}), on=["country", "population_year"], how="left")
        rows.append(sub)
    panel = pd.concat(rows, ignore_index=True)
    ages = single[["country", "year", "age", "population"]].rename(columns={"year": "population_year"})
    span = panel[["country", "population_year", "age_lo", "age_hi_used"]].drop_duplicates().merge(ages, on=["country", "population_year"])
    span = span.loc[span["age"].between(span["age_lo"], span["age_hi_used"])].assign(w=lambda d: d["age"] * d["population"])
    span = span.groupby(["country", "population_year", "age_lo", "age_hi_used"], as_index=False)[["w", "population"]].sum()
    span["mean_age"] = span["w"] / span["population"]
    panel = panel.merge(span[["country", "population_year", "age_lo", "age_hi_used", "mean_age"]],
                        on=["country", "population_year", "age_lo", "age_hi_used"], how="left")
    panel["band"] = np.where(panel["open_ended"], panel["age_lo"].astype(str) + "+",
                             panel["age_lo"].astype(str) + "–" + panel["age_hi_used"].astype(str))
    panel = panel.sort_values(["country", "outcome_year", "age_lo"])
    panel["exposure_alignment"] = "same_birth_cohort_previous_year"
    panel.to_csv(OUT / "measles_case_age_panel.csv", index=False)
    return panel


def build_case_vaccination_panel(single: pd.DataFrame, case_age: pd.DataFrame) -> pd.DataFrame:
    """Case vaccination status with the national >=1-dose coverage of the same ages.

    ``national_mcv1`` is the reconstructed MCV1 coverage (preceding year) of
    the cohorts in which the cases occurred: band-specific coverage weighted
    by each band's reported cases. It is the coverage that the cases would
    have had if they had been drawn at random from the national cohorts.
    """
    ecdc = fetch_ecdc_case_distributions()
    ecdc = ecdc.loc[ecdc["kind"].eq("vaccination") & ecdc["canonical_iso3"].notna() & ecdc["N"].gt(0)]
    wide = ecdc.pivot_table(index=["canonical_iso3", "TimeCode", "N"], columns="XLabel", values="YValue", aggfunc="first").fillna(0).reset_index()
    panel = pd.DataFrame({
        "source": "ECDC", "region": "EU/EEA", "country": wide["canonical_iso3"], "outcome_year": wide["TimeCode"].astype(int),
        "cases": wide["N"], "known_status": (wide["N"] * (1 - wide["Unknown"] / 100)).round(),
        "vaccinated": (wide["N"] * (wide["1 dose"] + wide[">=2 doses"] + wide["UNKDOSE"]) / 100).round()})
    cov = single.dropna(subset=["MCV1_coverage"]).assign(w=lambda d: d["MCV1_coverage"] * d["population"])
    bands = case_age.loc[case_age["cases"].gt(0), ["country", "outcome_year", "exposure_year", "age_lo", "age_hi_used", "cases"]]
    rows = []
    for (lo, hi), sub in bands.groupby(["age_lo", "age_hi_used"]):
        c = cov.loc[cov["age"].between(max(lo, 1), hi)].groupby(["country", "year"], as_index=False)[["w", "population"]].sum()
        c["band_mcv1"] = c["w"] / c["population"]
        rows.append(sub.merge(c[["country", "year", "band_mcv1"]].rename(columns={"year": "exposure_year"}), on=["country", "exposure_year"], how="left"))
    bands = pd.concat(rows).dropna(subset=["band_mcv1"])
    national = (bands.assign(w=bands["band_mcv1"] * bands["cases"]).groupby(["country", "outcome_year"], as_index=False)
                .agg(w=("w", "sum"), age_cases=("cases", "sum")))
    national["national_mcv1"] = national.pop("w") / national["age_cases"]
    panel = panel.merge(national, on=["country", "outcome_year"], how="left")
    panel.to_csv(OUT / "measles_case_vaccination_panel.csv", index=False)
    return panel


def compute_re(p: pd.DataFrame, contact: pd.DataFrame, scenario: str) -> float:
    bands = [c for c in contact.columns if c != "Unnamed: 0"]
    C = contact[bands].to_numpy(float)
    gap = np.full(len(bands), 0.30)
    mapping = {"0_4": ["0_0", "1_4"], "5_9": ["5_9"], "10_14": ["10_14"], "15_19": ["15_19"], "20_24": ["20_29"], "25_29": ["20_29"],
               "30_34": ["30_39"], "35_39": ["30_39"]}
    for i, b in enumerate(bands):
        source = mapping.get(b)
        if source:
            vals = [p.get(f"gap_{x}", np.nan) for x in source]
            vals = [v for v in vals if pd.notna(v)]
            if vals: gap[i] = float(np.mean(vals))
    if scenario == "catchup_10_14": gap[2] = min(gap[2], .10)
    elif scenario == "catchup_15_19": gap[3] = min(gap[3], .10)
    elif scenario == "high_child_history": gap[:4] = np.minimum(gap[:4], .10)
    # ``gap`` is already the complement of programme-derived protection,
    # including the 0.97 completed-series clinical effectiveness. Reapplying
    # VE here would double-scale measles protection. MMR acquisition ratio is 1.
    susceptibility = gap
    return float(np.linalg.eigvals(np.diag(susceptibility) @ C).real.max() / np.linalg.eigvals(C).real.max())


def transmission_scenarios(profiles: pd.DataFrame) -> pd.DataFrame:
    # DEU is the only EU/EEA country with a local reciprocal contact matrix.
    countries = {"DEU": "DEU"}
    rows = []
    for country, label in countries.items():
        contact_path = CONTACT_DIR / f"{country}.csv"
        if not contact_path.exists():
            continue
        p = profiles.loc[profiles["country"].eq(country) & profiles["year"].eq(2025)].iloc[0].to_dict()
        contact = pd.read_csv(contact_path)
        for scenario in ["observed", "catchup_10_14", "catchup_15_19", "high_child_history"]:
            rows.append({"country": country, "country_name": label, "scenario": scenario, "relative_Re": compute_re(p, contact, scenario)})
    result = pd.DataFrame(rows)
    if not result.empty:
        result["relative_reduction"] = result.groupby("country")["relative_Re"].transform(lambda s: 100 * (1 - s / s.iloc[0]))
    result.to_csv(OUT / "measles_transmission_scenarios.csv", index=False)
    return result


def main() -> None:
    single = build_single_age_gaps()
    profiles = reconstruct_age_profiles(single)
    epi_all, _ = load_who_monthly()
    panel = add_lagged_exposures(epi_all.loc[epi_all["year"].isin(FORMAL_YEARS)].copy(), profiles)
    case_age = build_case_age_panel(single)
    vaccination = build_case_vaccination_panel(single, case_age)
    scenarios = transmission_scenarios(profiles)
    status = single.loc[single["year"].eq(2023) & single["age"].between(1, 39)]
    manifest = {
        "archive_years": ARCHIVE_YEARS, "formal_years": FORMAL_YEARS,
        "countries": "ECDC case data: 30; WHO/UNICEF monthly data: 29", "age_bands": AGE_BANDS,
        "min_band_completeness": MIN_BAND_COMPLETENESS, "open_band_cap": OPEN_BAND_CAP,
        "mcv2_status_share_2023_age_1_39": status.groupby("mcv2_status")["population"].sum().div(status["population"].sum()).round(4).to_dict(),
        "n_archive_month_rows": len(epi_all), "n_month_rows": len(panel),
        "n_month_countries": int(panel.country.nunique()), "n_countries": int(panel.country.nunique()),
        "case_age_sources": case_age.groupby("source")["country"].nunique().to_dict(),
        "case_age_primary": "closed 1-29-year bands; age a in outcome year t uses age a-1 in t-1 from the same birth cohort",
        "case_age_open_band": "30+ cases retained only as sensitivity with 30-39 denominator and exposure",
        "case_vaccination_rows": vaccination.groupby("source").size().to_dict(),
        "screening_quantity": "screening-equivalent coverage under assumed vaccine effectiveness; model-based, not a directly observed community denominator",
        "measles_contact_mapping": "s_i = G_i; 0.97 completed-series effectiveness is already included in P = 1 - G",
        "scenario_countries": int(scenarios.country.nunique()) if not scenarios.empty else 0,
        "external_pressure": "lagged EU/EEA incidence excluding focal country; proxy, not direct importation measurement",
        "2026_population": "2025 population used because 2026 single-age WPP values are not in the local archive",
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

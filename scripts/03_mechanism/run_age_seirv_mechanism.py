#!/usr/bin/env python3
"""Age-structured SEIRV mechanism layer for the selected-country panel.

This layer is a mechanism view over the canonical primary kernel. It combines
the reconstructed cohort protection index with WPP single-age demography and the
country-specific Prem et al. (2017) contact matrices.  Annual deterministic
epidemic potentials are simulated for each country, disease and observation year.

The cohort gap is a protection-state index from the reconstruction.  It is mapped
to infection susceptibility through an explicit transmission-blocking scenario;
the mapping is separated from clinical protection because toxoid and pertussis
vaccines can prevent disease more strongly than acquisition.  The resulting
quantities are mechanism scenarios, rather than fitted effective susceptibility
or causal intervention effects.
"""
from __future__ import annotations

import json
import math
import os
import re
import zipfile
from pathlib import Path
from typing import Dict, Iterable, Mapping
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "02_models"))
from select_transmission_countries import selected_countries
from model_spec import SPEC, protection_factors, DISEASES as SPEC_DISEASES, AGE_LABELS as SPEC_AGE_LABELS
from transmission_kernel import simulate_with_diagnostics

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
COHORT = Path(os.environ.get("VPDS_COHORT_DIR", ROOT / "output" / "five_disease_revision"))
OUT = Path(os.environ.get("VPDS_MECHANISM_OUT", ROOT / "output" / "revised" / "age_seirv"))
OUT.mkdir(parents=True, exist_ok=True)

COUNTRIES, COUNTRY_SELECTION = selected_countries()
COUNTRY_NAMES = COUNTRY_SELECTION.set_index("country").country_name.to_dict()
PREM_NAMES = COUNTRY_SELECTION.set_index("country").prem_sheet.to_dict()
DISEASES = SPEC_DISEASES
YEARS = np.arange(int(SPEC["time"]["analysis_years"][0]),
                 int(SPEC["time"]["analysis_years"][1]) + 1, dtype=int)
AGE_LABELS = SPEC_AGE_LABELS
AGE_START = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75])

# Reference values are intentionally exposed as scenarios.  They are not
# estimated from the aggregate panel.  R0 values are literature-compatible
# central values used to scale each contact matrix's dominant eigenvalue.
PARAMETERS = {
    disease: {"R0": values["R0"], "latent_days": values["latent_days"],
              "infectious_days": values["infectious_days"],
              "transmission_blocking": values["acquisition_ratio"],
              "source_reduction": values.get("source_reduction", 0.0),
              "blocking_range": [values["acquisition_ratio_grid"][0], 1.0]}
    for disease, values in SPEC["diseases"].items()
}
ADULT_SCENARIOS = SPEC["scenarios"]["adult_gap"]
CENTRAL_ADULT = SPEC["age"]["adult_gap"]


def _column_number(cell_ref: str) -> int:
    letters = re.match(r"([A-Z]+)", cell_ref).group(1)
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n - 1


def _read_xlsx_sheets(path: Path) -> Dict[str, np.ndarray]:
    """Read Prem workbooks while ignoring their incorrect A1 worksheet bounds."""
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
        sheets = []
        for s in wb.find(ns + "sheets"):
            target = relmap[s.attrib[relns + "id"]]
            target = target.lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
            sheets.append((s.attrib["name"], target))
        out = {}
        for sheet_name, target in sheets:
            root = ET.fromstring(z.read(target))
            rows = {}
            max_row = max_col = 0
            for row in root.iter(ns + "row"):
                rnum = int(row.attrib["r"]) - 1
                for cell in row.findall(ns + "c"):
                    ref = cell.attrib.get("r", "A1")
                    cnum = _column_number(ref)
                    v = cell.find(ns + "v")
                    value = v.text if v is not None else None
                    if value is not None and cell.attrib.get("t") == "s":
                        value = shared[int(value)]
                    try:
                        if value is not None and cell.attrib.get("t") != "s":
                            value = float(value)
                    except ValueError:
                        pass
                    rows.setdefault(rnum, {})[cnum] = value
                    max_row = max(max_row, rnum)
                    max_col = max(max_col, cnum)
            arr = np.full((max_row + 1, max_col + 1), np.nan, dtype=object)
            for i, cols in rows.items():
                for j, value in cols.items():
                    arr[i, j] = value
            out[sheet_name] = arr
        return out


def load_prem_matrices() -> Dict[str, np.ndarray]:
    path1 = RAW / "contact_matrices_prem2017" / "MUestimates_all_locations_1.xlsx"
    path2 = RAW / "contact_matrices_prem2017" / "MUestimates_all_locations_2.xlsx"
    all_sheets = {}
    for path in [path1, path2]:
        all_sheets.update(_read_xlsx_sheets(path))
    result = {}
    for iso in COUNTRIES:
        name = PREM_NAMES[iso]
        if name not in all_sheets:
            raise KeyError(f"Prem matrix missing {iso} ({name})")
        raw = all_sheets[name]
        # The workbook stores a 16-column header followed by 16 matrix rows.
        if raw.shape[0] < 16 or raw.shape[1] < 16:
            raise ValueError(f"Unexpected Prem matrix dimensions for {name}: {raw.shape}")
        # Some sheets carry X1..X16 column labels; later sheets omit that row.
        start = 1 if raw.shape[0] >= 17 and isinstance(raw[0, 0], str) else 0
        result[iso] = np.asarray(raw[start:start + 16, :16], dtype=float)
    return result


def load_wpp_age() -> pd.DataFrame:
    p = RAW / "population_by_single_age.csv"
    d = pd.read_csv(p, usecols=["country", "year", "age", "population"])
    d = d[d.country.isin(COUNTRIES) & d.year.isin(YEARS)].copy()
    d["age"] = pd.to_numeric(d.age, errors="coerce")
    d["population"] = pd.to_numeric(d.population, errors="coerce")
    d.loc[d.age >= 75, "age_band"] = "75_plus"
    d.loc[d.age.between(0, 4), "age_band"] = "0_4"
    for lo in range(5, 75, 5):
        d.loc[d.age.between(lo, lo + 4), "age_band"] = f"{lo}_{lo+4}"
    d = d.dropna(subset=["age_band"])
    # Preserve the infant share before collapsing to the 0–4 contact band.
    # Routine doses in the reconstructed profiles begin after infancy, so
    # applying the 1–4 gap to every 0–4 resident would overstate protection.
    infant = (d[d.age.eq(0)].groupby(["country", "year"], as_index=False)
              .population.sum().rename(columns={"population": "population_age0"}))
    early = (d[d.age.between(1, 4)].groupby(["country", "year"], as_index=False)
             .population.sum().rename(columns={"population": "population_age1_4"}))
    out = d.groupby(["country", "year", "age_band"], as_index=False).population.sum()
    wide = out.pivot_table(index=["country", "year"], columns="age_band", values="population", fill_value=0)
    wide = wide.reindex(columns=AGE_LABELS, fill_value=0).reset_index()
    wide = wide.merge(infant, on=["country", "year"], how="left", validate="one_to_one")
    wide = wide.merge(early, on=["country", "year"], how="left", validate="one_to_one")
    wide[["population_age0", "population_age1_4"]] = wide[["population_age0", "population_age1_4"]].fillna(0.0)
    if (wide[AGE_LABELS].sum(axis=1) <= 0).any():
        raise ValueError("WPP age population contains an empty country-year")
    return wide


def load_cohort_gaps() -> pd.DataFrame:
    files = {
        "Measles": COHORT / "cohort_level_all_ages.csv",
        "Rubella": COHORT / "cohort_level_all_ages.csv",
        "Pertussis": COHORT / "cohort_level_all_ages.csv",
        "Diphtheria": COHORT / "cohort_level_all_ages.csv",
        "Mumps": COHORT / "mumps_all_mcv_proxy_cohort_level.csv",
    }
    frames = []
    for disease, path in files.items():
        d = pd.read_csv(path, low_memory=False)
        d = d[(d.country.isin(COUNTRIES)) & (d.observation_year.isin(YEARS)) & d.age.between(1, 14)].copy()
        # The shared cohort file contains all component diseases.  Filter on
        # its original label before assigning the analysis label; assigning
        # first would duplicate every disease profile within each target.
        d = d[d.disease.eq(disease)].copy()
        d["disease"] = disease
        d["gap_waning"] = pd.to_numeric(d["gap_waning_projected"], errors="coerce")
        d["gap_unwaning"] = pd.to_numeric(d["gap_no_waning_projected"], errors="coerce")
        frames.append(d[["country", "observation_year", "age", "cohort_size", "gap_waning", "gap_unwaning", "disease"]])
    x = pd.concat(frames, ignore_index=True)
    # Keep profile gaps separate from infection-susceptibility mapping.
    x["cohort_size"] = pd.to_numeric(x.cohort_size, errors="coerce")
    x["gap_waning"] = x.gap_waning.clip(0, 1)
    x["gap_unwaning"] = x.gap_unwaning.clip(0, 1)
    # The same declared neutral childhood gap is used in both model paths;
    # absent profiles never become zero protection or reported case zeros.
    x["profile_fallback"] = x.gap_waning.isna() | x.gap_unwaning.isna()
    x.loc[x.profile_fallback, ["gap_waning", "gap_unwaning"]] = SPEC["age"]["missing_gap"]
    if x[["gap_waning", "gap_unwaning", "cohort_size"]].isna().any().any():
        missing = x[x[["gap_waning", "gap_unwaning", "cohort_size"]].isna().any(axis=1)]
        raise ValueError(f"Unresolved cohort profiles in mechanism layer: {missing[['country','disease','observation_year','age']].head().to_dict('records')}")
    return x


def age_gap_vector(g: pd.DataFrame, mode: str, adult: float,
                   infant_fraction: float = 0.0, infant_gap: float = 1.0) -> np.ndarray:
    """Map ages 1–14 to WPP bands and keep an explicit infant scenario."""
    col = "gap_waning" if mode == "waning" else "gap_unwaning"
    out = np.full(16, float(adult), dtype=float)
    for idx, ages in enumerate([range(1, 5), range(5, 10), range(10, 15)]):
        q = g[g.age.isin(list(ages))]
        out[idx] = np.average(q[col].to_numpy(float), weights=q.cohort_size.to_numpy(float))
    out[0] = ((1.0 - infant_fraction) * out[0] + infant_fraction * infant_gap)
    return out


def reciprocal_matrix(matrix: np.ndarray, pop: np.ndarray) -> np.ndarray:
    """Enforce N_i C_ij = N_j C_ji without altering the mean contact scale."""
    p = np.clip(pop, 1.0, None)
    corrected = 0.5 * (matrix + matrix.T * (p[None, :] / p[:, None]))
    corrected = np.clip(corrected, 0.0, None)
    return corrected


def spectral_radius(x: np.ndarray) -> float:
    vals = np.linalg.eigvals(x)
    return float(np.max(np.real(vals)))


def simulate_year(n: np.ndarray, s_ratio: np.ndarray, c: np.ndarray, disease: str,
                  source_ratio: np.ndarray | None = None) -> Mapping[str, float]:
    """Run the canonical daily kernel for one deterministic year.

    ``s_ratio`` is the effective susceptible fraction. The kernel accepts that
    fraction as its unprotected gap when acquisition VE is one.
    """
    p = PARAMETERS[disease]
    # Preserve the fully susceptible contact denominator when source protection
    # is present. Renormalizing the modified matrix would erase its effect.
    source_ratio = np.ones(16) if source_ratio is None else source_ratio
    source_contacts = c @ np.diag(source_ratio)
    rho_c = spectral_radius(c)
    s_ratio = np.clip(np.asarray(s_ratio, dtype=float), 0.0, 1.0)
    targets = s_ratio[None, :]
    initial_re = p["R0"] * spectral_radius(np.diag(s_ratio) @ source_contacts) / rho_c
    annual, diagnostics = simulate_with_diagnostics(
        n, np.zeros(1), source_contacts, targets,
        p["R0"] * spectral_radius(source_contacts) / rho_c, p["latent_days"],
        p["infectious_days"], acquisition_ratio=1.0,
        clinical_fraction=1.0, importation=0.0, reporting=1.0)
    cum = diagnostics["cumulative_infections"]
    return {
        "effective_R_initial": float(initial_re),
        "attack_rate_365d": float(cum / np.sum(n)),
        "peak_daily_infections_per_100k": float(diagnostics["peak_daily_infections"] / np.sum(n) * 1e5),
        "cumulative_infections": float(cum),
        "vaccine_protected_fraction_initial": float(np.sum(n * (1.0 - s_ratio)) / np.sum(n)),
        "initial_effective_R_above_one": bool(initial_re > 1.0),
    }


def build_gap_panel() -> pd.DataFrame:
    gaps = load_cohort_gaps()
    pop = load_wpp_age()
    matrices = load_prem_matrices()
    rows = []
    for country in COUNTRIES:
        for year in YEARS:
            prow = pop[(pop.country == country) & (pop.year == year)].iloc[0]
            n = prow[AGE_LABELS].to_numpy(float)
            infant_fraction = float(prow.population_age0 / max(prow.population_age0 + prow.population_age1_4, 1.0))
            c = reciprocal_matrix(matrices[country], n)
            # Matrix rows are recipient ages; columns are source ages. The
            # dominant eigenvalue is used for disease-specific R0 scaling.
            for disease in DISEASES:
                g = gaps[(gaps.country == country) & (gaps.observation_year == year) & (gaps.disease == disease)]
                for adult in [CENTRAL_ADULT]:
                    for mode in ["waning", "unwaning"]:
                        profile = age_gap_vector(g, mode, adult, infant_fraction=infant_fraction)
                        b = PARAMETERS[disease]["transmission_blocking"]
                        # gap is a protection-state index; B maps clinical
                        # protection to infection blocking in this scenario.
                        s_ratio = np.clip(1.0 - b * (1.0 - profile), 0.0, 1.0)
                        _, source_ratio, _ = protection_factors(profile, b, PARAMETERS[disease]["source_reduction"])
                        sim = simulate_year(n, s_ratio, c, disease, source_ratio)
                        rows.append({"country": country, "country_name": COUNTRY_NAMES[country],
                                     "year": int(year), "disease": disease, "mode": mode,
                                     "adult_susceptibility_scenario": adult,
                                     "transmission_blocking_scenario": b,
                                     "source_reduction_scenario": PARAMETERS[disease]["source_reduction"],
                                     "profile_fallback_fraction": float(g.profile_fallback.mean()),
                                     "mean_profile_gap": float(np.average(profile, weights=n)),
                                     "gap_0_4": profile[0], "gap_5_9": profile[1], "gap_10_14": profile[2],
                                     "adult_gap": adult, "age_gradient_10_14_minus_0_4": profile[2] - profile[0],
                                     "prem_matrix_spectral_radius": spectral_radius(c),
                                     "population": float(np.sum(n)), **sim})
    return pd.DataFrame(rows)


def sensitivity_panel() -> pd.DataFrame:
    gaps = load_cohort_gaps()
    pop = load_wpp_age()
    matrices = load_prem_matrices()
    rows = []
    for country in COUNTRIES:
        for year in YEARS:
            prow = pop[(pop.country == country) & (pop.year == year)].iloc[0]
            n = prow[AGE_LABELS].to_numpy(float)
            infant_fraction = float(prow.population_age0 / max(prow.population_age0 + prow.population_age1_4, 1.0))
            c = reciprocal_matrix(matrices[country], n)
            for disease in DISEASES:
                g = gaps[(gaps.country == country) & (gaps.observation_year == year) & (gaps.disease == disease)]
                for adult in ADULT_SCENARIOS:
                    for blocking in [PARAMETERS[disease]["blocking_range"][0], PARAMETERS[disease]["transmission_blocking"], 1.0]:
                        for mode in ["waning", "unwaning"]:
                            profile = age_gap_vector(g, mode, adult, infant_fraction=infant_fraction)
                            s_ratio = np.clip(1.0 - blocking * (1.0 - profile), 0.0, 1.0)
                            _, source_ratio, _ = protection_factors(profile, blocking, PARAMETERS[disease]["source_reduction"])
                            sim = simulate_year(n, s_ratio, c, disease, source_ratio)
                            rows.append({"country": country, "year": int(year), "disease": disease,
                                         "mode": mode, "adult_susceptibility_scenario": adult,
                                         "transmission_blocking_scenario": blocking,
                                         "source_reduction_scenario": PARAMETERS[disease]["source_reduction"],
                                         "profile_fallback_fraction": float(g.profile_fallback.mean()),
                                         "mean_profile_gap": float(np.average(profile, weights=n)),
                                         "age_gradient_10_14_minus_0_4": profile[2] - profile[0], **sim})
    return pd.DataFrame(rows)


def main() -> None:
    panel = build_gap_panel()
    sens = sensitivity_panel()
    panel.to_csv(OUT / "annual_age_seirv_results.csv", index=False)
    sens.to_csv(OUT / "annual_age_seirv_sensitivity.csv", index=False)
    matrix_dir = OUT / "contact_matrices_reciprocal"
    matrix_dir.mkdir(exist_ok=True)
    matrices = load_prem_matrices()
    p = load_wpp_age()
    matrix_audit = []
    for country, matrix in matrices.items():
        n = p[(p.country == country) & (p.year == 2019)].iloc[0][AGE_LABELS].to_numpy(float)
        corr = reciprocal_matrix(matrix, n)
        raw_err = float(np.max(np.abs(np.diag(n) @ matrix - (np.diag(n) @ matrix).T)))
        corr_err = float(np.max(np.abs(np.diag(n) @ corr - (np.diag(n) @ corr).T)))
        matrix_audit.append({"country": country, "raw_max_reciprocity_error": raw_err,
                             "corrected_max_reciprocity_error": corr_err,
                             "raw_spectral_radius": spectral_radius(matrix),
                             "corrected_spectral_radius": spectral_radius(corr)})
        pd.DataFrame(corr, index=AGE_LABELS, columns=AGE_LABELS).to_csv(matrix_dir / f"{country}.csv")
    pd.DataFrame(matrix_audit).to_csv(OUT / "contact_matrix_audit.csv", index=False)
    manifest = {
        "model": "annual deterministic 16-age SEIRV epidemic potential",
        "countries": COUNTRIES, "country_selection": "output/revised/transmission_country_selection/selection_manifest.json", "diseases": DISEASES, "years": YEARS.tolist(),
        "age_groups": AGE_LABELS,
        "demography": "data/raw/population_by_single_age.csv; UN WPP 2024 single-age estimates",
        "contact_matrix": "Prem et al. 2017 all-location country matrices from PLOS supplementary S2",
        "contact_source_url": "https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1005697.s002&type=supplementary",
        "reciprocity": "C_ij corrected so N_i C_ij equals N_j C_ji using each country-year WPP age vector",
        "primary_cohort_input": "projected cohort gap from output/five_disease_revision; mumps uses all-MCV proxy",
        "missing_profile_gap": SPEC["age"]["missing_gap"],
        "adult_scenarios": ADULT_SCENARIOS,
        "parameter_status": "literature-compatible central R0 and natural-history scenarios; no aggregate-case calibration",
        "susceptibility_definition": ("clinical protection mass mapped separately to recipient acquisition and source infectiousness; operator diag(s) C diag(t)"
            if SPEC.get("protection_mapping") == "clinical_protection_mass_v2" else
            "clinical protection-state gap mapped to infection susceptibility by an explicit transmission-blocking scenario"),
        "natural_infection_immunity": "age-specific natural-infection immunity is unidentified from aggregate inputs; adult residual susceptibility and protection-state scenarios absorb this uncertainty",
        "missing_case_policy": "case data are not used for SEIRV fitting; missing annual cases remain missing in primary data products",
        "outputs": ["annual_age_seirv_results.csv", "annual_age_seirv_sensitivity.csv", "disease_mechanism_summary.csv", "disease_mechanism_contrasts.csv", "country_mechanism_contrasts.csv", "profile_fallback_audit.csv", "contact_matrix_audit.csv", "contact_matrices_reciprocal/"],
        "parameters": PARAMETERS,
        "protection_mapping": SPEC.get("protection_mapping", "legacy_index_mapping"),
        "cohort_input_directory": str(COHORT),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = ["# Age-structured SEIRV mechanism layer", "",
             f"This is an independent mechanism analysis for {len(COUNTRIES)} selected countries and five diseases.",
             "It uses WPP 2024 single-age populations, Prem et al. (2017) all-location contact matrices,",
             "and the reconstructed cohort protection index. Aggregate case counts are not used to fit beta",
             "or reporting parameters; missing case rows therefore cannot become zero observations.", "",
             "The cohort gap describes a protection-state scenario. Infection susceptibility is obtained",
             "with `s = 1 - B * (1 - gap)`, where B is an explicit transmission-blocking scenario. This",
             "separates clinical protection from acquisition or onward-transmission protection. Adult",
             "susceptibility and B are sensitivity dimensions. Outputs are epidemic potentials and",
             "structural comparisons, not causal effects or observed incidence forecasts.", ""]
    if SPEC.get("protection_mapping") == "clinical_protection_mass_v2":
        lines += ["Review mapping: the cohort protection mass already includes clinical VE.",
                  "Use an acquisition-to-clinical ratio for recipients and a separate source",
                  "factor for protected people among infections. The operator is `diag(s) C diag(t)`.",
                  "Diphtheria can have zero acquisition contrast and a positive source contrast.",
                  "Source protection is fixed at the central value in the acquisition grid;",
                  "the separate operator sensitivity table varies the source effect.", ""]
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")
    summary = (panel[panel["mode"].eq("waning")].groupby("disease", as_index=False)
               .agg(mean_Re=("effective_R_initial", "mean"), mean_attack_rate=("attack_rate_365d", "mean"),
                    mean_peak_per_100k=("peak_daily_infections_per_100k", "mean"),
                    mean_age_gradient=("age_gradient_10_14_minus_0_4", "mean")))
    summary.to_csv(OUT / "disease_mechanism_summary.csv", index=False)
    paired = (panel.pivot_table(index=["country", "country_name", "year", "disease"],
                                 columns="mode", values=["effective_R_initial", "attack_rate_365d"],
                                 aggfunc="first").reset_index())
    paired.columns = ["_".join(x).strip("_") for x in paired.columns]
    paired["delta_R_waning_minus_unwaning"] = (paired["effective_R_initial_waning"] -
                                                paired["effective_R_initial_unwaning"])
    paired["delta_attack_rate_pct"] = 100.0 * (
        paired["attack_rate_365d_waning"] / paired["attack_rate_365d_unwaning"].clip(lower=1e-12) - 1.0)
    paired.to_csv(OUT / "country_mechanism_contrasts.csv", index=False)
    (paired.groupby("disease", as_index=False)
     .agg(mean_delta_R=("delta_R_waning_minus_unwaning", "mean"),
          mean_delta_attack_rate_pct=("delta_attack_rate_pct", "mean"))
     .to_csv(OUT / "disease_mechanism_contrasts.csv", index=False))
    gaps = load_cohort_gaps()
    fallback_audit = (gaps.groupby(["country", "observation_year", "disease"], as_index=False)
                      .agg(age_rows=("age", "size"), fallback_age_rows=("profile_fallback", "sum")))
    fallback_audit["fallback_fraction"] = fallback_audit.fallback_age_rows / fallback_audit.age_rows
    fallback_audit.to_csv(OUT / "profile_fallback_audit.csv", index=False)
    print(panel.shape, sens.shape)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

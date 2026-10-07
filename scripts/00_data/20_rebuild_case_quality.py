#!/usr/bin/env python3
"""Extract WHO case reports, retaining blank fields and exact-calendar lags."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data" / "raw" / "who_immunization_portal_reported_cases.xlsx"
OUT = ROOT / "data" / "processed" / "case_data_quality_preserved.csv"
MAPPING = {
    "MEASLES": "Measles", "RUBELLA": "Rubella", "PERTUSSIS": "Pertussis",
    "DIPHTHERIA": "Diphtheria", "MUMPS": "Mumps",
}
KEY = ["country", "disease", "year"]


def extract_cases(source: pd.DataFrame) -> pd.DataFrame:
    x = source.loc[source["GROUP"].eq("COUNTRIES")].copy()
    x["disease"] = x["DISEASE"].map(MAPPING)
    x["year"] = pd.to_numeric(x["YEAR"], errors="coerce").astype("Int64")
    x = x.loc[x["disease"].notna() & x["year"].between(2000, 2025)].copy()
    x["cases"] = pd.to_numeric(x["CASES"], errors="coerce")
    invalid_text = x["CASES"].notna() & x["cases"].isna()
    if invalid_text.any():
        raise ValueError("Non-numeric WHO case values require explicit review")
    finite = x["cases"].dropna()
    if ((finite < 0) | ~np.isfinite(finite) | (finite % 1 != 0)).any():
        raise ValueError("Case reports must be nonnegative integer counts")
    x["case_status"] = np.select(
        [x["cases"].eq(0), x["cases"].gt(0)],
        ["reported_zero", "reported_positive"], default="unreported",
    )
    out = x[["CODE", "year", "disease", "cases", "case_status"]].rename(
        columns={"CODE": "country"})
    if out[KEY].isna().any().any() or out.duplicated(KEY).any():
        raise ValueError("WHO country/disease/year keys must be present and unique")
    return add_calendar_lag(out.sort_values(KEY).reset_index(drop=True))


def add_calendar_lag(frame: pd.DataFrame) -> pd.DataFrame:
    """Join t-1 explicitly; skipping an unreported year never bridges its gap."""
    out = frame.drop(columns=["cases_lag1", "case_status_lag1"], errors="ignore")
    if out.duplicated(KEY).any():
        raise ValueError("Case table contains duplicate calendar keys")
    lag = out[KEY + ["cases", "case_status"]].rename(
        columns={"cases": "cases_lag1", "case_status": "case_status_lag1"})
    lag["year"] = lag["year"] + 1
    out = out.merge(lag, on=KEY, how="left", validate="one_to_one")
    out["case_status_lag1"] = out["case_status_lag1"].fillna("absent_from_source")
    return out


def summarize(frame: pd.DataFrame, period: str, start: int, end: int,
              countries: list[str] | None = None) -> pd.DataFrame:
    x = frame.loc[frame.year.between(start, end)].copy()
    if countries is not None:
        x = x.loc[x.country.isin(countries)]
    result = x.groupby(["disease", "case_status"], as_index=False).agg(
        size=("country", "size"), n_country=("country", "nunique"))
    expected = len(countries) * (end - start + 1) if countries is not None else None
    observed = x.groupby("disease").size()
    result.insert(0, "period", period)
    result["year_start"] = start
    result["year_end"] = end
    result["expected_country_years"] = expected
    result["observed_rows_all_statuses"] = result["disease"].map(observed)
    result["source_rows_absent"] = (expected - result["observed_rows_all_statuses"]
                                     if expected is not None else pd.NA)
    return result


def main() -> None:
    out = extract_cases(pd.read_excel(SOURCE, sheet_name="Data"))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    study = pd.read_csv(ROOT / "data" / "processed" / "data_with_vaccination_gaps.csv",
                        usecols=["country"])["country"].dropna().unique().tolist()
    summary = summarize(out, "analysis_window_study_universe", 2014, 2023, study)
    summary.to_csv(OUT.with_name("case_data_quality_summary.csv"), index=False)
    summarize(out, "source_extract_window", 2000, 2025).to_csv(
        OUT.with_name("case_data_quality_summary_full_window.csv"), index=False)
    manifest = {
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sheet": "Data", "source_group": "COUNTRIES",
        "extract_years": [2000, 2025], "analysis_years": [2014, 2023],
        "n_rows": len(out), "n_country_codes": out.country.nunique(),
        "study_universe_n": len(study), "study_universe_source": "data/processed/data_with_vaccination_gaps.csv",
        "lag_definition": "same country and disease, calendar year t-1",
        "blank_fields": "unreported; retained as missing",
        "zero_fields": "reported_zero; retained as observed counts",
        "country_scope": "all WHO source country codes; the model uses the harmonized study universe",
    }
    OUT.with_name("case_data_quality_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")
    print("Wrote", OUT, out.shape)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

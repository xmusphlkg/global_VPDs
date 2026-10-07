#!/usr/bin/env python3
"""Harmonize country codes and the analysis scope.

The WPP Country/Area ISO3 list is the canonical geography. WHO country files are
restricted to that list, while the original downloads remain unchanged in
data/raw and are still the source of record.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw"
OUT = ROOT / "data/processed/harmonized"
START, END = 2000, 2025

DATASETS = [
    "coverage_data", "coverage_core_data", "coverage_wuenic_all",
    "coverage_reported_all", "coverage_core_reported", "case_data",
    "incidence_rate_data", "population_data", "births_data",
    "population_by_single_age", "population_age_groups",
    "country_year_demography",
    "vaccine_introduction_data", "vaccine_introduction_core",
    "vaccine_schedule_data", "vaccine_schedule_core",
    "monitoring_system_indicators", "coverage_survey_data",
    "measles_rubella_epi_curve_data", "measles_rubella_reporting_data",
    "country_covariates", "country_covariates_long",
]


def country_codes(path: Path) -> set[str]:
    if not path.exists():
        return set()
    try:
        frame = pd.read_csv(path, usecols=["country"], low_memory=False)
    except (ValueError, pd.errors.EmptyDataError):
        return set()
    return set(frame["country"].dropna().astype(str).str.upper().unique())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    wpp = pd.read_csv(RAW / "wpp_country_metadata.csv", dtype={"country": str})
    wpp["country"] = wpp["country"].str.upper()
    wpp_codes = set(wpp.country)

    observed: dict[str, set[str]] = {}
    all_codes: set[str] = set()
    for name in DATASETS:
        codes = country_codes(RAW / f"{name}.csv")
        observed[name] = codes
        all_codes |= codes

    rows = []
    for code in sorted(all_codes | wpp_codes):
        match = wpp.loc[wpp.country.eq(code)]
        if not match.empty:
            row = match.iloc[0]
            canonical, status, reason = code, "included_wpp_country_area", ""
            name = row.country_name
            iso2, locid = row.iso2, row.un_location_id
        else:
            canonical, status = pd.NA, "excluded_not_in_wpp_country_area"
            reason = "WHO code has no matching WPP Country/Area; preserved in raw source"
            name, iso2, locid = pd.NA, pd.NA, pd.NA
            for dataset, codes in observed.items():
                if code in codes:
                    source = pd.read_csv(RAW / f"{dataset}.csv", usecols=["country"], nrows=0)
                    del source
                    break
        rows.append({
            "native_code": code, "canonical_iso3": canonical, "iso2": iso2,
            "un_location_id": locid, "country_name_wpp": name,
            "code_status": status, "exclusion_reason": reason,
            "observed_in": ";".join(sorted(d for d, codes in observed.items() if code in codes)),
        })
    crosswalk = pd.DataFrame(rows)
    crosswalk.to_csv(OUT.parent / "country_code_crosswalk.csv", index=False)

    coverage_codes = observed.get("coverage_data", set()) & wpp_codes
    case_codes = observed.get("case_data", set()) & wpp_codes
    population_codes = observed.get("population_data", set()) & wpp_codes
    universe_codes = sorted(coverage_codes & case_codes & population_codes)
    universe = wpp[wpp.country.isin(universe_codes)].copy()
    universe["analysis_country"] = True
    universe["analysis_start_year"] = START
    universe["analysis_end_year"] = END
    universe["scope_definition"] = "WPP Country/Area ∩ WHO WUENIC coverage ∩ WHO four-disease cases ∩ WPP population"
    universe = universe.rename(columns={"country": "canonical_iso3"})
    universe.to_csv(OUT.parent / "analysis_country_universe.csv", index=False)
    scope = pd.DataFrame([
        {"field": "canonical_code", "value": "WPP 2024 ISO3_code"},
        {"field": "period_start", "value": START},
        {"field": "period_end", "value": END},
        {"field": "geography", "value": "Country/Area only; no WHO regions or aggregates"},
        {"field": "country_count", "value": len(universe_codes)},
        {"field": "excluded_codes", "value": ";".join(sorted(all_codes - wpp_codes))},
        {"field": "created_at_utc", "value": datetime.now(timezone.utc).isoformat()},
    ])
    scope.to_csv(OUT.parent / "analysis_scope.csv", index=False)

    # Create a consistent analysis-ready copy; raw downloads are not changed.
    for name in DATASETS:
        source = RAW / f"{name}.csv"
        if not source.exists():
            continue
        frame = pd.read_csv(source, low_memory=False)
        if "country" not in frame.columns:
            continue
        frame["country"] = frame["country"].astype(str).str.upper()
        frame = frame[frame.country.isin(universe_codes)].copy()
        if "year" in frame.columns:
            year = pd.to_numeric(frame.year, errors="coerce")
            frame = frame[year.between(START, END)].copy()
        if "cohort_year" in frame.columns:
            cohort = pd.to_numeric(frame.cohort_year, errors="coerce")
            frame = frame[cohort.between(START, END) | cohort.isna()].copy()
        frame.to_csv(OUT / f"{name}.csv", index=False)

    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "canonical_source": "data/raw/wpp_country_metadata.csv",
        "period": [START, END], "n_analysis_countries": len(universe_codes),
        "excluded_native_codes": sorted(all_codes - wpp_codes),
        "raw_data_preserved": True,
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(f"Canonical country universe: {len(universe_codes)} countries ({START}-{END})")
    print(f"Excluded native codes: {', '.join(sorted(all_codes - wpp_codes)) or 'none'}")
    print(f"Harmonized files: {len(list(OUT.glob('*.csv')))}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Select population anchors and data-complete extensions for transmission models."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
YEARS = tuple(range(2014, 2024))
DISEASES = ("Measles", "Rubella", "Pertussis", "Diphtheria", "Mumps")
N_POPULATION_ANCHORS = 5
N_EXTENSIONS = 5
MIN_CASE_YEARS = 7
MIN_POSITIVE_YEARS = 1
MIN_PROFILE_FRACTION = 0.95
EXPECTED_PROFILE_ROWS = len(YEARS) * 14
PREM_ALIASES = {"GBR": "United Kingdom of Great Britain", "TUR": "Turkey"}
OUT = ROOT / "output/revised/transmission_country_selection"


def _prem_sheet_names() -> set[str]:
    namespace = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    names = set()
    for suffix in (1, 2):
        path = ROOT / f"data/raw/contact_matrices_prem2017/MUestimates_all_locations_{suffix}.xlsx"
        with zipfile.ZipFile(path) as workbook:
            tree = ET.fromstring(workbook.read("xl/workbook.xml"))
        names.update(sheet.attrib["name"] for sheet in tree.iter(namespace + "sheet"))
    return names


def selection_audit() -> pd.DataFrame:
    population = pd.read_csv(
        ROOT / "data/raw/population_by_single_age.csv",
        usecols=["country", "year", "age", "population"],
    )
    population = population[population.year.isin(YEARS) & population.population.notna()]
    age_counts = population.groupby(["country", "year"]).age.nunique()
    pop_years = age_counts.eq(101).groupby("country").sum()
    pop_2023 = population[population.year.eq(2023)].groupby("country").population.sum()
    births = pd.read_csv(ROOT / "data/raw/births_data.csv", usecols=["country", "year", "births"])
    birth_years = births[births.year.isin(YEARS) & births.births.notna()].groupby("country").year.nunique()
    names = pd.read_csv(ROOT / "data/raw/wpp_country_metadata.csv", usecols=["country", "country_name"])
    audit = names.drop_duplicates("country").set_index("country")
    audit["population_2023"] = pop_2023
    audit["population_years"] = pop_years
    audit["birth_years"] = birth_years
    audit = audit.dropna(subset=["population_2023"]).copy()
    audit["prem_sheet"] = [PREM_ALIASES.get(code, name) for code, name in zip(audit.index, audit.country_name)]
    audit["has_prem_matrix"] = audit.prem_sheet.isin(_prem_sheet_names())

    cases = pd.read_csv(
        ROOT / "data/processed/case_data_quality_preserved.csv",
        usecols=["country", "year", "disease", "cases", "case_status"],
    )
    cases = cases[cases.year.isin(YEARS) & cases.disease.isin(DISEASES)]
    cases = cases[cases.case_status.ne("absent_from_source") & cases.cases.notna()]
    if cases.duplicated(["country", "year", "disease"]).any():
        raise ValueError("Case country-year-disease keys are not unique")
    case_summary = cases.groupby(["country", "disease"]).cases.agg(
        reported_years="size", positive_years=lambda values: int((values > 0).sum())
    )

    cohort_dir = ROOT / "output/five_disease_revision"
    columns = ["country", "observation_year", "age", "disease", "gap_waning_projected", "gap_no_waning_projected"]
    strict = pd.read_csv(cohort_dir / "cohort_level_all_ages.csv", usecols=columns)
    proxy = pd.read_csv(cohort_dir / "mumps_all_mcv_proxy_cohort_level.csv", usecols=columns)
    profiles = pd.concat([strict[strict.disease.ne("Mumps")], proxy], ignore_index=True)
    profiles = profiles[profiles.observation_year.isin(YEARS) & profiles.age.between(1, 14)]
    if profiles.duplicated(["country", "observation_year", "age", "disease"]).any():
        raise ValueError("Cohort country-year-age-disease keys are not unique")
    profiles["complete"] = profiles.gap_waning_projected.notna() & profiles.gap_no_waning_projected.notna()
    profile_summary = profiles.groupby(["country", "disease"]).complete.agg(
        profile_rows="size", complete_profile_rows="sum"
    )

    for disease in DISEASES:
        key = disease.lower()
        for column in ("reported_years", "positive_years"):
            audit[f"{column}_{key}"] = case_summary[column].xs(disease, level="disease").reindex(audit.index).fillna(0).astype(int)
        for column in ("profile_rows", "complete_profile_rows"):
            audit[f"{column}_{key}"] = profile_summary[column].xs(disease, level="disease").reindex(audit.index).fillna(0).astype(int)
        audit[f"profile_fraction_{key}"] = audit[f"complete_profile_rows_{key}"] / EXPECTED_PROFILE_ROWS

    audit["data_eligible"] = (
        audit.has_prem_matrix
        & audit.population_years.eq(len(YEARS))
        & audit.birth_years.eq(len(YEARS))
    )
    for disease in DISEASES:
        key = disease.lower()
        audit["data_eligible"] &= (
            audit[f"reported_years_{key}"].ge(MIN_CASE_YEARS)
            & audit[f"positive_years_{key}"].ge(MIN_POSITIVE_YEARS)
            & audit[f"profile_rows_{key}"].eq(EXPECTED_PROFILE_ROWS)
            & audit[f"profile_fraction_{key}"].ge(MIN_PROFILE_FRACTION)
        )
    audit = audit.sort_values(["population_2023", "country_name"], ascending=[False, True])
    audit["population_rank"] = range(1, len(audit) + 1)
    anchors = list(audit.head(N_POPULATION_ANCHORS).index)
    if not (audit.loc[anchors, ["has_prem_matrix", "population_years", "birth_years"]].has_prem_matrix.all()
            and audit.loc[anchors, "population_years"].eq(len(YEARS)).all()
            and audit.loc[anchors, "birth_years"].eq(len(YEARS)).all()):
        raise ValueError("A population anchor lacks a contact matrix or complete demography")
    extensions = list(audit.loc[audit.data_eligible & ~audit.index.isin(anchors)].head(N_EXTENSIONS).index)
    if len(extensions) != N_EXTENSIONS:
        raise ValueError(f"Only {len(extensions)} countries meet the extension criteria")
    selected = anchors + extensions
    audit["selection_tier"] = "excluded"
    audit.loc[anchors, "selection_tier"] = "population_anchor"
    audit.loc[extensions, "selection_tier"] = "data_complete_extension"
    audit["selection_order"] = pd.NA
    for order, country in enumerate(selected, start=1):
        audit.loc[country, "selection_order"] = order
    return audit.reset_index().rename(columns={"index": "country"})


def selected_countries() -> tuple[list[str], pd.DataFrame]:
    audit = selection_audit()
    selected = audit[audit.selection_order.notna()].sort_values("selection_order")
    return selected.country.tolist(), audit


def main() -> None:
    selected, audit = selected_countries()
    OUT.mkdir(parents=True, exist_ok=True)
    audit.to_csv(OUT / "selection_audit.csv", index=False)
    selected_rows = audit[audit.selection_order.notna()].sort_values("selection_order")
    selected_rows.to_csv(OUT / "selected_countries.csv", index=False)
    table_dir = ROOT / "output/revised/supplementary_tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    selected_rows.to_csv(table_dir / "TableS41_transmission_country_selection.csv", index=False)
    (OUT / "selection_manifest.json").write_text(json.dumps({
        "years": list(YEARS),
        "diseases": list(DISEASES),
        "population_anchor_rule": "Five highest 2023 WPP country populations; anchors retain their observed data quality",
        "extension_rule": "Five highest 2023 populations among non-anchors meeting every disease and input threshold",
        "min_reported_years_per_disease": MIN_CASE_YEARS,
        "min_positive_years_per_disease": MIN_POSITIVE_YEARS,
        "min_projected_profile_fraction_per_disease": MIN_PROFILE_FRACTION,
        "required_profile_rows_per_disease": EXPECTED_PROFILE_ROWS,
        "required_population_and_birth_years": len(YEARS),
        "contact_matrix_required": True,
        "mumps_profile": "labelled all-MCV proxy",
        "selected_countries": selected,
        "population_anchors": selected[:N_POPULATION_ANCHORS],
        "data_complete_extensions": selected[N_POPULATION_ANCHORS:],
        "eligible_extension_candidates": int((audit.data_eligible & ~audit.country.isin(selected[:N_POPULATION_ANCHORS])).sum()),
    }, indent=2) + "\n", encoding="utf-8")
    print("selected", ", ".join(selected))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Summarise aggregate external serology without creating pseudo-individual data.

The source files contain published country-by-age or survey-by-age estimates.
This script keeps each survey and assay separate and reports within-survey
age contrasts. No pooled prevalence, synthetic records, or cross-era effect is
calculated.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "data" / "external"
OUT = ROOT / "output" / "revised" / "external_validation"
OUT.mkdir(parents=True, exist_ok=True)


def percentile(values: pd.Series, q: float) -> float:
    return float(np.percentile(values.dropna().to_numpy(float), q))


def contrast_row(
    source: str,
    disease: str,
    contrast: str,
    values: pd.Series,
    unit_label: str,
    source_url: str,
    notes: str,
) -> dict[str, object]:
    values = pd.to_numeric(values, errors="coerce").dropna()
    return {
        "source": source,
        "disease": disease,
        "contrast": contrast,
        "unit": unit_label,
        "n_units": int(len(values)),
        "estimate_pp": round(float(values.median()), 3) if len(values) else np.nan,
        "p25_pp": round(percentile(values, 25), 3) if len(values) else np.nan,
        "p75_pp": round(percentile(values, 75), 3) if len(values) else np.nan,
        "min_pp": round(float(values.min()), 3) if len(values) else np.nan,
        "max_pp": round(float(values.max()), 3) if len(values) else np.nan,
        "n_positive_direction": int((values > 0).sum()),
        "n_negative_direction": int((values < 0).sum()),
        "n_zero_direction": int((values == 0).sum()),
        "source_url": source_url,
        "notes": notes,
    }


def audit_inputs() -> dict[str, object]:
    measles = pd.read_csv(EXT / "esen2_2000s" / "esen2_measles_age_bands.csv")
    rubella = pd.read_csv(EXT / "esen2_2000s" / "esen2_rubella_age_bands.csv")
    mumps = pd.read_csv(EXT / "esen2_2000s" / "esen2_mumps_age_bands.csv")
    kiggs = pd.read_csv(EXT / "kiggs_2003_2006" / "kiggs_mmr_age_bands.csv")
    pienter = pd.read_csv(EXT / "pienter2_2006_2007" / "pienter2_mumps_age_bands.csv")

    # These checks protect the synthesis from accidental table/row drift.
    assert measles.country.nunique() == 18 and len(measles) == 54
    assert rubella.country.nunique() == 17 and len(rubella) == 51
    assert mumps.country.nunique() == 18 and len(mumps) == 71
    assert set(measles.source_table) == {"Table 3"}
    assert set(rubella.source_table) == {"Table 3"}
    assert set(mumps.source_table) == {"Table 4"}
    romania = rubella[rubella.country.eq("Romania")].sort_values("age_lower")
    assert romania.seronegative_pct.tolist() == [82.2, 62.6, 36.7]
    assert set(kiggs.disease) == {"Measles", "Mumps", "Rubella"}
    assert kiggs.groupby("disease").size().eq(5).all()
    pienter_years = pienter[pienter.age_note.eq("years") & pienter.age_lower.between(5, 14)]
    assert int(pienter_years["sample_n"].sum()) > 0
    return {
        "esen2_measles_countries": int(measles.country.nunique()),
        "esen2_mumps_countries": int(mumps.country.nunique()),
        "esen2_rubella_countries": int(rubella.country.nunique()),
        "kiggs_serology_records": int(len(kiggs)),
        "pienter2_records_age_5_14": int(len(pienter_years)),
        "pienter2_sample_n_age_5_14": int(pienter_years["sample_n"].sum()),
        "checks": "passed",
    }


def esen2_rows() -> list[dict[str, object]]:
    base = EXT / "esen2_2000s"
    measles = pd.read_csv(base / "esen2_measles_age_bands.csv")
    rubella = pd.read_csv(base / "esen2_rubella_age_bands.csv")
    mumps = pd.read_csv(base / "esen2_mumps_age_bands.csv")
    rows: list[dict[str, object]] = []
    urls = {
        "Measles": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2647410/",
        "Rubella": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2647372/",
        "Mumps": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9152080/",
    }
    specs = {
        "Measles": (measles, (2, 4), (10, 19), "seronegative_pct", "seropositive"),
        "Rubella": (rubella, (2, 4), (10, 14), "seronegative_pct", "seropositive"),
        "Mumps": (mumps, (1, 4), (10, 14), "seropositive_pct", "seropositive"),
    }
    for disease, (data, young_band, old_band, value_col, _) in specs.items():
        young = data[data.age_lower.eq(young_band[0]) & data.age_upper.eq(young_band[1])]
        old = data[data.age_lower.eq(old_band[0]) & data.age_upper.eq(old_band[1])]
        x = young[["country", value_col]].rename(columns={value_col: "young"})
        y = old[["country", value_col]].rename(columns={value_col: "old"})
        z = x.merge(y, on="country", how="inner")
        # Seronegative source values are converted to seropositive percentages.
        if value_col == "seronegative_pct":
            z["young"] = 100.0 - z.young
            z["old"] = 100.0 - z.old
        values = z.old - z.young
        selected_note = "Sweden is retained with the source-selected ages and flagged in the input age_note column." if disease == "Mumps" else ""
        rows.append(contrast_row(
            "ESEN2_1996_2004", disease,
            f"{old_band[0]}-{old_band[1]} minus {young_band[0]}-{young_band[1]} age-band seropositivity",
            values, "country", urls[disease],
            "Country-level median and quartiles; no individual records or design weights. " + selected_note,
        ))
        if disease == "Mumps":
            # Sweden contributes selected ages rather than complete age bands.
            # Report a standard-band sensitivity alongside the full source set.
            standard = data[data.age_note.isna()]
            young_std = standard[standard.age_lower.eq(young_band[0]) & standard.age_upper.eq(young_band[1])]
            old_std = standard[standard.age_lower.eq(old_band[0]) & standard.age_upper.eq(old_band[1])]
            z_std = young_std[["country", value_col]].rename(columns={value_col: "young"}).merge(
                old_std[["country", value_col]].rename(columns={value_col: "old"}), on="country", how="inner"
            )
            rows.append(contrast_row(
                "ESEN2_1996_2004", disease,
                "10-14 minus 1-4 age-band seropositivity (standard age bands only)",
                z_std.old - z_std.young, "country", urls[disease],
                "Sensitivity excluding selected-age Sweden and other annotated youngest bands; no individual records or design weights.",
            ))
    return rows


def kiggs_rows() -> list[dict[str, object]]:
    data = pd.read_csv(EXT / "kiggs_2003_2006" / "kiggs_mmr_age_bands.csv")
    rows: list[dict[str, object]] = []
    url = "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0042867"
    for disease, z in data.groupby("disease", sort=False):
        young = z[z.age_lower.eq(1) & z.age_upper.eq(2)].iloc[0]
        old = z[z.age_lower.eq(14) & z.age_upper.eq(17)].iloc[0]
        for col, label in (("positive_pct", "positive"), ("negative_pct", "negative"), ("borderline_pct", "borderline")):
            rows.append({
                "source": "KiGGS_2003_2006", "disease": disease,
                "contrast": f"14-17 minus 1-2 age-band {label} percentage",
                "unit": "survey age bands", "n_units": 1,
                "estimate_pp": round(float(old[col] - young[col]), 3), "p25_pp": np.nan,
                "p75_pp": np.nan, "min_pp": np.nan, "max_pp": np.nan,
                "n_positive_direction": int(old[col] > young[col]),
                "n_negative_direction": int(old[col] < young[col]),
                "n_zero_direction": int(old[col] == young[col]),
                "source_url": url,
                "notes": "National weighted age-band percentages; positive, borderline and negative titre categories remain separate.",
            })
    return rows


def pienter_rows() -> list[dict[str, object]]:
    data = pd.read_csv(EXT / "pienter2_2006_2007" / "pienter2_mumps_age_bands.csv")
    age = data[(data.age_note.eq("years")) & data.age_lower.between(5, 14)].copy()
    # The paper places the second routine dose at age 9. Averages bracket this
    # dose and are descriptive across published age bands rather than pooled
    # population estimates.
    pre = age[age.age_lower.between(5, 8)].seropositive_pct
    post = age[age.age_lower.between(10, 14)].seropositive_pct
    age8 = float(age.loc[age.age_lower.eq(8), "seropositive_pct"].iloc[0])
    age10 = float(age.loc[age.age_lower.eq(10), "seropositive_pct"].iloc[0])
    url = "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0058234"
    rows = [contrast_row(
        "PIENTER2_2006_2007", "Mumps", "mean ages 10-14 minus mean ages 5-8 seropositivity",
        pd.Series([round(float(post.mean() - pre.mean()), 3)]), "age-band means", url,
        "Unweighted descriptive contrast of published national age-band estimates; second dose scheduled at age 9.",
    )]
    rows.append(contrast_row(
        "PIENTER2_2006_2007", "Mumps", "age 10 minus age 8 seropositivity",
        pd.Series([round(age10 - age8, 3)]), "single-year age bands", url,
        "Published national estimates; the rise brackets the routine second MMR dose.",
    ))
    return rows


def main() -> None:
    audit = audit_inputs()
    rows = esen2_rows() + kiggs_rows() + pienter_rows()
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "aggregate_serology_summary.csv", index=False)
    pd.DataFrame([audit]).to_csv(OUT / "aggregate_serology_audit.csv", index=False)
    print(out.to_string(index=False))
    print("input audit:", audit)


if __name__ == "__main__":
    main()

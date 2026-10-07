#!/usr/bin/env python3
"""Rebuild the NHANES serology construct-calibration layer."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import t


ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "data" / "external" / "nhanes_2017_2020"
COHORT = ROOT / "output" / "five_disease_revision" / "cohort_level_all_ages.csv"
MUMPS_PROXY = ROOT / "output" / "five_disease_revision" / "mumps_all_mcv_proxy_cohort_level.csv"
OUT = ROOT / "output" / "revised" / "external_validation"
OUT.mkdir(parents=True, exist_ok=True)

AGE_MIN, AGE_MAX = 6, 14
MODEL_YEAR = 2018
WEIGHT_EPSILON = 1e-6
DESIGN_AGE_MIN, DESIGN_AGE_MAX = 6, 59
VARIANCE_METHOD = "Taylor linearization of stratified PSU survey ratio; logit-t interval with domain design degrees of freedom"
ANALYTES = {
    "Measles": {"value": "SSMEASI", "quant": "SSMEASQ", "cutoff": "200 mIU/mL"},
    "Mumps": {"value": "SSMUMPI", "quant": "SSMUMPQ", "cutoff": "46 AU/mL"},
    "Rubella": {"value": "SSRUBI", "quant": "SSRUBQ", "cutoff": "10 IU/mL"},
}


def read_nhanes() -> pd.DataFrame:
    lab = pd.read_sas(EXT / "P_SSMMRV.xpt", format="xport", encoding="utf-8")
    demo = pd.read_sas(EXT / "P_DEMO.xpt", format="xport", encoding="utf-8")
    lab["SEQN"] = pd.to_numeric(lab["SEQN"], errors="coerce")
    demo["SEQN"] = pd.to_numeric(demo["SEQN"], errors="coerce")
    columns = ["SEQN", "RIDAGEYR", "SDMVPSU", "SDMVSTRA"]
    analyte_columns = [v for spec in ANALYTES.values() for v in (spec["value"], spec["quant"])]
    d = demo[columns].merge(
        lab[["SEQN", "WTSSMMPP", *analyte_columns]],
        on="SEQN", how="inner", validate="one_to_one",
    )
    d["age"] = pd.to_numeric(d["RIDAGEYR"], errors="coerce")
    d["weight_raw"] = pd.to_numeric(d["WTSSMMPP"], errors="coerce")
    # SAS zero-weight records can arrive as a tiny positive floating value.
    d["weight"] = d["weight_raw"].where(d["weight_raw"] > WEIGHT_EPSILON)
    # Preserve the entire serum-eligible frame for domain variance estimation.
    d = d[d["age"].between(DESIGN_AGE_MIN, DESIGN_AGE_MAX)].copy()
    for disease, spec in ANALYTES.items():
        d[f"{disease}_interpretation"] = pd.to_numeric(d[spec["value"]], errors="coerce")
        d[f"{disease}_positive"] = d[f"{disease}_interpretation"].map({1.0: 1.0, 2.0: 0.0})
        d[f"{disease}_quantitative"] = pd.to_numeric(d[spec["quant"]], errors="coerce")
    return d


def load_model_profile(year: int) -> pd.DataFrame:
    def load(path: Path, profile: str) -> pd.DataFrame:
        x = pd.read_csv(path, low_memory=False)
        x = x[
            (x.country == "USA")
            & (x.observation_year == year)
            & x.age.between(AGE_MIN, AGE_MAX)
        ].copy()
        x["model_profile"] = profile
        x["profile_status"] = x["state_projected"].astype(str)
        return x[
            [
                "country", "observation_year", "age", "disease",
                "gap_waning_projected", "gap_no_waning_projected",
                "cohort_size", "model_profile", "profile_status",
            ]
        ]

    strict = load(COHORT, "component_strict_projected")
    proxy = load(MUMPS_PROXY, "mumps_all_MCV_proxy")
    x = pd.concat([strict[strict.disease.ne("Mumps")], proxy], ignore_index=True)
    x["model_gap"] = pd.to_numeric(x["gap_waning_projected"], errors="coerce")
    x["model_gap_no_waning"] = pd.to_numeric(x["gap_no_waning_projected"], errors="coerce")
    x["model_protection"] = 1.0 - x["model_gap"]
    x["model_protection_no_waning"] = 1.0 - x["model_gap_no_waning"]
    return x


def weighted_estimate(frame: pd.DataFrame, outcome: str) -> dict[str, float]:
    z = frame.loc[
        frame[outcome].isin([0.0, 1.0]) & frame.weight.gt(WEIGHT_EPSILON),
        [outcome, "weight"],
    ].copy()
    if z.empty:
        return {
            "n": 0, "weighted_n": np.nan, "n_eff": np.nan,
            "estimate": np.nan, "se_kish": np.nan,
            "ci_low_kish": np.nan, "ci_high_kish": np.nan,
        }
    w = z.weight.to_numpy(float)
    y = z[outcome].to_numpy(float)
    total = w.sum()
    estimate = float(np.sum(w * y) / total)
    n_eff = float(total * total / np.sum(w * w))
    se = float(np.sqrt(max(estimate * (1.0 - estimate), 0.0) / max(n_eff, 1.0)))
    return {
        "n": int(len(z)), "weighted_n": float(total), "n_eff": n_eff,
        "estimate": estimate, "se_kish": se,
        "ci_low_kish": max(0.0, estimate - 1.96 * se),
        "ci_high_kish": min(1.0, estimate + 1.96 * se),
    }


def survey_ratio_estimate(
    frame: pd.DataFrame, outcome: str, domain: pd.Series
) -> dict[str, object]:
    """Taylor variance for a domain ratio, retaining zero-contribution PSUs.

    The full weight-positive serum frame supplies the masked strata/PSU index.
    Residual PSU totals are e_hj=sum_i w_i D_i (y_i-p)/sum_i w_i D_i.
    The with-replacement variance is sum_h m_h/(m_h-1)*sum_j(e_hj-e_hbar)^2.
    """
    design = frame.loc[frame.weight.gt(WEIGHT_EPSILON)].copy()
    if design[["SDMVSTRA", "SDMVPSU"]].isna().any().any():
        raise ValueError("Positive-weight records require complete strata and PSU codes")
    eligible = domain.reindex(design.index, fill_value=False).astype(bool)
    valid = eligible & design[outcome].isin([0.0, 1.0])
    sub = design.loc[valid]
    point = weighted_estimate(sub, outcome)
    psu_index = design[["SDMVSTRA", "SDMVPSU"]].drop_duplicates()
    design_counts = psu_index.groupby("SDMVSTRA").size()
    if design_counts.le(1).any():
        raise ValueError("Full serum design contains a singleton stratum")
    domain_index = sub[["SDMVSTRA", "SDMVPSU"]].drop_duplicates()
    domain_strata = int(domain_index.SDMVSTRA.nunique())
    domain_psu = int(len(domain_index))
    df = domain_psu - domain_strata
    p = point["estimate"]
    se, lo, hi = np.nan, np.nan, np.nan
    interval_status = "empty_domain"
    if len(sub):
        design["linearized"] = 0.0
        design.loc[valid, "linearized"] = (
            sub.weight * (sub[outcome] - p) / point["weighted_n"]
        )
        psu_totals = design.groupby(["SDMVSTRA", "SDMVPSU"])["linearized"].sum()
        variance = 0.0
        for _, values in psu_totals.groupby(level=0):
            m = len(values)
            variance += m / (m - 1) * float(np.sum((values - values.mean()) ** 2))
        se = float(np.sqrt(max(variance, 0.0)))
        if 0 < p < 1 and df > 0:
            radius = float(t.ppf(0.975, df) * se / (p * (1.0 - p)))
            lo, hi = [float(v) for v in expit([logit(p) - radius, logit(p) + radius])]
            interval_status = "logit_t"
        elif p in (0.0, 1.0):
            interval_status = "boundary_estimate_logit_interval_undefined"
        else:
            interval_status = "insufficient_domain_design_df"
    return {
        **point,
        "n_positive": int(sub[outcome].eq(1).sum()),
        "n_negative": int(sub[outcome].eq(0).sum()),
        "survey_se": se,
        "ci_low": lo,
        "ci_high": hi,
        "design_df": df,
        "domain_psu_count": domain_psu,
        "domain_strata_count": domain_strata,
        "full_design_psu_count": int(len(psu_index)),
        "full_design_strata_count": int(len(design_counts)),
        "full_design_df": int(len(psu_index) - len(design_counts)),
        "interval_status": interval_status,
    }


def serology_row(
    frame: pd.DataFrame, disease: str, age_min: int, age_max: int
) -> dict[str, object]:
    domain = frame.age.between(age_min, age_max)
    est = survey_ratio_estimate(frame, f"{disease}_positive", domain)
    return {
        "source": "CDC_NHANES_2017_March2020", "disease": disease,
        "age": age_min if age_min == age_max else np.nan,
        "age_min": age_min, "age_max": age_max,
        "age_band": f"{age_min} years" if age_min == age_max else f"{age_min}-{age_max} years",
        "endpoint": "IgG_seropositivity", "cutoff": ANALYTES[disease]["cutoff"],
        "seropositive": est["estimate"], "seronegative": 1.0 - est["estimate"],
        "seropositive_se": est["survey_se"], "survey_se": est["survey_se"],
        "seropositive_se_kish": est["se_kish"],
        "seropositive_ci_low": est["ci_low"], "seropositive_ci_high": est["ci_high"],
        "seronegative_ci_low": 1.0 - est["ci_high"],
        "seronegative_ci_high": 1.0 - est["ci_low"],
        "seropositive_ci_low_kish": est["ci_low_kish"],
        "seropositive_ci_high_kish": est["ci_high_kish"],
        "seronegative_ci_low_kish": 1.0 - est["ci_high_kish"],
        "seronegative_ci_high_kish": 1.0 - est["ci_low_kish"],
        "n": est["n"], "n_valid_assay": est["n"], "valid_weighted_assay": est["n"],
        "n_positive": est["n_positive"], "n_negative": est["n_negative"],
        "single_age_sparse_negative": age_min == age_max and est["n_negative"] < 10,
        "n_age_eligible": int(domain.sum()), "weighted_n": est["weighted_n"],
        "n_eff": est["n_eff"], "n_eff_kish": est["n_eff"],
        **{key: est[key] for key in ["design_df", "domain_psu_count", "domain_strata_count",
           "full_design_psu_count", "full_design_strata_count", "full_design_df", "interval_status"]},
        "survey_ci_method": VARIANCE_METHOD,
    }


def calibration_metrics(
    frame: pd.DataFrame, profile_scope: str, model_year: int
) -> dict[str, object]:
    z = frame.dropna(subset=["model_gap", "seronegative"]).sort_values("age")
    if profile_scope == "known_only":
        z = z[z.profile_status.eq("known")]
    if len(z) >= 3:
        rho = float(z.model_gap.corr(z.seronegative, method="spearman"))
        rmse = float(np.sqrt(np.mean((z.model_gap - z.seronegative) ** 2)))
        x, y = z.model_gap.to_numpy(float), z.seronegative.to_numpy(float)
        if np.ptp(x) > 1e-12:
            slope, intercept = np.polyfit(x, y, 1)
            predicted = intercept + slope * x
            r2 = float(1 - np.sum((y - predicted) ** 2) / np.sum((y - y.mean()) ** 2)) if np.ptp(y) > 0 else np.nan
        else:
            intercept, slope, r2 = float(y.mean()), np.nan, np.nan
        age_first, age_last = int(z.age.iloc[0]), int(z.age.iloc[-1])
        model_gradient = float(z.model_gap.iloc[-1] - z.model_gap.iloc[0])
        sero_gradient = float(z.seronegative.iloc[-1] - z.seronegative.iloc[0])
    else:
        rho = rmse = intercept = slope = r2 = model_gradient = sero_gradient = np.nan
        age_first = age_last = np.nan
    direction = "increasing" if model_gradient > 0 else "decreasing" if model_gradient < 0 else "flat"
    sero_direction = "increasing" if sero_gradient > 0 else "decreasing" if sero_gradient < 0 else "flat"
    return {
        "disease": str(frame.disease.iloc[0]), "model_year": model_year,
        "profile_scope": profile_scope, "n_age_groups": int(len(z)),
        "n_valid_assay": int(z.n_valid_assay.sum()) if len(z) else 0,
        "spearman_model_gap_vs_seronegative": rho, "rmse_gap_scale_descriptive": rmse,
        "calibration_intercept": intercept, "calibration_slope": slope, "calibration_r2": r2,
        "age_first": age_first, "age_last": age_last,
        "model_gap_last_minus_first": model_gradient, "seronegative_last_minus_first": sero_gradient,
        "model_gradient_direction": direction, "seronegative_gradient_direction": sero_direction,
        "direction_agreement": bool(direction == sero_direction and direction != "flat"),
        "survey_ci_method": VARIANCE_METHOD,
        "interpretation": "external serological construct calibration; endpoints remain distinct and agreement is descriptive",
    }


def build_age_table(d: pd.DataFrame, profile: pd.DataFrame) -> pd.DataFrame:
    sero = pd.DataFrame([
        serology_row(d, disease, age, age)
        for disease in ANALYTES for age in range(AGE_MIN, AGE_MAX + 1)
    ])
    out = sero.merge(
        profile[MODEL_COLUMNS], on=["disease", "age"], how="left", validate="one_to_one",
    )
    out["gap_minus_seronegative"] = out.model_gap - out.seronegative
    out["protection_minus_seropositive"] = out.model_protection - out.seropositive
    return out


MODEL_COLUMNS = [
    "disease", "age", "model_gap", "model_gap_no_waning", "model_protection",
    "model_protection_no_waning", "cohort_size", "model_profile", "profile_status",
]


def main() -> None:
    d = read_nhanes()
    profile = load_model_profile(MODEL_YEAR)
    age_table = build_age_table(d, profile)
    age_table.to_csv(OUT / "nhanes_age_validation.csv", index=False)

    metrics_rows = []
    for disease in ANALYTES:
        disease_table = age_table[age_table.disease.eq(disease)]
        for scope in ("all_profile", "known_only"):
            metrics_rows.append(calibration_metrics(disease_table, scope, MODEL_YEAR))
    metrics = pd.DataFrame(metrics_rows)
    metrics[metrics.profile_scope.eq("known_only")].to_csv(OUT / "nhanes_validation_metrics.csv", index=False)
    metrics.to_csv(OUT / "nhanes_validation_metrics_profile_scope.csv", index=False)

    # Broad age bands stabilize the descriptive serology summaries.
    bands = [(6, 14), (6, 8), (9, 11), (12, 14)]
    band_rows = [serology_row(d, disease, lo, hi) for disease in ANALYTES for lo, hi in bands]
    pd.DataFrame(band_rows).to_csv(OUT / "nhanes_age_band_serology.csv", index=False)

    sensitivity_rows = []
    for year in range(2017, 2021):
        prof = load_model_profile(year)
        tmp = age_table.drop(columns=[c for c in MODEL_COLUMNS if c in age_table and c not in ("disease", "age")], errors="ignore").merge(
            prof[MODEL_COLUMNS], on=["disease", "age"], how="left", validate="one_to_one"
        )
        for disease in ANALYTES:
            for scope in ("all_profile", "known_only"):
                sensitivity_rows.append(calibration_metrics(tmp[tmp.disease.eq(disease)], scope, year))
    pd.DataFrame(sensitivity_rows).to_csv(OUT / "nhanes_model_year_sensitivity.csv", index=False)

    waning_rows = []
    for disease in ANALYTES:
        z = age_table[age_table.disease.eq(disease)].copy()
        z = z[z.profile_status.eq("known")]
        for scenario, col in (("waning", "model_gap"), ("no_waning", "model_gap_no_waning")):
            valid = z.dropna(subset=[col, "seronegative"])
            waning_rows.append({
                "disease": disease, "model_year": MODEL_YEAR, "profile_scope": "known_only",
                "scenario": scenario, "n_age_groups": len(valid),
                "mean_model_gap": float(valid[col].mean()),
                "mean_gap_minus_seronegative": float((valid[col] - valid.seronegative).mean()),
                "spearman_model_gap_vs_seronegative": float(valid[col].corr(valid.seronegative, method="spearman")),
                "interpretation": "descriptive scenario contrast; serology remains an external construct",
            })
    pd.DataFrame(waning_rows).to_csv(OUT / "nhanes_waning_sensitivity.csv", index=False)

    missing = []
    for disease in ANALYTES:
        outcome = f"{disease}_positive"
        for age, group in d[d.age.between(AGE_MIN, AGE_MAX)].groupby("age", sort=True):
            valid = group[outcome].isin([0.0, 1.0]) & group.weight.gt(WEIGHT_EPSILON)
            missing.append({
                "disease": disease, "age": int(age), "age_eligible_rows": int(len(group)),
                "positive_weight_rows": int(group.weight.gt(WEIGHT_EPSILON).sum()),
                "valid_weighted_assay_rows": int(valid.sum()),
                "valid_assay_rows": int(group[outcome].isin([0.0, 1.0]).sum()),
                "zero_or_missing_weight_rows": int(group.weight.le(WEIGHT_EPSILON).sum() + group.weight.isna().sum()),
                "assay_missing_rows": int(group[outcome].isna().sum()),
            })
    pd.DataFrame(missing).to_csv(OUT / "nhanes_missingness_audit.csv", index=False)

    profile_2018 = load_model_profile(MODEL_YEAR)
    manifest = {
        "analysis": "Independent serological construct calibration against CDC/NCHS NHANES",
        "source": "US CDC/NCHS NHANES 2017-March 2020 pre-pandemic MMRV surplus serum file",
        "source_urls": {
            "serology_codebook": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_SSMMRV.htm",
            "variance_guidance": "https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx",
        },
        "age_range": [AGE_MIN, AGE_MAX], "design_frame_age_range": [DESIGN_AGE_MIN, DESIGN_AGE_MAX],
        "model_observation_year": MODEL_YEAR,
        "weight_rule": f"WTSSMMPP > {WEIGHT_EPSILON}; tiny SAS zero-weight values treated as zero",
        "assay_rule": "interpretation code 1 or 2 jointly with positive weight; missing assays excluded, never recoded as negative",
        "variance": VARIANCE_METHOD,
        "kish_status": "sensitivity only",
        "age_eligible_rows": int(len(d[d.age.between(AGE_MIN, AGE_MAX)])),
        "positive_weight_rows_full_design": int(d.weight.gt(WEIGHT_EPSILON).sum()),
        "positive_weight_rows_age_domain": int(d[d.age.between(AGE_MIN, AGE_MAX)].weight.gt(WEIGHT_EPSILON).sum()),
        "valid_weighted_assay_rows_per_antigen": {disease: int((d[d.age.between(AGE_MIN, AGE_MAX)][f"{disease}_positive"].isin([0.0, 1.0]) & d[d.age.between(AGE_MIN, AGE_MAX)].weight.gt(WEIGHT_EPSILON)).sum()) for disease in ANALYTES},
        "valid_weighted_assay_rows_full_design_per_antigen": {disease: int((d[f"{disease}_positive"].isin([0.0, 1.0]) & d.weight.gt(WEIGHT_EPSILON)).sum()) for disease in ANALYTES},
        "projected_profile_ages_2018": {disease: profile_2018.loc[(profile_2018.disease == disease) & profile_2018.profile_status.ne("known"), "age"].astype(int).tolist() for disease in ANALYTES},
        "interpretation": "calibration boundary; endpoints remain distinct and results are descriptive",
        "input_files": [str(EXT / "P_SSMMRV.xpt"), str(EXT / "P_DEMO.xpt")],
        "outputs": ["nhanes_age_validation.csv", "nhanes_validation_metrics.csv", "nhanes_validation_metrics_profile_scope.csv", "nhanes_model_year_sensitivity.csv", "nhanes_age_band_serology.csv", "nhanes_waning_sensitivity.csv", "nhanes_missingness_audit.csv"],
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    audit = {
        "age_eligible_rows": manifest["age_eligible_rows"], "positive_weight_rows_full_design": manifest["positive_weight_rows_full_design"],
        "positive_weight_rows_age_domain": manifest["positive_weight_rows_age_domain"], "valid_weighted_assay_rows_per_antigen": manifest["valid_weighted_assay_rows_per_antigen"],
        "age_disease_rows": int(len(age_table)), "projected_profile_ages_2018": manifest["projected_profile_ages_2018"],
        "survey_variance_method": VARIANCE_METHOD, "construct_comparison_decision": "retain as calibration boundary; no direct protection validation claim",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text("""# Independent NHANES serological construct calibration

This analysis compares age-specific NHANES MMRV IgG serology with a programme-derived cohort gap. These endpoints remain distinct: the model gap is a vaccine-history protection scenario, while IgG reflects detected antibody after vaccination or infection.

The design frame retains all WTSSMMPP-positive serum records aged 6–59 years. Age domains are represented by indicators so strata and PSUs outside a domain contribute zero residuals, following CDC Taylor-linearization guidance. The primary interval is a stratified PSU Taylor ratio variance with logit-t confidence limits and domain degrees of freedom equal to represented PSUs minus represented strata. Kish intervals are retained as sensitivity columns.

Aged 6–14 years includes 2,541 age-eligible records and 1,592 records with both valid positive weight and valid assay interpretation. Missing assays are excluded rather than assigned a negative result. Single-year estimates with fewer than ten negatives are flagged. Broad age bands are in `nhanes_age_band_serology.csv`; model years 2017–2020 are in `nhanes_model_year_sensitivity.csv`; waning and no-waning scenarios are descriptive in `nhanes_waning_sensitivity.csv`.

USA 2018 projected cohort ages are flagged in `profile_status`; `known_only` metrics exclude them. Mumps uses the labelled all-MCV proxy and remains exploratory because the assay threshold is a seropositivity threshold without an established correlate of protection.
""", encoding="utf-8")
    print(f"NHANES calibration rebuilt: age-domain rows={len(d[d.age.between(AGE_MIN, AGE_MAX)])}, full design rows={int(d.weight.gt(WEIGHT_EPSILON).sum())}, age-disease rows={len(age_table)}")
    print(metrics[metrics.profile_scope.eq("known_only")].to_string(index=False))


if __name__ == "__main__":
    main()

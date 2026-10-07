#!/usr/bin/env python3
"""Leakage-safe, age-complete five-disease cohort reconstruction.

The unit of reconstruction is country, observation year, and age 1--14.  Dose
coverage is joined by country and scheduled vaccination year.  A dose scheduled
after the observation year is ``not_due`` with coverage zero; a missing value or
unknown vaccine component stays ``unknown``.  Mumps component evidence is kept
in a small curated table, while the all-MCV interpretation is a labelled
sensitivity proxy.
"""
from __future__ import annotations
import json, math, os
from pathlib import Path
from typing import Iterable
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
COHORT_DIR = ROOT / "data" / "processed" / "cohort_immunity_profiles"
CURATED_DIR = Path(os.environ.get("VPDS_CURATED_DIR", ROOT / "data" / "curated"))
OUT = Path(os.environ.get("VPDS_COHORT_DIR", ROOT / "output" / "five_disease_revision"))
OUT.mkdir(parents=True, exist_ok=True)
CURATED_DIR.mkdir(parents=True, exist_ok=True)
OBS_YEARS = list(range(2014, 2024))
AGES = list(range(1, 15))
COHORT_SPEC = json.loads(Path(os.environ.get("VPDS_COHORT_SPEC", ROOT / "data/curated/five_disease_parameters.json")).read_text())
PARAMETERS = COHORT_SPEC["disease_parameters"]

def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists(): raise FileNotFoundError(path)
    return pd.read_csv(path)

def num(s): return pd.to_numeric(s,errors="coerce")

def load_curated_inputs() -> pd.DataFrame:
    """Read manually preserved component evidence without rewriting it."""
    path = CURATED_DIR / "mumps_program_evidence.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing manually curated input: {path}")
    evidence = pd.read_csv(path)
    required = {"country", "dose", "birth_start", "birth_end", "component_status",
                "target_age_months", "timing_status", "historical_source_urls"}
    missing = required.difference(evidence.columns)
    if missing:
        raise ValueError(f"Curated mumps evidence is missing columns: {sorted(missing)}")
    return evidence

def load_cohorts() -> pd.DataFrame:
    c=read_csv(COHORT_DIR/"age_stratified_cohorts.csv")
    for col in ["observation_year","age","birth_year"]: c[col]=num(c[col]).astype("Int64")
    c["cohort_size"]=num(c["cohort_size"])
    c=c[c.observation_year.isin(OBS_YEARS)&c.age.isin(AGES)].copy()
    countries=c[["country","country_name"]].drop_duplicates()
    full=pd.MultiIndex.from_product([countries.country.unique(),OBS_YEARS,AGES],names=["country","observation_year","age"]).to_frame(index=False)
    full=full.merge(countries,on="country",how="left")
    full["birth_year"]=full.observation_year-full.age
    return full.merge(c,on=["country","country_name","observation_year","age","birth_year"],how="left")

def load_coverage(filename: str, columns: Iterable[str]) -> pd.DataFrame:
    d=read_csv(COHORT_DIR/filename); d["year"]=num(d["year"]).astype("Int64")
    for col in columns:
        if col not in d: d[col]=np.nan
        d[col]=num(d[col])/100.
    return d[["country","year",*columns]].drop_duplicates(["country","year"])

def schedule_months(disease: str) -> dict[str,float]:
    """Fixed analysis scenario, independent of current schedule snapshots."""
    return {dose:months for dose,months,_ in PARAMETERS[disease]["doses"]}

def evidence_lookup(frame: pd.DataFrame, dose: int, evidence: pd.DataFrame) -> pd.DataFrame:
    """Attach curated component and target-age evidence without rowwise filters."""
    rows=[]
    for _,r in evidence[evidence.dose.eq(dose)].iterrows():
        for year in range(int(r.birth_start),int(r.birth_end)+1):
            rows.append({"country":r.country,"birth_year":year,
                         "component_status":r.component_status,
                         "target_age_months":r.target_age_months,
                         "timing_status":r.timing_status,
                         "source_url":r.historical_source_urls})
    ev=pd.DataFrame(rows)
    if ev.duplicated(["country","birth_year"]).any():
        raise ValueError("Overlapping curated country-birth-cohort evidence")
    return frame[["country","birth_year"]].merge(ev,on=["country","birth_year"],how="left",validate="many_to_one")

def merge_dose_coverage(frame: pd.DataFrame,cov: pd.DataFrame,dose: str,target_months)->pd.DataFrame:
    x=frame.copy()
    target=pd.Series(target_months,index=x.index,dtype=float) if np.isscalar(target_months) else pd.Series(target_months.to_numpy(),index=x.index,dtype=float)
    x[f"{dose}_target_age_months"]=target
    x[f"{dose}_year"]=x.birth_year.astype(float)+np.floor(target/12.)
    x=x.merge(cov[["country","year",dose]].rename(columns={"year":f"{dose}_year",dose:f"{dose}_raw_coverage"}),on=["country",f"{dose}_year"],how="left",validate="many_to_one")
    future=x[f"{dose}_year"].gt(x.observation_year)
    x[f"{dose}_coverage"]=x[f"{dose}_raw_coverage"]
    x.loc[future,f"{dose}_coverage"]=0.
    invalid=x[f"{dose}_coverage"].notna()&~x[f"{dose}_coverage"].between(0,1)
    x.loc[invalid,f"{dose}_coverage"]=np.nan
    x[f"{dose}_state"]=np.where(future,"not_due",np.where(x[f"{dose}_coverage"].notna(),"observed","unknown_coverage"))
    x.loc[x[f"{dose}_year"].isna(),f"{dose}_state"]="unknown_timing"
    return x

def component_status(country: str,birth_year: int,dose: int,e: pd.DataFrame)->str:
    r=e[e.country.eq(country)&e.dose.eq(dose)&e.birth_start.le(birth_year)&e.birth_end.ge(birth_year)]
    return str(r.iloc[0].component_status) if len(r) else "unknown_component"

def project_states(c1,c2,s1,s2,projection):
    unknown={"unknown_coverage","unknown_component","unknown_timing","transition_unknown"}
    if s1 in unknown: return np.nan,np.nan,np.nan,"unknown"
    if c2 is None:
        if not np.isfinite(c1): return np.nan,np.nan,np.nan,"unknown"
        return 1-c1,c1,0.,"known"
    if s2 in unknown or not np.isfinite(c1) or not np.isfinite(c2): return np.nan,np.nan,np.nan,"unknown"
    if s2=="not_due": c2=0.
    bad=c2>c1+1e-12
    if bad and not projection: return np.nan,np.nan,np.nan,"inconsistent"
    status="projected" if bad else "known"
    if bad: c2=min(c1,c2)
    return 1-c1,c1-c2,c2,status

def gap(p0,p1,p2,ve1,ve2,y1,y2,decay_years):
    if not np.isfinite(p0+p1+p2): return np.nan
    d1=1. if decay_years is None or not np.isfinite(y1) else math.exp(-max(y1,0.)/decay_years)
    d2=1. if decay_years is None or not np.isfinite(y2) else math.exp(-max(y2,0.)/decay_years)
    return float(np.clip(p0+p1*(1-ve1*d1)+p2*(1-(ve1 if ve2 is None else ve2)*d2),0.,1.))

def build_rows(base: pd.DataFrame,disease: str,e: pd.DataFrame,all_mcv_proxy=False)->pd.DataFrame:
    par=PARAMETERS[disease]; specs=par["doses"]; cov=load_coverage(par["coverage_file"],[x[0] for x in specs]); x=base.copy()
    d1,_,ve1=specs[0]; d2,_,ve2=specs[1] if len(specs)>1 else (None,None,None)
    for i,(dose,months,_) in enumerate(specs,1):
        if disease=="Mumps" and not all_mcv_proxy:
            ev=evidence_lookup(x,i,e)
            x[f"{dose}_component_state"]=ev.component_status.fillna("unknown_component").to_numpy()
            x[f"{dose}_timing_status"]=ev.timing_status.fillna("unknown_timing").to_numpy()
            x[f"{dose}_evidence_sources"]=ev.source_url.fillna("").to_numpy()
            target=ev.target_age_months
        else:
            target=months
            x[f"{dose}_timing_status"]="fixed_standardized_scenario"
            if disease=="Mumps":
                x[f"{dose}_component_state"]="proxy_all_MCV"
                x[f"{dose}_evidence_sources"]="analysis_assumption"
        x=merge_dose_coverage(x,cov,dose,target)
        if disease=="Mumps" and not all_mcv_proxy:
            due=x[f"{dose}_state"].ne("not_due")
            absent=x[f"{dose}_component_state"].eq("confirmed_absent")
            unknown=x[f"{dose}_component_state"].isin(["transition_unknown","unknown_component"])
            x.loc[due&absent,f"{dose}_state"]="component_absent"
            x.loc[due&absent,f"{dose}_coverage"]=0.
            x.loc[due&unknown,f"{dose}_state"]="unknown_component"
            x.loc[due&unknown,f"{dose}_coverage"]=np.nan
    c1=x[f"{d1}_coverage"].copy(); s1=x[f"{d1}_state"].copy()
    c2=x[f"{d2}_coverage"].copy() if d2 else pd.Series(0.,index=x.index)
    s2=x[f"{d2}_state"].copy() if d2 else pd.Series("component_absent",index=x.index)
    y1=x[f"{d1}_year"].where(s1.eq("observed")&c1.gt(0))
    y2=x[f"{d2}_year"].where(s2.eq("observed")&c2.gt(0)) if d2 else pd.Series(np.nan,index=x.index)
    if disease=="Mumps" and not all_mcv_proxy:
        # An absent first component followed by an observed MMR appointment is
        # a single mumps-containing dose at the second MCV appointment.
        remap=s1.eq("component_absent")&s2.eq("observed")
        c1.loc[remap]=c2.loc[remap]; y1.loc[remap]=y2.loc[remap]
        c2.loc[remap]=0.; y2.loc[remap]=np.nan
        s1.loc[remap]="observed"; s2.loc[remap]="component_absent"
    unknown=s1.str.startswith("unknown")|s2.str.startswith("unknown")|c1.isna()|c2.isna()
    bad=(c2>c1+1e-12)&~unknown
    for projection,label in [(False,"strict"),(True,"projected")]:
        c2_used=c2.clip(upper=c1) if projection else c2
        valid=~unknown if projection else (~unknown&~bad)
        x[f"state_{label}"]=np.where(unknown,"unknown",np.where(bad,"projected" if projection else "inconsistent","known"))
        x[f"p0_{label}"]=(1-c1).where(valid)
        x[f"p1_{label}"]=(c1-c2_used).where(valid)
        x[f"p2_{label}"]=c2_used.where(valid)
        t1=(x.observation_year-y1).clip(lower=0).fillna(0)
        t2=(x.observation_year-y2).clip(lower=0).fillna(0)
        p0=x[f"p0_{label}"]; p1=x[f"p1_{label}"]; p2=x[f"p2_{label}"]
        def calc_gap(tau):
            decay1=np.exp(-t1/tau) if tau else 1.
            decay2=np.exp(-t2/tau) if tau else 1.
            return p0+p1*(1-ve1*decay1)+p2*(1-(ve1 if ve2 is None else ve2)*decay2)
        x[f"gap_no_waning_{label}"]=calc_gap(None)
        x[f"gap_waning_{label}"]=calc_gap(par["waning_years"])
        if disease=="Mumps" and label=="strict":
            for tau in [16.,27.,51.]: x[f"gap_waning_{int(tau)}y_strict"]=calc_gap(tau)
    x["effective_dose1_coverage"]=c1; x["effective_dose2_coverage"]=c2
    x["effective_dose1_year"]=y1; x["effective_dose2_year"]=y2
    if disease=="Mumps":
        x["mumps_effective_dose1_year"]=y1; x["mumps_effective_dose2_year"]=y2
    for dose,_,_ in specs:
        x[f"years_since_{dose}"]=(x.observation_year-x[f"{dose}_year"]).where(x[f"{dose}_state"].eq("observed"))
    x["disease"]=disease
    x["scenario"]="all_mcv_proxy" if all_mcv_proxy else "component_strict"
    x["reconstruction_assumption"]="documented_component_with_selected_timing" if disease=="Mumps" and not all_mcv_proxy else "standardized_dose_timing_scenario"
    x["parameter_status"]="selected_index_coefficients_and_decay_scenario"
    x["n_age_target"]=14
    x["weight_available"]=np.isfinite(x.cohort_size)&x.cohort_size.gt(0)
    x["cohort_gap_total"]=x.gap_waning_strict; x["cohort_gap_no_waning"]=x.gap_no_waning_strict
    # Preserve all provenance fields, including scheduled ages and effective
    # mumps dose years.  The scheduled future year is metadata; its contribution
    # to effective vaccination coverage is exactly zero.
    return x

def wmean(g,col):
    v=num(g[col]); w=num(g.cohort_size); ok=v.notna()&w.notna()&w.gt(0)
    return float(np.average(v[ok],weights=w[ok])) if ok.any() else np.nan

def aggregate(frame):
    keys=["country","country_name","disease","scenario","observation_year"]; rows=[]
    for key,g in frame.groupby(keys,dropna=False):
        r=dict(zip(keys,key)); r["n_cohorts"]=int(g.age.nunique()); r["n_age_rows"]=int(len(g)); r["total_pop_1_14"]=float(num(g.cohort_size).sum(min_count=14)); r["n_valid_weights"]=int((np.isfinite(num(g.cohort_size))&num(g.cohort_size).gt(0)).sum()); r["weight_complete"]=bool(r["n_valid_weights"]==14 and len(g)==14)
        g2=g.copy(); g2["age_group"]=pd.cut(g2.age,[0,4,9,14],labels=["1_4","5_9","10_14"])
        for label in ["strict","projected"]:
            known=g[f"state_{label}"].isin(["known","projected"]); r[f"n_known_{label}"]=int(known.sum()); r[f"coverage_fraction_{label}"]=float(known.mean()) if len(g) else np.nan
            denom=num(g.cohort_size).sum(min_count=14); r[f"population_coverage_fraction_{label}"]=float(num(g.loc[known,"cohort_size"]).sum()/denom) if denom>0 else np.nan
            r[f"gap_no_waning_{label}"]=wmean(g,f"gap_no_waning_{label}"); r[f"gap_waning_{label}"]=wmean(g,f"gap_waning_{label}")
            for ag,sub in g2.groupby("age_group",observed=True): r[f"gap_waning_{label}_{ag}"]=wmean(sub,f"gap_waning_{label}")
            r[f"eligible_14_{label}"]=bool(r["n_cohorts"]==14 and r["weight_complete"] and r[f"n_known_{label}"]==14); r[f"eligible_80pct_{label}"]=bool(r["n_cohorts"]==14 and r["weight_complete"] and r[f"coverage_fraction_{label}"]>=.8 and r[f"population_coverage_fraction_{label}"]>=.8)
        if r["disease"]=="Mumps":
            for col in ["gap_waning_16y_strict","gap_waning_27y_strict","gap_waning_51y_strict"]: r[col]=wmean(g,col)
        # Compatibility names used by downstream model code.
        r["cohort_gap_total"]=r["gap_waning_strict"]; r["cohort_gap_no_waning"]=r["gap_no_waning_strict"]
        rows.append(r)
    return pd.DataFrame(rows)

def main():
    e=load_curated_inputs(); base=load_cohorts(); frames=[build_rows(base,d,e,False) for d in PARAMETERS]
    strict=pd.concat(frames,ignore_index=True); proxy=build_rows(base,"Mumps",e,True); levels=pd.concat([strict,proxy],ignore_index=True); cy=aggregate(levels)
    strict.to_csv(OUT/"cohort_level_all_ages.csv",index=False); proxy.to_csv(OUT/"mumps_all_mcv_proxy_cohort_level.csv",index=False); cy.to_csv(OUT/"country_year_14cohort.csv",index=False)
    cy[cy.eligible_14_strict].to_csv(OUT/"country_year_14cohort_primary.csv",index=False); cy[cy.eligible_80pct_strict].to_csv(OUT/"country_year_14cohort_80pct.csv",index=False)
    cy[(cy.disease=="Mumps")&(cy.scenario=="component_strict")].to_csv(OUT/"mumps_component_country_year.csv",index=False); cy[(cy.disease=="Mumps")&(cy.scenario=="all_mcv_proxy")].to_csv(OUT/"mumps_all_mcv_proxy_country_year.csv",index=False)
    strict.groupby(["disease","state_strict"],dropna=False).size().reset_index(name="n_rows").to_csv(OUT/"audit_state_counts.csv",index=False)
    strict[strict.state_strict.eq("inconsistent")].to_csv(OUT/"audit_inconsistent_two_dose_rows.csv",index=False)
    summary=cy.groupby(["disease","scenario"],as_index=False).agg(country_years=("country","size"),countries=("country","nunique"),complete_14=("eligible_14_strict","sum"),covered_80pct=("eligible_80pct_strict","sum"),mean_gap=("gap_waning_strict","mean"),mean_gap_projected=("gap_waning_projected","mean"))
    summary.to_csv(OUT/"revision_cohort_summary.csv",index=False)
    inputs=[COHORT_DIR/"age_stratified_cohorts.csv"]+[COHORT_DIR/PARAMETERS[d]["coverage_file"] for d in PARAMETERS]
    manifest={"analysis_observation_years":OBS_YEARS,"analysis_ages":AGES,"primary_definition":"14 age rows, all cohort weights positive, and all dose states known under strict nested coverage","sensitivity_definition":"14 age rows, all cohort weights positive, at least 80% age rows and at least 80% population weight known; projected two-dose coverage is separate","unknown_handling":"missing coverage, missing timing or unknown vaccine component remains unknown","not_due_handling":"scheduled dose after observation year is zero coverage and state not_due","parameter_status":"selected_index_coefficients_and_decay_scenarios; no coefficients fitted to reported cases","timing_status":"fixed standardized timing for global scenarios; curated timing for mumps component subset","mumps_component_sources":"USA ACIP 1998; UK Green Book/1996 two-dose schedule; Finland national two-dose programme since 1982; China CDC historical transition","input_files":[str(p.relative_to(ROOT)) for p in inputs if p.exists()],"curated_input_files":[str(p.relative_to(ROOT)) for p in [CURATED_DIR/"mumps_program_evidence.csv",CURATED_DIR/"five_disease_parameters.json"] if p.exists()],"summary":summary.to_dict(orient="records")}
    (OUT/"reconstruction_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    (OUT/"README_CN.md").write_text(f"""# 五病队列重建（revision）

观察年份为 {OBS_YEARS[0]}–{OBS_YEARS[-1]}，每个国家—年份保留 1–14 岁完整队列。剂次覆盖率按国家和计划接种年份显式合并。

- 计划年份晚于观察年份：`not_due`，覆盖率为 0。
- 覆盖率缺失或腮腺炎成分无法确认：`unknown`，不填充为 0。
- 两剂覆盖率违反嵌套关系时，`strict` 行标为不一致并排除严格主样本；`projected` 行把第二剂投影至第一剂以内。
- 主样本要求 14 个年龄行、严格剂次状态全部已知；敏感性分析要求至少 80% 年龄行已知。
- 腮腺炎严格成分证据限于 `data/curated/mumps_program_evidence.csv` 的 USA、GBR、FIN、CHN；历史程序来源和出生队列边界均保存在该表。全 MCV 结果另列 `all_mcv_proxy`。
- 腮腺炎初始应答为0.964，指数持续时间均值27.4年；16/27/51年为舍入敏感性网格。白喉部分/完整临床保护指数0.71/0.87，30年为参考衰减情景，检验10/30/100年及不衰减。百日咳5年为产品混合下的指数情景；所有定义见补充参数表。
- 全球疾病程序使用固定标准化接种年龄；腮腺炎严格子集使用证据表中的年龄场景。个人联合剂次、补种和迁移没有观测，出生年份加整年年龄映射写入每剂 `*_target_age_months` 与 `*_year`。
- `cohort_gap_total` 与 `cohort_gap_no_waning` 兼容下游模型；效力只应用一次，时间衰减只应用于效力，年份差不再扣除年龄。

主要表：`country_year_14cohort_primary.csv`、`country_year_14cohort_80pct.csv`、`cohort_level_all_ages.csv`。输入文件和定义见 `reconstruction_manifest.json`。
""",encoding="utf-8")
    print(summary.to_string(index=False))

if __name__=="__main__": main()

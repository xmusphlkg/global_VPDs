#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(tidyverse);library(glmmTMB);library(broom.mixed)})
args <- commandArgs(trailingOnly=TRUE)
IN <- if(length(args)>=1) args[[1]] else "output/five_disease_revision"
OUT <- if(length(args)>=2) args[[2]] else "output/revised"
CASE_FILE <- if(length(args)>=3) args[[3]] else NA_character_
dir.create(OUT,recursive=TRUE,showWarnings=FALSE);dir.create(file.path(OUT,"models"),showWarnings=FALSE)
raw <- read_csv("data/processed/data_with_vaccination_gaps.csv",show_col_types=FALSE)
case_status_from_count <- function(x) if_else(is.na(x), "unreported",
  if_else(x == 0, "reported_zero", "reported_positive"))
if(!is.na(CASE_FILE)) {
  q <- read_csv(CASE_FILE,show_col_types=FALSE) %>%
    select(country,disease,year,cases_quality=cases,case_status,
           cases_lag1_quality=cases_lag1,case_status_lag1)
  # Replace both the response and its one-year history from the same quality
  # table. This prevents blank WHO fields from becoming structural zeros in
  # the prior covariate, and the lag join cannot bridge an unreported year.
  raw <- raw %>%
    select(-any_of(c("cases","cases_lag1","case_status","case_status_lag1"))) %>%
    left_join(q,by=c("country","disease","year")) %>%
    rename(cases=cases_quality,cases_lag1=cases_lag1_quality)
} else {
  # Compatibility path for archived zero-filled inputs; primary analyses pass
  # CASE_FILE and therefore use the quality-preserved branch above.
  raw <- raw %>% mutate(
    case_status=case_status_from_count(cases),
    case_status_lag1=case_status_from_count(cases_lag1))
}
cy <- read_csv(file.path(IN,"country_year_14cohort.csv"),show_col_types=FALSE)
stopifnot(!anyDuplicated(raw[c("country","disease","year")]),
          !anyDuplicated(cy[c("country","disease","scenario","observation_year")]))
prepare <- function(elig="eligible_14_strict", lab="strict") {
  x <- cy %>% filter(.data[[elig]],
    (disease!="Mumps"&scenario=="component_strict")|
    (disease=="Mumps"&scenario=="all_mcv_proxy")) %>%
    transmute(country,disease,year=observation_year,scenario,n_cohorts,
      known_fraction=.data[[paste0("coverage_fraction_",lab)]],
      gap0=.data[[paste0("gap_no_waning_",lab)]],
      gap=.data[[paste0("gap_waning_",lab)]])
  raw %>% inner_join(x,by=c("country","disease","year")) %>%
    mutate(coverage_10=coverage/10,hist_10=historical_vulnerability_index/10,
      gap0_10=gap0*10,gap_10=gap*10,prior=log1p(cases_lag1),
      gdp=log(gdp_per_capita_constant_ppp),uhc=uhc_service_coverage_index/10,
      urban=urban_population_pct/10,health=health_expenditure_pct_gdp,
      time=(year-2018)/5) %>%
    filter(year>=2014,year<=2023,population>0,cases>=0,
      if_all(all_of(c("coverage_10","hist_10","gap0_10","gap_10","prior","gdp","uhc","urban","health")),is.finite))
}
# Keep an explicit missingness audit alongside the complete-case analysis. The
# model requires a current report and an observed t-1 report because prior is a
# prespecified covariate; all other complete-case losses are shown separately.
joined_audit <- function(elig="eligible_14_strict", lab="strict") {
  x <- cy %>% filter(.data[[elig]],
    (disease!="Mumps"&scenario=="component_strict")|
    (disease=="Mumps"&scenario=="all_mcv_proxy")) %>%
    transmute(country,disease,year=observation_year,scenario,
      known_fraction=.data[[paste0("coverage_fraction_",lab)]],
      gap0=.data[[paste0("gap_no_waning_",lab)]],gap=.data[[paste0("gap_waning_",lab)]])
  raw %>% inner_join(x,by=c("country","disease","year")) %>%
    mutate(coverage_10=coverage/10,hist_10=historical_vulnerability_index/10,
      gap0_10=gap0*10,gap_10=gap*10,prior=log1p(cases_lag1),
      gdp=log(gdp_per_capita_constant_ppp),uhc=uhc_service_coverage_index/10,
      urban=urban_population_pct/10,health=health_expenditure_pct_gdp,
      time=(year-2018)/5,
      current_report=case_status %in% c("reported_zero","reported_positive"),
      lag_report=case_status_lag1 %in% c("reported_zero","reported_positive"),
      covariates_complete=if_all(all_of(c("coverage_10","hist_10","gap0_10","gap_10","gdp","uhc","urban","health")),is.finite),
      model_complete=current_report & lag_report & covariates_complete & is.finite(prior))
}
joined <- joined_audit()
d <- prepare()
write_csv(d,file.path(OUT,"analysis_data.csv"))
write_csv(d%>%count(disease,scenario,name="n"),file.path(OUT,"analysis_sample_by_disease.csv"))
write_csv(tibble(n_complete=nrow(d),n_country=n_distinct(d$country),n_disease=n_distinct(d$disease),year_min=min(d$year),year_max=max(d$year),mumps_proxy=sum(d$disease=="Mumps")),file.path(OUT,"analysis_sample_summary.csv"))
write_csv(joined %>% group_by(disease) %>% summarise(
  n_joined=n(), n_country_joined=n_distinct(country),
  n_current_missing=sum(!current_report), n_lag_missing=sum(!lag_report),
  n_covariate_missing=sum(!covariates_complete), n_model_complete=sum(model_complete),
  n_country_model=n_distinct(country[model_complete]), .groups="drop"),
  file.path(OUT,"case_quality_observation_audit.csv"))

base <- cases~coverage_10+prior+gdp+uhc+urban+health+factor(year)+offset(log(population))+(1|country)
diagnostics <- list()
fit <- function(form,dat,id) {
 f <- glmmTMB(form,data=dat,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=2000,eval.max=2000)))
 if(f$fit$convergence!=0||!isTRUE(f$sdr$pdHess)) f <- glmmTMB(form,data=dat,family=nbinom2,control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"),optCtrl=list(maxit=4000)))
 ok <- f$fit$convergence==0&&isTRUE(f$sdr$pdHess)&&is.finite(AIC(f))
 diagnostics[[length(diagnostics)+1]] <<- tibble(id,n=nrow(dat),convergence=f$fit$convergence,positive_hessian=isTRUE(f$sdr$pdHess),AIC=AIC(f),ok)
 if(!ok) stop("Model failed numerical checks: ",id)
 f
}
compare <- function(a,b,dis,contrast,nested=FALSE,analysis="primary") {
 k <- attr(logLik(b),"df")-attr(logLik(a),"df"); stat <- 2*as.numeric(logLik(b)-logLik(a))
 tibble(disease=dis,comparison=contrast,analysis,n=nobs(a),aic_base=AIC(a),aic_extended=AIC(b),delta_aic=AIC(b)-AIC(a),lrt_statistic=if(nested)stat else NA_real_,df=if(nested)k else NA_real_,lrt_p_value=if(nested&&k>0)pchisq(max(stat,0),k,lower.tail=FALSE) else NA_real_)
}
# Report the estimand explicitly as a 0.10 absolute gap increase instead of
# the ambiguous pp shorthand used in earlier tables.
coefrow <- function(f,model,dis,analysis) broom.mixed::tidy(f,effects="fixed",conf.int=TRUE)%>%filter(term%in%c("hist_10","gap0_10","gap_10"))%>%mutate(gap_term=grepl("gap",term),effect_scale=if_else(gap_term,"per 0.10 absolute gap increase","per 0.10 historical exposure increase"))%>%transmute(disease=dis,model,analysis,n=nobs(f),term,estimate,std_error=std.error,p_value=p.value,rate_ratio=exp(estimate),rr_low=exp(conf.low),rr_high=exp(conf.high),effect_scale)
ladder <- function(dat,dis,tag,save=FALSE) {
 lab <- if(dis=="Mumps")paste0(tag,"_MCV_proxy")else tag
 message(lab," ",dis," n=",nrow(dat))
 f0<-fit(base,dat,paste(lab,dis,"M0"));f1<-fit(update(base,.~.+hist_10),dat,paste(lab,dis,"M1"));f2<-fit(update(base,.~.+gap0_10),dat,paste(lab,dis,"M2"))
 f3<-if(isTRUE(all.equal(dat$gap0_10,dat$gap_10)))f2 else fit(update(base,.~.+gap_10),dat,paste(lab,dis,"M3"))
 if(save)saveRDS(list(M0=f0,M1=f1,M2=f2,M3=f3),file.path(OUT,"models",paste0(tolower(dis),".rds")))
 co3<-coefrow(f3,"M3_Cohort_time",dis,lab)%>%mutate(term="gap_10")
 list(coef=bind_rows(coefrow(f1,"M1_Historical",dis,lab),coefrow(f2,"M2_Cohort_no_waning",dis,lab),co3),
  cmp=bind_rows(compare(f0,f1,dis,"M1 vs M0",TRUE,lab),compare(f0,f2,dis,"M2 vs M0",TRUE,lab),compare(f0,f3,dis,"M3 vs M0",TRUE,lab),compare(f1,f2,dis,"M2 vs M1",FALSE,lab),compare(f1,f3,dis,"M3 vs M1",FALSE,lab),compare(f2,f3,dis,"M3 vs M2",FALSE,lab)))
}
run_ladders <- function(dat,tag,save=FALSE)lapply(split(dat,dat$disease),function(x)ladder(x,x$disease[1],tag,save))
res<-run_ladders(d,"complete14",TRUE)
co<-bind_rows(lapply(res,function(x)x$coef));co$q_value<-NA_real_;ix<-co$model=="M3_Cohort_time"&co$disease!="Mumps";co$q_value[ix]<-p.adjust(co$p_value[ix],"BH")
cm<-bind_rows(lapply(res,function(x)x$cmp))
write_csv(co,file.path(OUT,"disease_model_coefficients.csv"));write_csv(cm,file.path(OUT,"disease_model_comparisons.csv"));write_csv(filter(cm,comparison%in%c("M3 vs M1","M3 vs M2")),file.path(OUT,"results_summary.csv"))
pool_formula<-cases~coverage_10+prior+gdp+uhc+urban+health+factor(year)+disease+offset(log(population))+(1|country_disease)
pools<-list()
for(scope in c("four_disease_primary","five_disease_proxy")){
 z<-if(scope=="four_disease_primary")filter(d,disease!="Mumps")else d;z$country_disease<-interaction(z$country,z$disease,drop=TRUE)
 f0<-fit(pool_formula,z,paste(scope,"M0"));f3<-fit(update(pool_formula,.~.+gap_10),z,paste(scope,"M3"));f4<-fit(update(pool_formula,.~.+gap_10*disease),z,paste(scope,"M4"))
 saveRDS(list(M0=f0,M3=f3,M4=f4),file.path(OUT,"models",paste0(scope,".rds")))
 pools[[scope]]<-bind_rows(compare(f0,f3,"All","M3 vs pooled baseline",TRUE,scope),compare(f3,f4,"All","M4 interaction vs M3",TRUE,scope))
}
write_csv(bind_rows(pools),file.path(OUT,"pooled_model_comparisons.csv"))
sens<-list()
for(setting in c("coverage80","projection14")){
 z<-if(setting=="coverage80")prepare("eligible_80pct_strict","strict")else prepare("eligible_14_projected","projected");rr<-run_ladders(z,setting);sens[[setting]]<-bind_rows(lapply(rr,function(x)x$cmp))
}
write_csv(bind_rows(sens),file.path(OUT,"sample_sensitivity_comparisons.csv"))
mz<-filter(d,disease=="Mumps");mcols<-cy%>%filter(disease=="Mumps",scenario=="all_mcv_proxy")%>%select(country,year=observation_year,starts_with("gap_waning_"));mz<-left_join(mz,mcols,by=c("country","year"))
dur<-list()
for(tau in c(16,27,51)){z<-mz%>%mutate(gap_10=.data[[paste0("gap_waning_",tau,"y_strict")]]*10);f<-fit(update(base,.~.+gap_10),z,paste("mumps_duration",tau));dur[[as.character(tau)]]<-coefrow(f,"M3_Cohort_time","Mumps","duration_MCV_proxy")%>%mutate(tau,AIC=AIC(f),mean_gap=mean(z$gap_10)/10)}
write_csv(bind_rows(dur),file.path(OUT,"mumps_duration_sensitivity.csv"))
preds<-list();hold_formula<-cases~coverage_10+prior+gdp+uhc+urban+health+time+disease+offset(log(population))+(1|country_disease)
for(scope in c("four_disease_primary","five_disease_proxy"))for(cut in c(2016,2019)){
 z<-if(scope=="four_disease_primary")filter(d,disease!="Mumps")else d;z$country_disease<-interaction(z$country,z$disease,drop=TRUE);tr<-filter(z,year<=cut);te<-filter(z,year>cut,year<=if(cut==2016)2019 else 2023)
 for(m in c("M0","M1","M2","M3")){fm<-switch(m,M0=hold_formula,M1=update(hold_formula,.~.+hist_10),M2=update(hold_formula,.~.+gap0_10),M3=update(hold_formula,.~.+gap_10));f<-fit(fm,tr,paste(scope,"holdout",cut,m));pr<-predict(f,newdata=te,type="response",allow.new.levels=TRUE);preds[[length(preds)+1]]<-tibble(scope,cutoff=cut,model=m,country=te$country,disease=te$disease,year=te$year,actual=te$cases,prediction=pr,dispersion=sigma(f),n_train=nrow(tr))}}
pred<-bind_rows(preds);write_csv(pred,file.path(OUT,"holdout_predictions.csv"))
metrics<-function(x)tibble(n_test=nrow(x),n_train=first(x$n_train),MAE=mean(abs(x$actual-x$prediction)),RMSE=sqrt(mean((x$actual-x$prediction)^2)),R2=1-sum((x$actual-x$prediction)^2)/sum((x$actual-mean(x$actual))^2),log_MSE=mean((log1p(x$actual)-log1p(x$prediction))^2),mean_log_predictive_density=mean(dnbinom(x$actual,mu=x$prediction,size=x$dispersion,log=TRUE)))
write_csv(pred%>%group_by(scope,cutoff,model)%>%group_modify(~metrics(.x))%>%ungroup(),file.path(OUT,"time_holdout_metrics.csv"))
write_csv(pred%>%group_by(scope,cutoff,disease,model)%>%group_modify(~metrics(.x))%>%ungroup(),file.path(OUT,"time_holdout_by_disease.csv"))
write_csv(bind_rows(diagnostics),file.path(OUT,"model_diagnostics.csv"));capture.output(sessionInfo(),file=file.path(OUT,"sessionInfo.txt"));message("All model fits passed numerical checks; outputs: ",OUT)

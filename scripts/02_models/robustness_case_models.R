#!/usr/bin/env Rscript
# Robustness analyses for the cohort-gap/case association.
# Adds (i) country fixed effects, (ii) pre-pandemic 2014-2019 restriction,
# and (iii) nonlinear and piecewise gap specifications. It intentionally
# leaves the primary manuscript outputs untouched.
suppressPackageStartupMessages({library(tidyverse); library(glmmTMB); library(broom.mixed)})
args <- commandArgs(trailingOnly=TRUE)
IN <- if(length(args)>=1) args[[1]] else "output/five_disease_revision"
OUT <- if(length(args)>=2) args[[2]] else "output/robustness_case_models"
CASE_FILE <- if(length(args)>=3) args[[3]] else "data/processed/case_data_quality_preserved.csv"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE); dir.create(file.path(OUT,"models"),showWarnings=FALSE)
raw <- read_csv("data/processed/data_with_vaccination_gaps.csv", show_col_types=FALSE)
q <- read_csv(CASE_FILE, show_col_types=FALSE) %>%
  select(country,disease,year,cases_quality=cases,case_status,cases_lag1_quality=cases_lag1,case_status_lag1) 
raw <- raw %>% select(-any_of(c("cases","cases_lag1","case_status","case_status_lag1"))) %>%
  left_join(q,by=c("country","disease","year")) %>% rename(cases=cases_quality,cases_lag1=cases_lag1_quality)
cy <- read_csv(file.path(IN,"country_year_14cohort.csv"),show_col_types=FALSE)
x <- cy %>% filter(eligible_14_strict,
  (disease!="Mumps" & scenario=="component_strict") | (disease=="Mumps" & scenario=="all_mcv_proxy")) %>%
  transmute(country,disease,year=observation_year,scenario,n_cohorts,
    known_fraction=coverage_fraction_strict,gap0=gap_no_waning_strict,gap=gap_waning_strict)
d <- raw %>% inner_join(x,by=c("country","disease","year")) %>%
 mutate(coverage_10=coverage/10,hist_10=historical_vulnerability_index/10,
   gap0_10=gap0*10,gap_10=gap*10,prior=log1p(cases_lag1),
   gdp=log(gdp_per_capita_constant_ppp),uhc=uhc_service_coverage_index/10,
   urban=urban_population_pct/10,health=health_expenditure_pct_gdp,time=(year-2018)/5,
   current_report=case_status %in% c("reported_zero","reported_positive"),
   lag_report=case_status_lag1 %in% c("reported_zero","reported_positive")) %>%
 filter(year>=2014,year<=2023,population>0,cases>=0,current_report,lag_report,
   if_all(all_of(c("coverage_10","hist_10","gap0_10","gap_10","prior","gdp","uhc","urban","health")),is.finite))
write_csv(d,file.path(OUT,"analysis_data.csv"))
# Formula uses the same covariates, year effects and population offset as primary.
base_rhs <- "coverage_10 + prior + gdp + uhc + urban + health + factor(year) + offset(log(population)) + (1|country)"
fit_one <- function(dat, rhs, id) {
  f <- as.formula(paste("cases ~",rhs))
  z <- tryCatch(glmmTMB(f,data=dat,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=3000,eval.max=3000))),error=function(e)e)
  if(inherits(z,"error") || z$fit$convergence!=0 || !isTRUE(z$sdr$pdHess)) {
    z <- tryCatch(glmmTMB(f,data=dat,family=nbinom2,control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"),optCtrl=list(maxit=5000))),error=function(e)e)
  }
  if(inherits(z,"error")) stop(id,": ",z$message)
  ok <- z$fit$convergence==0 && isTRUE(z$sdr$pdHess) && is.finite(AIC(z))
  if(!ok) stop("Numerical failure: ",id)
  saveRDS(z,file.path(OUT,"models",paste0(gsub("[^A-Za-z0-9_]","_",id),".rds")))
  z
}
extract <- function(f, dat, disease, analysis, term_pattern="gap") {
  broom.mixed::tidy(f,effects="fixed",conf.int=TRUE) %>% filter(grepl(term_pattern,term)) %>%
    transmute(disease,analysis,n=nobs(f),term,estimate,std_error=std.error,p_value=p.value,
      rate_ratio=exp(estimate),rr_low=exp(conf.low),rr_high=exp(conf.high),AIC=AIC(f),n_country=n_distinct(dat$country),year_min=min(dat$year),year_max=max(dat$year))
}
all_out <- list(); diag <- list()
pre_window <- list()
for(dis in sort(unique(d$disease))) {
 z <- filter(d,disease==dis); z$country_fe <- factor(z$country)
 specs <- list(
  primary=list(rhs=paste(base_rhs,"+ gap_10"),dat=z),
  country_fixed_effect=list(rhs=paste("coverage_10 + prior + gdp + uhc + urban + health + factor(year) + offset(log(population)) + gap_10 + country_fe"),dat=z),
  pre_pandemic=list(rhs=paste(base_rhs,"+ gap_10"),dat=filter(z,year<=2019)),
  nonlinear_spline=list(rhs=paste(base_rhs,"+ splines::ns(gap_10, df=3)"),dat=z),
  piecewise_median=list(rhs=paste(base_rhs,"+ gap_low + gap_high"),dat=mutate(z,gap_low=pmin(gap_10,median(gap_10)),gap_high=pmax(gap_10-median(gap_10),0)))
 )
 for(nm in names(specs)) {
  zz <- specs[[nm]]$dat; fit <- fit_one(zz,specs[[nm]]$rhs,paste(dis,nm,sep="_"));
  # Extract all gap-related terms for spline/piecewise; primary terms are comparable per 0.1 gap.
  all_out[[length(all_out)+1]] <- extract(fit,zz,dis,nm,"gap")
  diag[[length(diag)+1]] <- tibble(disease=dis,analysis=nm,n=nrow(zz),n_country=n_distinct(zz$country),year_min=min(zz$year),year_max=max(zz$year),AIC=AIC(fit),convergence=fit$fit$convergence,positive_hessian=isTRUE(fit$sdr$pdHess),formula=deparse(formula(fit)))
 }
 # Directly compare the historical-coverage and cohort-time representations
 # on identical pre-pandemic rows, using the same covariates as the primary
 # model. This is the relevant test for the pathogen contrast.
 zp <- filter(z, year <= 2019)
 if (nrow(zp) > 0) {
   m1 <- fit_one(zp, paste(base_rhs, "+ hist_10"), paste(dis, "pre_pandemic_history", sep="_"))
   m3 <- readRDS(file.path(OUT, "models", paste0(gsub("[^A-Za-z0-9_]", "_", paste(dis, "pre_pandemic", sep="_")), ".rds")))
   pre_window[[length(pre_window)+1]] <- tibble(
     disease=dis, n=nobs(m1), n_country=n_distinct(zp$country),
     AIC_M1=AIC(m1), AIC_M3=AIC(m3), delta_AIC_M3_minus_M1=AIC(m3)-AIC(m1),
     RR_M1=exp(fixef(m1)$cond[["hist_10"]]),
     RR_M3=exp(fixef(m3)$cond[["gap_10"]]))
 }
}
coef <- bind_rows(all_out); write_csv(coef,file.path(OUT,"robustness_coefficients.csv")); write_csv(bind_rows(diag),file.path(OUT,"robustness_diagnostics.csv"))
write_csv(bind_rows(pre_window), file.path(OUT, "pre_pandemic_model_comparison.csv"))
# A compact comparison for the prespecified M3 gap effect, where available.
comp <- coef %>% filter(analysis %in% c("primary","country_fixed_effect","pre_pandemic")) %>%
 mutate(effect=case_when(term=="gap_10"~"gap_10",TRUE~term))
write_csv(comp,file.path(OUT,"robustness_comparison.csv"))
# Save per-disease sample audit and reproducibility metadata.
write_csv(d %>% count(disease,year,name="n") %>% arrange(disease,year),file.path(OUT,"sample_by_disease_year.csv"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo.txt"))
writeLines(c("Input cohort directory: ",IN,paste("Input case quality file:",CASE_FILE),paste("Rows analysed:",nrow(d)),"Primary: negative-binomial GLMM, year fixed effects, country random intercept; robustness adds country fixed effects, 2014-2019 restriction, cubic natural spline and median piecewise gap."),file.path(OUT,"README.txt"))
message("Wrote robustness outputs to ",OUT)

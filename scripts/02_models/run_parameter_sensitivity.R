#!/usr/bin/env Rscript
# Robustness to the uncalibrated duration/partial-series assumptions on fixed rows.
suppressPackageStartupMessages({library(tidyverse);library(glmmTMB);library(broom.mixed)})
out <- Sys.getenv("VPDS_ANALYSIS_OUT", "output/revised")
dat <- read_csv(file.path(out,"analysis_data.csv"),show_col_types=FALSE)
ages <- read_csv(file.path(Sys.getenv("VPDS_COHORT_DIR", "output/five_disease_revision"), "cohort_level_all_ages.csv"),show_col_types=FALSE)
base <- cases~coverage_10+prior+gdp+uhc+urban+health+factor(year)+offset(log(population))+(1|country)
rows <- list()
for(dis in c("Pertussis","Diphtheria")){
 z <- filter(dat,disease==dis)
 x <- filter(ages,disease==dis,state_strict=="known")
 old <- readRDS(file.path(out,"models",paste0(tolower(dis),".rds")))
 review <- TRUE
 for(tau in if(dis=="Pertussis") c(2.5,5,10) else if(review) c(10,30,100) else c(5,10,20))
 for(e1 in if(dis=="Pertussis") c(.25,.50,.75) else if(review) c(.17,.71,.85) else c(.50,.80,.90)){
   e2 <- if(dis=="Pertussis") .85 else if(review) .87 else .95
   a <- x %>% mutate(t1=pmax(0,observation_year-coalesce(effective_dose1_year,as.numeric(observation_year))),
         t2=pmax(0,observation_year-coalesce(effective_dose2_year,as.numeric(observation_year))),
         gap_new=p0_strict+p1_strict*(1-e1*exp(-t1/tau))+p2_strict*(1-e2*exp(-t2/tau))) %>%
     group_by(country,year=observation_year) %>% summarise(gap_new=weighted.mean(gap_new,cohort_size),.groups="drop")
   z2 <- left_join(z,a,by=c("country","year")) %>% mutate(gap_10=10*gap_new)
   stopifnot(nrow(z2)==nrow(z),all(is.finite(z2$gap_10)))
   f <- glmmTMB(update(base,.~.+gap_10),data=z2,family=nbinom2,
          control=glmmTMBControl(optCtrl=list(iter.max=2000,eval.max=2000)))
   stopifnot(f$fit$convergence==0,isTRUE(f$sdr$pdHess))
   cc <- broom.mixed::tidy(f,effects="fixed",conf.int=TRUE)%>%filter(term=="gap_10")
   rows[[length(rows)+1]] <- tibble(disease=dis,n=nrow(z2),tau,partial_series_coefficient=e1,
       full_series_coefficient=e2,mean_gap=mean(z2$gap_new),
       estimate=cc$estimate,rr=exp(cc$estimate),rr_low=exp(cc$conf.low),rr_high=exp(cc$conf.high),
       p_value=cc$p.value,delta_aic_vs_historical=AIC(f)-AIC(old$M1),converged=TRUE)
 }
}
write_csv(bind_rows(rows),file.path(out,"dtp_parameter_sensitivity.csv"))
# Common historical window on exactly the primary comparison sample.
h <- list()
for(dis in c("Pertussis","Diphtheria","Measles","Rubella","Mumps")){
 z <- filter(dat,disease==dis) %>% mutate(hist_10=gap_window_10/10)
 stopifnot(all(is.finite(z$hist_10)))
 f<-glmmTMB(update(base,.~.+hist_10),data=z,family=nbinom2,control=glmmTMBControl(optCtrl=list(iter.max=2000,eval.max=2000)))
 stopifnot(f$fit$convergence==0,isTRUE(f$sdr$pdHess))
 old<-readRDS(file.path(out,"models",paste0(tolower(dis),".rds")))
 h[[dis]]<-tibble(disease=dis,n=nrow(z),historical_window=10,
     delta_aic_M3_vs_M1=AIC(old$M3)-AIC(f),converged=TRUE)
}
write_csv(bind_rows(h),file.path(out,"common_window_sensitivity.csv"))
message("Parameter and historical-window robustness fits passed.")

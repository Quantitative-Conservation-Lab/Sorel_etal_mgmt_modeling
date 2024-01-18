

#libraries
library(here)
library(TMB)
library(TMBhelper)
library(tidyverse)
library(viridisLite)
library(readxl)

#load r objects
load(here("data/dat_IPM_hatch.Rdata")) # data
load(here("data/par_IPM_hatch.Rdata")) # initial parameters
load(here("data/rand_par_IPM_hatch.Rdata")) # vector of random parameter names (not necessary given we will use TMBstan)
load(here("data/map_IPM_hatch.Rdata")) # map of parameters not to fit

setwd(here("src"))
TMB::compile("IPM_non_centered_hatchery_scenarios.cpp")
dyn.load(dynlib("IPM_non_centered_hatchery_scenarios"))

mod<-TMB::MakeADFun(data=dat_IPM,parameters = par_IPM,
                    random=(rand_par_IPM),
                    map=map_IPM,DLL ="IPM_non_centered_hatchery_scenarios",silent = FALSE,
                    inner.control = list(maxit = 1000))

save(mod,file=here("results","ipm_fit_non_centered_4_13_hatch_scen.Rdata"))


set.seed(1234)
par_IPM1<-mod$env$parList(par=rnorm(length(mod$env$last.par.best),mod$env$last.par.best,.05))
par_IPM2<-mod$env$parList(par=rnorm(length(mod$env$last.par.best),mod$env$last.par.best,.05))
par_IPM3<-mod$env$parList(par=rnorm(length(mod$env$last.par.best),mod$env$last.par.best,.05))

#The posterior sampling is computationally and time intensive. It will take several hours to complete. The results are available for download at  https://zenodo.org/records/10526151
mod$env$data$do_tmbstan<-1
stan_run_non_centered_hatch_4_12<- tmbstan::tmbstan(mod,init=list(par_IPM1,par_IPM2,par_IPM3),
                                         chains=3,cores=3,iter=2000,warmup=500,thin=3,
                                 control=list(max_treedepth =20,adapt_delta=0.999))

#Save posterior samples. This file is available for download at  https://zenodo.org/records/10526151
save(stan_run_non_centered_hatch_4_12,
     file=here("results","stan_run_non_centered_hatch_4_12.Rdata"))


#### subsample posterior samples for population projection simulations and save
list_of_draws<-as.array(stan_run_non_centered_hatch_4_12)
list_of_draws<-rbind(list_of_draws[,1,],list_of_draws[,2,],list_of_draws[,3,])
list_of_draws<-list_of_draws[,-dim(list_of_draws)[2]]
list_of_draws<-list_of_draws[seq(1,dim(list_of_draws)[1],by=180),]

save(list_of_draws,file="list_of_draws.rda")


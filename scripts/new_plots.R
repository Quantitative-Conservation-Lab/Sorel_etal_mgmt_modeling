library(here)
library(TMB)
library(TMBhelper)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(forcats)
library(viridisLite)
library(readxl)

setwd(here("src"))
TMB::compile("IPM_non_centered_hatchery_scenarios.cpp")
dyn.load(dynlib("IPM_non_centered_hatchery_scenarios"))

load(here("results","ipm_fit_non_centered_5_02_hatch_scen.Rdata"))



##load posterior samples. 
load(here("results","list_of_draws.rda"))



#*************************************
# model is fit to redds in years 1996-2019
# projection is years 2020-2069
#*************************************

# set some parameters for simulations
## tiers for allowed PNI levels
mod$env$data$NOR_tier<-matrix(c(372,278,208,176,
                                350,259,176,80)*1.2,nrow=2,byrow=TRUE)
mod$env$data$proj_years<-50 #number of projections years
dat_IPM<-mod$env$data
dat_IPM$proj_years<-50
mod$env$data$bs_prop<-c(.3,.3)# maximum proporiton of natural origin returns to take for broodstock
mod$env$data$BS_Max<-c(74,64) #hatchery broodstock 
mod$env$data$NS_restoration_scalar<-matrix(1,4,3)#rep(1,4) # start at baseline habitat
mod$env$data$DS_restoration_scalar<-matrix(1,4,3)#rep(1,4) # start at baseline habitat
mod$env$data$restoration_scalar<-matrix(1,4,3)#rep(1,4) # start at baseline habitat

sim<-mod$simulate()

# number of simulations to do for each posterior samples
n_sims_per_par<-50

sim_sum_yrs<-seq(to=(dat_IPM$last_t+dat_IPM$proj_years),by=1,length.out=10
)
## objects to hold outputs of simulations
total_years<-dat_IPM$last_t+mod$env$data$proj_years

sim_out<-sim_pHOS<-sim_p_female<-array(NA,dim=c(dim(sim$S_hat),dim(list_of_draws)[1]*n_sims_per_par))

hatchery_returns<-matrix(NA,50,dim(list_of_draws)[1]*n_sims_per_par)


sim_PNI<-sim_NOB<-sim_pNOB<-sim_BS_tot<-array(NA,dim=c(dim(sim$broodstock_proj),dim(list_of_draws)[1]*n_sims_per_par))


#indices of parameters to save
p_fem_loc<-which(names(mod$env$last.par.best)=="mu_fem")

PSM_loc<-which(names(mod$env$last.par.best)=="mu_pss")


# load environmental covariate projections for ismulations
load(file=here("data/proj_arrays_hatch_8_11_2022.Rdata"))




multipliers<-c(1.1,1.25,1.5,1.75,2)
# multipliers<-c(1.5,2)
n_mlts<-length(multipliers)
n_rest<-n_mlts*6+1
{
  rest_scalar<-array(1,c(4,3,n_rest))
  #1 = baseline
  cnt<-1
  rest_scalar[4,1,((cnt+1):(cnt+n_mlts))]<-rest_scalar[4,1,((cnt+1):(cnt+n_mlts))]*multipliers # restore Chiwawa
  cnt<-cnt+n_mlts
  rest_scalar[4,2, ((cnt+1):(cnt+n_mlts))]<-rest_scalar[4,2, ((cnt+1):(cnt+n_mlts))]*multipliers #restore Nason
  cnt<-cnt+n_mlts
  rest_scalar[4,3, ((cnt+1):(cnt+n_mlts))]<-rest_scalar[4,3, ((cnt+1):(cnt+n_mlts))]*multipliers# restore white
  cnt<-cnt+n_mlts
  rest_scalar[4,, ((cnt+1):(cnt+n_mlts))]<-rest_scalar[4,, ((cnt+1):(cnt+n_mlts))]*rep(multipliers,each=3)# restore all natal
  cnt<-cnt+n_mlts
  rest_scalar[1:3,, ((cnt+1):(cnt+n_mlts))]<-rest_scalar[1:3,, ((cnt+1):(cnt+n_mlts))]*rep(multipliers,each=9)# restore downstream
  cnt<-cnt+n_mlts
  rest_scalar[,, ((cnt+1):(cnt+n_mlts))]<-rest_scalar[,, ((cnt+1):(cnt+n_mlts))]*rep(multipliers,each=12)# restore everything
  cnt<-cnt+n_mlts
}



# object with baseline and reduced hatchery broodstock size targets
BS_size<-matrix(c(74,64),7,2,byrow = T)*seq(0,1.5,by=.25)

gc()
#main projection 

library(doParallel)
library(foreach)
cl <- makeCluster(6)
registerDoParallel(cl)
start<-Sys.time()

sim_out <- foreach::foreach(RS =1:n_rest, .combine = 'rbind',.packages = c("tidyverse","TMB")) %dopar% {
  
  out<-dplyr::tibble()
  dyn.load(dynlib("IPM_non_centered_hatchery_scenarios"))
  # for(RS in select){
  mod$env$data$restoration_scalar<-rest_scalar[,,RS]
  
  RS_i<-rest_scalar[,,RS]
  for(BS in 1:7){
    mod$env$data$BS_Max<-BS_size[BS,]#[[BS]]
    
    n_env<-dim(proj_arrays$x_proj_array)[3]      
    
    set.seed(1234)
    count<-1
    for ( i in 1:(dim(list_of_draws)[1])){
      for(j in 1:n_sims_per_par){
        
        mod$env$data$X_proj <-( proj_arrays$x_proj_array[,,j] %>% as( "dMatrix") |> as("generalMatrix") |> as("TsparseMatrix"))
        mod$env$data$X_phi_proj <- (proj_arrays$phi_proj_array[,,j] %>% as( "dMatrix") |> as("generalMatrix") |> as("TsparseMatrix"))
        set.seed(j)
        sim<-mod$simulate(par=list_of_draws[i,])
        sim_out[,,count]<-sim$S_hat   #natural-origin female spawners
        sim_NOB[,,count]<-sim$broodstock_proj
        
        sim_BS_tot[,,count]<-sim$total_broodstock[25:74,]
        
        
        hatchery_returns[,count] <- (sim$hatch_ret_a[25:74,c(2,3,5,6)]*plogis(list_of_draws[i,PSM_loc]) )|> apply(1,sum)

        
        #adult life history proportions
        # sim_LH_props[,count]<-(sim$A_tum_y[25:74,,,] %>% apply(c(3),sum)) %>% proportions
        
        
        count_ts<-0
        count_ts2<-0
        count_tsl<-0
        for(t in 1:(total_years)){
          for(s in 1:dat_IPM$n_s){
            if(t<(dat_IPM$first_t[s]+1)){next}
            
            count_ts<-count_ts+1
            sim_pHOS[t,s,count]<-sim$pHOS[count_ts]               #pHOS
            
            if(t>=(dat_IPM$first_t[s]+6)){
              count_ts2<-count_ts2+1
              sim_p_female[t,s,count]<-sim$p_female[count_ts2]       #pFemale
            }
          }
        } 
        
        
        
        # sim_pars_array[1,,count]<-plogis(list_of_draws[i,PSM_loc])
        # # 
        # sim_pars_array[2,,count]<-plogis(list_of_draws[i,RRS_loc])
        
        count<-count+1
      }
      
    }
    
    

    
    
    Strategy = ifelse(
      
      all(RS_i[,]>1), "All restored",
      ifelse(
        all(RS_i[4,]>1), "All natal only",
        ifelse(
          all(RS_i[1:3,]>1), "Downstream only",
          ifelse(
            all(RS_i[4,1]>1) ,"Chiwawa only",
            ifelse(
              all(RS_i[4,2]>1) , "Nason only",
              ifelse(
                all(RS_i[4,3]>1) ,"White only",
                "Baseline"
              ))))))   
    
    multiplier<-ifelse(Strategy=="Baseline",1,
                       mean(RS_i[RS_i>1])
    )
    
    
    BS_scenario <-  mod$env$data$BS_Max[1]/74
    
    
    sim_S_wild<-sim_out*(1-sim_pHOS)
    
    sim_S_hatch_all<-(sim_out-sim_S_wild)/sim_p_female
    
    sim_S_wild_all<-sim_S_wild/sim_p_female #add males
    
    
    sim_NOB_w_tot<-sim_NOB %>%abind::abind(apply(sim_NOB,c(1,3),sum),along=2)
    
    sim_BS_w_tot<-sim_BS_tot %>%abind::abind(apply(sim_BS_tot,c(1,3),sum),along=2)
    
    sim_pNOB_w_tots<-sim_NOB_w_tot/sim_BS_w_tot
    
    
    
    
    
    
    
    
    dat_w_tots<-sim_S_wild_all %>%abind::abind(apply(sim_S_wild_all,c(1,3),sum),along=2)
    
    hatch_dat_w_tots<-sim_S_hatch_all %>%abind::abind(apply(sim_S_hatch_all,c(1,3),sum),along=2)
    
    
    pHOS_w_tots<-hatch_dat_w_tots/(dat_w_tots+hatch_dat_w_tots)
    
    
    sim_PNI_w_tots<-(sim_pNOB_w_tots)/(pHOS_w_tots[25:74,c(1:2,4),]+sim_pNOB_w_tots)
    
    
    # geometric mean wild spawners in year 20-50 of projection
    geo_mean<-apply(dat_w_tots[sim_sum_yrs,,],2:3,function(x){
      exp(sum(log(x)/length(x)))-.01
    })
    
    # geometric mean hatchery spawners in year 20-50 of projection
    geo_mean_hatch<-apply(hatch_dat_w_tots[sim_sum_yrs,,],2:3,function(x){
      exp(sum(log(x+.01)/length(x)))-.01
    })

    # geometric mean hatchery returns in year 20-50 of projection
    geo_mean_hatch_ret<-apply(hatchery_returns [41:50,],2,function(x){
      exp(sum(log(x+.01)/length(x)))-.01
    })

    #  mean pHOS 20-50 of projection
    mean_pHOS<-apply(pHOS_w_tots[sim_sum_yrs,,],2:3,function(x){
      mean(x)
    })
    
    #  mean pHOS 20-50 of projection
    mean_pNOB<-apply(sim_pNOB_w_tots[41:50,,],2:3,function(x){
      mean(x)
    })
    
    #  mean pHOS 20-50 of projection
    mean_PNI<-apply(sim_PNI_w_tots[41:50,,],2:3,function(x){
      mean(x)
    })
    
    
    
    prop_NA<-sum(is.na(geo_mean))/length(geo_mean)
    
    #summarized  by stream
    geo_mean_sum<-apply(geo_mean,1,quantile,probs=c(.05,.25,.5,.75,.95),na.rm=T) %>% t() %>% as_tibble()  %>% `colnames<-` (1:5)%>% mutate(stream=c("Chiwawa","Nason","White","Total"),                                                                                  stream=fct_relevel(stream,c("Chiwawa","Nason","White","Total"))) %>% mutate(var="abund",Hatchery=BS_scenario,Strategy=Strategy,multiplier=multiplier,propNA=prop_NA)
    
    out<-bind_rows(out,geo_mean_sum)
    
    
    #hatchery spawners
    geo_mean_sum_hatch<-apply(geo_mean_hatch,1,quantile,probs=c(.05,.25,.5,.75,.95),na.rm=T) %>% t() %>% as_tibble()  %>% `colnames<-` (1:5)%>% mutate(stream=c("Chiwawa","Nason","White","Total"),                                                                                  stream=fct_relevel(stream,c("Chiwawa","Nason","White","Total"))) %>% mutate(var="hatch_abund",Hatchery=BS_scenario,Strategy=Strategy,multiplier=multiplier,propNA=prop_NA)
    
    out<-bind_rows(out,geo_mean_sum_hatch)
    
    
    
    # hatchery return
    geo_mean_sum_ret<-quantile(geo_mean_hatch_ret,probs=c(.05,.25,.5,.75,.95),na.rm=T) %>% t() %>% as_tibble()  %>% `colnames<-` (1:5) %>% mutate(var="hatch_ret",Hatchery=BS_scenario,Strategy=Strategy,multiplier=multiplier,propNA=prop_NA)
    
    out<-bind_rows(out,geo_mean_sum_ret)
    
    # pHOS
    mean_pHOS2<-apply(mean_pHOS,1,quantile,probs=c(.05,.25,.5,.75,.95),na.rm=T) %>% t() %>% as_tibble()  %>% `colnames<-` (1:5)%>% mutate(stream=c("Chiwawa","Nason","White","Total"),                                                                                  stream=fct_relevel(stream,c("Chiwawa","Nason","White","Total"))) %>% mutate(var="pHOS",Hatchery=BS_scenario,Strategy=Strategy,multiplier=multiplier,propNA=prop_NA)
    
    out<-bind_rows(out,mean_pHOS2)
    
    #pNOB
    mean_pNOB2<-apply(mean_pNOB,1,quantile,probs=c(.05,.25,.5,.75,.95),na.rm=T) %>% t() %>% as_tibble()  %>% `colnames<-` (1:5)%>% mutate(stream=c("Chiwawa","Nason","Total"),                                                                                  stream=fct_relevel(stream,c("Chiwawa","Nason","Total"))) %>% mutate(var="pNOB",Hatchery=BS_scenario,Strategy=Strategy,multiplier=multiplier,propNA=prop_NA)
    
    
    out<-bind_rows(out,mean_pNOB2)
    # PNI
    mean_PNI2<-apply(mean_PNI,1,quantile,probs=c(.05,.25,.5,.75,.95),na.rm=T) %>% t() %>% as_tibble()  %>% `colnames<-` (1:5)%>% mutate(stream=c("Chiwawa","Nason","Total"),                                                                                  stream=fct_relevel(stream,c("Chiwawa","Nason","Total"))) %>% mutate(var="PNI",Hatchery=BS_scenario,Strategy=Strategy,multiplier=multiplier,propNA=prop_NA)
    
    
    out<-bind_rows(out,mean_PNI2)
    
  
  }
  out
}
# }
end<-Sys.time()
end-start
stopCluster(cl)

save(sim_out,file=here("results",
                       "sim_out_5_03_25_lots_hatch.rda"))



rm(list=ls())
gc()
load(file=here("results","sim_out_5_03_25_lots_hatch.rda"))

# Get unique scenarios excluding "baseline"
non_baseline <- sim_out %>% filter(Strategy != "Baseline") 

baseline <- sim_out %>% filter(Strategy == "Baseline") |> select(-Strategy)



new_summarized_sims<-non_baseline |> 
  bind_rows(crossing(baseline,tibble(Strategy=unique(non_baseline$Strategy)))) |> 
  mutate(Strategy=fct_relevel(Strategy,c("All restored",
                                         "All natal only",
                                         "Downstream only",
                                         "Chiwawa only",
                                         "Nason only",
                                         "White only"
  )))

new_summarized_sims|> filter(var=="abund",stream=="Total",multiplier!=1.1,Hatchery<=1.5) |> 
  ggplot(aes(x=multiplier,fill=factor(Hatchery)))+
  geom_bar(aes(y=`3`/1000),size=3,position = "dodge",stat="identity")+facet_wrap(~Strategy)+geom_linerange(aes(ymax=`4`/1000,ymin=`2`/1000,color=factor(Hatchery)),lwd=1,color="black",position = position_dodge(width = 0.22))+scale_x_continuous(breaks=c(1,1.5,2))+scale_fill_manual(values = c(colorRampPalette(c("lightblue","black","darkgreen"))(7)))+ylab("Geomean natural-origin spawners")+xlab("Survival increase %")


new_summarized_sims|> filter(var=="hatch_abund",stream=="Total",multiplier!=1.1,Hatchery<=1.5) |> 
  ggplot(aes(x=multiplier,fill=factor(Hatchery)))+
  geom_bar(aes(y=`3`/1000),size=3,position = "dodge",stat="identity")+facet_wrap(~Strategy)+geom_linerange(aes(ymax=`4`/1000,ymin=`2`/1000,color=factor(Hatchery)),lwd=1,color="black",position = position_dodge(width = 0.22))+scale_x_continuous(breaks=c(1,1.5,2))+scale_fill_manual(values = c(colorRampPalette(c("lightblue","black","darkgreen"))(7)))+ylab("Geomean hatchery-origin spawners")+xlab("Survival increase %")


new_summarized_sims|> filter(var=="pHOS",stream=="Total",multiplier!=1.1,Hatchery<=1.5) |> 
  ggplot(aes(x=multiplier,fill=factor(Hatchery)))+
  geom_bar(aes(y=`3`*100),size=3,position = "dodge",stat="identity")+facet_wrap(~Strategy)+geom_linerange(aes(ymax=`4`*100,ymin=`2`*100,color=factor(Hatchery)),lwd=1,color="black",position = position_dodge(width = 0.22))+scale_x_continuous(breaks=c(1,1.5,2))+scale_fill_manual(values = c(colorRampPalette(c("lightblue","black","darkgreen"))(7)))+ylab("pHOS")+xlab("Survival increase %")

new_summarized_sims|> filter(var=="pNOB",stream=="Total",multiplier!=1.1,Hatchery<=1.25) |> 
  ggplot(aes(x=multiplier,fill=factor(Hatchery)))+
  geom_bar(aes(y=`3`*100),size=3,position = "dodge",stat="identity")+facet_wrap(~Strategy)+geom_linerange(aes(ymax=`4`*100,ymin=`2`*100,color=factor(Hatchery)),lwd=1,color="black",position = position_dodge(width = 0.22))+scale_x_continuous(breaks=c(1,1.5,2))+scale_fill_manual(values = c(colorRampPalette(c("lightblue","black","darkgreen"))(6)))+ylab("pNOB")+xlab("Survival increase %")

new_summarized_sims|> filter(var=="PNI",stream=="Total",multiplier!=1.1,Hatchery<=1.25) |> 
  ggplot(aes(x=multiplier,fill=factor(Hatchery)))+
  geom_bar(aes(y=`3`*100),size=3,position = "dodge",stat="identity")+facet_wrap(~Strategy)+geom_linerange(aes(ymax=`4`*100,ymin=`2`*100,color=factor(Hatchery)),lwd=1,color="black",position = position_dodge(width = 0.22))+scale_x_continuous(breaks=c(1,1.5,2))+scale_fill_manual(values = c(colorRampPalette(c("lightblue","black","darkgreen"))(6)))+ylab("Geomean hatchery-origin spawners")+xlab("Survival increase %")+geom_hline(aes(yintercept=67))


new_summarized_sims|> filter(Hatchery==1,var=="abund",stream=="Total") |> 
  ggplot(aes(x=(multiplier-1)*100,color=Strategy))+#geom_ribbon(aes(ymax=`4`,ymin=`2`),fill="gray")+
  geom_point(aes(y=`3`))+#facet_wrap(~Strategy)+
  geom_line(aes(y=`3`))+ylab("Geomean natural-origin spawners")+xlab("Survival increase %")+xlim(0,200)+ylim(0,5000)


facet_labels <- data.frame(
  Strategy = factor(levels(new_summarized_sims$Strategy),levels=levels(new_summarized_sims$Strategy)),
  label = LETTERS[1:length(levels(new_summarized_sims$Strategy))] # A, B, C, ...
) 


top<-new_summarized_sims|> filter(Hatchery==1,var=="abund",stream=="Total",Strategy%in%c("All restored","All natal only","Downstream only")) |> 
  ggplot(aes(x=(multiplier)))+geom_linerange(aes(ymax=`4`/1000,ymin=`2`/1000),lwd=1.1)+
  geom_point(aes(y=`3`/1000),size=3)+facet_wrap(~Strategy)+
  #geom_line(aes(y=`3`/1000))+
  ylab("")+theme(text = element_text(size = 16))+ 
  geom_text(data = facet_labels[1:3,], aes(x = 1.1, y = 5, label = label), inherit.aes = FALSE,size=6)+
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+scale_x_continuous(breaks=c(1,1.5,2))

bottom<-new_summarized_sims|> filter(Hatchery==1,var=="abund",stream=="Total",!Strategy%in%c("All restored","All natal only","Downstream only")) |> 
  ggplot(aes(x=(multiplier)))+geom_linerange(aes(ymax=`4`/1000,ymin=`2`/1000),lwd=1.1)+
  geom_point(aes(y=`3`/1000),size=3)+facet_wrap(~Strategy)+
  #geom_line(aes(y=`3`/1000))+
  ylab("")+xlab("")+theme(text = element_text(size = 16))+ 
  geom_text(data = facet_labels[4:6,], aes(x = 1.1, y = 1.75, label = label), inherit.aes = FALSE,size=6)+scale_x_continuous(breaks=c(1,1.5,2))


library(cowplot)
combined <- plot_grid(
  top, bottom,
  ncol = 1,
  align = "v",
  axis = "lr",   # align left and right axes
  rel_heights = c(.85, 1)
)

ggdraw() +
  draw_plot(combined, 0, 0, 1, 1) +
  draw_label("Survival multiplier", x = 0.5, y = 0, vjust = -.65, size = 16) +
  draw_label("Spawners (thousands)", x = 0, y = 0.5, angle = 90, vjust = 1.2, size = 16)



ggsave("Figure 3.png")


facet_labels2 <- data.frame(
  Hab_mult = as.factor(sort(unique(new_summarized_sims$multiplier))),
  label = LETTERS[1:length(unique(new_summarized_sims$multiplier))] # A, B, C, ...
) 


top4<-new_summarized_sims|> filter(Strategy=="All restored",var=="abund",stream=="Total",multiplier<1.5) |> mutate(Hab_mult=as.factor(multiplier)) |> ggplot(aes(x=Hatchery))+geom_linerange(aes(ymax=`4`/1000,ymin=`2`/1000),lwd=1.05)+
  geom_point(aes(y=`3`/1000),size=2)+facet_wrap(~Hab_mult)+
  # geom_line(aes(y=`3`/1000))+
  ylab("")+xlab("")+theme(text = element_text(size = 16))+
  theme(
    axis.title.x =   element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+ 
  geom_text(data = facet_labels2[1:3,], aes(x = 0.1, y = 1.55, label = label), inherit.aes = FALSE,size=6)
  
  


  
bottom4<-  new_summarized_sims|> filter(Strategy=="All restored",var=="abund",stream=="Total",multiplier>=1.5) |> mutate(Hab_mult=as.factor(multiplier)) |> ggplot(aes(x=Hatchery))+geom_linerange(aes(ymax=`4`/1000,ymin=`2`/1000),lwd=1.05)+
  geom_point(aes(y=`3`/1000),size=2)+facet_wrap(~Hab_mult)+
  # geom_line(aes(y=`3`/1000))+
  ylab("")+xlab("Hatchery size proportion")+theme(text = element_text(size = 16))+ 
  geom_text(data = facet_labels2[4:6,], aes(x = 0.1, y = 5.25, label = label), inherit.aes = FALSE,size=6)

combined4 <- plot_grid(
  top4, bottom4,
  ncol = 1,
  align = "v",
  axis = "lr",   # align left and right axes
  rel_heights = c(.85, 1)
)

ggdraw() +
  draw_plot(combined4, 0, 0, 1, 1) +
  # draw_label("Survival multiplier", x = 0.5, y = 0, vjust = -.5, size = 14) +
  draw_label("Spawners (thousands)", x = 0, y = 0.5, angle = 90, vjust = 1.25, size = 16)


  
new_summarized_sims|> filter(Strategy=="All restored",var=="abund",stream=="Total") |> mutate(`Habitat\nsuvival\nmultiplier`=as.factor(multiplier)) |> ggplot(aes(x=Hatchery,color=`Habitat\nsuvival\nmultiplier`))+#geom_ribbon(aes(ymax=`4`,ymin=`2`),fill="gray")+
  geom_point(aes(y=`3`))+#facet_wrap(~Strategy)+
  geom_line(aes(y=`3`))+ylab("Geomean natural-origin spawners")+xlab("Hatchery multiplier")

new_summarized_sims|> filter(Strategy=="All restored",var=="pHOS",stream=="Total",multiplier<=2.75) |> mutate(`Habitat\nsuvival\nmultiplier`=as.factor(multiplier)) |> ggplot(aes(x=Hatchery))+geom_linerange(aes(ymax=`4`,ymin=`2`),lwd=1.05)+
  geom_point(aes(y=`3`),size=2)+facet_wrap(~`Habitat\nsuvival\nmultiplier`)+
  # geom_line(aes(y=`3`/1000))+
  ylab("pHOS")+xlab("Hatchery size proportion")+theme(text = element_text(size = 16))

new_summarized_sims|> filter(Strategy=="All restored",var=="pNOB",stream=="Total",multiplier<=2.75) |> mutate(`Habitat\nsuvival\nmultiplier`=as.factor(multiplier)) |> ggplot(aes(x=Hatchery))+geom_linerange(aes(ymax=`4`,ymin=`2`),lwd=1.05)+
  geom_point(aes(y=`3`),size=2)+facet_wrap(~`Habitat\nsuvival\nmultiplier`,scales="free_y")+
  # geom_line(aes(y=`3`/1000))+
  ylab("pNOB")+xlab("Hatchery size proportion")+theme(text = element_text(size = 16))


new_summarized_sims|> filter(Strategy=="All restored",var=="PNI",stream=="Total",multiplier<=2.75) |> mutate(`Habitat\nsuvival\nmultiplier`=as.factor(multiplier)) |> ggplot(aes(x=Hatchery))+geom_linerange(aes(ymax=`4`,ymin=`2`),lwd=1.05)+
  geom_point(aes(y=`3`),size=2)+facet_wrap(~`Habitat\nsuvival\nmultiplier`,scales="free_y")+
  # geom_line(aes(y=`3`/1000))+
  ylab("PNI")+xlab("Hatchery size proportion")+theme(text = element_text(size = 16))





new_summarized_sims |> filter(Hatchery==1,var=="abund",stream=="Total") |> select(x=`3`,multiplier,Strategy) |> group_by(Strategy) |> mutate(diff_from_baseline=x-min(x)) |>
  View()

new_summarized_sims |> filter(Strategy=="All restored",var=="abund",stream=="Total") |> select(x=`3`,multiplier,Hatchery) |> group_by(multiplier) |> arrange(multiplier,Hatchery) |> 
  mutate(diff_from_baseline=x-min(x),
         prop_inc=round((diff_from_baseline/min(x))*100),
         diff2=c(NA,diff(diff_from_baseline))) |> 
  View()

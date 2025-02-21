library(RColorBrewer)
library(mizer)
library('mizerHowTo')
library('mizerExperimental')
library(foreach)
library(ggplot2)
library(ggplotFL)
library(reshape2)
library(assertthat)
library(foreach)
library(snowfall)
library(ggrepel)
library(dplyr) 
library(tidyverse) 
library(viridis)
library(cowplot)
library(ggmap)
library(ggspatial)
library(sf)

#install.packages("ggplotFL", repos="http://flr-project.org/R")
# remotes::install_github("sizespectrum/mizerExperimental",force = TRUE)
# remotes::install_github("sizespectrum/mizerHowTo",force = TRUE)
#setwd("~/Desktop/R/R3")
#setwd("E:/Wjia/R3")
#setwd("~/JiaW/Rsun")
load("tunereprovalue.Rdata")#params_cali,calibration
source("summarySE.R")

# new base model------------
paramsbase<-upgradeParams(setInitialValues(params_tunerepro,mizer_tunerepro))
initial_effort <- array(data=rep(1.585998,200),dim=c(200,1),
                        dimnames =list(time=1:200,gear='trawl'))
mizerbase<- project(paramsbase, effort = initial_effort, dt = 0.1, t_save = 1,rm=FALSE)
spnames <- species_params(paramsbase)$species
mizer::plot(mizerbase)

# read information of names and biomass -----------------------------------
scitocomnames <- read_csv("target nontarget fish names.csv") %>% 
  rename("species"="Abbreviation") %>% 
  select(-c('Type','Scientific name')) %>% 
  rename("commonnames"="Common name")

biomass <- read_csv("fish biomass observation.csv") %>% 
  left_join(as_tibble(colMeans(tail(getBiomass(mizerbase),10)/10^6),rownames= 'Species')) %>% 
  rename('Calibrated biomass (t)' = 'value') 

#random recruitment--------
BevertonHoltRDD<-function (rdi, species_params,rm, ...) 
{
  if (missing(rm)) {
    rm <-FALSE }
  if (!("R_max" %in% names(species_params))) {
    stop("The R_max column is missing in species_params.")
  }
  if(isTRUE(rm))
  {return(rdi/(1 + rdi/species_params$R_max)*rlnorm(21,log(1)-0.2^2/2,0.2))}
  return(rdi/(1 + rdi/species_params$R_max))
}

#setInitialvalues2 function------------
setInitialValues2<-function (params, sim,no_t) 
{
  assertthat::assert_that(is(params, "MizerParams"), is(sim, "MizerSim"))
  no_t <- dim(sim@n)[1]
  if (missing(no_t)) 
    no_t <- dim(sim@n)[1]
  if (!identical(dim(sim@n)[2:3], dim(params@initial_n))) {
    stop("The consumer size spectrum of the simulation in `sim` has a ", 
         "different size from that in `params`.")
  }
  if (!identical(length(sim@n_pp[no_t, ]), length(params@initial_n_pp))) {
    stop("The resource size spectrum of the simulation in `sim` has a ", 
         "different size from that in `params`.")
  }
  if (!identical(length(sim@n_other[no_t, ]), length(params@initial_n_other))) {
    stop("The number of other components in the simulation in `sim` is ", 
         "different from that in `params`.")
  }
  if (!identical(length(sim@effort[no_t, ]), length(params@initial_effort))) {
    stop("The number of gears in the simulation in `sim` is ", 
         "different from that in `params`.")
  }
  if (!identical(dimnames(sim@effort)[[2]], names(params@initial_effort))) {
    stop("The gears in the simulation in `sim` have different names ", 
         "from those in `params`.")
  }
  params@initial_n[] <- sim@n[no_t, , ]
  params@initial_n_pp[] <- sim@n_pp[no_t, ]
  params@initial_n_other[] <- sim@n_other[no_t, ]
  params@initial_effort[] <- sim@effort[no_t, ]
  params@time_modified <- lubridate::now()
  params
}

#getYield2 function----------
getYield2<-function (sim,step) 
{
  if (missing(step)) 
    step <- 10
  assertthat::assert_that(is(sim, "MizerSim"))
  biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw, 
                   "*")
  f <- getFMort(sim, drop = FALSE)/step#each row =year-f
  yield <- apply(f * biomass, c(1, 2), sum)
  return(yield)
}
#dynamic tac by species mode --------
eisstac<-function(run){
  tac<-as.array(colMeans(tail(getYield(mizerbase),10))*0.8 )*rlnorm(21,log(1)-0.2^2/2,0.2)#initial tac for the first year
  step=50;year=50 #step is t_save/dt*2 in each year; year is length of simulation
  targetsp<-c("P.polyactis","O.ochellatus","S.niphonius","O.oratoria")
  nontargetsp<-spnames[-match(targetsp,spnames)]
  tacnon<-sum(colMeans(tail(getYield(mizerbase),10))[nontargetsp]) #total nontarget catch of base model
  eisstacrun<-foreach(t=1:year,.combine = 'rbind',.errorhandling = "remove") %do%{
    startt<-t-1
    if (t==1){lastmizer<-paramsbase;effortini=3}
    mizeryear1<-project(lastmizer,effort =effortini,t_start=startt,t_save=1/step,t_max = 1,dt=1/step/2,rm=FALSE)
    yield_step1 <-getYield2(mizeryear1,step)
    #identical(yield,getYield2(mizeryear1))
    yieldsum1<-apply(tail(yield_step1,step),2,cumsum)#cumsum from time t+1/step
    
    #decided the stop time in a year
    nnf<-foreach(i=1:length(targetsp),.combine = c) %do% {
      target<-targetsp[i]
      nof<-max(which(yieldsum1[,target]<tac[target]))+(t-1)*step #from 1/step
      if(nof==-Inf) {nof<-(t-1)*step+1}
      names(nof)<-target
      return(nof)
    }
    nnf2<-max(which(rowSums(yieldsum1[,nontargetsp])<tacnon))+(t-1)*step 
    names(nnf2)<-"nontargetsp"
    nf<-min(c(nnf,nnf2))
    chokesp<-names(which.min(c(nnf,nnf2)))
    
    #change the mizeryear1 to mizeryear1.5 
    mizeryear1.5<-mizeryear1
    mizeryear1.5@params<-setInitialValues2(paramsbase,mizeryear1,no_t=nf+1)#change initial_(n,n_pp,n_other,effort)
    mizeryear1.5@effort<-head(mizeryear1@effort,nf+1)
    mizeryear1.5@n<-head(mizeryear1@n,nf+1)
    mizeryear1.5@n_pp<-head(mizeryear1@n_pp,nf+1)
    mizeryear1.5@n_other<-head(mizeryear1@n_other,nf+1)
    fint1<-getTimes(mizeryear1.5)[idxFinalT(mizeryear1.5)]
    t_max <- t-nf/step
    if(t-nf/step==0){fint1<-fint1-1/step;t_max<-1/step}
    #append mizeryear1.5 to mizeryear2
    mizeryear2 <- project(mizeryear1.5,effort =0,t_start=fint1, t_max = t_max,t_save =1/step,dt=1/step/2,rm=FALSE)#tmax=1-(nf-(t-1)*step)/step
    #extract ei
    yield_step2<-getYield2(mizeryear2,step)
    yieldsum2<-apply(tail(yield_step2,step),2,cumsum)[step,]# yield of t year by species
    prop=nf/step-(t-1)
    Y0<-data.frame(time=t,Yield=yieldsum2,species=names(yieldsum2))
    Y<-melt(Y0,id=c("time","species"))
    SSB0<-tail(getSSB(mizeryear2),1)
    SSB1<-data.frame(time=t,SSB=SSB0[1,],species=colnames(SSB0))
    SSB<-melt(SSB1,id=c("time","species"))
    TY<-sum(yieldsum2)
    Ytar<-sum(yieldsum2[targetsp])
    Ynont<-sum(yieldsum2[nontargetsp])
    TSSB<-sum(SSB0)
    SSBtar<-sum(SSB0[,targetsp])
    SSBnont<-sum(SSB0[,nontargetsp])
    Eraw<-effortini*prop
    eiall0<-data.frame(time=t,effortini=effortini,prop=prop,Eraw=Eraw,chokesp=chokesp,
                       TY=TY,Ytar=Ytar,Ynont=Ynont,TSSB=TSSB,SSBtar=SSBtar,SSBnont=SSBnont)
    eiall<-melt(eiall0,id='time')
    eiall$species<-"all"
    ei<-rbind(eiall,Y,SSB)
    
    #set appropriate effort and initial mizer for next year
    if(chokesp=="nontargetsp"){effortini<-effortini*prop/sum(yieldsum2[nontargetsp])*tacnon*2}
    if(!chokesp=="nontargetsp"){effortini<-effortini*prop/yieldsum2[chokesp]*tac[chokesp]*2}
    effortini<-ifelse(prop==1,2*effortini,effortini)#all step lower than tac, then double the effortini
    effortini<-ifelse(effortini>20,20,effortini)#restrict unreal efforrt
    lastmizer<-mizeryear2
    tac<-0.5*(yieldsum2+0.2*as.vector(SSB0))*rlnorm(21,log(1)-0.2^2/2,0.2)
    return(ei)
  }
  eisstacrun$run<-run
  return(eisstacrun)
}
#mstac of target species mode --------
eimstac<-function(run){
  tac<-as.array(colMeans(tail(getYield(mizerbase),10))*0.8 )*rlnorm(21,log(1)-0.2^2/2,0.2)#initial tac for the first year
  step=50;year=50 #step is t_save/dt*2 in each year; year is length of simulation
  targetsp<-c("P.polyactis","O.ochellatus","S.niphonius","O.oratoria")
  nontargetsp<-spnames[-match(targetsp,spnames)]
  tactar<-sum(tac[targetsp])
  tacnon<-sum(colMeans(tail(getYield(mizerbase),10))[nontargetsp]) #total nontarget catch of base model
  eimstacrun<-foreach(t=1:year,.combine = 'rbind',.errorhandling = "remove") %do%{
    startt<-t-1
    if (t==1){lastmizer<-paramsbase;effortini=3}
    mizeryear1<-project(lastmizer,effort =effortini,t_start=startt,t_save=1/step,t_max = 1,dt=1/step/2,rm=FALSE)
    yield_step1 <-getYield2(mizeryear1,step)
    #identical(yield,getYield2(mizeryear1))
    yieldsum1<-apply(tail(yield_step1,step),2,cumsum)#cumsum from time t+1/step
    
    #decided the stop time in a year
    nnf1<-max(which(rowSums(yieldsum1[,targetsp])<tactar))+(t-1)*step 
    names(nnf1)<-"targetsp"
    nnf2<-max(which(rowSums(yieldsum1[,nontargetsp])<tacnon))+(t-1)*step 
    names(nnf2)<-"nontargetsp"
    nf<-min(c(nnf1,nnf2))
    chokesp<-names(which.min(c(nnf1,nnf2)))
    
    #change the mizeryear1 to mizeryear1.5 
    mizeryear1.5<-mizeryear1
    mizeryear1.5@params<-setInitialValues2(paramsbase,mizeryear1,no_t=nf+1)#change initial_(n,n_pp,n_other,effort)
    mizeryear1.5@effort<-head(mizeryear1@effort,nf+1)
    mizeryear1.5@n<-head(mizeryear1@n,nf+1)
    mizeryear1.5@n_pp<-head(mizeryear1@n_pp,nf+1)
    mizeryear1.5@n_other<-head(mizeryear1@n_other,nf+1)
    fint1<-getTimes(mizeryear1.5)[idxFinalT(mizeryear1.5)]
    t_max <- t-nf/step
    if(t-nf/step==0){fint1<-fint1-1/step;t_max<-1/step}
    #append mizeryear1.5 to mizeryear2
    mizeryear2 <- project(mizeryear1.5,effort =0,t_start=fint1, t_max = t_max,t_save =1/step,dt=1/step/2,rm=FALSE)#tmax=1-(nf-(t-1)*step)/step
    #extract ei
    yield_step2<-getYield2(mizeryear2,step)
    yieldsum2<-apply(tail(yield_step2,step),2,cumsum)[step,]# yield of t year by species
    prop=nf/step-(t-1)
    Y0<-data.frame(time=t,Yield=yieldsum2,species=names(yieldsum2))
    Y<-melt(Y0,id=c("time","species"))
    SSB0<-tail(getSSB(mizeryear2),1)
    SSB1<-data.frame(time=t,SSB=SSB0[1,],species=colnames(SSB0))
    SSB<-melt(SSB1,id=c("time","species"))
    TY<-sum(yieldsum2)
    Ytar<-sum(yieldsum2[targetsp])
    Ynont<-sum(yieldsum2[nontargetsp])
    TSSB<-sum(SSB0)
    SSBtar<-sum(SSB0[,targetsp])
    SSBnont<-sum(SSB0[,nontargetsp])
    Eraw<-effortini*prop
    eiall0<-data.frame(time=t,effortini=effortini,prop=prop,Eraw=Eraw,chokesp=as.character(chokesp),
                       TY=TY,Ytar=Ytar,Ynont=Ynont,TSSB=TSSB,SSBtar=SSBtar,SSBnont=SSBnont)
    eiall<-melt(eiall0,id='time')
    eiall$species<-"all"
    ei<-rbind(eiall,Y,SSB)
    
    #set appropriate effort and initial mizer for next year
    if(chokesp=="targetsp"){effortini<-effortini*prop/sum(yieldsum2[targetsp])*tactar*2}
    if(chokesp=="nontargetsp"){effortini<-effortini*prop/sum(yieldsum2[nontargetsp])*tacnon*2}
    effortini<-ifelse(prop==1,2*effortini,effortini)#all step lower than tac, then double the effortini
    effortini<-ifelse(effortini>20,20,effortini)#restrict unreal efforrt
    lastmizer<-mizeryear2
    return(ei)
  }
  eimstacrun$run<-run
  return(eimstacrun)
}

#mstacall of all species mode --------
eimstacall<-function(run){
  tac<-as.array(colMeans(tail(getYield(mizerbase),10))*0.8 )*rlnorm(21,log(1)-0.2^2/2,0.2)#initial tac for the first year
  step=50;year=50 #step is t_save/dt*2 in each year; year is length of simulation
  targetsp<-c("P.polyactis","O.ochellatus","S.niphonius","O.oratoria")
  nontargetsp<-spnames[-match(targetsp,spnames)]
  tacall <- sum(tac)
  # tactar<-sum(tac[targetsp])
  # tacnon<-sum(colMeans(tail(getYield(mizerbase),10))[nontargetsp]) #total nontarget catch of base model
  eimstacallrun<-foreach(t=1:year,.combine = 'rbind',.errorhandling = "remove") %do%{
    startt<-t-1
    if (t==1){lastmizer<-paramsbase;effortini=3}
    mizeryear1<-project(lastmizer,effort =effortini,t_start=startt,t_save=1/step,t_max = 1,dt=1/step/2,rm=FALSE)
    yield_step1 <-getYield2(mizeryear1,step)
    #identical(yield,getYield2(mizeryear1))
    yieldsum1<-apply(tail(yield_step1,step),2,cumsum)#cumsum from time t+1/step
    
    #decided the stop time in a year
    nf <- max(which(rowSums(yieldsum1)<tacall))+(t-1)*step 
    # nnf1<-max(which(rowSums(yieldsum1[,targetsp])<tactar))+(t-1)*step 
    # names(nnf1)<-"targetsp"
    # nnf2<-max(which(rowSums(yieldsum1[,nontargetsp])<tacnon))+(t-1)*step 
    # names(nnf2)<-"nontargetsp"
    # nf<-min(c(nnf1,nnf2))
    chokesp<-"allsp"
    #chokesp<-names(which.min(c(nnf1,nnf2)))
    
    #change the mizeryear1 to mizeryear1.5 
    mizeryear1.5<-mizeryear1
    mizeryear1.5@params<-setInitialValues2(paramsbase,mizeryear1,no_t=nf+1)#change initial_(n,n_pp,n_other,effort)
    mizeryear1.5@effort<-head(mizeryear1@effort,nf+1)
    mizeryear1.5@n<-head(mizeryear1@n,nf+1)
    mizeryear1.5@n_pp<-head(mizeryear1@n_pp,nf+1)
    mizeryear1.5@n_other<-head(mizeryear1@n_other,nf+1)
    fint1<-getTimes(mizeryear1.5)[idxFinalT(mizeryear1.5)]
    t_max <- t-nf/step
    if(t-nf/step==0){fint1<-fint1-1/step;t_max<-1/step}
    #append mizeryear1.5 to mizeryear2
    mizeryear2 <- project(mizeryear1.5,effort =0,t_start=fint1, t_max = t_max,t_save =1/step,dt=1/step/2,rm=FALSE)#tmax=1-(nf-(t-1)*step)/step
    #extract ei
    yield_step2<-getYield2(mizeryear2,step)
    yieldsum2<-apply(tail(yield_step2,step),2,cumsum)[step,]# yield of t year by species
    prop=nf/step-(t-1)
    Y0<-data.frame(time=t,Yield=yieldsum2,species=names(yieldsum2))
    Y<-melt(Y0,id=c("time","species"))
    SSB0<-tail(getSSB(mizeryear2),1)
    SSB1<-data.frame(time=t,SSB=SSB0[1,],species=colnames(SSB0))
    SSB<-melt(SSB1,id=c("time","species"))
    TY<-sum(yieldsum2)
    Ytar<-sum(yieldsum2[targetsp])
    Ynont<-sum(yieldsum2[nontargetsp])
    TSSB<-sum(SSB0)
    SSBtar<-sum(SSB0[,targetsp])
    SSBnont<-sum(SSB0[,nontargetsp])
    Eraw<-effortini*prop
    eiall0<-data.frame(time=t,effortini=effortini,prop=prop,Eraw=Eraw,chokesp=as.character(chokesp),
                       TY=TY,Ytar=Ytar,Ynont=Ynont,TSSB=TSSB,SSBtar=SSBtar,SSBnont=SSBnont)
    eiall<-melt(eiall0,id='time')
    eiall$species<-"all"
    ei<-rbind(eiall,Y,SSB)
    
    #set appropriate effort and initial mizer for next year
    # if(chokesp=="targetsp"){effortini<-effortini*prop/sum(yieldsum2[targetsp])*tactar*2}
    # if(chokesp=="nontargetsp"){effortini<-effortini*prop/sum(yieldsum2[nontargetsp])*tacnon*2}
    effortini<-effortini*prop/sum(yieldsum2)*tacall*2
    effortini<-ifelse(prop==1,2*effortini,effortini)#all step lower than tac, then double the effortini
    effortini<-ifelse(effortini>20,20,effortini)#restrict unreal efforrt
    lastmizer<-mizeryear2
    return(ei)
  }
  eimstacallrun$run<-run
  return(eimstacallrun)
}
#parallel snowfall---------------------


#sfInit(parallel=TRUE,type='SOCK',socketHosts=  rep(c("node01", "node02", "node03","node04"),each=25)) ##slower in creating cluster and faster in computing
sfInit(parallel=TRUE, cpus=7)
sfExport('mizerbase','paramsbase','spnames') 
sfExport('eisstac','eimstac','getYield2','setInitialValues2','BevertonHoltRDD') 
sfLibrary(mizer)
sfLibrary(mizerHowTo)
sfLibrary(mizerExperimental)
sfLibrary(reshape2)
sfLibrary(assertthat)
sfLibrary(foreach)

#### call function
system.time({res_eisstac= sfLapply(1:1000,eisstac)})
system.time({res_eimstac= sfLapply(1:1000,eimstac)})
system.time({res_eimstacall= sfLapply(1:1000,eimstacall)})#takds 10 hours

res_eisstac<-do.call(rbind,res_eisstac)
res_eimstac<-do.call(rbind,res_eimstac)
res_eimstacall<-do.call(rbind,res_eimstacall)

save(res_eisstac,file = "res_eisstac.Rdata")
save(res_eimstac,file = "res_eimstac.Rdata")
save(res_eimstacall,file = "res_eimstacall.Rdata") 
sfStop()

#plot---------------------------
load("res_eisstac.Rdata")
load("res_eimstac.Rdata")
load("res_eimstacall.Rdata")
res_eisstac$time <- round(res_eisstac$time)
res_eimstac$time <- round(res_eimstac$time)
res_eimstacall$time <- round(res_eimstacall$time)
# table(subset(res_eisstac,variable=="chokesp")$value)
# table(subset(res_eimstac,variable=="chokesp")$value)
targetsp<-c("P.polyactis","O.ochellatus","S.niphonius","O.oratoria")
SSBini<-colMeans(tail(getSSB(mizerbase),10));SSBini<-data.frame(species=names(SSBini),value=SSBini)
dataset <- bind_rows(mutate(res_eisstac,mp='ss'),mutate(res_eimstac,mp='ms'),mutate(res_eimstacall,mp='msall')) %>% 
  mutate(mp=factor(mp,levels=c('ss','ms','msall'))) 

create_p1 <- function(data, variable_label, y_label) {
  data %>% 
    mutate(value=as.numeric(value)) %>% 
    mutate(value = ifelse(variable %in% c("Ytar","TSSB"),value/10^9,value)) %>% 
    ggplot(aes(x = time, y = value, group = mp)) +
    stat_flquantiles(aes(fill = mp), probs = c(0.25, 0.75), geom = "ribbon", alpha = 0.20) +
    geom_smooth(aes(group = mp), color = "black", se = FALSE, linewidth = 0.4) +
    labs(x = "", y = y_label,title =variable_label) +
    scale_color_manual(aesthetics = c("colour", "fill"), values = c("ss" = "#377eb8", "ms" = "#e41a1c", "msall" = "grey40")) +
    theme_bw() +
    theme(legend.position = "none")
}

p1_1 <- dataset %>% filter(variable == "Eraw") %>% create_p1("(a) Relative fishing effort", "E")
p1_2 <- dataset %>% filter(variable == "Ytar") %>% create_p1("(b) Total fishery yield", "TY [kt]")
p1_3 <- dataset %>% filter(variable == "TSSB") %>% create_p1("(c) Total biomass", "B [kt]")

p1 <- plot_grid(p1_1, p1_2, p1_3, nrow = 1, align = "h")+ draw_label("Years", x = 0.5, y = 0, vjust = -0.5,size=12)
#ggsave("p1.png",plot = p1, height=4, width=8, units="in")

chokecolors <- c("O.ochellatus" = "#fdb462", "O.oratoria" = "#80b1d3", 
                 "S.niphonius" = "#bebada", "P.polyactis" = "#fb8072", 
                 "nontargetsp" = "#ffffb3", "targetsp" = "#8dd3c7")
chokenames <- c("O.ochellatus" = "Octopus", "O.oratoria" = "Mantis shrimp", 
                "S.niphonius" = "Spanish mackerel", "P.polyactis" = "Yellow croaker", 
                "nontargetsp" = "nontarget", "targetsp" = "target")

chokedata <- dataset %>% filter(!mp=="msall") %>% 
  subset(variable=="chokesp") %>% 
  dplyr::count(value,mp) %>%
  mutate(mp = factor(mp,levels=c("ss","ms"))) %>% 
  group_by(mp) %>% 
  mutate(percentage = n / sum(n) * 100) %>% 
  mutate(commonsp=chokenames[value]) %>% 
  as.data.frame()

p2 <- ggplot(chokedata, aes(x = "", y = percentage, fill = value)) +
  geom_bar(stat = "identity", width = 1,color="white") +
  coord_polar(theta = "y") +
  geom_label_repel(aes(label = paste0(commonsp,":",percentage, "%")), 
                   nudge_y = 0.5) +
  scale_fill_manual(values = chokecolors) +
  facet_wrap(vars(mp), labeller = labeller(mp = c(ss = "Hybrid SSTAC MP", ms = "Hybrid MSTAC MP"))) + 
  theme_void()+
  labs(fill = "") +
  theme(legend.position = "none")
p2
#ggsave("p2.png", height=3, width=5, units="in")
targetsp<-c("P.polyactis","O.ochellatus","S.niphonius","O.oratoria")
nontargetsp <- spnames[-match(targetsp,spnames)]
classcolors <-c('target'="#8dd3c7",'nontarget'="yellow3")
p3 <- subset(dataset,variable=="SSB") %>% 
  mutate(class=ifelse(species %in% targetsp,'target','nontarget')) %>% 
  mutate(class=factor(class,levels=c('nontarget','target'))) %>% 
  mutate(species=factor(species,levels=c(nontargetsp,targetsp))) %>% 
  mutate(mp=factor(mp,levels=c("ss",'ms','msall'))) %>% 
  mutate(value=as.numeric(value)/10^9) %>% 
  group_by(time,variable, species, mp, class)  %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot()+geom_col(aes(x=time,y=mean_value,group=species,fill=class),
                    color='white', position = position_stack())+
  labs(x='Years',y='Biomass [kt]')+
  scale_fill_manual(values=classcolors)+
  theme_bw()+  guides(fill = guide_legend(title = NULL))+
  facet_wrap(vars(mp),labeller =labeller(mp = c(ss = "Hybrid SSTAC MP", ms = "Hybrid MSTAC MP", msall="Reference MP"))) 
p3
#ggsave("p3.png", height=4, width=8, units="in")

p4 <- subset(dataset,variable=="SSB") %>% 
  mutate(class=ifelse(species %in% targetsp,'target','nontarget')) %>% 
  dplyr::filter(class=='nontarget') %>% 
  mutate(mp=factor(mp,levels=c("ss",'ms','msall'))) %>% 
  mutate(value=as.numeric(value)/10^9) %>% 
  group_by(time,variable, species, mp, class)  %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  left_join(scitocomnames) %>% 
  ggplot()+geom_col(aes(x=time,y=mean_value,fill=reorder(commonnames,mean_value,decreasing=TRUE)),
                    color='white', position = "fill")+
  labs(x='Years',y='Biomass percentage',fill='Species')+
  scale_fill_viridis(discrete=TRUE,direction=-1,begin=0,end=0.85)+
  theme_bw()+
  facet_wrap(vars(mp),labeller = labeller(mp=c(ss = "Hybrid SSTAC MP", ms = "Hybrid MSTAC MP", msall="Reference MP"))) 
p4
#ggsave("p4.png", height=4.1, width=9, units="in")

withpropSSB<-function(x){
  joinSSB<-left_join(subset(x,variable=="SSB"),SSBini,by="species",suffix=c(".SSB",".SSBini"))
  joinSSB$propSSB<-apply(joinSSB,1,function(x) as.numeric(x["value.SSB"])/as.numeric(x["value.SSBini"]))
  propSSB<-melt(joinSSB[,c("time","species","run","propSSB")],id=c("time","species","run"))
  x0<-rbind(x,propSSB)
  x0<-subset(x0,!variable=="chokesp")
  x0$value<-as.numeric(x0$value)
  return(x0)
}
res_eisstac0<-withpropSSB(res_eisstac)
res_eimstac0<-withpropSSB(res_eimstac)
res_eimstacall0<-withpropSSB(res_eimstacall)
# risk with time
lbfun <- function(x,lim){
  n<-dcast(data=subset(x,variable=='propSSB'&value>=0), time~.,length)
  m<-dcast(data=subset(x,variable=='propSSB'&value<lim),time~.,length,fill=0,drop=FALSE)
  mn<-merge(m,n,by='time',suffixes = c("m","n"),all.y = TRUE)
  mn[is.na(mn$.m),'.m']<-0
  cbind(time=mn[,1],value=mn$.m/mn$.n)
}
# risk with time * species
lbfunsp <- function(x,lim){
  n<-dcast(data=subset(x,variable=='propSSB'&value>=0), time+species~.,length)
  m<-dcast(data=subset(x,variable=='propSSB'&value<lim),time+species~.,length,fill=0,drop=FALSE)
  mn<-merge(m,n,by=c('time','species'),suffixes = c("m","n"),all.y = TRUE)
  mn[is.na(mn$.m),'.m']<-0
  cbind(mn[,1:2],value=mn$.m/mn$.n)
}
# risk with time * run
lbfunrun <- function(x,lim){
  n<-dcast(data=subset(x,variable=='propSSB'&value>=0), time+run~.,length)
  m<-dcast(data=subset(x,variable=='propSSB'&value<lim),time+run~.,length,fill=0,drop=FALSE)
  mn<-merge(m,n,by=c('time','run'),suffixes = c("m","n"),all.y = TRUE)
  mn[is.na(mn$.m),'.m']<-0
  cbind(mn[,1:2],value=mn$.m/mn$.n)
}

ssrisk<-lbfun(res_eisstac0,0.2)
ssrisksp<-lbfunsp(res_eisstac0,0.2)
ssriskrun<-lbfunrun(res_eisstac0,0.2)
msrisk<-lbfun(res_eimstac0,0.2)
msrisksp<-lbfunsp(res_eimstac0,0.2)
msriskrun<-lbfunrun(res_eimstac0,0.2)
msallrisk<-lbfun(res_eimstacall0,0.2)
msallrisksp<-lbfunsp(res_eimstacall0,0.2)
msallriskrun<-lbfunrun(res_eimstacall0,0.2)
library(plyr); library(dplyr)
sssumrisk<- summarySE(ssriskrun, measurevar = "value", groupvars="time",na.rm = FALSE)
mssumrisk<- summarySE(msriskrun, measurevar = "value", groupvars="time",na.rm = FALSE)
msallsumrisk<- summarySE(msallriskrun, measurevar = "value", groupvars="time",na.rm = FALSE)
Bstatus<-rbind(data.frame(ssrisksp,type="ss"),data.frame(msrisksp,type="ms"),data.frame(msallrisksp,type="msall"))
res_ei<-rbind(data.frame(res_eisstac0,type="ss"),data.frame(res_eimstac0,type="ms"),data.frame(res_eimstacall0,type="msall"))
Bstatus_median<-melt(acast(subset(res_ei,variable=="propSSB"),time~species~type,median))
Bstatus_median[which(Bstatus_median$value<0.01),'value']<-0.01
colnames(Bstatus_median)<-c("time",'species','type','value')
Bstatus$target<-"2";Bstatus$target[which(Bstatus$species %in% targetsp)]<-"1"
Bstatus_median$target<-"2";Bstatus_median$target[which(Bstatus_median$species %in% targetsp)]<-"1"
table(Bstatus_median)
p5<-Bstatus %>% 
    mutate(type=factor(type,levels=c("ss","ms","msall"))) %>% 
    left_join(scitocomnames) %>% 
    ggplot() +
    geom_tile(aes(x = time, y = commonnames,fill = value),colour="black") +
    scale_fill_gradientn(colors=brewer.pal(9,'YlGnBu')[1:6], breaks=c(0,0.25,0.5,0.75,1),
                         labels=c("0%","25%","50%","75%","100%"),
                         guide = guide_colorbar(title = "Pdep(B/Bini<0.2)",
                                                title.vjust = 1,
                                                direction = "horizontal",
                                                title.theme=element_text(size=9),
                                                barheight=0.7,barwidth=25,
                                                frame.colour = "grey",frame.linewidth = 0.5,
                                                draw.ulim = TRUE, draw.llim = TRUE ))+
    #  geom_text(data=Bstatus_median,size=3,colour="white",aes(x = time, y = species,fontface=2,label=ifelse(value<0.2,round(value,2),'')))+
    facet_grid(target~type,scales = "free_y",space  = "free_y",
               labeller = labeller(target=c("1"='target species',"2"='nontarget species')
                                   ,type=c(ss = "Hybrid SSTAC MP", ms = "Hybrid MSTAC MP", msall="Reference MP"))) +
    theme_bw()+labs(x='Years',y='')+
    theme(panel.grid =element_blank(),panel.spacing.y =unit(-0.14, "lines"),
          legend.position = "bottom",legend.margin=margin(0,10,0,10),
          panel.border=element_blank(),legend.box.margin=margin(-10,0,0,0),
          axis.text.x = element_text(angle=30,hjust=1,colour='black'),
          axis.text.y = element_text(angle=0,hjust=1,colour='black'))
  p5
#ggsave("p5.png", height=8, width=10, units="in")

risk<-rbind(data.frame(ssrisk,type="ss"),data.frame(msrisk,type="ms"),data.frame(msallrisk,type="msall")) 
p6<-risk %>% 
    mutate(type=factor(type,levels=c("ss","ms","msall"))) %>% 
    ggplot(aes(x = time, y = value,color = type,fill=type))+
    geom_point(size=0.5)+geom_smooth(size=0.5,alpha=0.2)+
    #   geom_jitter(data=riskrun,aes(x = time,colour=type, y = value),size=0.04)+
    theme_bw()+ylim(0,1)+
    scale_color_manual(aesthetics = c("colour","fill"),values =  c("ss"="#377eb8", "ms"="#e41a1c","msall"="grey40")) +
    labs(y="Depletion Risk (Pdep, threshold=0.2)",x="Years")+
    theme(axis.text = element_text(size=9),
          axis.title.y = element_text(size=9),
          legend.position = "none")
p6
#ggsave("p6.png", height=3, width=4, units="in")


# SA plot-----------------------------------------------------
# plot study area map---

load("Shandong_basemap.Rdata")
arrow <- data.frame(x1 = 11.5, x2 = 13, y1 = 10, y2 = 10)
ps1 <- ggdraw(xlim = c(0, 20), ylim = c(0, 20)) +
  draw_plot(shandong, x = 0, y = 0, width = 12, height = 20) +
  draw_plot(site, x = 12.8, y = 3, width = 7, height =15) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y1+2), data = arrow, lineend = "round") +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y1-2), data = arrow, lineend = "round") 

ps1
#ggsave("ps1.png", plot = ps1, width = 9, height = 7)

#my plot function---
plot<-function(simplot, 
               ...) {
  library(ggpubr)
  library(ggplot2)
  source('plotbiomass.R')
  source('plotfeedinglevel.R')
  source('plotfmort.R')
  source('plotm2.R')
  source('plotspectra.R')
  df <- data.frame(x = 0:9, y = 0:9)
  base <- ggplot(df, aes(x, y)) +
    geom_blank() +
    theme_bw() + 
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p1<-plotfeedinglevel(simplot,print_it = F, ...)+theme_bw()+theme(panel.grid.major = element_blank(),
                                                                   panel.grid.minor = element_blank(),
                                                                   legend.position="none")
  p2<-plotbiomass(simplot,print_it = F, ...)+theme_bw()+theme(panel.grid.major = element_blank(),
                                                              panel.grid.minor = element_blank(),
                                                              legend.position="none")
  p3<-plotm2(simplot,print_it = F, ...)+theme_bw()+theme(panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank(),
                                                         legend.position="none")
  p4<-plotfmort(simplot,print_it = F, ...)+theme_bw()+theme(panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank(),
                                                            legend.position="none")
  p5<-plotspectra(simplot,ylim=c(NA,1e-12),print_it = F, ...)+theme_bw()+theme(panel.grid.major = element_blank(),
                                                                               panel.grid.minor = element_blank(),
                                                                               legend.position="none")
  p6<-as_ggplot(get_legend(plotspectra(simplot,ylim=c(NA,1e-12),total = F, plankton = T, 
                                       background = F, print_it = F, ...)+theme_bw()+
                             theme(legend.key.size=unit(1.5,'mm'),legend.key.width=unit(8.5,'mm'),
                                   legend.text =element_text(size=9),legend.key.height=unit(5,'mm'),
                                   legend.title =element_text(size=9))))
  base +
    annotation_custom(grob = ggplotGrob(p1), xmin = 0, xmax =4.5,ymin =6, ymax =9) +
    annotation_custom(grob = ggplotGrob(p2), xmin = 4.5, xmax =9,ymin =6, ymax =9) +
    annotation_custom(grob = ggplotGrob(p3), xmin = 0, xmax =4.5,ymin =3, ymax =6) +
    annotation_custom(grob = ggplotGrob(p4), xmin = 4.5, xmax =9,ymin =3, ymax =6) +
    annotation_custom(grob = ggplotGrob(p5), xmin = 0, xmax =4.5, ymin =0, ymax =3) +
    annotation_custom(grob = ggplotGrob(p6), xmin = 4.8, xmax =9, ymin =0.3, ymax =2.5)+ 
    geom_text(x=0.8,y=9,label="(a)",size=3.5)+geom_text(x=5.3,y=9,label="(b)",size=3.5)+
    geom_text(x=0.8,y=6,label="(c)",size=3.5)+geom_text(x=5.3,y=6,label="(d)",size=3.5)+
    geom_text(x=0.8,y=3,label="(e)",size=3.5)
}

# plot base model---

ps2 <- plot(mizerbase)
ps2
#ggsave("ps2.png", plot = ps2, height=10, width=7.5, units="in")

# relative yield--

iniyield<-colMeans(tail(getYield(mizerbase), 10)) 
iniyieldtable <-   tibble( species = names(iniyield),iniyield = iniyield)
propyield<- dataset %>% 
  filter(variable == "Yield") %>% 
  filter(time == 50) %>% 
  mutate(class = ifelse(species %in% targetsp, 'target', 'nontarget')) %>% 
  mutate(class=factor(class,levels=c('target','nontarget'))) %>% 
  left_join(iniyieldtable, by = join_by(species)) %>% 
  left_join(scitocomnames, by = join_by(species)) %>% 
  mutate(propyield = as.numeric(value) / iniyield) %>% 
  group_by(commonnames,class,species, mp) %>% 
  dplyr::summarise(mean_propyield = mean(propyield, na.rm = TRUE)) %>% 
  ungroup()

mean_propyield_table <-ggplot(propyield, aes(x = mp, y = commonnames,color=mp)) +
  geom_tile(color = "grey60", linewidth = 0.5, na.rm = FALSE, fill = "white") +
  geom_text(aes(label = sprintf("%.4f", mean_propyield)), vjust = 0.5) + 
  theme_minimal(base_size = 12) +
  scale_x_discrete(expand = c(0, 0),
                   labels = c("ss" = "Hybrid SSTAC MP", "ms" = "Hybrid MSTAC MP", "msall" = "Reference MP")) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_color_manual(aesthetics = c("colour", "fill"), 
                     values = c("ss" = "#377eb8", "ms" = "#e41a1c", "msall" = "grey40"))+
  theme_bw() + facet_grid(class~.,scales = "free_y",space  = "free_y")+
  labs(y='Species',x=" ")+theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.spacing = unit(0, "lines"),
        legend.position = "none") 

tac <- colMeans(tail(getYield(mizerbase), 10)) * 0.8 
tactable <- tibble( species = names(tac),tac = tac)
propyield_plot <- dataset %>% 
  filter(variable == "Yield") %>% 
  filter(time == 50) %>% 
  mutate(class = ifelse(species %in% targetsp, 'target', 'nontarget')) %>% 
  mutate(class = factor(class, levels = c('target', 'nontarget'))) %>% 
  left_join(tactable, by = join_by(species)) %>% 
  left_join(scitocomnames, by = join_by(species)) %>% 
  mutate(propyield = as.numeric(value) / tac) %>% 
  ggplot() +
  stat_summary(aes(x = propyield, y = commonnames, colour = mp), fun = "mean",
               geom = "point", size = 1.1, orientation = "y",
               position = position_dodge(0.60)) +
  stat_summary(aes(x = propyield, y = commonnames, colour = mp),
               geom = "linerange", size = 0.6, orientation = "y",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               position = position_dodge(0.60)) +
  geom_vline(aes(xintercept = 1), colour = 'grey30', linetype = 'dotted') +
  scale_color_manual(aesthetics = c("colour", "fill"), 
                     values = c("ss" = "#377eb8", "ms" = "#e41a1c", "msall" = "grey40"),
                     labels = c("ss" = "Hybrid SSTAC MP", "ms" = "Hybrid MSTAC MP", "msall" = "Reference MP")) +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  labs(y = 'Species', x = "Relative Yield") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.spacing = unit(0, "lines"),
        legend.position = "bottom") +
  guides(colour = guide_legend(title = " "))

ps4 <- plot_grid(
  propyield_plot + theme(plot.margin = margin(5, 5, 5, 5)),
  mean_propyield_table + theme(plot.margin = margin(5, 5, 5, 5)),
  labels = c("(a)", "(b)"),
  nrow = 1, rel_widths = c(1, 1) ,align = "hv",
  axis = "tb")
ps4
#ggsave("ps4.png", ps4, width = 11, height = 8)


# end ---------------------------------------------------------------------



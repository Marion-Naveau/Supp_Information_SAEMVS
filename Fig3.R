## >>>>>>>>>>>>>>>>> Script for Figure 3 <<<<<<<<<<<<<<<<<<<<<<

rm(list=ls())

library(microbenchmark)
library(ggplot2)
library(ggpubr)

## We recover the calculation times from MCMC-SAEM

file_name_MCMCSAEM=c("Saves/Time_MCMCSAEM_p500.Rdata","Saves/Time_MCMCSAEM_p700.Rdata","Saves/Time_MCMCSAEM_p1000.Rdata",
            "Saves/Time_MCMCSAEM_p1500.Rdata","Saves/Time_MCMCSAEM_p2000.Rdata","Saves/Time_MCMCSAEM_p2500.Rdata")

S=50
time_MCMCSAEM=rep(0,6)
for (k in c(1:length(file_name_MCMCSAEM))){
  load(file_name_MCMCSAEM[k])
  time=rep(0,S)
  for (s in 1:S){
    time[s]=summary(Total_Time[[s]])$min
  }
  time
  time_MCMCSAEM[k]=min(time) #in seconds
}

## We recover the calculation times from MCMC

file_name_MCMC=c("Saves/Time_MCMC_p500.Rdata","Saves/Time_MCMC_p700.Rdata","Saves/Time_MCMC_p1000.Rdata",
                  "Saves/Time_MCMC_p1500.Rdata","Saves/Time_MCMC_p2000.Rdata","Saves/Time_MCMC_p2500.Rdata")

S=50
time_MCMC=rep(0,6)
for (k in c(1:length(file_name_MCMC))){
  load(file_name_MCMC[k])
  time=rep(0,S)
  for (s in 1:S){
    time[s]=Total_Time[[s]]
  }
  time
  time_MCMC[k]=min(time) #in seconds
}

## Script for Figure 3:

data_time=data.frame(Method=rep(c("MCMC-SAEM","MCMC"),each=6),p=rep(c(500,700,1000,1500,2000,2500),2),Time=c(time_MCMCSAEM,time_MCMC))

head(data_time)
data_time$Method<-as.factor(data_time$Method)
Figure3=ggplot(data_time,aes(x=p,y=Time,color=Method))+geom_point()+labs(y="Time (in seconds)")+ theme_minimal()+scale_x_continuous(trans = 'log10')+scale_y_continuous(trans = 'log10')+
  scale_color_manual(values = c("#DDA0DD", "#00a3a6"))+
  geom_smooth(method="lm",se=FALSE)+stat_regline_equation() + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"),title=element_text(size=20,face="bold"),legend.title = element_text(size=20,face="bold"),legend.text = element_text(size=20))
Figure3

ggsave("Figure3.pdf")

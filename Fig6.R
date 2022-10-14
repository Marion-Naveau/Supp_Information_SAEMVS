## >>>>>>>>>>>> Script for Figure 6 <<<<<<<<<<<

rm(list=ls())

library(ggplot2)
library(stringr)
library(dplyr)

####################### Observed data #######################

data <- read.csv2("13INRmon_NUE_LN.csv")
head(data)
dates <- which(str_detect(colnames(data),"date") == T)

#Formatting the data set:
newdata <- data[,1:(dates[1])]
nbcol <- dim(newdata)[2]
colnames(newdata)[nbcol-1] <- "Y"
colnames(newdata)[nbcol] <- "Date"

for (i in dates[-1]){
  datatemp <- data[,c(1:(dates[1]-2),i-1,i)]
  colnames(datatemp)[nbcol-1] <- "Y"
  colnames(datatemp)[nbcol] <- "Date"
  newdata <- rbind(newdata,datatemp)
}


newdata <- newdata[,-c(2,3,6,7)] #One year, one location, one constraint and one condition: unnecessary variables

newdata$Date <- as.Date(newdata$Date,"%d/%m/%y")

newdata$Line <- paste(newdata$GENOTYPE,newdata$REPLICATE,sep="_")

# Only block 1 is kept when there are several blocks:
for (i in unique(newdata$GENOTYPE)){
  if (length(unique(newdata[newdata$GENOTYPE==i,]$BLOCK))!=1){
    ind_block=which((newdata$GENOTYPE==i)&(newdata$BLOCK!=1))
    newdata=newdata[-ind_block,]
  }
}

#We keep only replicate 1:
newdata1 <- filter(newdata,REPLICATE=="1")

minD <- min(julian(newdata1$Date))
newdata1$Day <- julian(newdata1$Date) - minD

data_obs=newdata1
save(data_obs,file="data_obs.Rdata")

par(mfrow=c(3,3))
for (i in unique(data_obs$GENOTYPE)){
  ind <- which((data_obs$GENOTYPE==i))
  plot(y=data_obs[ind,]$Y,x=data_obs[ind,]$Date,col=1,type='l',main=substitute(exp = paste("Genotype : ", gen_i),env = list(gen_i = i)))
}
genotype_plot=unique(data_obs$GENOTYPE)[200:209]
data_plot=data_obs[which(data_obs$GENOTYPE%in%genotype_plot),]
ggplot(data_plot, aes(x=Date, y=Y, group=GENOTYPE,color=GENOTYPE)) + geom_point(size=2) +
  geom_line(linetype = "dashed") + labs(color="VARIETIES") + ylab("Senescence (in %)") +theme_bw() +
  theme(axis.text=element_text(size=30), axis.title=element_text(size=30,face="bold"),legend.title = element_text(size=15), legend.text = element_text(size=15))



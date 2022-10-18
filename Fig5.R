## >>>>>>>>>>>> Script for Figure 5 <<<<<<<<<<<

rm(list=ls())

library(ggplot2)
library(cowplot)
library(stringr)
library(dplyr)
library(mvnfast)
library(doParallel)

####################### Observed data #######################

data <- read.csv2("Data/13INRmon_NUE_LN.csv")
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

###################### Data-set with the SNPs ######################

SNPs <- read.csv2("Data/genotypeFilterPlink_100_5_0.8.csv",sep="",row.names = 1)
SNPs <- as.matrix(SNPs)
SNPs <- t(SNPs)
SNPs <- as.data.frame(SNPs)
SNPs <- SNPs%>%mutate_at(colnames(SNPs), as.character)
SNPs <- SNPs%>%mutate_at(colnames(SNPs), as.numeric)
SNPs <- as.matrix(SNPs)

data_chromosome <- read.csv2("Data/carte_Axiom-TABW420k_WGAv1.csv",header=TRUE, sep = " ") #dataframe containing the position on the genome of each SNP
#data_chromosome$V1: chromosome, data_chromosome$V2: position on the chromosome, data_chromosome$name: name of the SNP

# Missing values were imputed as the marker observed frequency by the biologists,
# and then we replaced these imputed values by 0 or 1 using a threshold of 0.5:
for (i in 1:(dim(SNPs)[1])){
  for (j in 1:(dim(SNPs)[2]))
    if (SNPs[i,j]<0.5){
      SNPs[i,j]=0
    }
  else {
    SNPs[i,j]=1
  }
}

#Monomorphic and unmapped markers were removed from the data-set:
variability=rep(0,dim(SNPs)[2])
for (j in 1:dim(SNPs)[2]){
  variability[j]=length(unique(SNPs[,j]))
}
remove=which(variability==1)
SNPs=SNPs[,-remove]
dim(SNPs)

sum(colnames(SNPs)%in%data_chromosome$name) #Only 26000 SNPs are in data_chromosome
data_chromosome=data_chromosome[which(data_chromosome$name%in%colnames(SNPs)),] #only SNPs are kept in this data

#chrUn is a "dummy" chromosome used to mean that the position of the SNP in question on the genome is not known, they are deleted:
data_chromosome=data_chromosome[which(data_chromosome$V1!="chrUn"),]
dim(data_chromosome)[1] #Finally, we have 26189 SNPs in our data-set


heading_QTL <- read.csv("Data/marker_HD_INRmon13LN.csv",sep="")[,1]  #list of heading QTLs names
heading_QTL <- as.character(heading_QTL)


###################### Example for Chromosome 6A ######################

#### Firt we construct the matrix V_tilde associated with the chromosome 6A:

markers_chr6A=data_chromosome$name[(data_chromosome$V1=="chr6A")]
V_chr6A=SNPs[,colnames(SNPs)%in%markers_chr6A]
n=dim(V_chr6A)[1]
markers_noQTL=setdiff(colnames(SNPs),heading_QTL) #list of SNPs names that are not heading QTLs

# To limit collinearity between covariates, if several SNPs in a same chromosome are collinear for all varieties,
# only one of them is retained for analysis:
## Fisrt for the heading QTLs:
V_chr6A_QTL=as.matrix(V_chr6A[,colnames(V_chr6A)%in%heading_QTL])
nb_QTL=dim(V_chr6A_QTL)[2] #number of heading QTLs in chr6A
if (nb_QTL!=1){
  couples1 <- c()
  couples2 <- c()
  similarities <- matrix(NA,nb_QTL,nb_QTL)
  similarities2 <- matrix(NA,nb_QTL,nb_QTL)
  for (i in 1:(nb_QTL-1)){
    for (j in (i+1):nb_QTL){
      similarities[i,j] <- sum(V_chr6A_QTL[,i]==(V_chr6A_QTL[,j]))
      similarities2[i,j] <- sum(V_chr6A_QTL[,i]==(1-V_chr6A_QTL[,j]))
      if (similarities[i,j] == n){
        couples1 <- c(couples1,i,j)
      }
      if (similarities2[i,j] == n){
        couples2 <- c(couples2,i,j)
      }
    }
  }
  indremove1=c()
  if (is.null(couples1)==FALSE){
    couples1 <- matrix(couples1,ncol=2,byrow=T)
    v1 <- unique(couples1[,1])
    v2 <- unique(couples1[,2])
    indkeep1 <- setdiff(v1,v2)
    indremove1 <- v2
  }

  indremove2=c()
  if (is.null(couples2)==FALSE){
    couples2 <- matrix(couples2,ncol=2,byrow=T)
    w1 <- unique(couples2[,1])
    w2 <- unique(couples2[,2])
    indkeep2 <- setdiff(w1,w2)
    indremove2 <- w2
  }

  remove <- union(indremove1,indremove2)
  if (is.null(remove)==FALSE){
    keep=colnames(V_chr6A_QTL)[-remove]
    V_chr6A=V_chr6A[,colnames(V_chr6A)%in%c(markers_noQTL,keep)]
    markers_chr6A=colnames(V_chr6A)
  }

}


nb_QTL=sum(markers_chr6A%in%heading_QTL)
nb_QTL

## And then for the others:
couples1 <- c()
nb <- dim(V_chr6A)[2]
similarities <- matrix(NA,nb,nb)
for (i in 1:(nb-1)){
  for (j in (i+1):nb){
    similarities[i,j] <- sum(V_chr6A[,i]==(V_chr6A[,j]))
    if (similarities[i,j] == n){
      couples1 <- c(couples1,i,j)
    }
  }
}

couples1 <- matrix(couples1,ncol=2,byrow=T)

v1 <- unique(couples1[,1])
v2 <- unique(couples1[,2])
indkeep1 <- setdiff(v1,v2)
#We check that we are not deleting any heading QTLs:
keep=intersect(which(markers_chr6A%in%heading_QTL),v2)
# keep=integer(0): it's ok, we are not deleting heading QTL
indremove1 <- v2

couples2 <- c()
nb <- dim(V_chr6A)[2]
similarities <- matrix(NA,nb,nb)
for (i in 1:(nb-1)){
  for (j in (i+1):nb){
    similarities[i,j] <- sum(V_chr6A[,i]==(1-V_chr6A[,j]))
    if (similarities[i,j] == n){
      couples2 <- c(couples2,i,j)
    }
  }
}

couples2 <- matrix(couples2,ncol=2,byrow=T)
w1 <- unique(couples2[,1])
w2 <- unique(couples2[,2])
#We check that we are not deleting any heading QTLs:
keep=intersect(which(markers_chr6A%in%heading_QTL),w2)
# keep=integer(0): it's ok, we are not deleting heading QTL
indremove2 <- w2

remove <- union(indremove1,indremove2)

V_chr6A <- V_chr6A[,-remove]
nb_QTL=sum(heading_QTL%in%colnames(V_chr6A))
ind_QTL=which(colnames(V_chr6A)%in%heading_QTL) #indices of the columns of V_chr6A that correspond to heading QTLs
dim(V_chr6A)

sspop <- readRDS("Data/resPCOAdf.Rds") #5 new covariates to control for subpopulation structure, obtained by PCA
V_tilde_chr6A= as.matrix(cbind(Intercept=rep(1,n),sspop,V_chr6A))

V_tilde_chr6A[,-1] <- scale(V_tilde_chr6A[,-1],center=T,scale=T)
data_chr6A=list(V_tilde=V_tilde_chr6A,nb_QTL=nb_QTL,ind_QTL=ind_QTL+6,markers=colnames(V_chr6A))
save(data_chr6A,file="data_real_chr6A.Rdata")



## Then we applied SAEMVS:

load("Saves/data_real_chr6A.Rdata")
load("Saves/data_obs.Rdata")
sspop <- readRDS("Data/resPCOAdf.Rds")
data_chromosome <- read.csv2("Data/carte_Axiom-TABW420k_WGAv1.csv",header=TRUE, sep = " ")

varieties=unique(data_obs$GENOTYPE)
obs_date=unique(data_obs$Day)
n=length(unique(data_obs$GENOTYPE))
J=length(obs_date)
Id=rep(c(1:n),each=J)
id=as.matrix(Id)
Y=c()
for (i in varieties){
  for (j in obs_date){
    Y=c(Y,data_obs[(data_obs$GENOTYPE==i)&(data_obs$Day==j),]$Y)
  }
}
t=rep(unique(data_obs$Day),n)

V_tilde=data_chr6A$V_tilde
nb_QTL=data_chr6A$nb_QTL
markers=data_chr6A$markers
p=length(markers)
K=dim(sspop)[2]

source('R/Functions_SAEMVS_real_data.R')

# Initial values for the parameters to be estimated

ind_QTL=data_chr6A$ind_QTL
pos_QTLs=data_chromosome$V2[data_chromosome$name%in%data_chr6A$markers[ind_QTL-6]]
names_init_QTL=c()
for (k in 1:length(pos_QTLs)){
  names_init_QTL=c(names_init_QTL,data_chromosome$name[which(data_chromosome$V1=="chr6A" & data_chromosome$V2<=pos_QTLs[k]+10^6 & data_chromosome$V2>=pos_QTLs[k]-10^6)])
}
ind_QTL=which(colnames(V_tilde)%in%names_init_QTL)

mu_init <- 20
beta_tilde_init <- c(mu_init,rep(0.25,K),rep(0.1,p))
beta_tilde_init[ind_QTL]=rep(0.25,length(ind_QTL))
alpha_init=0.5
sigma2_init <- 80
Gamma2_init <- 50
eta_init <- 5
Omega2_init <- 50

# Hyperparameter
nu0 <- 0.00006        #spike parameter
nu1 <- 10             #slab parameter
nu_Gamma <- 1         #prior parameter of Gamma2
lb_Gamma <- 1         #prior parameter of Gamma2
nu_sigma <- 1         #prior parameter of sigma2
lb_sigma <- 1         #prior parameter of sigma2
a <- 1                #prior parameter of alpha
b <- p                #prior parameter of alpha
sigma2_mu <- 100      #prior parameter of mu
sigma2_lb <- 100      #prior parameter of lambda
sigma2_eta <- 100     #prior parameter of eta
nu_Omega <- 1         #prior parameter of Omega2
lb_Omega <- 1         #prior parameter of Omega2

tau <- 0.95           #annealing parameter
param_init <- list(beta_tilde=beta_tilde_init,alpha=alpha_init,Gamma2=Gamma2_init,sigma2=sigma2_init,eta=eta_init,Omega2=Omega2_init)
hyperparam <- list(nu0=nu0,nu1=nu1,nu_Gamma=nu_Gamma,lb_Gamma=lb_Gamma,nu_sigma=nu_sigma,lb_sigma=lb_sigma,a=a,b=b,sigma2_mu=sigma2_mu,sigma2_lb=sigma2_lb,sigma2_eta=sigma2_eta,nu_Omega=nu_Omega,lb_Omega=lb_Omega,tau=tau)

niter <- 400       #number of iterations of the algorithm
nburnin <- 250     #number of burn-in iterations
niterMH_phi <- 1   #number of iterations of the Metropolis Hastings algorithm for phi at S-step
niterMH_psi<- 1    #number of iterations of the Metropolis Hastings algorithm for psi at S-step

## To be uncommented to test the MCMC-SAEM algorithm for MAP calculation:
# res <- SAEM_MAP(niter, nburnin, niterMH_phi,niterMH_psi,Y,t,id,V_tilde, param_init,hyperparam)
# par(mfrow=c(2,2))
# for (k in 1:24){
#   if (k==1){
#     plot(res$beta_tilde[k,],xlab="Iterations",ylab="",main=expression(mu),type='l')
#     abline(v=nburnin,col='red')
#   }
#   if (k%in%c(2:6)){
#     plot(res$beta_tilde[k,],xlab="Iterations",ylab="",main=substitute(exp = lambda[valeur_k],env = list(valeur_k = k-1)),type='l')
#     abline(v=nburnin,col='red')
#   }
#   if (k%in%c(7:24)){
#     plot(res$beta_tilde[k,],xlab="Iterations",ylab="",main=substitute(exp = beta[valeur_k],env = list(valeur_k = k-6)),type='l')
#     abline(v=nburnin,col='red')
#     abline(h=0,col='blue')
#   }
# }
# plot(res$Gamma2,xlab="Iterations",ylab="",main=expression(Gamma^2),type='l')
# plot(res$sigma2,xlab="Iterations",ylab="",main=expression(sigma^2),type='l')
# plot(res$alpha,xlab="Iterations",ylab="",main=expression(alpha),type='l')
# plot(res$eta,xlab="Iterations",ylab="",main=expression(eta),type='l')
# plot(res$Omega2,xlab="Iterations",ylab="",main=expression(Omega^2),type='l')


### SAEMVS procedure

M <- 10
Delta <- 10^(seq(-4.5,-3,length.out = M))

res <- Model_selection(Delta,niter,nburnin,niterMH_phi,niterMH_psi,Y,t,id,V_tilde,param_init,hyperparam)

#save(res,file="Res_real_chr6A.Rdata")

load("Saves/Res_real_chr6A.Rdata")
res$graph

ggsave("Figure5.pdf")

## >>>>>>>>>>>>>>> Script for Figure 1 with example <<<<<<<<<<<<<<

rm(list=ls())

library(ggplot2)
library(cowplot)
library(mvnfast)
library(doParallel)
library(ggpattern)

######################### Example to obtain one bar from Figure 1 #########################

source('../R/Functions_SAEMVS_simu.R')

S=10                           #number of simulated data-sets tested for each combination of parameter
# In the paper S=100, but for reasons of computation time,
# in this example we put S=10 to have an execution time of about 12 min

n=100                          #number of individuals in {100,200}
J=10                           #number of repetitions
p=500                          #number of covariates in {500,2000,5000}
sigma2=30                      #residual variance of y
Gamma2=200                     #inter-individual variance of phi_i in {200,1000,2000}
mu=1200                        #intercept
beta=c(100,50,20,rep(0,p-3))   #covariate fixed effects vector
beta_tilde=c(mu,beta)          #useful reformulation
psi1=200                       #parameter of the function g
psi2=300                       #parameter of the function g
psi=c(psi1,psi2)

M=20
Delta=10^(seq(-2,2,length.out = M)) #grid of nu0 values

niter=500         #number of iterations of the algorithm
nburnin=350       #number of burn-in iterations
niterMH_phi=1     #number of iterations of the Metropolis Hastings algorithm for phi at S-step
niterMH_psi=10    #number of iterations of the Metropolis Hastings algorithm for psi at S-step

#We apply the SAEMVS method on the S data sets in parallel using ncore cores:
ncore = 10
cl = makeCluster(ncore)
registerDoParallel(cl)

debut=Sys.time()
resTot<-foreach(s = 1:S, .packages = c("nlme","mvnfast","doParallel","ggplot2","cowplot","glmnet")) %dopar% {

  print(s)
  ############# Model simulation #############
  set.seed(s)

  # Simulation of individual parameters phi_i
  phi=rep(0,n)
  V_tilde=matrix(NA,nrow=n,ncol=p+1)
  V_tilde[,1]=rep(1,n)
  for (i in 1:n) {
    V_tilde[i,2:(p+1)]=rmvn(1,mu=rep(0,p),sigma=diag(1,p))
  }
  V_tilde[,c(2:(p+1))]=scale(V_tilde[,c(2:(p+1))])

  for (i in 1:n){
    phi[i]=sum(beta_tilde*V_tilde[i,]) + rnorm(1,mean=0,sd=sqrt(Gamma2))
  }
  V=V_tilde[,-1]

  # Simulation of observed data y
  Y=rep(0,n*J)                    #vector of y_ij with y_ij=Y[(i-1)*J+j]
  t=rep(seq(150,3000,length=J),n) #vector of t_ij with x_ij=x[(i-1)*J+j]
  for (i in 1:n){
    for (j in 1:J) {
      Y[(i-1)*J+j]=g(phi_i=phi[i],psi=psi,t_ij=t[(i-1)*J+j])+rnorm(1,mean=0,sd=sqrt(sigma2))
    }
  }

  Id=rep(c(1:n),each=J)
  data=data.frame(Id,Y,t)
  id=as.matrix(Id)

  # Initial values for the parameters to be estimated
  betatilde_init=c(1400,rep(100,10),rep(1,p-10))
  sigma2_init=100
  Gamma2_init=5000
  eta_init=rep(400,2)
  alpha_init=0.5

  Omega_init=rep(20,2) #Omega=(w_1^2,w_2^2)

  # Fixed hyperparameters
  nu0=0.5           #spike parameter
  nu1=12000         #slab parameter
  nu_Gamma=1        #prior parameter of Gamma2
  lb_Gamma=1        #prior parameter of Gamma2
  nu_sigma=1        #prior parameter of sigma2
  lb_sigma=1        #prior parameter of sigma2
  a=1               #prior parameter of alpha
  b=p               #prior parameter of alpha
  sigma2_mu=3000^2  #prior parameter of mu
  rho2=rep(1200,2)  #prior parameter of eta

  tau2=c(rep(c(rep(1,39),0.9),12),rep(1,20)) #vector to be multiplied to Omega to limit the estimation error between the initial model and the extended model
  tau=c(0.98,tau2) #vector containing the annealing parameter and tau2

  #Initialisation
  param_init=list(beta_tilde=betatilde_init,alpha=alpha_init,Gamma2=Gamma2_init,sigma2=sigma2_init,eta=eta_init,Omega=Omega_init)
  hyperparam=list(nu0=nu0,nu1=nu1,nu_Gamma=nu_Gamma,lb_Gamma=lb_Gamma,nu_sigma=nu_sigma,lb_sigma=lb_sigma,a=a,b=b,sigma2_mu=sigma2_mu,rho2=rho2,tau=tau)

  ############# SAEMVS #############
  res=Model_selection(Delta,niter,nburnin,niterMH_phi,niterMH_psi,Y,t,id,V_tilde,param_init,hyperparam,s=s)
  nu0_select=res$nu0_select                                #nu0 selected by the eBIC criterion
  beta_tildehat_select=res$beta_tildehat_select            #beta_tilde MAP estimate associated with this nu0
  Gamma2hat_select=res$Gamma2hat_select                    #Gamma^2 MAP estimate associated with this nu0
  sigma2hat_select=res$sigma2hat_select                    #sigma^2 MAP estimate associated with this nu0
  psihat_select=res$psihat_select                          #psi estimate associated with this nu0

  result=list(nu0_select=nu0_select,beta_tildehat_select=beta_tildehat_select,Gamma2hat_select=Gamma2hat_select,sigma2hat_select=sigma2hat_select,psihat_select=psihat_select)

  result
}
stopCluster(cl)
save(resTot,file="Res_test_Fig1.Rdata")

load("Res_test_Fig1.Rdata")

nu0_select=rep(NA,S)
beta_tildehat_select=matrix(0,nrow=p+1,ncol=S)
Gamma2hat_select=rep(NA,S)
sigma2hat_select=rep(NA,S)
psihat_select=matrix(NA,nrow=2,ncol=S)

for (s in 1:S){
  nu0_select[s]=resTot[[s]]$nu0_select
  beta_tildehat_select[,s]=resTot[[s]]$beta_tildehat_select
  Gamma2hat_select[s]=resTot[[s]]$Gamma2hat_select
  sigma2hat_select[s]=resTot[[s]]$sigma2hat_select
  psihat_select[,s]=resTot[[s]]$psihat_select
}

betahat_select=beta_tildehat_select[-1,]

## On the S datasets, we look at the number of data-sets on which SAEMVS selects the correct model: nb_exact_model,
## a model that contains the correct model but not equal (there are false positives (FP) but not false negatives (FN)): nb_cont,
## a model that is included in the correct model (FN but not FP): nb_includ,
## or a model that contains both FP and FN: nb_FP_FN

Model_size=colSums(betahat_select!=0)      #Sizes of model selected for the S datasets
table(Model_size)

v=which(Model_size==3)  #Sets of datasets for which SAEMVS has selected a model of size 3

nb_exact_model=0
for (s in v){
  nb_exact_model=nb_exact_model+(all(betahat_select[1:3,s]!=rep(0,3)))
}
nb_exact_model

nb_cont=0
for (s in (setdiff(S,v))){
  nb_cont=nb_cont + (all(betahat_select[1:3,s]!=rep(0,3)))
}
nb_cont

nb_includ=0
for (s in 1:S){
  cond=(all(betahat_select[4:p,s]==rep(0,p-3)) & (betahat_select[1,s]==0 | betahat_select[2,s]==0 | betahat_select[3,s]==0))
  nb_includ=nb_includ+cond
}
nb_includ

nb_FP_FN=100-nb_exact_model-nb_cont-nb_includ


######################### Script for Figure 1 #########################

rm(list=ls())

library(ggplot2)
library(ggpattern)

# The saved results are retrieved:
file_name=c("Res_simuA1.Rdata","Res_simuA2.Rdata","Res_simuA3.Rdata","Res_simuA4.Rdata","Res_simuA5.Rdata",
            "Res_simuA6.Rdata","Res_simuA7.Rdata","Res_simuA8.Rdata","Res_simuA9.Rdata","Res_simuA10.Rdata",
            "Res_simuA11.Rdata","Res_simuA12.Rdata","Res_simuA13.Rdata","Res_simuA14.Rdata","Res_simuA15.Rdata",
            "Res_simuA16.Rdata","Res_simuA17.Rdata","Res_simuA18.Rdata")

N=length(file_name)
nb_exact_model=rep(0,N)
nb_cont=rep(0,N)
nb_includ=rep(0,N)
nb_FP_FN=rep(0,N)

S=100

for (k in 1:N){
  load(file_name[k])
  resTot=res$resTot

  p=length(resTot[[1]]$betachap_select)-1
  beta_tildehat_select=matrix(0,nrow=p+1,ncol=S)

  for (s in 1:S){
    beta_tildehat_select[,s]=resTot[[s]]$betachap_select
  }

  betahat_select=beta_tildehat_select[-1,]

  Model_size=colSums(betahat_select!=0)      #Sizes of model selected for the S datasets
  v=which(Model_size==3)  #Sets of datasets for which SAEMVS has selected a model of size 3

  for (s in v){
    nb_exact_model[k]=nb_exact_model[k]+(all(betahat_select[1:3,s]!=rep(0,3)))
  }

  for (s in (setdiff(c(1:100),v))){
    nb_cont[k]=nb_cont[k] + (all(betahat_select[1:3,s]!=rep(0,3)))
  }

  for (s in 1:S){
    cond=(all(betahat_select[4:p,s]==rep(0,p-3)) & (betahat_select[1,s]==0 | betahat_select[2,s]==0 | betahat_select[3,s]==0))
    nb_includ[k]=nb_includ[k]+cond
  }

  nb_FP_FN[k]=100-nb_exact_model[k]-nb_cont[k]-nb_includ[k]

}

data_n100 <- data.frame(Gamma2=rep(c(200,1000,2000),each=3),p=rep(c(500,2000,5000),3),Result=rep(c("Exact","FP but not FN","FN but not FP","FP and FN"),each=9),count=c(nb_exact_model[1:9],nb_cont[1:9],nb_includ[1:9],nb_FP_FN[1:9]))
data_n100
data_n100$Gamma2<-as.factor(data_n100$Gamma2)
data_n100$p<-as.factor(data_n100$p)
data_n100$Result<-factor(data_n100$Result,levels=c("FP and FN","FN but not FP","FP but not FN","Exact"))

ggplot(mapping=aes(x=p, y=count, fill=p,pattern=Result))+ geom_bar_pattern(data=data_n100,width=0.7,color="black",pattern_fill="black",pattern_density=0.05,pattern_spacing=0.025,stat='identity')+
  labs(x=expression(Gamma^2), y = "Proportion (in %)") +scale_fill_brewer(palette="Set2")+ scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100))+ scale_pattern_manual(values = c("Exact"="none", "FP but not FN"="stripe","FN but not FP"="circle","FP and FN"="crosshatch"),breaks = c("Exact","FP but not FN","FN but not FP","FP and FN"))+
  facet_grid(~Gamma2,switch = "x")+
  theme_bw()+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size=30),axis.text=element_text(size=30), axis.title=element_text(size=30,face="bold"),title=element_text(size=30,face="bold"),legend.title = element_text(size=30,face="bold"),legend.text = element_text(size=30))+guides(pattern=guide_legend(override.aes = list(fill="white")),fill=guide_legend(override.aes = list(pattern="none")))


data_n200 <- data.frame(Gamma2=rep(rep(c(200,1000,2000),each=3),2),p=rep(rep(c(500,2000,5000),3),2),Result=rep(c("Exact","FP but not FN","FN but not FP","FP and FN"),each=9),count=c(nb_exact_model[10:N],nb_cont[10:N],nb_includ[10:N],nb_FP_FN[10:N]))
data_n200
data_n200$Gamma2<-as.factor(data_n200$Gamma2)
data_n200$p<-as.factor(data_n200$p)
data_n200$Result<-factor(data_n200$Result,levels=c("FP and FN","FN but not FP","FP but not FN","Exact"))

ggplot(mapping=aes(x=p, y=count, fill=p,pattern=Result))+ geom_bar_pattern(data=data_n200,width=0.7,color="black",pattern_fill="black",pattern_density=0.05,pattern_spacing=0.025,stat='identity')+
  labs(x=expression(Gamma^2), y = "Proportion (in %)") +scale_fill_brewer(palette="Set2")+ scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100))+ scale_pattern_manual(values = c("Exact"="none", "FP but not FN"="stripe","FN but not FP"="circle","FP and FN"="crosshatch"),breaks = c("Exact","FP but not FN","FN but not FP","FP and FN"))+
  facet_grid(~Gamma2,switch = "x")+
  theme_bw()+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size=30),axis.text=element_text(size=30), axis.title=element_text(size=30,face="bold"),title=element_text(size=30,face="bold"),legend.title = element_text(size=30,face="bold"),legend.text = element_text(size=30))+guides(pattern=guide_legend(override.aes = list(fill="white")),fill=guide_legend(override.aes = list(pattern="none")))


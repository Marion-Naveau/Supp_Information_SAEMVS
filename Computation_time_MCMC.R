## >>>>>>>>>>>> Script to determine the computation time of an MCMC implementation <<<<<<<<<<<

rm(list=ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(doParallel)
library(microbenchmark)

source("R/nimble_fit_functions.R")

## As explained in subsection 5.4, to stabilise the MCMC procedure, the prior on sigma^2
## is modified to an uniform distribution on [0,200] for both methods MCMC and MCMC-SAEM.

############## Computation time ###########

S=50                           #number of simulated data-sets tested
n=200                          #number of individuals
J=10                           #number of repetitions
p=500                          #number of covariates in {500,700,1000,1500,2000,2500}
sigma2=30                      #residual variance of y
Gamma2=200                     #inter-individual variance of phi_i in {200,1000,2000}
mu=1200                        #intercept
beta=c(100,50,20,rep(0,p-3))   #covariate fixed effects vector
beta_tilde=c(mu,beta)          #useful reformulation
psi1=200                       #parameter of the function g
psi2=300                       #parameter of the function g
psi=c(psi1,psi2)

g <- function(phi_i,psi,t_ij){
  return(psi[1]/(1+exp(-(t_ij-phi_i)/psi[2])))
}

ncore = 10
cl = makeCluster(ncore)
registerDoParallel(cl)

Total_Time<-foreach(s = 1:S, .packages = c("nlme","mvnfast","doParallel","ggplot2","cowplot","glmnet","microbenchmark","tidyverse","nimble")) %dopar% {

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

  data_input = data %>%
    as_tibble() %>%
    left_join(V %>%
                as_tibble() %>%
                rowid_to_column("Id")) %>%
    rename(times = t, i = Id)

  nu0=0.04
  niter=3000

  nimbleMCMC_object = OneVariableParamSpikeAndSlabDelay_marginal_delta_setup(data_input,nu_0=nu0,niter=niter)

  ## Time of the execution for MCMC sampling (1 chain)
  time=TimeOneVariableParamSpikeAndSlabDelay_marginal_delta(Cmcmc = nimbleMCMC_object, nu_0 = nu0, niter = niter)

  time
}
stopCluster(cl)

save(Total_Time,file="Res_Test_comp_time_MCMC.Rdata")

load("Saves/Res_Test_comp_time_MCMC.Rdata")

time=rep(0,50)
for (s in 1:50){
  time[s]=Total_Time[[s]]
}
time
min(time)
mean(time)

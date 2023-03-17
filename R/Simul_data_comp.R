## This file contains a script that simulates the 100 data-sets used in the comparison study in section 5.1

rm(list=ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(doParallel)
library(psych)

S=100
n=200                                             #number of individuals
Ji=rep(12,n)                              #number of repetitions
p=500                                             #number of covariates
q=2                                               #number of random effects
sigma2=10^(-3)                                    #residual variance of y
Gamma=matrix(c(0.2,0.05,0.05,0.1),nrow=q,ncol=q)  #inter-individual variance of phi_i
mu=c(6,8)                                         #intercept
beta=matrix(c(3,2,1,rep(0,p-3),0,0,3,2,1,rep(0,p-5)),nrow=p,ncol=q,byrow=FALSE)  #covariate fixed effects matrix
beta_tilde=rbind(mu,beta)                         #useful reformulation
dose=100                                          #parameter of the function g


g <- function(dose,t_ij,phi_i){
  return((dose*phi_i[1]/(30*phi_i[1]-phi_i[2]))*(exp(-phi_i[2]*t_ij/30)-exp(-phi_i[1]*t_ij)))
}

data=list()
phi_simul=list()
V_simul=list()
for (s in 1:S){
  set.seed(s)
  phi=matrix(0,nrow=n,ncol=q)
  V_tilde=matrix(NA,nrow=n,ncol=p+1)
  V_tilde[,1]=rep(1,n)
  for (i in 1:n) {
    V_tilde[i,-1]=rbinom(p,size=1,prob=0.2)
  }
  V_tilde[,-1]=scale(V_tilde[,-1])
  
  for (i in 1:n){
    phi[i,]=c(t(beta_tilde)%*%V_tilde[i,]) + c(rmvn(1,mu=rep(0,q),sigma=Gamma))
  }
  V=V_tilde[,-1]
  
  # Simulation of observed data y
  
  Y=rep(0,sum(Ji))                    #vector of y_ij with y_ij=Y[sum(Ji[1:(i-1)])+j]
  t=rep(c(0.05,0.15,0.25,0.4,0.5,0.8,1,2,7,12,24,40),n)
  
  for (i in 1:n){
    for (j in 1:Ji[i]) {
      if (i==1){
        Y[j]=g(phi_i=phi[i,],dose=dose,t_ij=t[j])+rnorm(1,mean=0,sd=sqrt(sigma2))
      } else {
        Y[sum(Ji[1:(i-1)])+j]=g(phi_i=phi[i,],dose=dose,t_ij=t[sum(Ji[1:(i-1)])+j])+rnorm(1,mean=0,sd=sqrt(sigma2))
      }
    }
  }
  
  Id=c()
  for (i in 1:n){
    Id=c(Id,rep(i,Ji[i]))
  }
  
  data[[s]]=data.frame(Id,Y,t)
  phi_simul[[s]]=phi
  V_simul[[s]]=V
  
}

save(data,file="data_comp.Rdata")
save(phi_simul,file="phi_simul.Rdata")
save(V_simul,file="V_simul.Rdata")

load("data_comp.Rdata")
data=data[[1]]
Id=data$Id
Y=data$Y
t=data$t

par(mfrow=c(2,2))
for (i in 1:20){
  plot(data[Id==i,]$t, data[Id==i,]$Y,type="p")
}

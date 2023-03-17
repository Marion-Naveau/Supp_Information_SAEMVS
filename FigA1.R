rm(list=ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(doParallel)
library(psych)

source('R/Functions_SAEMVS_multi.R') 

################################### Model simulation ###################################

n=200                                             #number of individuals
Ji=rep(12,n)                                      #number of repetitions
p=500                                             #number of covariates
q=2                                               #number of random effects
sigma2=10^(-3)                                    #residual variance of y
Gamma=matrix(c(0.2,0.05,0.05,0.1),nrow=q,ncol=q)  #inter-individual variance of phi_i
mu=c(6,8)                                         #intercept
beta=matrix(c(3,2,1,rep(0,p-3),0,0,3,2,1,rep(0,p-5)),nrow=p,ncol=q,byrow=FALSE)  #fixed effects matrix
beta_tilde=rbind(mu,beta)                         #useful reformulation
dose=100                                          #parameter of the function g


g <- function(dose,t_ij,phi_i){
  return((dose*phi_i[1]/(30*phi_i[1]-phi_i[2]))*(exp(-phi_i[2]*t_ij/30)-exp(-phi_i[1]*t_ij)))
}


# Simulation of individual parameters phi_i
s=2
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
data=data.frame(Id,Y,t)
id=as.matrix(Id)

# Data representation
data2=data[id<=10,]
ggplot(data2, aes(x=t, y=Y, group=Id,color=Id)) + geom_point() + geom_line() + theme_bw()


################################### MCMC-SAEM ###################################

# Initial values for the parameters to be estimated
betatilde_init=matrix(c(10,rep(1,10),rep(0.1,p-10),10,rep(1,10),rep(0.1,p-10)),nrow=p+1,ncol=q,byrow=FALSE)
sigma2_init=10^(-2)
Gamma_init=matrix(c(0.5,0.1,0.1,0.5),nrow=q,ncol=q)
alpha_init=c(0.5,0.5)

# Fixed hyperparameters
nu0=0.005               #spike parameter
nu1=1000               #slab parameter
nu_sigma=1             #prior parameter of sigma2
lb_sigma=1             #prior parameter of sigma2
a=c(1,1)               #prior parameter of alpha
b=c(p,p)               #prior parameter of alpha
sigma2_mu=5^2          #prior parameter of mu
Q=0.2*diag(1,q)
d=q+2

tau=0.98
param_init=list(beta_tilde=betatilde_init,alpha=alpha_init,Gamma=Gamma_init,sigma2=sigma2_init)
hyperparam=list(dose=dose,nu0=nu0,nu1=nu1,nu_sigma=nu_sigma,lb_sigma=lb_sigma,a=a,b=b,sigma2_mu=sigma2_mu,q=q,Q=Q,d=d,tau=tau)

niter=300         #number of iterations of the algorithm
nburnin=150       #number of burn-in iterations
niterMH_phi=1     #number of iterations of the Metropolis Hastings algorithm for phi at S-step

## To be uncommented to test the MCMC-SAEM algorithm for MAP calculation:
res=SAEM_MAP(niter, nburnin, niterMH_phi,Y,t,id, V_tilde, param_init = param_init,hyperparam = hyperparam,s=s)

par(mfrow=c(2,2))
for (k in 1:12){
  if (k==1){
    plot(res$beta_tilde[k,1,],xlab="Iterations",ylab="",main=expression(mu[1]),cex.main=2,cex.lab=2,type='l')
    abline(h=beta_tilde[k,1],col="red",lty="dashed")
  }
  else {
    plot(res$beta_tilde[k,1,],xlab="Iterations",ylab="",main=substitute(exp = paste(beta[valeur_k], " for ", phi[1]),
                                                                        env = list(valeur_k = k-1)),cex.main=2,cex.lab=2,type='l') 
    abline(h=beta_tilde[k,1],col="red",lty="dashed")
  }
}
for (k in 1:12){
  if (k==1){
    plot(res$beta_tilde[k,2,],xlab="Iterations",ylab="",main=expression(mu[2]),cex.main=2,cex.lab=2,type='l')
    abline(h=beta_tilde[k,2],col="red",lty="dashed")
  }
  else {
    plot(res$beta_tilde[k,2,],xlab="Iterations",ylab="",main=substitute(exp = paste(beta[valeur_k], " for ", phi[2]),
                                                                        env = list(valeur_k = k-1)),cex.main=2,cex.lab=2,type='l')
    abline(h=beta_tilde[k,2],col="red",lty="dashed")
  }
}
plot(res$Gamma[1,1,],xlab="Iterations",ylab="",main=expression(Gamma[11]),cex.main=2,cex.lab=2,type='l')
abline(h=Gamma[1,1],col="red",lty="dashed")
plot(res$Gamma[1,2,],xlab="Iterations",ylab="",main=expression(Gamma[12]),cex.main=2,cex.lab=2,type='l')
abline(h=Gamma[1,2],col="red",lty="dashed")
plot(res$Gamma[2,1,],xlab="Iterations",ylab="",main=expression(Gamma[21]),cex.main=2,cex.lab=2,type='l')
abline(h=Gamma[2,1],col="red",lty="dashed")
plot(res$Gamma[2,2,],xlab="Iterations",ylab="",main=expression(Gamma[22]),cex.main=2,cex.lab=2,type='l')
abline(h=Gamma[2,2],col="red",lty="dashed")
plot(res$sigma2,xlab="Iterations",ylab="",main=expression(sigma^2),cex.main=2,cex.lab=2,type='l')
abline(h=sigma2,col="red",lty="dashed")
plot(res$alpha[1,],xlab="Iterations",ylab="",main=expression(alpha[1]),cex.main=2,cex.lab=2,type='l')
abline(h=3/p,col="red",lty="dashed")
plot(res$alpha[2,],xlab="Iterations",ylab="",main=expression(alpha[2]),cex.main=2,cex.lab=2,type='l')
abline(h=3/p,col="red",lty="dashed")


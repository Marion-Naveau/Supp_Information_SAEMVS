## The file "Test_Simu_SAEMVS.R" allows to combine the functions of the file "Functions_SAEMVS_simu.R"
## to launch the SAEMVS method on the model as presented in the simulation study Section 5 on one simulated data-set.

rm(list=ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(doParallel)

source('R/Functions_SAEMVS_simu.R')

################################### Model simulation ###################################

n=200                          #number of individuals
J=10                           #number of repetitions
p=500                          #number of covariates
sigma2=30                      #residual variance of y
Gamma2=200                     #inter-individual variance of phi_i
mu=1200                        #intercept
beta=c(100,50,20,rep(0,p-3))   #covariate fixed effects vector
beta_tilde=c(mu,beta)          #useful reformulation
psi1=200                       #parameter of the function g
psi2=300                       #parameter of the function g
psi=c(psi1,psi2)

# Simulation of individual parameters phi_i
s=2
set.seed(s)
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

# Data representation
data2=data[id<=30,]
ggplot(data2, aes(x=t, y=Y, group=Id)) + geom_point() + geom_line() + theme_bw()


################################### MCMC-SAEM ###################################

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

param_init=list(beta_tilde=betatilde_init,alpha=alpha_init,Gamma2=Gamma2_init,sigma2=sigma2_init,eta=eta_init,Omega=Omega_init)
hyperparam=list(nu0=nu0,nu1=nu1,nu_Gamma=nu_Gamma,lb_Gamma=lb_Gamma,nu_sigma=nu_sigma,lb_sigma=lb_sigma,a=a,b=b,sigma2_mu=sigma2_mu,rho2=rho2,tau=tau)

niter=500         #number of iterations of the algorithm
nburnin=350       #number of burn-in iterations
niterMH_phi=1     #number of iterations of the Metropolis Hastings algorithm for phi at S-step
niterMH_psi=10    #number of iterations of the Metropolis Hastings algorithm for psi at S-step

## To be uncommented to test the MCMC-SAEM algorithm for MAP calculation:
# res=SAEM_MAP(niter, nburnin, niterMH_phi,niterMH_psi,Y,t,id, V_tilde, param_init = param_init,hyperparam = hyperparam,s=s)
# par(mfrow=c(2,2))
# for (k in 1:21){
#   if (k==1){
#     plot(res$beta_tilde[k,],xlab="Iterations",ylab="",main=expression(mu),type='l')
#     abline(h=beta_tilde[k],col="red")
#   }
#   else {
#     plot(res$beta_tilde[k,],xlab="Iterations",ylab="",main=substitute(exp = beta[valeur_k],
#                                                                 env = list(valeur_k = k-1)),type='l')
#     abline(h=beta_tilde[k],col="red")
#   }
# }
# plot(res$Gamma2,xlab="Iterations",ylab="",main=expression(Gamma^2),type='l')
# abline(h=Gamma2,col="red")
# plot(res$sigma2,xlab="Iterations",ylab="",main=expression(sigma^2),type='l')
# abline(h=sigma2,col="red")
# plot(res$alpha,xlab="Iterations",ylab="",main=expression(alpha),type='l')
# abline(h=3/p,col="red")
# plot(res$eta[1,],xlab="Iterations",ylab="",main=expression(eta[1]),type='l')
# abline(h=psi[1],col="red")
# plot(res$eta[2,],xlab="Iterations",ylab="",main=expression(eta[2]),type='l')
# abline(h=psi[2],col="red")

################################### SAEMVS procedure ###################################

M=10
Delta=10^(seq(-2,2,length.out = M)) #grid of nu0 values

res <- Model_selection(Delta,niter,nburnin,niterMH_phi,niterMH_psi,Y,t,id,V_tilde,param_init,hyperparam,s=s)

#save(res,file="Res_test.Rdata")

load("Saves/Res_test.Rdata")
res$graph


## >>>>>>>>>>>> Script to determine the computation time of MCMC-SAEM <<<<<<<<<<<

rm(list=ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(doParallel)
library(microbenchmark)

source('R/Functions_SAEMVS_simu.R')

## As explained in subsection 5.4, to stabilise the MCMC procedure, the prior on sigma^2
## is modified to an uniform distribution on [0,200] for both methods MCMC and MCMC-SAEM.

## Here you have the MCMC-SAEM algorithm with this change of prior:

SAEM_MAP <- function(niter,nburnin,niterMH_phi,niterMH_psi,Y,t,id,V_tilde,param_init,hyperparam,s) {
  #niter = number of iterations of the algorithm
  #nburnin = number of burn-in iterations (nburnin < niter)
  #niterMH_phi = number of iterations of the Metropolis Hastings algorithm for phi at S-step
  #niterMH_psi = number of iterations of the Metropolis Hastings algorithm for psi at S-step
  #Y = observed data ((y_ij)_{1<=j<=J})_{1<=i<=n}
  #t = times ((t_ij)_{1<=j<=J})_{1<=i<=n}
  #id = individual identifiers
  #V_tilde = matrix n*(p+1) such as for 1<=i<=n, line i is the vector of Vtilde_i=(1,V_i) in R^(p+1)
  #param_init = initialisation of the parameters to be estimated and Omega: list(beta_tilde,alpha,Gamma2,sigma2,eta,Omega)
  #hyperparam = list of fixed hyperparameters list(nu0,nu1,nu_Gamma,lb_Gamma,nu_sigma,lb_sigma,a,b,sigma2_mu,rho2,tau), where
  #  tau=vector containing the annealing parameter and the vector to be multiplied to Omega to limit the estimation error between the initial model and the extended model
  #s = seed

  set.seed(s)
  n = length(unique(id))   #number of individuals
  J = length(Y[id==1])     #number of observations per individual
  yi <- matrix(Y,J,n)
  ti <- matrix(t,J,n)

  Ntot = n*J               #total number of observations in the sample
  tV_tilde <- t(V_tilde)
  tV_tildeV_tilde <- tV_tilde%*%V_tilde
  p = dim(V_tilde)[2]-1

  # Fixed hyperparameters
  nu0=hyperparam$nu0
  nu1=hyperparam$nu1
  nu_Gamma=hyperparam$nu_Gamma
  lb_Gamma=hyperparam$lb_Gamma
  nu_sigma=hyperparam$nu_sigma
  lb_sigma=hyperparam$lb_sigma
  a=hyperparam$a
  b=hyperparam$b
  rho2=hyperparam$rho2
  sigma2_mu=hyperparam$sigma2_mu
  tau=hyperparam$tau

  # Initialiqation of the estimation vectors over the iterations
  beta_tilde = matrix(NA,nrow=p+1,ncol=niter+1)
  alpha = rep(NA,niter+1)
  Gamma2 = rep(NA,niter+1)
  sigma2 = rep(NA,niter+1)
  eta=matrix(NA,nrow=2,ncol=niter+1)
  Omega=matrix(NA,nrow=2,ncol=niter+1) #each line corresponds to the updates of (omega_1)^2 and (omega_2)^2 respectively

  beta_tilde[,1] = param_init$beta_tilde
  alpha[1] <- param_init$alpha
  Gamma2[1] <- param_init$Gamma2
  sigma2[1] <- param_init$sigma2
  eta[,1] <- param_init$eta
  Omega[,1] <- param_init$Omega

  ## Initialisation of sufficient statistics useful for updating estimates
  s1 = rep(0,niter+1)                    #stat associated with y_ij
  s2 = rep(0,niter+1)                    #stat associated with the sum of phi^2
  s3 = matrix(0,nrow=n,ncol=niter+1)     #stats associated with phi_i
  s4 = matrix(0,nrow=2,ncol=niter+1)     #stat associated with psi^2
  s5 = matrix(0,nrow=2,ncol=niter+1)     #stat associated with psi
  mco=matrix(0,n,J)                      #will be useful for the s1 update

  ## Definition of the step sequence for stochastic approximation
  gamma = rep(1,niter)
  for (k in (nburnin+1):niter){
    gamma[k] = 1/(k-nburnin)^(2/3)
  }

  ## Initialisation of individual parameters phi and fixed effects psi

  phi=matrix(NA,nrow=n,ncol=niter+1) #matrix containing in column the (phi_i)_i simulated at each iteration k
  psi=matrix(NA,nrow=2,ncol=niter+1) #idem for psi
  phi[,1] <- V_tilde%*%beta_tilde[,1]
  psi[,1] <- eta[,1]

  ## MCMC-SAEM:
  for (k in 1:niter) {
    ## S-step: Metropolis Hastings within Gibbs
    phi_tilde=matrix(NA,nrow=n,ncol=niterMH_phi+1) #matrix with in column r the vector (phi_i^(k,r))_i created by M-H
    psi_tilde=matrix(NA,nrow=2,ncol=niterMH_psi+1) #idem for psi
    phi_tilde[,1]=phi[,k]
    psi_tilde[,1]=psi[,k]
    Vbetak=V_tilde%*%beta_tilde[,k]

    for (i in 1:n) {
      for (r in 1:niterMH_phi){
        phi_ci=rnorm(1,mean=Vbetak[i],sd=sqrt(Gamma2[k])) #candidate
        logratio= sum(dnorm(yi[,i],mean=g(phi_ci,psi[,k],ti[,i]),sd=sqrt(sigma2[k]),log=TRUE))-
          sum(dnorm(yi[,i],mean=g(phi_tilde[i,r],psi[,k],ti[,i]),sd=sqrt(sigma2[k]),log=TRUE))
        u = runif(1)
        cond = (log(u)<=logratio)
        if (cond) {
          phi_tilde[i,r+1]=phi_ci
        }
        else {
          phi_tilde[i,r+1]=phi_tilde[i,r]
        }
      }
    }
    phi[,k+1]=phi_tilde[,niterMH_phi+1]

    for (m in 1:2){
      for (r in 1:niterMH_psi){
        psi_mc=rnorm(1,mean=eta[m,k],sd=sqrt(Omega[m,k]))
        phi_k=rep(phi[,k+1],each=J)
        if (m==1){
          psi_c=c(psi_mc,psi[2,k])
          psi_r=c(psi_tilde[1,r],psi[2,k])
        }
        if (m==2){
          psi_c=c(psi_tilde[1,niterMH_psi+1],psi_mc)
          psi_r=c(psi_tilde[1,niterMH_psi+1],psi_tilde[2,r])
        }

        logratio = sum(dnorm(Y,mean=g(phi_k,psi_c,t),sd=sqrt(sigma2[k]),log=TRUE))-
          sum(dnorm(Y,mean=g(phi_k,psi_r,t),sd=sqrt(sigma2[k]),log=TRUE))

        u = runif(1)
        cond = (log(u)<=logratio)
        if (cond) {
          psi_tilde[m,r+1]=psi_mc
        }
        if (cond==FALSE) {
          psi_tilde[m,r+1]=psi_tilde[m,r]
        }
      }
    }

    psi[,k+1]=psi_tilde[,niterMH_psi+1]

    ## SA-step: updating of sufficient statistics by stochastic approximation
    for (i in 1:n){
      mco[i,] = (yi[,i]-g(phi[i,k+1],psi[,k+1],ti[,i]))^2
    }

    s1[k+1]=s1[k]+gamma[k]*(sum(mco)-s1[k])
    s2[k+1]=s2[k]+gamma[k]*(sum(phi[,k+1]^2)-s2[k])
    s3[,k+1]=s3[,k]+gamma[k]*(phi[,k+1]-s3[,k])
    p_star_k=p_star(beta_tilde[-1,k],alpha[k],nu0,nu1,p)
    d_tilde_star=c(1/sigma2_mu,(1-p_star_k)/nu0+p_star_k/nu1)
    s4[,k+1]=s4[,k]+gamma[k]*(psi[,k+1]^2-s4[,k])
    s5[,k+1]=s5[,k]+gamma[k]*(psi[,k+1]-s5[,k])

    ## M-step: update of (beta_tilde,alpha,Gamma2,sigma2,eta) and Omega
    D_star=diag(d_tilde_star)
    beta_tilde[,k+1]=solve(Gamma2[k]*D_star+tV_tildeV_tilde,tV_tilde%*%s3[,k+1])

    alpha[k+1]=(sum(p_star_k)+a-1)/(p+b+a-2)

    Vbeta=V_tilde%*%beta_tilde[,k+1]
    Gamma2[k+1]=max(tau[1]*Gamma2[k],(sum((Vbeta)^2)+nu_Gamma*lb_Gamma+s2[k+1]-2*sum((Vbeta)*s3[,k+1]))/(n+nu_Gamma+2))

    if ((s1[k+1])/(n*J)<=200){
      sigma2[k+1]=(s1[k+1])/(n*J)
    }
    else {
      sigma2[k+1]=200
    }

    eta[1,k+1]=s5[1,k+1]/(1+Omega[1,k]/rho2[1])
    eta[2,k+1]=s5[2,k+1]/(1+Omega[2,k]/rho2[2])

    Omega[,k+1]=tau[k+1]*Omega[,k]

  }

  return(list(beta_tilde=beta_tilde,alpha=alpha,Gamma2=Gamma2,sigma2=sigma2,eta=eta,Omega=Omega,phi=phi,psi=psi,s1=s1,s2=s2,s3=s3,s4=s4,s5=s5))
}


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

niter=500         #number of iterations of the algorithm
nburnin=350       #number of burn-in iterations
niterMH_phi=1     #number of iterations of the Metropolis Hastings algorithm for phi at S-step
niterMH_psi=10    #number of iterations of the Metropolis Hastings algorithm for psi at S-step

ncore = 10
cl = makeCluster(ncore)
registerDoParallel(cl)

Total_Time<-foreach(s = 1:S, .packages = c("nlme","mvnfast","doParallel","ggplot2","cowplot","glmnet","microbenchmark")) %dopar% {

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

  # Initial values for the parameters to be estimated with MCMC-SAEM
  betatilde_init=c(1400,rep(100,10),rep(1,p-10))
  sigma2_init=100
  Gamma2_init=500
  eta_init=rep(400,2)
  alpha_init=0.5

  Omega_init=rep(20,2) #Omega=(w_1^2,w_2^2,w_3^2)

  # Fixed hyperparameters
  nu0=0.04           #spike parameter
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

  # MCMC-SAEM
  time=microbenchmark(SAEM_MAP(niter, nburnin, niterMH_phi,niterMH_psi,Y,t,id, V_tilde, param_init,hyperparam,s),times=1,unit="s")

  time
}
stopCluster(cl)

save(Total_Time,file="Res_Test_comp_time_MCMCSAEM.Rdata")

load("Saves/Res_Test_comp_time_MCMCSAEM.Rdata")
time=rep(0,S)
for (s in 1:S){
  time[s]=summary(Total_Time[[s]])$min
}
time
min(time) #in seconds
mean(time)

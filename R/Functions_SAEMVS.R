## This file contains all the functions needed to run the SAEMVS procedure on the detailed example in one-dimension setting

g <- function(phi_i,psi,t_ij){
  return(psi[1]/(1+exp(-(t_ij-phi_i)/psi[2])))
}

p_star <- function(beta,alpha,nu0,nu1,p){
  beta=matrix(beta)
  norm1 <- apply(beta,1,function(x) dnorm(x,0,sqrt(nu1)))
  norm0 <- apply(beta,1,function(x) dnorm(x,0,sqrt(nu0)))
  p_star <- norm1*alpha/(norm1*alpha+norm0*(1-alpha))
  return(p_star)
}

SAEM_MAP <- function(niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init,hyperparam,s) {
  #niter = number of iterations of the algorithm
  #nburnin = number of burn-in iterations (nburnin < niter)
  #niterMH_phi = number of iterations of the Metropolis Hastings algorithm for phi at S-step
  #Y = observed data ((y_ij)_{1<=j<=J})_{1<=i<=n}
  #t = times ((t_ij)_{1<=j<=J})_{1<=i<=n}
  #id = individual identifiers
  #V_tilde = matrix n*(p+1) such as for 1<=i<=n, line i is the vector of Vtilde_i=(1,V_i) in R^(p+1)
  #param_init = initialisation of the parameters to be estimated: list(beta_tilde,alpha,Gamma2,sigma2)
  #hyperparam = list of fixed hyperparameters list(nu0,nu1,nu_Gamma,lb_Gamma,nu_sigma,lb_sigma,a,b,psi,sigma2_mu,tau), where tau=annealing parameter
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
  psi=hyperparam$psi
  sigma2_mu=hyperparam$sigma2_mu
  tau=hyperparam$tau

  # Initialiqation of the estimation vectors over the iterations
  beta_tilde = matrix(NA,nrow=p+1,ncol=niter+1)
  alpha = rep(NA,niter+1)
  Gamma2 = rep(NA,niter+1)
  sigma2 = rep(NA,niter+1) 

  beta_tilde[,1] = param_init$beta_tilde
  alpha[1] <- param_init$alpha
  Gamma2[1] <- param_init$Gamma2
  sigma2[1] <- param_init$sigma2

  ## Initialisation of sufficient statistics useful for updating estimates
  s1 = rep(0,niter+1)                    #stat associated with y_ij
  s2 = rep(0,niter+1)                    #stat associated with the sum of phi^2
  s3 = matrix(0,nrow=n,ncol=niter+1)     #stats associated with phi_i
  mco=matrix(0,n,J)                      #will be useful for the s1 update

  ## Definition of the step sequence for stochastic approximation
  gamma = rep(1,niter)
  for (k in (nburnin+1):niter){
    gamma[k] = 1/(k-nburnin)^(2/3)
  }

  ## Initialisation of individual parameters phi and fixed effects psi
  phi=matrix(NA,nrow=n,ncol=niter+1) #matrix containing in column the (phi_i)_i simulated at each iteration k
  phi[,1] <- V_tilde%*%beta_tilde[,1]

  for (k in 1:niter) {
    #print(k)

    ## S-step: Metropolis Hastings within Gibbs
    phi_tilde=matrix(NA,nrow=n,ncol=niterMH_phi+1) #matrix with in column r the vector (phi_i^(k,r))_i created by M-H
    phi_tilde[,1]=phi[,k]
    Vbetak=V_tilde%*%beta_tilde[,k]

    for (i in 1:n) {
      for (r in 1:niterMH_phi){
        phi_ci=rnorm(1,mean=Vbetak[i],sd=sqrt(Gamma2[k]))
        logratio= sum(dnorm(yi[,i],mean=g(phi_ci,psi,ti[,i]),sd=sqrt(sigma2[k]),log=TRUE))-
          sum(dnorm(yi[,i],mean=g(phi_tilde[i,r],psi,ti[,i]),sd=sqrt(sigma2[k]),log=TRUE))
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


    ## SA-step: updating of sufficient statistics by stochastic approximation
    for (i in 1:n){
      mco[i,] = (yi[,i]-g(phi[i,k+1],psi,ti[,i]))^2
    }

    s1[k+1]=s1[k]+gamma[k]*(sum(mco)-s1[k])
    s2[k+1]=s2[k]+gamma[k]*(sum(phi[,k+1]^2)-s2[k])
    s3[,k+1]=s3[,k]+gamma[k]*(phi[,k+1]-s3[,k])
    p_star_k=p_star(beta_tilde[-1,k],alpha[k],nu0,nu1,p)
    d_tilde_star=c(1/sigma2_mu,(1-p_star_k)/nu0+p_star_k/nu1)

    ## M-step: update of (beta_tilde,alpha,Gamma2,sigma2)
    D_etoile=diag(d_tilde_star)
    beta_tilde[,k+1]=solve(Gamma2[k]*D_etoile+tV_tildeV_tilde,tV_tilde%*%s3[,k+1])

    alpha[k+1]=(sum(p_star_k)+a-1)/(p+b+a-2)

    Vbeta=V_tilde%*%beta_tilde[,k+1]
    Gamma2[k+1]=max(tau*Gamma2[k],(sum((Vbeta)^2)+nu_Gamma*lb_Gamma+s2[k+1]-2*sum((Vbeta)*s3[,k+1]))/(n+nu_Gamma+2))

    sigma2[k+1]=max(tau*sigma2[k],(nu_sigma*lb_sigma+s1[k+1])/(n*J+nu_sigma+2))

  }

  return(list(beta_tilde=beta_tilde,alpha=alpha,Gamma2=Gamma2,sigma2=sigma2))
}

SAEM_EMV <- function(niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init,hyperparam,I,s){
  #Input notations identical to SAEM_MAP except:
  #hyperparam = list(psi,tau) where tau=annealing parameter
  #I = set of indices l such that beta_tilde_l !=0

  set.seed(s)
  n = length(unique(id))
  J = length(Y[id==1])
  yi <- matrix(Y,J,n)
  ti <- matrix(t,J,n)

  Ntot = n*J
  p = dim(V_tilde)[2]-1

  # Hyperparameter
  psi=hyperparam$psi
  tau=hyperparam$tau

  # Initialiqation of the estimation vectors over the iterations
  beta_tilde = matrix(0,nrow=p+1,ncol=niter+1)
  Gamma2 = rep(NA,niter+1)
  sigma2 = rep(NA,niter+1)

  beta_tilde[,1] <- param_init$beta_tilde
  Gamma2[1] <- param_init$Gamma2
  sigma2[1] <- param_init$sigma2

  # Initialisation of sufficient statistics useful for updating estimates
  s1 = rep(0,niter+1)                    #stat associated with y_ij
  s2 = rep(0,niter+1)                    #stat associated with the sum of phi^2
  s3 = matrix(0,nrow=n,ncol=niter+1)     #stats associated with phi_i
  mco=matrix(0,n,J)                      #will be useful for the s1 update

  # Definition of the step sequence for stochastic approximation
  gamma = rep(1,niter)
  for (k in (nburnin+1):niter){
    gamma[k] = 1/(k-nburnin)^(2/3)
  }

  ## Initialisation of individual parameters phi and fixed effects psi
  phi=matrix(NA,nrow=n,ncol=niter+1) #matrix containing in column the (phi_i)_i simulated at each iteration k
  phi[,1] <- V_tilde%*%beta_tilde[,1]

  for (k in 1:niter) {
    #print(k)

    ## S-step: Metropolis Hastings within Gibbs
    phi_tilde=matrix(NA,nrow=n,ncol=niterMH_phi+1) #matrix with in column r the vector (phi_i^(k,r))_i created by M-H
    phi_tilde[,1]=phi[,k]
    Vbetak=V_tilde%*%beta_tilde[,k]

    for (i in 1:n) {
      for (r in 1:niterMH_phi){
        phi_ci=rnorm(1,mean=Vbetak[i],sd=sqrt(Gamma2[k]))
        logratio= sum(dnorm(yi[,i],mean=g(phi_ci,psi,ti[,i]),sd=sqrt(sigma2[k]),log=TRUE))-
          sum(dnorm(yi[,i],mean=g(phi_tilde[i,r],psi,ti[,i]),sd=sqrt(sigma2[k]),log=TRUE))
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

    # SA-step: updating of sufficient statistics by stochastic approximation
    for (i in 1:n){
      mco[i,] = (yi[,i]-g(phi[i,k+1],psi,ti[,i]))^2
    }
    s1[k+1]=s1[k]+gamma[k]*(sum(mco)-s1[k])
    s2[k+1]=s2[k]+gamma[k]*(sum(phi[,k+1]^2)-s2[k])
    s3[,k+1]=s3[,k]+gamma[k]*(phi[,k+1]-s3[,k])

    ## M-step: update of (beta_tilde,Gamma2,sigma2)
    V_hat=V_tilde[,I]
    tV_hat=t(V_hat)
    betatilde=solve(tV_hat%*%V_hat,tV_hat%*%s3[,k+1])
    beta_tilde[I,k+1]=betatilde

    Vbeta=V_tilde%*%beta_tilde[,k+1]
    Gamma2[k+1]=max(tau*Gamma2[k],(sum((Vbeta)^2)+s2[k+1]-2*sum((Vbeta)*s3[,k+1]))/n)

    sigma2[k+1]=max(tau*sigma2[k],(s1[k+1])/(n*J))

  }

  return(list(beta_tildeEMV=beta_tilde[,niter+1],Gamma2EMV=Gamma2[niter+1],sigma2EMV=sigma2[niter+1]))
}

Model_selection <- function(Delta,niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init,hyperparam,s){
  #Input notations identical to SAEM_MAP except:
  #Delta= grid of nu0 values

  p=dim(V_tilde)[2]-1
  M=length(Delta)

  nu1=hyperparam$nu1
  psi=hyperparam$psi
  tau=hyperparam$tau

  n = length(unique(id))
  J = length(Y[id==1])
  yi <- matrix(Y,J,n)
  ti <- matrix(t,J,n)
  Ntot = n*J

  ncore = 10
  cl = makeCluster(ncore)
  registerDoParallel(cl)

  threshold_support<-foreach(m = 1:M, .packages = c("nlme","mvnfast","doParallel","ggplot2","cowplot","glmnet","microbenchmark")) %dopar% {

    source('R/Functions_SAEMVS.R')

    ## Calculation of the threshold
    hyperparam$nu0=Delta[m]
    res = SAEM_MAP(niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init,hyperparam,s)
    beta_tildehat=res$beta[,niter+1]
    Gamma2hat=res$Gamma2[niter+1]
    alphahat=res$alpha[niter+1]
    sigma2hat=res$sigma2[niter+1]
    nu0=Delta[m]
    threshold=sqrt(2*nu0*nu1*log(sqrt(nu1/nu0)*(1-alphahat)/alphahat)/(nu1-nu0))

    norm1 <- apply(matrix(beta_tildehat[-1]),1,function(x) dnorm(x,0,sqrt(nu1)))
    norm0 <- apply(matrix(beta_tildehat[-1]),1,function(x) dnorm(x,0,sqrt(nu0)))
    P_inclu=norm1*alphahat/(norm1*alphahat+norm0*(1-alphahat))

    ## Calculation of the support
    support=c(1,(abs(beta_tildehat[-1])>=threshold)) #indicators beta_l!=0

    result= list(threshold=threshold,support=support,beta_tildehat=beta_tildehat,Gamma2hat=Gamma2hat,alphahat=alphahat,sigma2hat=sigma2hat,P_inclu=P_inclu)
    result
  }
  stopCluster(cl)

  save(threshold_support,file="threshold_support.Rdata")

  beta_tildehat=matrix(NA,nrow=p+1,ncol=M) #each column is the beta_tilde MAP estimate for the k-th value of nu0
  Gamma2hat=rep(NA,M)                      #Gamma^2 MAP estimate for each value of nu0
  sigma2hat=rep(NA,M)                      #sigma^2 MAP estimate for each value of nu0
  alphahat=rep(NA,M)                       #alpha MAP estimate for each value of nu0
  threshold=rep(NA,M)                      #threshold value for each value of nu0
  eBIC = rep(Inf,M)                        #eBIC value for each value of nu0
  support=matrix(NA,nrow=p+1,ncol=M)       #selected support for each value of nu0 (including intercept)
  P_inclu=matrix(NA,nrow=p,ncol=M)         #A posteriori inclusion probability of covariates knowing Theta_hat for each value of nu0
  loglike=rep(NA,M)                        #log-likelihood value for each model

  load("threshold_support.Rdata")
  for (m in 1:M){
    threshold[m]=threshold_support[[m]]$threshold
    support[,m]=threshold_support[[m]]$support
    beta_tildehat[,m]=threshold_support[[m]]$beta_tildehat
    Gamma2hat[m]=threshold_support[[m]]$Gamma2hat
    sigma2hat[m]=threshold_support[[m]]$sigma2hat
    alphahat[m]=threshold_support[[m]]$alphahat
    P_inclu[,m]=threshold_support[[m]]$P_inclu
  }

  ### Computation of eBIC
  #We put in unique_support the distinct supports found along the grid
  unique_support=matrix(support[,1],nrow=p+1,ncol=1)
  for (m in 2:M){
    m2=1
    while (m2<=dim(unique_support)[2]){
      if (identical(support[,m],unique_support[,m2])==TRUE){
        m2=Inf
      }
      else {
        m2=m2+1
      }
    }
    if (m2!=Inf){
      unique_support=cbind(unique_support,support[,m])
    }
  }
  l=dim(unique_support)[2]
  unique_eBIC=rep(NA,l)     #vector containing the eBIC associated with these unique supports
  loglike_unique=rep(NA,l)  #idem for the log-likelihood
  for (ll in 1:l){
    I=which(unique_support[,ll]==1)
    d=length(I)
    if (d!=0){
      res=SAEM_EMV(niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init = param_init,hyperparam=list(tau=0.98,psi=psi),I=I,s=s)
      beta_tildeEMV=res$beta_tildeEMV
      Gamma2EMV=res$Gamma2EMV
      sigma2EMV=res$sigma2EMV
      loglike_unique[ll]=0
      TT=5000
      for (i in 1:n){
        int=rep(0,TT+1)
        phi_i=rmvn(1,mu=rep(sum(beta_tildeEMV*V_tilde[i,]),TT),sigma=diag(Gamma2EMV,TT))
        for (tt in 1:TT){
          mco= (yi[,i]-g(phi_i[tt],psi,ti[,i]))^2
          int[tt+1]=int[tt]+exp(-sum(mco)/(2*sigma2EMV))
        }
        loglike_unique[ll]=loglike_unique[ll] + log((2*pi*sigma2EMV)^(-J/2)*1/(TT)*int[TT+1])
      }
    }
    if (d==0){
      loglike_unique[ll]=-Inf
    }
    unique_eBIC[ll]=-2*loglike_unique[ll]+(d-1)*log(n)+2*log(choose(p,d-1))
  }
  for (m in 1:M){
    ll=1
    while (identical(support[,m],unique_support[,ll])==FALSE){
      ll=ll+1
    }
    eBIC[m]=unique_eBIC[ll]
    loglike[m]=loglike_unique[ll]
  }

  Id=rep(c(1:(p+2)),M)
  y=rbind(beta_tildehat[-1,],threshold)
  y=rbind(y,-threshold)
  y=c(y)
  x=rep(Delta,each=p+2)
  data=data.frame(Id,y,x)
  g1=ggplot(data, aes(x=log(x), y=y, group=Id,color=as.factor(Id))) + geom_point() + geom_line()+
    scale_color_manual(values = c(rep("blue",3),rep("black",p-3),rep("red",2))) +
    theme_bw() + xlab(expression(paste("log(",nu[0]," ) "))) +
    ylab(expression(hat(beta)))+ theme(legend.position = "none") +ggtitle("SAEMVS Regularisation plot")

  data2=data.frame(Delta,eBIC)
  g2=ggplot(data2,aes(x=log(Delta),y=eBIC))+geom_point()+ theme_bw() + xlab(expression(paste("log(",nu[0]," ) "))) +
    ylab("eBIC") +ggtitle("eBIC criterion")
  graph=plot_grid(g1, g2, labels=c("A","B"), ncol = 2, nrow = 1)

  indmin=which.min(eBIC)
  nu0_select=Delta[indmin]
  model_select=which(abs(beta_tildehat[-1,indmin])>=threshold[indmin])
  beta_tildehat_select=c(beta_tildehat[1,indmin],rep(0,p))
  beta_tildehat_select[model_select+1]=beta_tildehat[model_select+1,indmin]
  Gamma2hat_select=Gamma2hat[indmin]
  sigma2hat_select=sigma2hat[indmin]
  return(list(graph=graph,model_select=model_select,beta_tildehat_select=beta_tildehat_select,Gamma2hat_select=Gamma2hat_select,sigma2hat_select=sigma2hat_select,nu0_select=nu0_select,beta_tildehat=beta_tildehat,alphahat=alphahat,sigma2hat=sigma2hat,Gamma2hat=Gamma2hat,P_inclu=P_inclu))
}


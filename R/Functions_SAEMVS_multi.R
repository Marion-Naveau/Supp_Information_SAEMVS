## This file contains all the functions needed to run the SAEMVS procedure with multi random effects

g <- function(dose,t_ij,phi_i){
  return((dose*phi_i[1]/(30*phi_i[1]-phi_i[2]))*(exp(-phi_i[2]*t_ij/30)-exp(-phi_i[1]*t_ij)))
}

p_star <- function(beta,alpha,nu0,nu1){
  norm1 <- apply(beta,1,function(x) dnorm(x,0,sqrt(nu1)))
  norm0 <- apply(beta,1,function(x) dnorm(x,0,sqrt(nu0)))
  p_star <- matrix(norm1*alpha/(norm1*alpha+norm0*(1-alpha)),nrow=dim(beta)[1],ncol=dim(beta)[2],byrow=TRUE)
  return(p_star)
}

SAEM_MAP <- function(niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init,hyperparam,s) {
  #niter = number of iterations of the algorithm
  #nburnin = number of burn-in iterations (nburnin < niter)
  #niterMH_phi = number of iterations of the Metropolis Hastings algorithm for phi at S-step
  #Y = observed data ((y_ij)_{1<=j<=n_i})_{1<=i<=n}
  #t = times ((t_ij)_{1<=j<=n_i})_{1<=i<=n}
  #id = individual identifiers
  #V_tilde = matrix n*(p+1) such as for 1<=i<=n, line i is the vector of Vtilde_i=(1,V_i) in R^(p+1)
  #param_init = initialisation of the parameters to be estimated: list(beta_tilde,alpha,Gamma,sigma2)
  #hyperparam = list of fixed hyperparameters list(dose,nu0,nu1,nu_sigma,lb_sigma,a,b,sigma2_mu,q,Q,d,tau), where tau=annealing parameter
  #s = seed
  
  set.seed(s)
  n = length(unique(id))   #number of individuals
  N = rep(0,n)             #number of observations per individual
  yi=list()
  ti=list()
  for (i in 1:n) {
    N[i]=length(Y[id==i])
    if (i!=1){
      yi[[i]]=Y[(sum(N[1:(i-1)])+1):sum(N[1:i])]
      ti[[i]]=t[(sum(N[1:(i-1)])+1):sum(N[1:i])]
    }
  }
  yi[[1]]=Y[1:N[1]]
  ti[[1]]=t[1:N[1]]
  
  Ntot = sum(N)               #total number of observations in the sample
  p = dim(V_tilde)[2]-1
  tVtilde_Vtilde=t(V_tilde)%*%V_tilde
  
  # Fixed hyperparameters
  nu0=hyperparam$nu0
  nu1=hyperparam$nu1
  nu_sigma=hyperparam$nu_sigma
  lb_sigma=hyperparam$lb_sigma
  a=hyperparam$a
  b=hyperparam$b
  q=hyperparam$q
  sigma2_mu=hyperparam$sigma2_mu
  Q=hyperparam$Q
  d=hyperparam$d
  tau=hyperparam$tau
  dose=hyperparam$dose
  
  # Initialiqation of the estimation vectors over the iterations
  beta_tilde = array(NA,dim=c(p+1,q,niter+1))
  alpha = matrix(NA,nrow=q,ncol=niter+1)
  Gamma = array(NA,dim=c(q,q,niter+1))
  sigma2 = rep(NA,niter+1)
  
  beta_tilde[,,1] = param_init$beta_tilde
  alpha[,1] <- param_init$alpha
  Gamma[,,1] <- param_init$Gamma
  sigma2[1] <- param_init$sigma2
  
  ## Initialisation of sufficient statistics useful for updating estimates
  s1 = rep(0,niter+1)                                #stat associated with y_ij
  s2 = array(0,dim=c(q,q,niter+1))                   #stat associated with the sum of phi*t(phi)
  s3 = array(0,dim=c(n,q,niter+1))                 #stat associated with the sum of phi_i*t(V_i_tilde)
  
  ## Definition of the step sequence for stochastic approximation
  gamma = rep(1,niter)
  for (k in (nburnin+1):niter){
    gamma[k] = 1/(k-nburnin)^(2/3)
  }
  
  ## Initialisation of individual parameters phi and fixed effects psi
  phi=array(NA,dim=c(n,q,niter+1)) #array containing in each matrix k, the matrix containing in column the (phi_i)_i simulated at each iteration k
  phi[,,1] <- V_tilde%*%beta_tilde[,,1]
  
  for (k in 1:niter) {
    
    ## S-step: Metropolis Hastings within Gibbs
    phi_tilde=array(NA,dim=c(n,q,niterMH_phi+1)) #array containing in each matrix r, the matrix containing the (phi_i^(k,r))_i created by M-H
    phi_tilde[,,1]=phi[,,k]
    Vbetak=V_tilde%*%beta_tilde[,,k]
    
    for (i in 1:n) {
      for (r in 1:niterMH_phi){
        phi_ci=rmvn(1,mu=Vbetak[i,],sigma=Gamma[,,k])
        logratio= sum(dnorm(yi[[i]],mean=g(dose,ti[[i]],phi_ci),sd=sqrt(sigma2[k]),log=TRUE))-
          sum(dnorm(yi[[i]],mean=g(dose,ti[[i]],phi_tilde[i,,r]),sd=sqrt(sigma2[k]),log=TRUE))
        u = runif(1)
        cond = (log(u)<=logratio)
        if (cond) {
          phi_tilde[i,,r+1]=phi_ci
        }
        else {
          phi_tilde[i,,r+1]=phi_tilde[i,,r]
        }
      }
    }
    phi[,,k+1]=phi_tilde[,,niterMH_phi+1]
    
    
    ## SA-step: updating of sufficient statistics by stochastic approximation
    mco=0
    for (i in 1:n){
      mco = mco+sum((yi[[i]]-g(dose,ti[[i]],phi[i,,k+1]))^2)
    }
    
    s1[k+1]=s1[k]+gamma[k]*(mco-s1[k])
    s2[,,k+1]=s2[,,k]+gamma[k]*(t(phi[,,k+1])%*%phi[,,k+1]-s2[,,k])
    s3[,,k+1]=s3[,,k]+gamma[k]*(phi[,,k+1]-s3[,,k])
    p_star_k=p_star(beta_tilde[-1,,k],alpha[,k],nu0,nu1)
    d_tilde_star=rbind(rep(1/sigma2_mu,q),(1-p_star_k)/nu0+p_star_k/nu1)
    
    ## M-step: update of (beta_tilde,alpha,Gamma,sigma2)
    D_etoile=diag(c(d_tilde_star))
    cbeta_tilde=solve(diag(1,q)%x%tVtilde_Vtilde + (Gamma[,,k]%x%diag(1,p+1))%*%D_etoile,c(t(V_tilde)%*%s3[,,k+1]))
    beta_tilde[,,k+1]=matrix(cbeta_tilde,nrow=p+1,ncol=q)
    
    for (m in 1:q){
      alpha[m,k+1]=(sum(p_star_k[,m])+a[m]-1)/(p+b[m]+a[m]-2)
    }
    
    argmaxGamma=(Q+s2[,,k+1]-t(V_tilde%*%beta_tilde[,,k+1])%*%s3[,,k+1]-t(s3[,,k+1])%*%V_tilde%*%beta_tilde[,,k+1]+t(beta_tilde[,,k+1])%*%t(V_tilde)%*%V_tilde%*%beta_tilde[,,k+1])/(n+d+q+1)
    if (tr(t(tau*Gamma[,,k])%*%(tau*Gamma[,,k]))>tr(t(argmaxGamma)%*%argmaxGamma)){
      Gamma[,,k+1]=tau*Gamma[,,k]
    } else {
      Gamma[,,k+1]=argmaxGamma
    }
    
    sigma2[k+1]=max(tau*sigma2[k],(nu_sigma*lb_sigma+s1[k+1])/(Ntot+nu_sigma+2))
    
  }
  
  return(list(beta_tilde=beta_tilde,alpha=alpha,Gamma=Gamma,sigma2=sigma2,phi_hat=phi[,,niter+1]))
}


SAEM_EMV <- function(niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init,hyperparam,I,s){
  #Input notations identical to SAEM_MAP except:
  #hyperparam = list(tau,q,dose) where tau=annealing parameter
  #I = matrix qx(p+1) representing the support of (beta_tilde)^T
  
  set.seed(s)
  n = length(unique(id))   #number of individuals
  N = rep(0,n)             #number of observations per individual
  yi=list()
  ti=list()
  for (i in 1:n) {
    N[i]=length(Y[id==i])
    if (i!=1){
      yi[[i]]=Y[(sum(N[1:(i-1)])+1):sum(N[1:i])]
      ti[[i]]=t[(sum(N[1:(i-1)])+1):sum(N[1:i])]
    }
  }
  yi[[1]]=Y[1:N[1]]
  ti[[1]]=t[1:N[1]]
  
  Ntot = sum(N)               #total number of observations in the sample
  
  p = dim(V_tilde)[2]-1
  
  # Fixed hyperparameters
  tau=hyperparam$tau
  q=hyperparam$q
  dose=hyperparam$dose
  
  S=list()
  s=c()
  for (m in 1:q){
    S[[m]]=which(I[m,]==1)
    s[m]=length(S[[m]])
  }
  cumsum_s=c(0,cumsum(s))
  
  X=list()
  for (i in 1:n){
    X[[i]]=matrix(0,nrow=q,ncol=sum(s))
    for (m in 1:q){
      X[[i]][m,c((cumsum_s[m]+1):cumsum_s[m+1])]=V_tilde[i,S[[m]]]
    }
  }
  
  # Initialiqation of the estimation vectors over the iterations
  beta_tilde = array(0,dim=c(p+1,q,niter+1))
  Gamma = array(NA,dim=c(q,q,niter+1))
  sigma2 = rep(NA,niter+1)
  
  beta_tilde[,,1] = param_init$beta_tilde
  Gamma[,,1] <- param_init$Gamma
  sigma2[1] <- param_init$sigma2
  
  ## Initialisation of sufficient statistics useful for updating estimates
  s1 = rep(0,niter+1)                                #stat associated with y_ij
  s2 = array(0,dim=c(q,q,niter+1))                   #stat associated with the sum of phi*t(phi)
  s3 = array(0,dim=c(n,q,niter+1))                   #stat associated with phi_i
  
  ## Definition of the step sequence for stochastic approximation
  gamma = rep(1,niter)
  for (k in (nburnin+1):niter){
    gamma[k] = 1/(k-nburnin)^(2/3)
  }
  
  ## Initialisation of individual parameters phi and fixed effects psi
  phi=array(NA,dim=c(n,q,niter+1)) #array containing in each matrix k, the matrix containing in column the (phi_i)_i simulated at each iteration k
  phi[,,1] <- V_tilde%*%beta_tilde[,,1]
  
  for (k in 1:niter) {
    print(k)
    
    ## S-step: Metropolis Hastings within Gibbs
    phi_tilde=array(NA,dim=c(n,q,niterMH_phi+1)) #array containing in each matrix r, the matrix containing the (phi_i^(k,r))_i created by M-H
    phi_tilde[,,1]=phi[,,k]
    Vbetak=V_tilde%*%beta_tilde[,,k]
    
    for (i in 1:n) {
      for (r in 1:niterMH_phi){
        phi_ci=rmvn(1,mu=Vbetak[i,],sigma=Gamma[,,k])
        logratio= sum(dnorm(yi[[i]],mean=g(dose,ti[[i]],phi_ci),sd=sqrt(sigma2[k]),log=TRUE))-
          sum(dnorm(yi[[i]],mean=g(dose,ti[[i]],phi_tilde[i,,r]),sd=sqrt(sigma2[k]),log=TRUE))
        u = runif(1)
        cond = (log(u)<=logratio)
        if (cond) {
          phi_tilde[i,,r+1]=phi_ci
        }
        else {
          phi_tilde[i,,r+1]=phi_tilde[i,,r]
        }
      }
    }
    phi[,,k+1]=phi_tilde[,,niterMH_phi+1]
    
    
    ## SA-step: updating of sufficient statistics by stochastic approximation
    tphi_phi=t(phi[,,k+1])%*%phi[,,k+1]
    mco=0
    for (i in 1:n){
      mco=mco + sum((yi[[i]]-g(dose,ti[[i]],phi[i,,k+1]))^2)
    }
    
    s1[k+1]=s1[k]+gamma[k]*(mco-s1[k])
    s2[,,k+1]=s2[,,k]+gamma[k]*(tphi_phi-s2[,,k])
    s3[,,k+1]=s3[,,k]+gamma[k]*(phi[,,k+1]-s3[,,k])
    
    ## M-step: update of (beta_tilde,Gamma,sigma2)
    sum_tXi_Gam_Xi=matrix(0,nrow=sum(s),ncol=sum(s))
    for (i in 1:n){
      sum_tXi_Gam_Xi=sum_tXi_Gam_Xi+ t(X[[i]])%*%solve(Gamma[,,k])%*%X[[i]]
    }
    
    sum_tXi_Gam_s3i=rep(0,sum(s))
    for (i in 1:n){
      sum_tXi_Gam_s3i=sum_tXi_Gam_s3i+ t(X[[i]])%*%solve(Gamma[,,k])%*%s3[i,,k+1]
    }
    
    cbeta_tilde_S=solve(sum_tXi_Gam_Xi,sum_tXi_Gam_s3i)
    cbeta_tilde=rep(0,q*(p+1))
    suppI=which(c(t(I))==1)
    cbeta_tilde[suppI]=cbeta_tilde_S
    beta_tilde[,,k+1]=matrix(cbeta_tilde,nrow=p+1,ncol=q,byrow=FALSE)
    
    sum_beta_Xi=matrix(0,nrow=q,ncol=q)
    for (i in 1:n){
      sum_beta_Xi=sum_beta_Xi+X[[i]]%*%cbeta_tilde_S%*%t(cbeta_tilde_S)%*%t(X[[i]])
    }
    sum_Xi_beta_s3i=matrix(0,nrow=q,ncol=q)
    for (i in 1:n){
      sum_Xi_beta_s3i=sum_Xi_beta_s3i+X[[i]]%*%cbeta_tilde_S%*%t(s3[i,,k+1])
    }
    sum_s3i_beta_Xi=t(sum_Xi_beta_s3i)
    argmaxGamma=(s2[,,k+1]+ sum_beta_Xi - sum_Xi_beta_s3i - sum_s3i_beta_Xi)/(n)
    if (tr(t(tau*Gamma[,,k])%*%(tau*Gamma[,,k]))>tr(t(argmaxGamma)%*%argmaxGamma)){
      Gamma[,,k+1]=tau*Gamma[,,k]
    } else {
      Gamma[,,k+1]=argmaxGamma
    }
    
    sigma2[k+1]=max(tau*sigma2[k],(s1[k+1])/(Ntot))
    
  }
  
  return(list(beta_tildeEMV=beta_tilde[,,niter+1],GammaEMV=Gamma[,,niter+1],sigma2EMV=sigma2[niter+1]))
}


Model_selection <- function(Delta,niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init,hyperparam,s){
  #Input notations identical to SAEM_MAP except:
  #Delta= grid of nu0 values
  
  p=dim(V_tilde)[2]-1
  M=length(Delta)
  
  nu1=hyperparam$nu1
  tau=hyperparam$tau
  q=hyperparam$q
  dose=hyperparam$dose
  
  n = length(unique(id))
  N = rep(0,n)             #number of observations per individual
  yi=list()
  ti=list()
  for (i in 1:n) {
    N[i]=length(Y[id==i])
    if (i!=1){
      yi[[i]]=Y[(sum(N[1:(i-1)])+1):sum(N[1:i])]
      ti[[i]]=t[(sum(N[1:(i-1)])+1):sum(N[1:i])]
    }
  }
  yi[[1]]=Y[1:N[1]]
  ti[[1]]=t[1:N[1]]
  
  Ntot = sum(N)               #total number of observations in the sample
  
  
  ncore = 10
  cl = makeCluster(ncore)
  registerDoParallel(cl)

  threshold_support<-foreach(u = 1:M, .packages = c("nlme","mvnfast","doParallel","ggplot2","cowplot","glmnet","microbenchmark","psych")) %dopar% {

    source('R/Functions_SAEMVS_multi.R')
    
    ## Calculation of the threshold
    hyperparam$nu0=Delta[u]
    res = SAEM_MAP(niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init,hyperparam,s)
    beta_tildehat=res$beta_tilde[,,niter+1]
    Gammahat=res$Gamma[,,niter+1]
    alphahat=res$alpha[,niter+1]
    sigma2hat=res$sigma2[niter+1]
    nu0=Delta[u]
    threshold=sqrt(2*nu0*nu1*log(sqrt(nu1/nu0)*(1-alphahat)/alphahat)/(nu1-nu0))

    ## Calculation of the support
    support=matrix(NA,nrow=q,ncol=p+1)
    for (m in 1:q){
      support[m,]=c(1,(abs(beta_tildehat[-1,m])>=threshold[m])) 
    }


    result= list(threshold=threshold,support=support,beta_tildehat=beta_tildehat,Gammahat=Gammahat,alphahat=alphahat,sigma2hat=sigma2hat)
    result
  }
  stopCluster(cl)

  beta_tildehat=array(NA,dim=c(p+1,q,M))  #each matrix (p+1)xq is the beta_tilde MAP estimate for the k-th value of nu0
  Gammahat=array(NA,dim=c(q,q,M))         #Gamma MAP estimate for each value of nu0
  sigma2hat=rep(NA,M)                      #sigma^2 MAP estimate for each value of nu0
  alphahat=matrix(NA,nrow=q,ncol=M)        #alpha MAP estimate for each value of nu0
  threshold=matrix(NA,nrow=q,ncol=M)       #threshold value for each value of nu0
  eBIC = rep(Inf,M)                        #eBIC value for each value of nu0
  support=array(NA,dim=c(q,p+1,M))        #selected support for each value of nu0 (including intercept)
  
  for (m in 1:M){
    threshold[,m]=threshold_support[[m]]$threshold
    support[,,m]=threshold_support[[m]]$support
    beta_tildehat[,,m]=threshold_support[[m]]$beta_tildehat
    Gammahat[,,m]=threshold_support[[m]]$Gammahat
    sigma2hat[m]=threshold_support[[m]]$sigma2hat
    alphahat[,m]=threshold_support[[m]]$alphahat
  }
  
  ### Computation of eBIC
  #We put in unique_support the distinct supports found along the grid
  unique_support=list(support[,,1])
  for (m in 2:M){
    m2=1
    while (m2<=length(unique_support)){
      if (identical(support[,,m],unique_support[[m2]])==TRUE){
        m2=Inf
      }
      else {
        m2=m2+1
      }
    }
    if (m2!=Inf){
      unique_support=c(unique_support,list(support[,,m]))
    }
  }
  l=length(unique_support)
  unique_eBIC=rep(NA,l)     #vector containing the eBIC associated with these unique supports
  loglike_unique=rep(NA,l)  #idem for the log-likelihood
  for (ll in 1:l){
    I=unique_support[[ll]]
    res=SAEM_EMV(niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init = param_init,hyperparam,I=I,s=s)
    beta_tildeEMV=res$beta_tildeEMV
    GammaEMV=res$GammaEMV
    sigma2EMV=res$sigma2EMV
    loglike_unique[ll]=0
    TT=5000
    for (i in 1:n){
      int=rep(0,TT+1)
      for (tt in 1:TT){
        phi_i=rmvn(1,mu=t(beta_tildeEMV)%*%V_tilde[i,],sigma=GammaEMV)
        mco= (yi[[i]]-g(dose,ti[[i]],phi_i))^2
        int[tt+1]=int[tt]+exp(-sum(mco)/(2*sigma2EMV))
      }
      loglike_unique[ll]=loglike_unique[ll] + log((2*pi*sigma2EMV)^(-N[i]/2)*1/(TT)*int[TT+1])
    }
    r=sum(I!=0)
    unique_eBIC[ll]=-2*loglike_unique[ll]+(r-q)*log(n)+2*log(choose(p,r-q))
  }
  for (m in 1:M){
    ll=1
    while (identical(support[,,m],unique_support[[ll]])==FALSE){
      ll=ll+1
    }
    eBIC[m]=unique_eBIC[ll]
    loglike[m]=loglike_unique[ll]
  }
  
  Id=rep(c(1:(p+2)),M)
  g=list()
  for (m in 1:q){
    y=rbind(beta_tildehat[-1,m,],threshold[m,])
    y=rbind(y,-threshold[m,])
    y=c(y)
    x=rep(Delta,each=p+2)
    data=data.frame(Id,y,x)
    g[[m]]=ggplot(data, aes(x=log(x), y=y, group=Id,color=as.factor(Id))) + geom_point() + geom_line()+
      scale_color_manual(values = c(rep("blue",5),rep("black",p-5),rep("red",2))) +
      theme_bw() + xlab(expression(paste("log(",nu[0]," ) "))) +
      ylab(expression(hat(beta)))+ theme(legend.position = "none") +ggtitle("SAEMVS Regularisation plot")
    
  }
  data2=data.frame(Delta,eBIC)
  g2=ggplot(data2,aes(x=log(Delta),y=eBIC))+geom_point()+ theme_bw() + xlab(expression(paste("log(",nu[0]," ) "))) +
    ylab("eBIC") +ggtitle("eBIC criterion")
  
  indmin=which.min(eBIC)
  nu0_select=Delta[indmin]
  model_select1=which(abs(beta_tildehat[-1,1,indmin])>=threshold[indmin])
  model_select2=which(abs(beta_tildehat[-1,2,indmin])>=threshold[indmin])
  beta_tildehat_select=matrix(0,nrow=p+1,ncol=q)
  beta_tildehat_select[1,]=beta_tildehat[1,,indmin]
  beta_tildehat_select[model_select1+1,1]=beta_tildehat[model_select1+1,1,indmin]
  beta_tildehat_select[model_select2+1,2]=beta_tildehat[model_select2+1,2,indmin]
  Gammahat_select=Gammahat[,,indmin]
  sigma2hat_select=sigma2hat[indmin]
  return(list(graph=g,eBIC=g2,nu0_select=nu0_select,model_select1=model_select1,model_select2=model_select2,beta_tildehat_select=beta_tildehat_select,Gammahat_select=Gammahat_select,sigma2hat_select=sigma2hat_select,beta_tildehat=beta_tildehat,alphahat=alphahat,sigma2hat=sigma2hat,Gammahat=Gammahat))
}


## >>>>>>>>>>>> Script for Figure 1 with examples <<<<<<<<<<<

######################### Example to obtain one bar with SAEMVS from Figure 1 #########################

rm(list=ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(doParallel)

source('R/Functions_SAEMVS_multi.R')

S=1                           #number of simulated data-sets tested 
# In the paper S=100, but for reasons of computation time,
# in this example we put S=10 to have an execution time of about 10 min

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

M=10
Delta=10^(seq(-3,0,length.out = M)) #grid of nu0 values

niter=300         #number of iterations of the algorithm
nburnin=150       #number of burn-in iterations
niterMH_phi=1     #number of iterations of the Metropolis Hastings algorithm for phi at S-step

#We apply the SAEMVS method on the S data sets in parallel using ncore cores:
ncore = 10
cl = makeCluster(ncore)
registerDoParallel(cl)

resTot<-foreach(s = 1:S, .packages = c("nlme","mvnfast","doParallel","ggplot2","cowplot","glmnet","psych")) %dopar% {
  ############# Model simulation #############
  set.seed(s)
  
  # The saved data-sets are retrieved:
  load("Saves/data_comp.Rdata")      #list of the 100 data-frame containing Id, Y and t (for all data-sets)
  load("Saves/V_simul.Rdata")        #list of the 100 matrices V (for all data-sets)
  V=V_simul[[s]]
  V_tilde=cbind(rep(1,n),V)
  data=data[[s]]
  
  p_partial=0                                                     #proportion of partially observed individuals in {0,0.1,0.2,0.3,0.4}
  n_partial=p_partial*n                                           #number of partially observed individuals
  Ji=c(rep(3,n_partial),rep(12,n-n_partial))                      #number of repetitions
  data=data[(data$Id%in%c(1:n_partial) & data$t%in%c(0.05,0.15,0.25)) | data$Id%in%c((n_partial+1):n),]
  Id=data$Id
  Y=data$Y
  t=data$t
  id=as.matrix(Id)
  
  # Initial values for the parameters to be estimated
  betatilde_init=matrix(c(10,rep(1,10),rep(0.1,p-10),10,rep(1,10),rep(0.1,p-10)),nrow=p+1,ncol=q,byrow=FALSE)
  sigma2_init=10^(-2)
  Gamma_init=matrix(c(0.5,0.1,0.1,0.5),nrow=q,ncol=q)
  alpha_init=c(0.5,0.5)
  
  # Fixed hyperparameters
  nu0=0.01                         #spike parameter
  nu1=1000                         #slab parameter
  nu_sigma=1                       #prior parameter of sigma2
  lb_sigma=1                       #prior parameter of sigma2
  a=c(1,1)                         #prior parameter of alpha
  b=c(p,p)                         #prior parameter of alpha
  sigma2_mu=5^2                    #prior parameter of mu
  Sigma_Gamma=0.2*diag(1,q)        #prior parameter of Gamma
  d=q+2                            #prior parameter of Gamma
  
  tau=0.98               #annealing parameter
  
  #Initialisation
  param_init=list(beta_tilde=betatilde_init,alpha=alpha_init,Gamma=Gamma_init,sigma2=sigma2_init)
  hyperparam=list(dose=dose,nu0=nu0,nu1=nu1,nu_sigma=nu_sigma,lb_sigma=lb_sigma,a=a,b=b,sigma2_mu=sigma2_mu,q=q,Sigma_Gamma=Sigma_Gamma,d=d,tau=tau)
  
  ############# SAEMVS #############
  res=Model_selection(Delta,niter,nburnin,niterMH_phi,Y,t,id,V_tilde,param_init,hyperparam,s=s)
  nu0_select=res$nu0_select                           #nu0 selected by the eBIC criterion
  beta_tildehat_select=res$beta_tildehat_select       #beta_tilde MAP estimate associated with this nu0
  Gammahat_select=res$Gammahat_select                 #Gamma MAP estimate associated with this nu0
  sigma2hat_select=res$sigma2hat_select               #sigma^2 MAP estimate associated with this nu0
  
  result=list(nu0_select=nu0_select,beta_tildehat_select=beta_tildehat_select,Gammahat_select=Gammahat_select,sigma2hat_select=sigma2hat_select)
  
  result
}
stopCluster(cl)

nu0_select=rep(NA,S)
beta_tildehat_select=array(NA,dim=c(p+1,q,S))
support_select1=list()                          #for the first individual parameter
support_select2=list()                          #for the second individual parameter

for (s in 1:S){
  nu0_select[s]=resTot[[s]]$nu0_select
  beta_tildehat_select[,,s]=resTot[[s]]$beta_tildehat_select
  support_select1[[s]]=resTot[[s]]$model_select1
  support_select2[[s]]=resTot[[s]]$model_select2
}

beta_select=beta_tildehat_select[-1,,]

## On the S datasets, we look at the number of data-sets on which SAEMVS selects the correct model for the first and 
## the second individual parameter: nb_exact_model1, nb_exact_model2
## or a model that contains the correct model for 
## the first and the second individual parameter: nb_contain1, nb_contain2

model_size1=colSums(beta_select[,1,]!=0)
table(model_size1) 
mean(model_size1)

model_size2=colSums(beta_select[,2,]!=0)
table(model_size2) 
mean(model_size2)

v1=which(model_size1==3)
nb_exact_model1=0
s_exact1=c()
for (s in v1){
  nb_exact_model1=nb_exact_model1+(all(beta_select[1:3,1,s]!=rep(0,3)))
  if (all(beta_select[1:3,1,s]!=rep(0,3))){
    s_exact1=c(s_exact1,s)
  }
}
nb_exact_model1

v2=which(model_size2==3)
nb_exact_model2=0
s_exact2=c()
for (s in v2){
  nb_exact_model2=nb_exact_model2+(all(beta_select[3:5,2,s]!=rep(0,3)))
  if (all(beta_select[3:5,2,s]!=rep(0,3))){
    s_exact2=c(s_exact2,s)
  }
}
nb_exact_model2

nb_contain1=0 
s_contain1=c()
for (s in 1:S){
  nb_contain1=nb_contain1 + (all(beta_select[1:3,1,s]!=rep(0,3)))
  if (all(beta_select[1:3,1,s]!=rep(0,3))){
    s_contain1=c(s_contain1,s)
  }
}
nb_contain1

nb_contain2=0 
s_contain2=c()
for (s in 1:S){
  nb_contain2=nb_contain2 + (all(beta_select[3:5,2,s]!=rep(0,3)))
  if (all(beta_select[3:5,2,s]!=rep(0,3))){
    s_contain2=c(s_contain2,s)
  }
}
nb_contain2


######################### Example to obtain one bar with mgaussian from Figure 1 #########################

rm(list=ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(doParallel)

S=100                           #number of simulated data-sets tested for each combination of parameter

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

result=list()
for (s in 1:S){ 
  
  ############# Model simulation #############
  set.seed(s)
  
  # The saved data-sets are retrieved:
  load("Saves/data_comp.Rdata")      #list of the 100 data-frame containing Id, Y and t (for all data-sets)
  load("Saves/V_simul.Rdata")        #list of the 100 matrices V (for all data-sets)
  V=V_simul[[s]]
  V_tilde=cbind(rep(1,n),V)
  data=data[[s]]
  
  p_partial=0                                                     #proportion of partially observed individuals in {0,0.1,0.2,0.3,0.4}
  n_partial=p_partial*n                                           #number of partially observed individuals
  Ji=c(rep(3,n_partial),rep(12,n-n_partial))                      #number of repetitions
  data=data[(data$Id%in%c(1:n_partial) & data$t%in%c(0.05,0.15,0.25)) | data$Id%in%c((n_partial+1):n),]
  Id=data$Id
  Y=data$Y
  t=data$t
  id=as.matrix(Id)
  
  ######### Estimation of phi ##############
  phi_hat=matrix(NA,nrow=n,ncol=q)
  for (i in 1:n){
    if (i==1){
      L <- function(phi_i,Y,t){
        sum=0
        for (j in 1:Ji[i]){
          sum=sum+(Y[j]-g(phi_i,dose=dose,t_ij=t[j]))^2
        }
        return(sum)
      }
      
      phi_hat[i,]=nlm(L,c(5,5),Y=Y,t=t)$estimate
    } else {
      L <- function(phi_i,Y,t){
        sum=0
        for (j in 1:Ji[i]){
          sum=sum+(Y[sum(Ji[1:(i-1)])+j]-g(phi_i,dose=dose,t_ij=t[sum(Ji[1:(i-1)])+j]))^2
        }
        return(sum)
      }
      
      phi_hat[i,]=nlm(L,c(5,5),Y=Y,t=t)$estimate
    }
  }
  
  phi_hat_c=phi_hat  #phi_hat - intercept
  for (i in 1:n){
    phi_hat_c[i,]=phi_hat[i,]-mu
  }
  
  ############### LASSO #######################
  mfit <- glmnet(V, phi_hat_c, family = "mgaussian",alpha=1)
  
  reg.cvlasso <- cv.glmnet(V,phi_hat_c,family="mgaussian",alpha=1)
  bestlam <- reg.cvlasso$lambda.1se
  beta_select <- glmnet(V, phi_hat_c, family = "mgaussian",alpha=1,lambda = bestlam)$beta
  
  result[[s]]=list(bestlam=bestlam, beta_select=beta_select, phi_hat=phi_hat)
  
}

beta_select=array(NA,dim=c(p,q,S))

for (s in 1:S){
  beta_select[,1,s]=as.vector(result[[s]]$beta_select$y1)
  beta_select[,2,s]=as.vector(result[[s]]$beta_select$y2)
}

## On the S datasets, we look at the number of data-sets on which "mgaussian" selects the correct model for the first and 
## the second individual parameter: nb_exact_model1, nb_exact_model2
## or a model that contains the correct model for 
## the first and the second individual parameter: nb_contain1, nb_contain2

model_size1=colSums(beta_select[,1,]!=0)
table(model_size1) 
mean(model_size1)

model_size2=colSums(beta_select[,2,]!=0)
table(model_size2) 
mean(model_size2)

v1=which(model_size1==3)
nb_exact_model1=0
s_exact1=c()
for (s in v1){
  nb_exact_model1=nb_exact_model1+(all(beta_select[1:3,1,s]!=rep(0,3)))
  if (all(beta_select[1:3,1,s]!=rep(0,3))){
    s_exact1=c(s_exact1,s)
  }
}
nb_exact_model1

v2=which(model_size2==3)
nb_exact_model2=0
s_exact2=c()
for (s in v2){
  nb_exact_model2=nb_exact_model2+(all(beta_select[3:5,2,s]!=rep(0,3)))
  if (all(beta_select[3:5,2,s]!=rep(0,3))){
    s_exact2=c(s_exact2,s)
  }
}
nb_exact_model2

nb_contain1=0 
s_contain1=c()
for (s in 1:S){
  nb_contain1=nb_contain1 + (all(beta_select[1:3,1,s]!=rep(0,3)))
  if (all(beta_select[1:3,1,s]!=rep(0,3))){
    s_contain1=c(s_contain1,s)
  }
}
nb_contain1

nb_contain2=0 
s_contain2=c()
for (s in 1:S){
  nb_contain2=nb_contain2 + (all(beta_select[3:5,2,s]!=rep(0,3)))
  if (all(beta_select[3:5,2,s]!=rep(0,3))){
    s_contain2=c(s_contain2,s)
  }
}
nb_contain2



######################### Example to obtain one bar with gaussian from Figure 1 #########################

rm(list=ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(doParallel)

S=100                           #number of simulated data-sets tested 

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

result=list()
for (s in 1:S){ 
  
  ############# Model simulation #############
  set.seed(s)
  
  # The saved data-sets are retrieved:
  load("Saves/data_comp.Rdata")      #list of the 100 data-frame containing Id, Y and t (for all data-sets)
  load("Saves/V_simul.Rdata")        #list of the 100 matrices V (for all data-sets)
  V=V_simul[[s]]
  V_tilde=cbind(rep(1,n),V)
  data=data[[s]]
  
  p_partial=0                                                     #proportion of partially observed individuals in {0,0.1,0.2,0.3,0.4}
  n_partial=p_partial*n                                           #number of partially observed individuals
  Ji=c(rep(3,n_partial),rep(12,n-n_partial))                      #number of repetitions
  data=data[(data$Id%in%c(1:n_partial) & data$t%in%c(0.05,0.15,0.25)) | data$Id%in%c((n_partial+1):n),]
  Id=data$Id
  Y=data$Y
  t=data$t
  id=as.matrix(Id)
  
  ######### Estimation of phi ##############
  phi_hat=matrix(NA,nrow=n,ncol=q)
  for (i in 1:n){
    if (i==1){
      L <- function(phi_i,Y,t){
        sum=0
        for (j in 1:Ji[i]){
          sum=sum+(Y[j]-g(phi_i,dose=dose,t_ij=t[j]))^2
        }
        return(sum)
      }
      
      phi_hat[i,]=nlm(L,c(5,5),Y=Y,t=t)$estimate
    } else {
      L <- function(phi_i,Y,t){
        sum=0
        for (j in 1:Ji[i]){
          sum=sum+(Y[sum(Ji[1:(i-1)])+j]-g(phi_i,dose=dose,t_ij=t[sum(Ji[1:(i-1)])+j]))^2
        }
        return(sum)
      }
      
      phi_hat[i,]=nlm(L,c(5,5),Y=Y,t=t)$estimate
    }
  }
  
  phi_hat_c=phi_hat  #phi_hat - intercept
  for (i in 1:n){
    phi_hat_c[i,]=phi_hat[i,]-mu
  }
  
  ############### LASSO #######################
  mfit1 <- glmnet(V, phi_hat_c[,1], family = "gaussian",alpha=1)
  
  reg.cvlasso1 <- cv.glmnet(V,phi_hat_c[,1],family="gaussian",alpha=1)
  bestlam1 <- reg.cvlasso1$lambda.1se
  beta_select1 <- glmnet(V, phi_hat_c[,1], family = "gaussian",alpha=1,lambda = bestlam1)$beta
  
  mfit2 <- glmnet(V, phi_hat_c[,2], family = "gaussian",alpha=1)
  
  reg.cvlasso2 <- cv.glmnet(V,phi_hat_c[,2],family="gaussian",alpha=1)
  bestlam2 <- reg.cvlasso2$lambda.1se
  beta_select2 <- glmnet(V, phi_hat_c[,2], family = "gaussian",alpha=1,lambda = bestlam2)$beta
  
  bestlam=c(bestlam1,bestlam2)
  result[[s]]=list(bestlam=bestlam, beta_select1=beta_select1,beta_select2=beta_select2, phi_hat=phi_hat)
  
}

beta_select=array(NA,dim=c(p,q,S))

for (s in 1:S){
  beta_select[,1,s]=as.vector(result[[s]]$beta_select1)
  beta_select[,2,s]=as.vector(result[[s]]$beta_select2)
}

## On the S datasets, we look at the number of data-sets on which "gaussian" selects the correct model for the first and 
## the second individual parameter: nb_exact_model1, nb_exact_model2
## or a model that contains the correct model for 
## the first and the second individual parameter: nb_contain1, nb_contain2

model_size1=colSums(beta_select[,1,]!=0)
table(model_size1) 
mean(model_size1)

model_size2=colSums(beta_select[,2,]!=0)
table(model_size2) 
mean(model_size2)

v1=which(model_size1==3)
nb_exact_model1=0
s_exact1=c()
for (s in v1){
  nb_exact_model1=nb_exact_model1+(all(beta_select[1:3,1,s]!=rep(0,3)))
  if (all(beta_select[1:3,1,s]!=rep(0,3))){
    s_exact1=c(s_exact1,s)
  }
}
nb_exact_model1

v2=which(model_size2==3)
nb_exact_model2=0
s_exact2=c()
for (s in v2){
  nb_exact_model2=nb_exact_model2+(all(beta_select[3:5,2,s]!=rep(0,3)))
  if (all(beta_select[3:5,2,s]!=rep(0,3))){
    s_exact2=c(s_exact2,s)
  }
}
nb_exact_model2

nb_contain1=0 
s_contain1=c()
for (s in 1:S){
  nb_contain1=nb_contain1 + (all(beta_select[1:3,1,s]!=rep(0,3)))
  if (all(beta_select[1:3,1,s]!=rep(0,3))){
    s_contain1=c(s_contain1,s)
  }
}
nb_contain1

nb_contain2=0 
s_contain2=c()
for (s in 1:S){
  nb_contain2=nb_contain2 + (all(beta_select[3:5,2,s]!=rep(0,3)))
  if (all(beta_select[3:5,2,s]!=rep(0,3))){
    s_contain2=c(s_contain2,s)
  }
}
nb_contain2


######################### Script for Figure 1 #########################

########## Figure 1A: for the first individual parameter ##########

rm(list=ls())

library(ggplot2)
library(ggpattern)
library(cowplot)

# The saved results are retrieved:
file_name_SAEMVS=c("Saves/betahat1_simuSA1.Rdata","Saves/betahat1_simuSA2.Rdata","Saves/betahat1_simuSA3.Rdata",
            "Saves/betahat1_simuSA4.Rdata","Saves/betahat1_simuSA5.Rdata")

file_name_mg=c("Saves/betahat1_simuMG1_mg.Rdata","Saves/betahat1_simuMG2_mg.Rdata","Saves/betahat1_simuMG3_mg.Rdata",
               "Saves/betahat1_simuMG4_mg.Rdata","Saves/betahat1_simuMG5_mg.Rdata")

file_name_g=c("Saves/betahat1_simuMG1_g.Rdata","Saves/betahat1_simuMG2_g.Rdata","Saves/betahat1_simuMG3_g.Rdata",
              "Saves/betahat1_simuMG4_g.Rdata","Saves/betahat1_simuMG5_g.Rdata")

## For each proportion of partially observed individuals, on the S datasets, we look at the number of data-sets on which SAEMVS, "gaussian", and "mgaussian" select 
## the correct model for the first individual parameter: nb_exact_model_SAEMVS, nb_exact_model_mg, nb_exact_model_g
## a model that contains the correct model for the first individual parameter: nb_contain_SAEMVS, nb_contain_mg, nb_contain_g

N_SAEMVS=length(file_name_SAEMVS)
N_mg=length(file_name_mg)
N_g=length(file_name_g)

nb_exact_model_SAEMVS=rep(0,N_SAEMVS)
nb_exact_model_mg=rep(0,N_mg)
nb_exact_model_g=rep(0,N_g)

nb_contain_SAEMVS=rep(0,N_SAEMVS)
nb_contain_mg=rep(0,N_mg)
nb_contain_g=rep(0,N_g)

S=100

for (k in 1:N_SAEMVS){
  load(file_name_SAEMVS[k])   #beta_select in R^(p x S): each column is the value of the MAP estimator of beta for the first individual parameter
  #obtain by SAEMVS for the selected nu0 for one proportion of partially observed individuals, for one data-set among the S
  p=length(beta_select[,1])
  
  Model_size=colSums(beta_select!=0)      #Sizes of model selected for the S datasets
  v=which(Model_size==3)  #Sets of datasets for which SAEMVS has selected a model of size 3
  
  for (s in v){
    nb_exact_model_SAEMVS[k]=nb_exact_model_SAEMVS[k]+(all(beta_select[1:3,s]!=rep(0,3)))
  }
  
  for (s in 1:S){
    nb_contain_SAEMVS[k]=nb_contain_SAEMVS[k] + (all(beta_select[1:3,s]!=rep(0,3)))
  }
}

for (k in 1:N_mg){
  load(file_name_mg[k])   #beta_select in R^(p x S): each column is the value of the MAP estimator of beta for the first individual parameter
  #obtain by "mgaussian" for lambda1se for one proportion of partially observed individuals, for one data-set among the S
  p=length(beta_select[,1])
  
  Model_size=colSums(beta_select!=0)      #Sizes of model selected for the S datasets
  v=which(Model_size==3)  #Sets of datasets for which "mgaussian" has selected a model of size 3
  
  for (s in v){
    nb_exact_model_mg[k]=nb_exact_model_mg[k]+(all(beta_select[1:3,s]!=rep(0,3)))
  }
  
  for (s in 1:S){
    nb_contain_mg[k]=nb_contain_mg[k] + (all(beta_select[1:3,s]!=rep(0,3)))
  }
}

for (k in 1:N_g){
  load(file_name_g[k])   #beta_select in R^(p x S): each column is the value of the MAP estimator of beta for the first individual parameter
  #obtain by "gaussian" for lambda1se for one proportion of partially observed individuals, for one data-set among the S
  p=length(beta_select[,1])
  
  Model_size=colSums(beta_select!=0)      #Sizes of model selected for the S datasets
  v=which(Model_size==3)  #Sets of datasets for which "gaussian" has selected a model of size 3
  
  for (s in v){
    nb_exact_model_g[k]=nb_exact_model_g[k]+(all(beta_select[1:3,s]!=rep(0,3)))
  }
  
  for (s in 1:S){
    nb_contain_g[k]=nb_contain_g[k] + (all(beta_select[1:3,s]!=rep(0,3)))
  }
}

nb_exact_model=c(rbind(nb_exact_model_SAEMVS,nb_exact_model_mg,nb_exact_model_g))
nb_contain=c(rbind(nb_contain_SAEMVS,nb_contain_mg,nb_contain_g))

data_phi1 <- data.frame(p_partial=rep(c(0,10,20,30,40),each=3),Method=rep(c("SAEMVS","mgaussian","gaussian"),5),Result=rep(c("Exact","Strictly included"),each=15),count=c(nb_exact_model, nb_contain))
data_phi1
data_phi1$p_partial<-as.factor(data_phi1$p_partial)
data_phi1$Method<-as.factor(data_phi1$Method)

g1=ggplot(mapping=aes(x=p_partial, y=count, fill=Method,pattern=Result)) +geom_col(data=data_phi1[data_phi1$Result=="Strictly included",], aes(group = Method), width = 0.7, position = position_dodge(width = 0.7),color="black")+ geom_bar_pattern(data=data_phi1[data_phi1$Result=="Exact",], position="dodge",width=0.7,color="black",pattern_fill="black",pattern_density=0.05,pattern_spacing=0.025,stat='identity')+
  labs(x="Percentage of partially observed individuals (in %)", y = "Proportion (in %)") +scale_fill_brewer(palette="Set2")+ scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100))+ scale_pattern_manual(values = c(Exact="stripe","Strictly included"="none"))+
  theme_bw()+theme(axis.text=element_text(size=30), axis.title=element_text(size=30,face="bold"),title=element_text(size=30,face="bold"),legend.position='none')+guides(pattern=guide_legend(override.aes = list(fill="white")),fill=guide_legend(override.aes = list(pattern="none")))
g1

########## Figure 1B: for the second individual parameter ##########

# The saved results are retrieved:
file_name_SAEMVS=c("Saves/betahat2_simuSA1.Rdata","Saves/betahat2_simuSA2.Rdata","Saves/betahat2_simuSA3.Rdata",
                   "Saves/betahat2_simuSA4.Rdata","Saves/betahat2_simuSA5.Rdata")

file_name_mg=c("Saves/betahat2_simuMG1_mg.Rdata","Saves/betahat2_simuMG2_mg.Rdata","Saves/betahat2_simuMG3_mg.Rdata",
               "Saves/betahat2_simuMG4_mg.Rdata","Saves/betahat2_simuMG5_mg.Rdata")

file_name_g=c("Saves/betahat2_simuMG1_g.Rdata","Saves/betahat2_simuMG2_g.Rdata","Saves/betahat2_simuMG3_g.Rdata",
              "Saves/betahat2_simuMG4_g.Rdata","Saves/betahat2_simuMG5_g.Rdata")

## For each proportion of partially observed individuals, on the S datasets, we look at the number of data-sets on which SAEMVS, "gaussian", and "mgaussian" select 
## the correct model for the first individual parameter: nb_exact_model_SAEMVS, nb_exact_model_mg, nb_exact_model_g
## a model that contains the correct model for the first individual parameter: nb_contain_SAEMVS, nb_contain_mg, nb_contain_g

N_SAEMVS=length(file_name_SAEMVS)
N_mg=length(file_name_mg)
N_g=length(file_name_g)

nb_exact_model_SAEMVS=rep(0,N_SAEMVS)
nb_exact_model_mg=rep(0,N_mg)
nb_exact_model_g=rep(0,N_g)

nb_contain_SAEMVS=rep(0,N_SAEMVS)
nb_contain_mg=rep(0,N_mg)
nb_contain_g=rep(0,N_g)

S=100

for (k in 1:N_SAEMVS){
  load(file_name_SAEMVS[k])   #beta_select in R^(p x S): each column is the value of the MAP estimator of beta for the first individual parameter
  #obtain by SAEMVS for the selected nu0 for one proportion of partially observed individuals, for one data-set among the S
  p=length(beta_select[,1])
  
  Model_size=colSums(beta_select!=0)      #Sizes of model selected for the S datasets
  v=which(Model_size==3)  #Sets of datasets for which SAEMVS has selected a model of size 3
  
  for (s in v){
    nb_exact_model_SAEMVS[k]=nb_exact_model_SAEMVS[k]+(all(beta_select[3:5,s]!=rep(0,3)))
  }
  
  for (s in 1:S){
    nb_contain_SAEMVS[k]=nb_contain_SAEMVS[k] + (all(beta_select[3:5,s]!=rep(0,3)))
  }
}

for (k in 1:N_mg){
  load(file_name_mg[k])   #beta_select in R^(p x S): each column is the value of the MAP estimator of beta for the first individual parameter
  #obtain by "mgaussian" for lambda1se for one proportion of partially observed individuals, for one data-set among the S
  p=length(beta_select[,1])
  
  Model_size=colSums(beta_select!=0)      #Sizes of model selected for the S datasets
  v=which(Model_size==3)  #Sets of datasets for which "mgaussian" has selected a model of size 3
  
  for (s in v){
    nb_exact_model_mg[k]=nb_exact_model_mg[k]+(all(beta_select[3:5,s]!=rep(0,3)))
  }
  
  for (s in 1:S){
    nb_contain_mg[k]=nb_contain_mg[k] + (all(beta_select[3:5,s]!=rep(0,3)))
  }
}

for (k in 1:N_g){
  load(file_name_g[k])   #beta_select in R^(p x S): each column is the value of the MAP estimator of beta for the first individual parameter
  #obtain by "gaussian" for lambda1se for one proportion of partially observed individuals, for one data-set among the S
  p=length(beta_select[,1])
  
  Model_size=colSums(beta_select!=0)      #Sizes of model selected for the S datasets
  v=which(Model_size==3)  #Sets of datasets for which "gaussian" has selected a model of size 3
  
  for (s in v){
    nb_exact_model_g[k]=nb_exact_model_g[k]+(all(beta_select[3:5,s]!=rep(0,3)))
  }
  
  for (s in 1:S){
    nb_contain_g[k]=nb_contain_g[k] + (all(beta_select[3:5,s]!=rep(0,3)))
  }
}

nb_exact_model=c(rbind(nb_exact_model_SAEMVS,nb_exact_model_mg,nb_exact_model_g))
nb_contain=c(rbind(nb_contain_SAEMVS,nb_contain_mg,nb_contain_g))

data_phi1 <- data.frame(p_partial=rep(c(0,10,20,30,40),each=3),Method=rep(c("SAEMVS","mgaussian","gaussian"),5),Result=rep(c("Exact","Strictly included"),each=15),count=c(nb_exact_model, nb_contain))
data_phi1
data_phi1$p_partial<-as.factor(data_phi1$p_partial)
data_phi1$Method<-as.factor(data_phi1$Method)

g2=ggplot(mapping=aes(x=p_partial, y=count, fill=Method,pattern=Result)) +geom_col(data=data_phi1[data_phi1$Result=="Strictly included",], aes(group = Method), width = 0.7, position = position_dodge(width = 0.7),color="black")+ geom_bar_pattern(data=data_phi1[data_phi1$Result=="Exact",], position="dodge",width=0.7,color="black",pattern_fill="black",pattern_density=0.05,pattern_spacing=0.025,stat='identity')+
  labs(x="Percentage of partially observed individuals (in %)", y = "Proportion (in %)") +scale_fill_brewer(palette="Set2")+ scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100))+ scale_pattern_manual(values = c(Exact="stripe","Strictly included"="none"))+
  theme_bw()+theme(axis.text=element_text(size=30), axis.title=element_text(size=30,face="bold"),title=element_text(size=30,face="bold"),legend.title = element_text(size=30,face="bold"),legend.text = element_text(size=30))+guides(pattern=guide_legend(override.aes = list(fill="white")),fill=guide_legend(override.aes = list(pattern="none")))
g2

Figure1=plot_grid(g1, g2, labels=c("a","b"), ncol = 2, nrow = 1,rel_widths=c(1.5,2),label_size = 30)
Figure1

ggsave("Fig1.eps",device="eps", width = 30,height = 10)

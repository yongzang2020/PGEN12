###Please note that the code is designed for 3 subgroups. If this is not the case, the functions will need to be modified accordingly.
###The R packages needed
library(mvtnorm);library(msm);library(dlm);library(truncnorm);library(rjags)


#########################################################################################################
##### Some functions utilized by the PGen I-II design ###################################################
#########################################################################################################

uti_benchmark_PGen12=function(score,targetT,targetE){
  ## targetT: highest acceptable toomegacity rate for dose-finding
  ## targetE: lowest acceptable efficacy rate  for dose-finding
  
  # Assume independence between toxicity and efficacy
  targetP<-matrix(c(0,(1-targetT)*(1-targetE),(1-targetT)*targetE,0,targetT*(1-targetE),targetT*targetE),2,3, byrow=TRUE)
  
  # Calculate the benchmark utility
  uu = sum(targetP*score) #highest unacceptable utility
  uu = uu+(100-uu)/2        # benchmark utility (i.e., desirable utility)
  return(uu)
}

get.boundary <- function(target, targetE, ncohort, cohortsize, p.saf=NA, p.tox=NA,  cutoff.eli, cutoff.eli.E){#Calculate the boundaries of BOIN-12 design
  # if the user does not provide p.saf and p.tox, use the default values
  if(is.na(p.saf)) p.saf=0.6*target;
  if(is.na(p.tox)) p.tox=1.4*target;
  
  ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
  npts = ncohort*cohortsize;
  ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;elimE=NULL;
  for(n in (1:ncohort)*cohortsize)
  {
    error.min=3;
    for(m1 in 0:(n-1))
    {
      for(m2 in (m1+1):n)
      {
        
        error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
        error2 = 1-pbinom(m1, n, p.saf);
        error3 = pbinom(m2-1, n, p.tox);
        
        error=error1+error2+error3;
        if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
      }
    }
    ntrt = c(ntrt, n);
    b.e = c(b.e, cutoff1);
    b.d = c(b.d, cutoff2);
    
    elimineed=0; # indicating whether elimination is needed
    elimineedE=0
    if(n<3) { elim = c(elim, NA); elimE = c(elimE,NA)}  # require treating at least 3 patients before eliminating a dose
    else
    {
      for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
      {
        if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
      }
      if(elimineed==1) { elim = c(elim, ntox); }
      else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
      
      for(neff in n:0){
        if(pbeta(targetE,neff+1,n-neff+1)>cutoff.eli.E){elimineedE=1; break;}
      }
      if(elimineedE==1){elimE=c(elimE,neff)} else {elimE=c(elimE,NA)}
    }
  }
  for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
  boundaries = rbind(ntrt, elim, b.d, b.e,elimE);
  rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=", 
                           "Deescalate if # of DLT >=",  "Escalate if # of DLT <=", "Eliminate if # of Eff <=");
  colnames(boundaries) = rep("", ncohort);
  
  return(boundaries);
}

adm=function(y,targetT,targetE,cutoff.eli.T,cutoff.eli.E){#Function to obtain the admissible set for toxicity and efficacy
  y=as.matrix(y)
  n=apply(y,2,sum)
  n.A=ifelse(n>0,1,0)
  yT=y[2,]+y[4,]+y[6,]
  yE=y[5,]+y[6,]
  postT=1-pbeta(targetT,1+yT,1+n-yT)
  T.A=ifelse(postT<cutoff.eli.T,1,0)
  postE=pbeta(targetE,1+yE,1+n-yE)
  E.A=ifelse(postE<cutoff.eli.E,1,0)
  re=n.A*T.A*E.A
  re1=which(re==1)
  return(re1)
}

adm_tox=function(y,targetT,cutoff.eli.T){#Function to obtain the admissible set only for toxicity
  y=as.matrix(y)
  n=apply(y,2,sum)
  n.A=ifelse(n>0,1,0)
  yT=y[2,]+y[4,]+y[6,]
  postT=1-pbeta(targetT,1+yT,1+n-yT)
  T.A=ifelse(postT<cutoff.eli.T,1,0)
  re=n.A*T.A
  re1=which(re==1)
  return(re1)
}

log.likelihood.X_E <- function(X_E,X_T,dose,alpha0,alpha1,alpha2,alpha3,beta0,beta1,sigma12) {
  return(-sum((X_E-alpha0-alpha1*dose^alpha3/(alpha2^alpha3+dose^alpha3)-sigma12*(X_T-beta0-beta1*dose))^2)/(2*(1-sigma12^2)))
}

log.likelihood.sigma12 <- function(X_E,X_T,dose,alpha0,alpha1,alpha2,alpha3,beta0,beta1,sigma12) {
  top.part1 <- sum((X_T-beta0-beta1*dose)^2)
  top.part2 <- sum((X_E-alpha0-alpha1*dose^alpha3/(alpha2^alpha3+dose^alpha3))^2)
  top.part3 <- -2*sigma12*sum((X_T-beta0-beta1*dose)*(X_E-alpha0-alpha1*dose^alpha3/(alpha2^alpha3+dose^alpha3)))
  return(-(length(X_E)/2)*log(1-sigma12^2)-(top.part1+top.part2+top.part3)/(2*(1-sigma12^2)))
}  

log.likelihood.Y_ET <- function(parameters,data,dose.level) {
  alpha0 <- parameters[1]
  alpha1 <- parameters[2]
  alpha2 <- parameters[3]
  alpha3 <- parameters[4]
  beta0 <- parameters[5]
  beta1 <- parameters[6]
  eta2 <- parameters[7]
  sigma12 <- parameters[8]
  
  if(eta2<=0|alpha2<=0|alpha3<=0|beta1<=0) return(-Inf)
  Y_T<-data[[1]]
  Y_E<-data[[2]]
  d<-data[[3]]
  re=array(0,c(2,3,4))
  cor=matrix( c(1,sigma12,sigma12,1),nrow=2     )
  for(k in 1:4){
    re[1,1,k]=max(0,pmvnorm(lower=c(-Inf,-Inf),upper=c(0,0), 
                            mean=c(beta0+beta1*dose.level[k], alpha0+alpha1*dose.level[k]^alpha3/(alpha2^alpha3+dose.level[k]^alpha3)), sigma=cor      )[1] )
    re[2,1,k]=max(0,pmvnorm(lower=c(0,-Inf),upper=c(Inf,0), 
                            mean=c(beta0+beta1*dose.level[k], alpha0+alpha1*dose.level[k]^alpha3/(alpha2^alpha3+dose.level[k]^alpha3)), sigma=cor      )[1] )
    re[1,2,k]=max(0,pmvnorm(lower=c(-Inf,0),upper=c(0,eta2), 
                            mean=c(beta0+beta1*dose.level[k], alpha0+alpha1*dose.level[k]^alpha3/(alpha2^alpha3+dose.level[k]^alpha3)), sigma=cor      )[1] )
    re[2,2,k]=max(0,pmvnorm(lower=c(0,0),upper=c(Inf,eta2), 
                            mean=c(beta0+beta1*dose.level[k], alpha0+alpha1*dose.level[k]^alpha3/(alpha2^alpha3+dose.level[k]^alpha3)), sigma=cor      )[1] )
    re[1,3,k]=max(0,pmvnorm(lower=c(-Inf,eta2),upper=c(0,Inf), 
                            mean=c(beta0+beta1*dose.level[k], alpha0+alpha1*dose.level[k]^alpha3/(alpha2^alpha3+dose.level[k]^alpha3)), sigma=cor      )[1] )
    re[2,3,k]=max(0,pmvnorm(lower=c(0,eta2),upper=c(Inf,Inf), 
                            mean=c(beta0+beta1*dose.level[k], alpha0+alpha1*dose.level[k]^alpha3/(alpha2^alpha3+dose.level[k]^alpha3)), sigma=cor      )[1] )
  }
  ll<-0
  for(i in 1:length(d))
    ll<-ll+log(re[Y_T[i]+1,Y_E[i]+1,d[i]])
  return(ll)
}

negative.log.likelihood.Y_ET <- function(parameters,data,dose.level) {
  return(-log.likelihood.Y_ET(parameters,data,dose.level))
}

log.likelihood.Y_S <- function(parameters,data,dose.level) {
  Y_E<-data[[2]]
  id<-which(Y_E!=0)
  Y_E<-Y_E[id]
  d<-data[[3]][id]
  d.level<-sort(unique(d))
  gamma_T <- parameters[1]
  gamma_E <- parameters[2]
  lambda_1 <- parameters[3]
  lambda_2 <- parameters[4]
  lambda_3 <- parameters[5]
  gamma_d <- c(0,parameters[6:(4+length(d.level))])
  
  if(lambda_1<0|lambda_2<0|lambda_3<0) return(-Inf)
  Y_T<-data[[1]][id]
  Y_S.obs<-data[[4]][id]
  status<-data[[5]][id]
  ll<-0
  T_n<-max(Y_S.obs)
  for(i in 1:length(d)){
    if(status[i]==1){
      if(Y_S.obs[i]<T_n/3) ll<-ll+log(lambda_1*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])*exp(-Y_S.obs[i]*lambda_1*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])))
      if(Y_S.obs[i]<T_n*2/3&Y_S.obs[i]>=T_n/3) ll<-ll+log(lambda_2*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])*exp(-(T_n/3*lambda_1+(Y_S.obs[i]-T_n/3)*lambda_2)*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])))
      if(Y_S.obs[i]<T_n&Y_S.obs[i]>=T_n*2/3) ll<-ll+log(lambda_3*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])*exp(-(T_n/3*lambda_1+T_n/3*lambda_2+(Y_S.obs[i]-T_n*2/3)*lambda_3)*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])))
    } else{
      chazard<-0
      if(Y_S.obs[i]<T_n/3) chazard<-Y_S.obs[i]*lambda_1
      if(Y_S.obs[i]<T_n*2/3&Y_S.obs[i]>=T_n/3) chazard<-T_n/3*lambda_1+(Y_S.obs[i]-T_n/3)*lambda_2
      if(Y_S.obs[i]>=T_n*2/3) chazard<-T_n/3*lambda_1+T_n/3*lambda_2+(Y_S.obs[i]-T_n*2/3)*lambda_3
      ll<-ll-chazard*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])
    }
  }
  return(ll)
}

negative.log.likelihood.Y_S <- function(parameters,data,dose.level) {
  return(-log.likelihood.Y_S(parameters,data,dose.level))
}

log.likelihood.Y_S_withoutconstraint <- function(parameters,data,dose.level) {
  Y_E<-data[[2]]
  id<-which(Y_E!=0)
  Y_E<-Y_E[id]
  d<-data[[3]][id]
  d.level<-sort(unique(d))
  gamma_T <- parameters[1]
  gamma_E <- parameters[2]
  lambda_1 <- parameters[3]
  lambda_2 <- parameters[4]
  lambda_3 <- parameters[5]
  gamma_d <- c(0,parameters[6:(4+length(d.level))])
  
  Y_T<-data[[1]][id]
  Y_S.obs<-data[[4]][id]
  status<-data[[5]][id]
  ll<-0
  T_n<-max(Y_S.obs)
  for(i in 1:length(d)){
    if(status[i]==1){
      if(Y_S.obs[i]<T_n/3) ll<-ll+log(lambda_1*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])*exp(-Y_S.obs[i]*lambda_1*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])))
      if(Y_S.obs[i]<T_n*2/3&Y_S.obs[i]>=T_n/3) ll<-ll+log(lambda_2*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])*exp(-(T_n/3*lambda_1+(Y_S.obs[i]-T_n/3)*lambda_2)*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])))
      if(Y_S.obs[i]<T_n&Y_S.obs[i]>=T_n*2/3) ll<-ll+log(lambda_3*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])*exp(-(T_n/3*lambda_1+T_n/3*lambda_2+(Y_S.obs[i]-T_n*2/3)*lambda_3)*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])))
    } else{
      chazard<-0
      if(Y_S.obs[i]<T_n/3) chazard<-Y_S.obs[i]*lambda_1
      if(Y_S.obs[i]<T_n*2/3&Y_S.obs[i]>=T_n/3) chazard<-T_n/3*lambda_1+(Y_S.obs[i]-T_n/3)*lambda_2
      if(Y_S.obs[i]>=T_n*2/3) chazard<-T_n/3*lambda_1+T_n/3*lambda_2+(Y_S.obs[i]-T_n*2/3)*lambda_3
      ll<-ll-chazard*exp(gamma_T*(Y_T[i]==1)+gamma_E*(Y_E[i]==2)+gamma_d[which(d.level==d[i])])
    }
  }
  return(ll)
}

rjmcmc.ETS<- function(data,dose.level,score,phi_S.lb,thin=1, N.post=2000,N.burnin=1000){#
  Y_T<-data$Y_T;Y_E<-data$Y_E;d<-data$d;dose<-data$dose;Y_S.obs<-data$Y_S.ob;status<-data$status;subgroup<-data$g           
  n<-length(d)
  d.all<-sort(unique(d))
  d.number<-length(d.all)
  T_n<-max(Y_S.obs)
  
  # priors
  mu.p.beta0<-0;sigma2.p.beta0<-10
  mu.p.beta1<-0;sigma2.p.beta1<-10
  mu.p.alpha0<-0;sigma2.p.alpha0<-10
  mu.p.alpha1<-0;sigma2.p.alpha1<-10
  mu.p.alpha2<-0;sigma2.p.alpha2<-1
  mu.p.alpha3<-0;sigma2.p.alpha3<-1
  mu.p.eta2<-0;sigma2.p.eta2<-10
  mu.p.gamma_T<-0;sigma2.p.gamma_T<-10
  mu.p.gamma_E<-0;sigma2.p.gamma_E<-10
  alpha.p.lambda_1<-0.1;beta.p.lambda_1<-0.1
  alpha.p.lambda_2<-0.1;beta.p.lambda_2<-0.1
  alpha.p.lambda_3<-0.1;beta.p.lambda_3<-0.1
  mu.p.gamma_d<-0;sigma2.p.gamma_d<-10
  
  ### "alpha0","alpha1","alpha2","alpha3","beta0","beta1","eta2","sigma12"  
  alpha0_t <- matrix(0,N.post*thin,3)
  alpha1_t <- matrix(0,N.post*thin,3)
  alpha2_t <- matrix(0,N.post*thin,3)
  alpha3_t <- matrix(0,N.post*thin,3)
  beta0_t <- matrix(0,N.post*thin,3)
  beta1_t <- matrix(0,N.post*thin,3)
  eta2_t <- matrix(0,N.post*thin,3)
  sigma12_t <- rep(0,N.post*thin)
  X_T_t <- matrix(0,N.post*thin,n)
  X_E_t <- matrix(0,N.post*thin,n)
  gamma_T_t<-matrix(0,N.post*thin,3)
  gamma_E_t<-matrix(0,N.post*thin,3)
  lambda_1_t<-matrix(0,N.post*thin,3)
  lambda_2_t<-matrix(0,N.post*thin,3)
  lambda_3_t<-matrix(0,N.post*thin,3)
  gamma_d_t<-array(0,c(N.post*thin,d.number-1,3))
  z_g_t<-matrix(1:3,N.post*thin,3,byrow=TRUE)
  omega_g_t<-matrix(1,N.post*thin,3)
  pi_t <- array(0,dim=c(2,3,4,N.post*thin,3))
  sur.p_t<-array(0,c(4,N.post*thin,3))
  
  ###set initial values
  alpha0<-rep(0,3);alpha1<-rep(0,3);alpha2<-rep(0,3);alpha3<-rep(0,3)
  beta0<-rep(0,3);beta1<-rep(0,3);eta2<-rep(0,3);sigma12<-0
  X_T <- rep(0,n);X_E <- rep(0,n)
  gamma_T<-rep(0,3);gamma_E<-rep(0,3)
  lambda_1<-rep(0,3);lambda_2<-rep(0,3);lambda_3<-rep(0,3)
  gamma_d <- matrix(0,d.number-1,3)
  z_g<-c(1,2,3);omega_g<-c(1,1,1)
  for(i in 1:3){
    id=(subgroup==i)
    data.ET=list(Y_T[id],Y_E[id],d[id])
    data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
    d.subgroup<-sort(unique(d[id&Y_E!=0]))
    MLE_ET<-nlminb(c(0,0,1,1,0,1,1,0),negative.log.likelihood.Y_ET,data=data.ET,dose.level=dose.level,
                   lower = c(-5, -5, 0.01, 0, -5,   0,  0, -0.99),
                   upper = c(5,   5, 1,    5,  5, 5, 10,  0.99))
    MLE_S<-nlminb(c(0,0,0.2,0.2,0.2,rep(0,length(d.subgroup)-1)),negative.log.likelihood.Y_S,data=data.S,dose.level=dose.level,
                  lower = c(-5, -5, 0, 0, 0, rep(-5,length(d.subgroup)-1)),
                  upper = c( 5,  5, 1, 1, 1, rep(5,length(d.subgroup)-1)))
    
    alpha0[i]<-MLE_ET[[1]][1]
    alpha1[i] <-MLE_ET[[1]][2]
    alpha2[i] <-MLE_ET[[1]][3]
    alpha3[i] <-MLE_ET[[1]][4]
    beta0[i] <-MLE_ET[[1]][5]
    beta1[i]  <-MLE_ET[[1]][6]
    eta2[i]   <-MLE_ET[[1]][7]
    sigma12<-sigma12+MLE_ET[[1]][8]
    gamma_T[i]<-MLE_S[[1]][1]
    gamma_E[i]<-MLE_S[[1]][2]
    lambda_1[i]<-MLE_S[[1]][3]
    lambda_2[i]<-MLE_S[[1]][4]
    lambda_3[i]<-MLE_S[[1]][5]
    gamma_d[match(d.subgroup[-1],d.all[-1]),i]<-MLE_S[[1]][-(1:5)]
    
    X_E_1 <- ifelse(Y_E[id]==0,-1,Y_E[id])
    X_E_2 <- ifelse(X_E_1==1,1,X_E_1)
    X_E[id] <- ifelse(X_E_2==2,100,X_E_2)
    X_T_1 <- ifelse(Y_T[id]==0,-1,Y_T[id])
    X_T[id] <- ifelse(X_T_1==1,1,X_T_1)
    # print(c(MLE_ET[[1]],MLE_S[[1]]))
  }
  sigma12<-sigma12/3
  alpha0_MLE<-alpha0;alpha1_MLE<-alpha1;alpha2_MLE<-alpha2;alpha3_MLE<-alpha3;beta0_MLE<-beta0;beta1_MLE<-beta1;eta2_MLE<-eta2;sigma12_MLE<-sigma12
  gamma_T_MLE<-gamma_T;gamma_E_MLE<-gamma_E;lambda_1_MLE<-lambda_1;lambda_2_MLE<-lambda_2;lambda_3_MLE<-lambda_3;gamma_d_MLE<-gamma_d
  
  for (ite in 1:(thin*N.post)) {
    #X_E and X_T
    for(i in 1:3){
      id=(subgroup==i)
      X_E_cut_lower <- c(-Inf,0,max(eta2[which(z_g==z_g[i])]))
      X_E_cut_upper <- c(0,min(eta2[which(z_g==z_g[i])]),Inf)
      X_T_cut <- c(-Inf,0,Inf)
      var.X <- 1-sigma12^2
      mean.X_E <- alpha0[i]+alpha1[i]*dose[id]^alpha3[i]/(alpha2[i]^alpha3[i]+dose[id]^alpha3[i])+sigma12*(X_T[id]-beta0[i]-beta1[i]*dose[id])
      X_E[id] <- rtnorm(sum(id),mean.X_E,sqrt(var.X),lower=X_E_cut_lower[Y_E[id]+1],upper=X_E_cut_upper[Y_E[id]+1])
      mean.X_T <-beta0[i]+beta1[i]*dose[id]+sigma12*(X_E[id]-alpha0[i]-alpha1[i]*dose[id]^alpha3[i]/(alpha2[i]^alpha3[i]+dose[id]^alpha3[i]))
      X_T[id] <- rtnorm(sum(id),mean.X_T,sqrt(var.X),lower=X_T_cut[Y_T[id]+1],upper=X_T_cut[Y_T[id]+2])
    }
    X_T_t[ite,] <- X_T
    X_E_t[ite,] <- X_E
    
    #beta0
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      if(omega_g[i]==1){
        mu.beta0.1 <- (sum(X_T[id]-beta1[i]*dose[id])-sigma12*sum(X_E[id]-alpha0[i]-alpha1[i]*dose[id]^alpha3[i]/(alpha2[i]^alpha3[i]+dose[id]^alpha3[i])))/sum(id)
        sigma2.beta0.1 <- (1-sigma12^2)/sum(id)
        
        sigma2.beta0 <- 1/(1/sigma2.beta0.1+1/sigma2.p.beta0)
        mu.beta0 <- sigma2.beta0*(mu.beta0.1/sigma2.beta0.1+mu.p.beta0/sigma2.p.beta0)
        beta0[i] <- rnorm(1,mu.beta0,sqrt(sigma2.beta0))
      } else{
        beta0[i] <- beta0[z_g[i]]
      }
    }
    beta0_t[ite,] <- beta0
    
    #beta1
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      if(omega_g[i]==1){
        var.bottom <- sum(dose[id]^2)
        mu.beta1.1 <- (sum(dose[id]*(X_T[id]-beta0[i]))-sigma12*sum(dose[id]*(X_E[id]-alpha0[i]-alpha1[i]*dose[id]^alpha3[i]/(alpha2[i]^alpha3[i]+dose[id]^alpha3[i]))))/var.bottom
        sigma2.beta1.1 <- (1-sigma12^2)/var.bottom
        
        sigma2.beta1 <- 1/(1/sigma2.beta1.1+1/sigma2.p.beta1)
        mu.beta1 <- sigma2.beta1*(mu.beta1.1/sigma2.beta1.1+mu.p.beta1/sigma2.p.beta1)
        beta1[i] <- rtnorm(1,mu.beta1,sqrt(sigma2.beta1),lower=0,upper=Inf)
      } else{
        beta1[i] <- beta1[z_g[i]]
      }
    }
    beta1_t[ite,] <- beta1
    
    #alpha0
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      if(omega_g[i]==1){
        mu.alpha0.1 <- (sum(X_E[id]-alpha1[i]*dose[id]^alpha3[i]/(alpha2[i]^alpha3[i]+dose[id]^alpha3[i]))-sigma12*sum(X_T[id]-beta0[i]-beta1[i]*dose[id]))/sum(id)
        sigma2.alpha0.1 <- (1-sigma12^2)/sum(id)
        
        sigma2.alpha0 <- 1/(1/sigma2.alpha0.1+1/sigma2.p.alpha0)
        mu.alpha0 <- sigma2.alpha0*(mu.alpha0.1/sigma2.alpha0.1+mu.p.alpha0/sigma2.p.alpha0)
        
        alpha0[i] <- rnorm(1,mu.alpha0,sqrt(sigma2.alpha0))
      } else{
        alpha0[i] <- alpha0[z_g[i]]
      }
    }
    alpha0_t[ite,] <- alpha0
    
    #alpha1
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      if(omega_g[i]==1){
        var.bottom <- sum((dose[id]^alpha3[i]/(alpha2[i]^alpha3[i]+dose[id]^alpha3[i]))^2)
        mu.alpha1.1 <- (sum((dose[id]^alpha3[i]/(alpha2[i]^alpha3[i]+dose[id]^alpha3[i]))*(X_E[id]-alpha0[i]))-
                          sigma12*sum((dose[id]^alpha3[i]/(alpha2[i]^alpha3[i]+dose[id]^alpha3[i]))*(X_T[id]-beta0[i]-beta1[i]*dose[id])))/var.bottom
        sigma2.alpha1.1 <- (1-sigma12^2)/var.bottom
        
        sigma2.alpha1 <- 1/(1/sigma2.alpha1.1+1/sigma2.p.alpha1)
        mu.alpha1 <- sigma2.alpha1*(mu.alpha1.1/sigma2.alpha1.1+mu.p.alpha1/sigma2.p.alpha1)
        alpha1[i] <- rnorm(1,mu.alpha1,sqrt(sigma2.alpha1))
      } else{
        alpha1[i] <- alpha1[z_g[i]]
      }
    }
    alpha1_t[ite,] <- alpha1
    
    
    #alpha2
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      if(omega_g[i]==1){
        logden <- function(x) {
          return(log.likelihood.X_E(X_E[id],X_T[id],dose[id],alpha0[i],alpha1[i],x,alpha3[i],beta0[i],beta1[i],sigma12)-(x-mu.p.alpha2)^2/(2*sigma2.p.alpha2))
        }
        alpha2[i] <- arms(alpha2[i],logden,function(x) ((x>=0)*(x<=1)),1)
      } else{
        alpha2[i] <- alpha2[z_g[i]]
      }
    }
    alpha2_t[ite,] <- alpha2
    
    #alpha3
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      if(omega_g[i]==1){
        logden <- function(x) {
          return(log.likelihood.X_E(X_E[id],X_T[id],dose[id],alpha0[i],alpha1[i],alpha2[i],x,beta0[i],beta1[i],sigma12)-(x-mu.p.alpha3)^2/(2*sigma2.p.alpha3))
        }
        alpha3[i] <- arms(alpha3[i],logden,function(x) ((x>=0)*(x<=5)),1)
      } else{
        alpha3[i] <- alpha3[z_g[i]]
      }
    }
    alpha3_t[ite,] <- alpha3
    
    # eta2
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      if(omega_g[i]==1){
        min_eta2 <- ifelse(sum(Y_E[id]==1)==0,0,max(max(subset(X_E[id],Y_E[id]==1)),0))
        max_eta2 <- min(10,ifelse(sum(Y_E[id]==2)==0,10,min(subset(X_E[id],Y_E[id]==2))))
        eta2[i]<- runif(1,min_eta2,max_eta2)
      }else{
        eta2[i] <- eta2[z_g[i]]
      }
    }
    eta2_t[ite,] <- eta2
    
    #sigma12
    logden <- function(x) {
      return(log.likelihood.sigma12(X_E,X_T,dose,alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],x))
    }
    sigma12 <- arms(sigma12,logden,function(x) ((x>-1)*(x<1)),1)
    sigma12_t[ite] <- sigma12
    
    
    #gamma_T
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      if(omega_g[i]==1){
        data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
        d.subgroup<-sort(unique(d[id&Y_E!=0]))
        logden <- function(x) {
          return(log.likelihood.Y_S_withoutconstraint(c(x,gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-(x-mu.p.gamma_T)^2/(2*sigma2.p.gamma_T))
        }
        gamma_T[i] <- arms(gamma_T[i],logden,function(x) ((x>=-5)*(x<=5)),1)
      }else{
        gamma_T[i] <- gamma_T[z_g[i]]
      }
    }
    gamma_T_t[ite,] <- gamma_T
    
    
    #gamma_E
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      if(omega_g[i]==1){
        data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
        d.subgroup<-sort(unique(d[id&Y_E!=0]))
        logden <- function(x) {
          return(log.likelihood.Y_S_withoutconstraint(c(gamma_T[i],x,lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-(x-mu.p.gamma_E)^2/(2*sigma2.p.gamma_E))
        }
        gamma_E[i] <- arms(gamma_E[i],logden,function(x) ((x>=-5)*(x<=5)),1)
      }else{
        gamma_E[i] <- gamma_E[z_g[i]]
      }
    }
    gamma_E_t[ite,] <- gamma_E
    
    
    #lambda_1
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      if(omega_g[i]==1){
        data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
        logden <- function(x) {
          return(log.likelihood.Y_S_withoutconstraint(c(gamma_T[i],gamma_E[i],x,lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)+(alpha.p.lambda_1-1)*log(x)-beta.p.lambda_1*x)
        }
        lambda_1[i] <- arms(lambda_1[i],logden,function(x) ((x>=0)*(x<=3)),1)
      }else{
        lambda_1[i] <- lambda_1[z_g[i]]
      }
    }
    lambda_1_t[ite,] <- lambda_1
    
    #lambda_2
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      if(omega_g[i]==1){
        data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
        logden <- function(x) {
          return(log.likelihood.Y_S_withoutconstraint(c(gamma_T[i],gamma_E[i],lambda_1[i],x,lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)+(alpha.p.lambda_2-1)*log(x)-beta.p.lambda_2*x)
        }
        lambda_2[i] <- arms(lambda_2[i],logden,function(x) ((x>=0)*(x<=3)),1)
      }else{
        lambda_2[i] <- lambda_2[z_g[i]]
      }
    }
    lambda_2_t[ite,] <- lambda_2
    
    #lambda_3
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      if(omega_g[i]==1){
        data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
        logden <- function(x) {
          return(log.likelihood.Y_S_withoutconstraint(c(gamma_T[i],gamma_E[i],lambda_2[i],x,lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)+(alpha.p.lambda_3-1)*log(x)-beta.p.lambda_3*x)
        }
        lambda_3[i] <- arms(lambda_3[i],logden,function(x) ((x>=0)*(x<=3)),1)
      }else{
        lambda_3[i] <- lambda_3[z_g[i]]
      }
    }
    lambda_3_t[ite,] <- lambda_3
    
    #gamma_d
    for(i in 1:3){
      id=(subgroup %in% which(z_g==i))
      if(omega_g[i]==1){
        data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
        d.subgroup<-sort(unique(d[id&Y_E!=0]))
        if(length(d.subgroup)>1){
          gamma_d_cluster<-gamma_d[match(d.subgroup[-1],d.all[-1]),i]
          for(j in 1:length(gamma_d_cluster)){
            logden <- function(x) {
              gamma_d_cluster[j]<-x
              return(log.likelihood.Y_S_withoutconstraint(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d_cluster),data=data.S,dose.level = dose.level)-(x-mu.p.gamma_d)^2/(2*sigma2.p.gamma_d))
            }
            gamma_d_cluster[j]<- arms(gamma_d_cluster[j],logden,function(x) ((x>=-5)*(x<=5)),1)
          }
          gamma_d[,i] <-0;gamma_d[match(d.subgroup[-1],d.all[-1]),i]<-gamma_d_cluster
        }
      }else{
        gamma_d[,i] <- gamma_d[,z_g[i]]
      }
    }
    gamma_d_t[ite,,] <- gamma_d
    
    for(i in 1:3){
      cor=matrix( c(1,sigma12,sigma12,1),nrow=2)
      for(k in 1:4){
        pi_t[1,1,k,ite,i]<-pmvnorm(lower=c(-Inf,-Inf),upper=c(0,0), mean=c( beta0[i]+beta1[i]*dose.level[k], alpha0[i]+alpha1[i]*dose.level[k]^alpha3[i]/(alpha2[i]^alpha3[i]+dose.level[k]^alpha3[i])), sigma=cor)[1]
        pi_t[2,1,k,ite,i]<-pmvnorm(lower=c(0,-Inf),upper=c(Inf,0), mean=c( beta0[i]+beta1[i]*dose.level[k], alpha0[i]+alpha1[i]*dose.level[k]^alpha3[i]/(alpha2[i]^alpha3[i]+dose.level[k]^alpha3[i])), sigma=cor)[1]
        pi_t[1,2,k,ite,i]<-pmvnorm(lower=c(-Inf,0),upper=c(0,eta2[i]), mean=c( beta0[i]+beta1[i]*dose.level[k], alpha0[i]+alpha1[i]*dose.level[k]^alpha3[i]/(alpha2[i]^alpha3[i]+dose.level[k]^alpha3[i])), sigma=cor)[1]
        pi_t[2,2,k,ite,i]<-pmvnorm(lower=c(0,0),upper=c(Inf,eta2[i]), mean=c( beta0[i]+beta1[i]*dose.level[k], alpha0[i]+alpha1[i]*dose.level[k]^alpha3[i]/(alpha2[i]^alpha3[i]+dose.level[k]^alpha3[i])), sigma=cor)[1]
        pi_t[1,3,k,ite,i]<-pmvnorm(lower=c(-Inf,eta2[i]),upper=c(0,Inf), mean=c( beta0[i]+beta1[i]*dose.level[k], alpha0[i]+alpha1[i]*dose.level[k]^alpha3[i]/(alpha2[i]^alpha3[i]+dose.level[k]^alpha3[i])), sigma=cor)[1]
        pi_t[2,3,k,ite,i]<-pmvnorm(lower=c(0,eta2[i]),upper=c(Inf,Inf), mean=c( beta0[i]+beta1[i]*dose.level[k], alpha0[i]+alpha1[i]*dose.level[k]^alpha3[i]/(alpha2[i]^alpha3[i]+dose.level[k]^alpha3[i])), sigma=cor)[1]
      }
    }
    # pi_t[,,,ite,]
    for(i in 1:3){
      id=(subgroup %in% which(z_g==z_g[i]))
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      gamma_d_all<-rep(0,4)
      gamma_d_all[d.all[-1]]<-gamma_d[,i]
      for(j in 1:4){
        pi.m<-pi_t[,,j,ite,i]
        sur.p_t[j,ite,i]<-pi.m[1,2]*exp(-T_n/3*sum(lambda_1[i],lambda_2[i],lambda_3[i])*exp(gamma_d_all[j]))+
          pi.m[1,3]*exp(-T_n/3*sum(lambda_1[i],lambda_2[i],lambda_3[i])*exp(gamma_E[i]+gamma_d_all[j]))+
          pi.m[2,2]*exp(-T_n/3*sum(lambda_1[i],lambda_2[i],lambda_3[i])*exp(gamma_T[i]+gamma_d_all[j]))+
          pi.m[2,3]*exp(-T_n/3*sum(lambda_1[i],lambda_2[i],lambda_3[i])*exp(gamma_T[i]+gamma_E[i]+gamma_d_all[j]))
      }
    }
    t(sur.p_t[,ite,])
    
    
    
    sigma2_star<-1
    
    if(sum(omega_g==1)==2){
      beta0_star<-beta0;beta1_star<-beta1;alpha0_star<-alpha0;alpha1_star<-alpha1
      alpha2_star<-alpha2;alpha3_star<-alpha3;eta2_star<-eta2
      gamma_T_star<-gamma_T;gamma_E_star<-gamma_E;lambda_1_star<-lambda_1
      lambda_2_star<-lambda_2;lambda_3_star<-lambda_3;gamma_d_star<-gamma_d
      i=which(omega_g==0)
      id=(subgroup==i)
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      beta0_star[i] <- rnorm(1,beta0_MLE[i],sqrt(sigma2_star))
      beta1_star[i] <- rtnorm(1,beta1_MLE[i],sqrt(sigma2_star),lower=0,upper=Inf)
      alpha0_star[i] <- rnorm(1,alpha0_MLE[i],sqrt(sigma2_star))
      alpha1_star[i] <- rnorm(1,alpha1_MLE[i],sqrt(sigma2_star))
      alpha2_star[i] <- rtnorm(1,alpha2_MLE[i],sqrt(0.1),lower=0,upper=1)
      alpha3_star[i] <- rtnorm(1,alpha3_MLE[i],sqrt(sigma2_star),lower=0,upper=5)
      eta2_star[i] <- rtnorm(1,eta2_MLE[i],sqrt(sigma2_star),lower=0,upper=10)
      gamma_T_star[i] <- rtnorm(1,gamma_T_MLE[i],sqrt(sigma2_star),lower=-5,upper=5)
      gamma_E_star[i] <- rtnorm(1,gamma_E_MLE[i],sqrt(sigma2_star),lower=-5,upper=5)
      lambda_1_star[i] <-rtnorm(1,lambda_1_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      lambda_2_star[i] <- rtnorm(1,lambda_2_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      lambda_3_star[i] <- rtnorm(1,lambda_3_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      gamma_d_star[,i]<-0
      for(j in match(d.subgroup[-1],d.all[-1])){
        gamma_d_star[j,i]<-rtnorm(1,gamma_d_MLE[j,i],sqrt(sigma2_star),lower=-5,upper=5)
      }
      data.ET=list(Y_T[id],Y_E[id],d[id])
      data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
      
      accept_pr<-log.likelihood.Y_ET(c(alpha0_star[i],alpha1_star[i],alpha2_star[i],alpha3_star[i],beta0_star[i],beta1_star[i],eta2_star[i],sigma12),data=data.ET,dose.level = dose.level)+
        log.likelihood.Y_S(c(gamma_T_star[i],gamma_E_star[i],lambda_1_star[i],lambda_2_star[i],lambda_3_star[i],gamma_d_star[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log.likelihood.Y_ET(c(alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],eta2[i],sigma12),data=data.ET,dose.level = dose.level)-
        log.likelihood.Y_S(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)+
        log(3/4*dnorm(beta0_star[i], mu.p.beta0, sqrt(sigma2.p.beta0))/dnorm(beta0_star[i], beta0_MLE[i], sqrt(sigma2_star)) *
              dtnorm(beta1_star[i], mu.p.beta1, sqrt(sigma2.p.beta1), lower = 0, upper = Inf)/dtnorm(beta1_star[i], beta1_MLE[i], sqrt(sigma2_star), lower = 0, upper = Inf)*
              dnorm(alpha0_star[i], mu.p.alpha0, sqrt(sigma2.p.alpha0))/dnorm(alpha0_star[i], alpha0_MLE[i], sqrt(sigma2_star))*
              dnorm(alpha1_star[i], mu.p.alpha1, sqrt(sigma2.p.alpha1))/dnorm(alpha1_star[i], alpha1_MLE[i], sqrt(sigma2_star))*
              dtnorm(alpha2_star[i], mu.p.alpha2, sqrt(sigma2.p.alpha2), lower = 0, upper = 1)/dtnorm(alpha2_star[i], alpha2_MLE[i], sqrt(0.1), lower = 0, upper = 1)*
              dtnorm(alpha3_star[i], mu.p.alpha3, sqrt(sigma2.p.alpha3), lower = 0, upper = 5)/dtnorm(alpha3_star[i], alpha3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 5)*
              dunif(eta2_star[i], min = 0, max = 10)/dtnorm(eta2_star[i], eta2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 10)*
              dtnorm(gamma_T_star[i], mu.p.gamma_T, sqrt(sigma2.p.gamma_T), lower = -5, upper = 5)/dtnorm(gamma_T_star[i], gamma_T_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dtnorm(gamma_E_star[i], mu.p.gamma_E, sqrt(sigma2.p.gamma_E), lower = -5, upper = 5)/dtnorm(gamma_E_star[i], gamma_E_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dgamma(lambda_1_star[i], alpha.p.lambda_1, beta.p.lambda_1)/(pgamma(3, alpha.p.lambda_1, beta.p.lambda_1)*dtnorm(lambda_1_star[i], lambda_1_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_2_star[i], alpha.p.lambda_2, beta.p.lambda_2)/(pgamma(3, alpha.p.lambda_2, beta.p.lambda_2)*dtnorm(lambda_2_star[i], lambda_2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_3_star[i], alpha.p.lambda_3, beta.p.lambda_3)/(pgamma(3, alpha.p.lambda_3, beta.p.lambda_3)*dtnorm(lambda_3_star[i], lambda_3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              ifelse(length(match(d.subgroup[-1],d.all[-1]))>0,prod(dtnorm(gamma_d_star[match(d.subgroup[-1],d.all[-1]), i], mu.p.gamma_d, sqrt(sigma2.p.gamma_d), lower = -5, upper = 5)/dtnorm(gamma_d_star[match(d.subgroup[-1],d.all[-1]), i],
                                                                                                                                                                                                 gamma_d_MLE[match(d.subgroup[-1],d.all[-1]), i], sqrt(sigma2_star), lower = -5, upper = 5)),1))
      accept_pr<-min(exp(accept_pr),1)
      # print(accept_pr)
      dec<-rbinom(1, size = 1, prob =accept_pr)
      if(dec==1){
        beta0 <- beta0_star;beta1 <- beta1_star;alpha0 <- alpha0_star;alpha1 <- alpha1_star;alpha2 <- alpha2_star;alpha3 <- alpha3_star
        eta2 <- eta2_star;gamma_T <- gamma_T_star;gamma_E <- gamma_E_star;lambda_1 <- lambda_1_star;lambda_2 <- lambda_2_star;lambda_3 <- lambda_3_star
        gamma_d <- gamma_d_star;z_g<-c(1,2,3);omega_g<-rep(1,3)
      }
    }
    if(sum(omega_g==1)==2){
      beta0_star<-beta0;beta1_star<-beta1;alpha0_star<-alpha0;alpha1_star<-alpha1
      alpha2_star<-alpha2;alpha3_star<-alpha3;eta2_star<-eta2
      gamma_T_star<-gamma_T;gamma_E_star<-gamma_E;lambda_1_star<-lambda_1
      lambda_2_star<-lambda_2;lambda_3_star<-lambda_3;gamma_d_star<-gamma_d
      i=which(!(duplicated(z_g)|duplicated(z_g, fromLast = TRUE)))  
      id=(subgroup==i)
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      beta0_star[i] <-  beta0[-i][1]
      beta1_star[i] <-  beta1[-i][1]
      alpha0_star[i] <- alpha0[-i][1]
      alpha1_star[i] <- alpha1[-i][1]
      alpha2_star[i] <- alpha2[-i][1]
      alpha3_star[i] <- alpha3[-i][1]
      eta2_star[i] <-   eta2[-i][1]
      gamma_T_star[i] <- gamma_T_MLE[-i][1]
      gamma_E_star[i] <- gamma_E[-i][1]
      lambda_1_star[i] <-lambda_1[-i][1]
      lambda_2_star[i] <- lambda_2[-i][1]
      lambda_3_star[i] <- lambda_3[-i][1]
      gamma_d_star[,i]<-0
      gamma_d_star[,i]<-gamma_d[,setdiff(1:3,i)[1]]
      
      data.ET=list(Y_T[id],Y_E[id],d[id])
      data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
      
      accept_pr<-log.likelihood.Y_ET(c(alpha0_star[i],alpha1_star[i],alpha2_star[i],alpha3_star[i],beta0_star[i],beta1_star[i],eta2_star[i],sigma12),data=data.ET,dose.level = dose.level)+
        log.likelihood.Y_S(c(gamma_T_star[i],gamma_E_star[i],lambda_1_star[i],lambda_2_star[i],lambda_3_star[i],gamma_d_star[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log.likelihood.Y_ET(c(alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],eta2[i],sigma12),data=data.ET,dose.level = dose.level)-
        log.likelihood.Y_S(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log(3/4*dnorm(beta0[i], mu.p.beta0, sqrt(sigma2.p.beta0))/dnorm(beta0[i], beta0_MLE[i], sqrt(sigma2_star)) *
              dtnorm(beta1[i], mu.p.beta1, sqrt(sigma2.p.beta1), lower = 0, upper = Inf)/dtnorm(beta1[i], beta1_MLE[i], sqrt(sigma2_star), lower = 0, upper = Inf)*
              dnorm(alpha0[i], mu.p.alpha0, sqrt(sigma2.p.alpha0))/dnorm(alpha0[i], alpha0_MLE[i], sqrt(sigma2_star))*
              dnorm(alpha1[i], mu.p.alpha1, sqrt(sigma2.p.alpha1))/dnorm(alpha1[i], alpha1_MLE[i], sqrt(sigma2_star))*
              dtnorm(alpha2[i], mu.p.alpha2, sqrt(sigma2.p.alpha2), lower = 0, upper = 1)/dtnorm(alpha2[i], alpha2_MLE[i], sqrt(0.1), lower = 0, upper = 1)*
              dtnorm(alpha3[i], mu.p.alpha3, sqrt(sigma2.p.alpha3), lower = 0, upper = 5)/dtnorm(alpha3[i], alpha3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 5)*
              dunif(eta2[i], min = 0, max = 10)/dtnorm(eta2[i], eta2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 10)*
              dtnorm(gamma_T[i], mu.p.gamma_T, sqrt(sigma2.p.gamma_T), lower = -5, upper = 5)/dtnorm(gamma_T[i], gamma_T_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dtnorm(gamma_E[i], mu.p.gamma_E, sqrt(sigma2.p.gamma_E), lower = -5, upper = 5)/dtnorm(gamma_E[i], gamma_E_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dgamma(lambda_1[i], alpha.p.lambda_1, beta.p.lambda_1)/(pgamma(3, alpha.p.lambda_1, beta.p.lambda_1)*dtnorm(lambda_1[i], lambda_1_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_2[i], alpha.p.lambda_2, beta.p.lambda_2)/(pgamma(3, alpha.p.lambda_2, beta.p.lambda_2)*dtnorm(lambda_2[i], lambda_2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_3[i], alpha.p.lambda_3, beta.p.lambda_3)/(pgamma(3, alpha.p.lambda_3, beta.p.lambda_3)*dtnorm(lambda_3[i], lambda_3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              ifelse(length(match(d.subgroup[-1],d.all[-1]))>0,prod(dtnorm(gamma_d[match(d.subgroup[-1],d.all[-1]), i], mu.p.gamma_d, sqrt(sigma2.p.gamma_d), lower = -5, upper = 5)/dtnorm(gamma_d[match(d.subgroup[-1],d.all[-1]), i],
                                                                                                                                                                                            gamma_d_MLE[match(d.subgroup[-1],d.all[-1]), i], sqrt(sigma2_star), lower = -5, upper = 5)),1))
      accept_pr<-min(exp(accept_pr),1)
      # print(accept_pr)
      dec<-rbinom(1, size = 1, prob =accept_pr)
      if(dec==1){
        beta0 <- beta0_star;beta1 <- beta1_star;alpha0 <- alpha0_star;alpha1 <- alpha1_star;alpha2 <- alpha2_star;alpha3 <- alpha3_star
        eta2 <- eta2_star;gamma_T <- gamma_T_star;gamma_E <- gamma_E_star;lambda_1 <- lambda_1_star;lambda_2 <- lambda_2_star;lambda_3 <- lambda_3_star
        gamma_d <- gamma_d_star;z_g<-c(1,1,1);omega_g<-c(1,0,0)
      }
    }
    if(sum(omega_g==1)==2){
      beta0_star<-beta0;beta1_star<-beta1;alpha0_star<-alpha0;alpha1_star<-alpha1
      alpha2_star<-alpha2;alpha3_star<-alpha3;eta2_star<-eta2
      gamma_T_star<-gamma_T;gamma_E_star<-gamma_E;lambda_1_star<-lambda_1
      lambda_2_star<-lambda_2;lambda_3_star<-lambda_3;gamma_d_star<-gamma_d
      z_g_star<-z_g;omega_g_star<-omega_g
      s.number<-rep(0,3)
      for(i in 1:3) s.number[i]=sum(z_g==i)
      index_switch<-which(z_g==which(s.number==2))[1]
      index_combine<-which(z_g==which(s.number==1))
      z_g_star[c(index_combine,index_switch)]=min(index_combine,index_switch)
      z_g_star[setdiff(1:3,c(index_combine,index_switch))]=setdiff(1:3,c(index_combine,index_switch))
      for(i in 1:3) omega_g_star[i]<-as.numeric(z_g_star[i]==i)
      z_g_star;omega_g_star
      beta0_star[index_switch]<-beta0[index_combine]
      beta1_star[index_switch] <-  beta1[index_combine]
      alpha0_star[index_switch] <- alpha0[index_combine]
      alpha1_star[index_switch] <- alpha1[index_combine]
      alpha2_star[index_switch] <- alpha2[index_combine]
      alpha3_star[index_switch] <- alpha3[index_combine]
      eta2_star[index_switch] <-   eta2[index_combine]
      gamma_T_star[index_switch] <- gamma_T_MLE[index_combine]
      gamma_E_star[index_switch] <- gamma_E[index_combine]
      lambda_1_star[index_switch] <-lambda_1[index_combine]
      lambda_2_star[index_switch] <- lambda_2[index_combine]
      lambda_3_star[index_switch] <- lambda_3[index_combine]
      gamma_d_star[,index_switch]<-gamma_d[,index_combine]
      
      accept_pr<-0
      for(i in 1:3){
        id=(subgroup==i)
        d.subgroup<-sort(unique(d[id&Y_E!=0]))
        data.ET=list(Y_T[id],Y_E[id],d[id])
        data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
        accept_pr<-accept_pr+log.likelihood.Y_ET(c(alpha0_star[i],alpha1_star[i],alpha2_star[i],alpha3_star[i],beta0_star[i],beta1_star[i],eta2_star[i],sigma12),data=data.ET,dose.level = dose.level)+
          log.likelihood.Y_S(c(gamma_T_star[i],gamma_E_star[i],lambda_1_star[i],lambda_2_star[i],lambda_3_star[i],gamma_d_star[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
          log.likelihood.Y_ET(c(alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],eta2[i],sigma12),data=data.ET,dose.level = dose.level)-
          log.likelihood.Y_S(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)
      }
      accept_pr<-min(exp(accept_pr),1)
      # print(accept_pr)
      dec<-rbinom(1, size = 1, prob =accept_pr)
      if(dec==1){
        beta0 <- beta0_star;beta1 <- beta1_star;alpha0 <- alpha0_star;alpha1 <- alpha1_star;alpha2 <- alpha2_star;alpha3 <- alpha3_star
        eta2 <- eta2_star;gamma_T <- gamma_T_star;gamma_E <- gamma_E_star;lambda_1 <- lambda_1_star;lambda_2 <- lambda_2_star;lambda_3 <- lambda_3_star
        gamma_d <- gamma_d_star;z_g<-z_g_star;omega_g<-omega_g_star
      }
    }
    if(sum(omega_g==1)==2){
      beta0_star<-beta0;beta1_star<-beta1;alpha0_star<-alpha0;alpha1_star<-alpha1
      alpha2_star<-alpha2;alpha3_star<-alpha3;eta2_star<-eta2
      gamma_T_star<-gamma_T;gamma_E_star<-gamma_E;lambda_1_star<-lambda_1
      lambda_2_star<-lambda_2;lambda_3_star<-lambda_3;gamma_d_star<-gamma_d
      z_g_star<-z_g;omega_g_star<-omega_g
      s.number<-rep(0,3)
      for(i in 1:3) s.number[i]=sum(z_g==i)
      index_switch<-which(z_g==which(s.number==2))[2]
      index_combine<-which(z_g==which(s.number==1))
      z_g_star[c(index_combine,index_switch)]=min(index_combine,index_switch)
      z_g_star[setdiff(1:3,c(index_combine,index_switch))]=setdiff(1:3,c(index_combine,index_switch))
      for(i in 1:3) omega_g_star[i]<-as.numeric(z_g_star[i]==i)
      z_g_star;omega_g_star
      beta0_star[index_switch]<-beta0[index_combine]
      beta1_star[index_switch] <-  beta1[index_combine]
      alpha0_star[index_switch] <- alpha0[index_combine]
      alpha1_star[index_switch] <- alpha1[index_combine]
      alpha2_star[index_switch] <- alpha2[index_combine]
      alpha3_star[index_switch] <- alpha3[index_combine]
      eta2_star[index_switch] <-   eta2[index_combine]
      gamma_T_star[index_switch] <- gamma_T_MLE[index_combine]
      gamma_E_star[index_switch] <- gamma_E[index_combine]
      lambda_1_star[index_switch] <-lambda_1[index_combine]
      lambda_2_star[index_switch] <- lambda_2[index_combine]
      lambda_3_star[index_switch] <- lambda_3[index_combine]
      gamma_d_star[,index_switch]<-gamma_d[,index_combine]
      
      accept_pr<-0
      for(i in 1:3){
        id=(subgroup==i)
        d.subgroup<-sort(unique(d[id&Y_E!=0]))
        data.ET=list(Y_T[id],Y_E[id],d[id])
        data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
        accept_pr<-accept_pr+log.likelihood.Y_ET(c(alpha0_star[i],alpha1_star[i],alpha2_star[i],alpha3_star[i],beta0_star[i],beta1_star[i],eta2_star[i],sigma12),data=data.ET,dose.level = dose.level)+
          log.likelihood.Y_S(c(gamma_T_star[i],gamma_E_star[i],lambda_1_star[i],lambda_2_star[i],lambda_3_star[i],gamma_d_star[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
          log.likelihood.Y_ET(c(alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],eta2[i],sigma12),data=data.ET,dose.level = dose.level)-
          log.likelihood.Y_S(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)
      }
      accept_pr<-min(exp(accept_pr),1)
      # print(accept_pr)
      dec<-rbinom(1, size = 1, prob =accept_pr)
      if(dec==1){
        beta0 <- beta0_star;beta1 <- beta1_star;alpha0 <- alpha0_star;alpha1 <- alpha1_star;alpha2 <- alpha2_star;alpha3 <- alpha3_star
        eta2 <- eta2_star;gamma_T <- gamma_T_star;gamma_E <- gamma_E_star;lambda_1 <- lambda_1_star;lambda_2 <- lambda_2_star;lambda_3 <- lambda_3_star
        gamma_d <- gamma_d_star;z_g<-z_g_star;omega_g<-omega_g_star
      }
    }
    
    if(sum(omega_g==1)==3){
      beta0_star<-beta0;beta1_star<-beta1;alpha0_star<-alpha0;alpha1_star<-alpha1
      alpha2_star<-alpha2;alpha3_star<-alpha3;eta2_star<-eta2
      gamma_T_star<-gamma_T;gamma_E_star<-gamma_E;lambda_1_star<-lambda_1
      lambda_2_star<-lambda_2;lambda_3_star<-lambda_3;gamma_d_star<-gamma_d
      index_cluster<-c(1,2)
      z_g_star<-z_g;omega_g_star<-omega_g
      index_1<-min(index_cluster);index_2<-max(index_cluster)
      z_g_star[index_2]<-index_1;omega_g_star[index_2]<-0
      
      beta0_star[index_2]<-beta0[index_1]
      beta1_star[index_2] <-  beta1[index_1]
      alpha0_star[index_2] <- alpha0[index_1]
      alpha1_star[index_2] <- alpha1[index_1]
      alpha2_star[index_2] <- alpha2[index_1]
      alpha3_star[index_2] <- alpha3[index_1]
      eta2_star[index_2] <-   eta2[index_1]
      gamma_T_star[index_2] <- gamma_T_MLE[index_1]
      gamma_E_star[index_2] <- gamma_E[index_1]
      lambda_1_star[index_2] <-lambda_1[index_1]
      lambda_2_star[index_2] <- lambda_2[index_1]
      lambda_3_star[index_2] <- lambda_3[index_1]
      gamma_d_star[,index_2]<-gamma_d[,index_1]
      
      i=index_2
      id=subgroup==i
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      data.ET=list(Y_T[id],Y_E[id],d[id])
      data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
      
      accept_pr<-log.likelihood.Y_ET(c(alpha0_star[i],alpha1_star[i],alpha2_star[i],alpha3_star[i],beta0_star[i],beta1_star[i],eta2_star[i],sigma12),data=data.ET,dose.level = dose.level)+
        log.likelihood.Y_S(c(gamma_T_star[i],gamma_E_star[i],lambda_1_star[i],lambda_2_star[i],lambda_3_star[i],gamma_d_star[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log.likelihood.Y_ET(c(alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],eta2[i],sigma12),data=data.ET,dose.level = dose.level)-
        log.likelihood.Y_S(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log(4/3*dnorm(beta0[i], mu.p.beta0, sqrt(sigma2.p.beta0))/dnorm(beta0[i], beta0_MLE[i], sqrt(sigma2_star)) *
              dtnorm(beta1[i], mu.p.beta1, sqrt(sigma2.p.beta1), lower = 0, upper = Inf)/dtnorm(beta1[i], beta1_MLE[i], sqrt(sigma2_star), lower = 0, upper = Inf)*
              dnorm(alpha0[i], mu.p.alpha0, sqrt(sigma2.p.alpha0))/dnorm(alpha0[i], alpha0_MLE[i], sqrt(sigma2_star))*
              dnorm(alpha1[i], mu.p.alpha1, sqrt(sigma2.p.alpha1))/dnorm(alpha1[i], alpha1_MLE[i], sqrt(sigma2_star))*
              dtnorm(alpha2[i], mu.p.alpha2, sqrt(sigma2.p.alpha2), lower = 0, upper = 1)/dtnorm(alpha2[i], alpha2_MLE[i], sqrt(0.1), lower = 0, upper = 1)*
              dtnorm(alpha3[i], mu.p.alpha3, sqrt(sigma2.p.alpha3), lower = 0, upper = 5)/dtnorm(alpha3[i], alpha3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 5)*
              dunif(eta2[i], min = 0, max = 10)/dtnorm(eta2[i], eta2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 10)*
              dtnorm(gamma_T[i], mu.p.gamma_T, sqrt(sigma2.p.gamma_T), lower = -5, upper = 5)/dtnorm(gamma_T[i], gamma_T_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dtnorm(gamma_E[i], mu.p.gamma_E, sqrt(sigma2.p.gamma_E), lower = -5, upper = 5)/dtnorm(gamma_E[i], gamma_E_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dgamma(lambda_1[i], alpha.p.lambda_1, beta.p.lambda_1)/(pgamma(3, alpha.p.lambda_1, beta.p.lambda_1)*dtnorm(lambda_1[i], lambda_1_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_2[i], alpha.p.lambda_2, beta.p.lambda_2)/(pgamma(3, alpha.p.lambda_2, beta.p.lambda_2)*dtnorm(lambda_2[i], lambda_2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_3[i], alpha.p.lambda_3, beta.p.lambda_3)/(pgamma(3, alpha.p.lambda_3, beta.p.lambda_3)*dtnorm(lambda_3[i], lambda_3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              ifelse(length(match(d.subgroup[-1],d.all[-1]))>0,prod(dtnorm(gamma_d[match(d.subgroup[-1],d.all[-1]), i], mu.p.gamma_d, sqrt(sigma2.p.gamma_d), lower = -5, upper = 5)/dtnorm(gamma_d[match(d.subgroup[-1],d.all[-1]), i],
                                                                                                                                                                                            gamma_d_MLE[match(d.subgroup[-1],d.all[-1]), i], sqrt(sigma2_star), lower = -5, upper = 5)),1))
      accept_pr<-min(exp(accept_pr),1)
      # print(accept_pr)
      dec<-rbinom(1, size = 1, prob =accept_pr)
      if(dec==1){
        beta0 <- beta0_star;beta1 <- beta1_star;alpha0 <- alpha0_star;alpha1 <- alpha1_star;alpha2 <- alpha2_star;alpha3 <- alpha3_star
        eta2 <- eta2_star;gamma_T <- gamma_T_star;gamma_E <- gamma_E_star;lambda_1 <- lambda_1_star;lambda_2 <- lambda_2_star;lambda_3 <- lambda_3_star
        gamma_d <- gamma_d_star;z_g<-z_g_star;omega_g<-omega_g_star
      }
    }
    if(sum(omega_g==1)==3){
      beta0_star<-beta0;beta1_star<-beta1;alpha0_star<-alpha0;alpha1_star<-alpha1
      alpha2_star<-alpha2;alpha3_star<-alpha3;eta2_star<-eta2
      gamma_T_star<-gamma_T;gamma_E_star<-gamma_E;lambda_1_star<-lambda_1
      lambda_2_star<-lambda_2;lambda_3_star<-lambda_3;gamma_d_star<-gamma_d
      index_cluster<-c(1,3)
      z_g_star<-z_g;omega_g_star<-omega_g
      index_1<-min(index_cluster);index_2<-max(index_cluster)
      z_g_star[index_2]<-index_1;omega_g_star[index_2]<-0
      
      beta0_star[index_2]<-beta0[index_1]
      beta1_star[index_2] <-  beta1[index_1]
      alpha0_star[index_2] <- alpha0[index_1]
      alpha1_star[index_2] <- alpha1[index_1]
      alpha2_star[index_2] <- alpha2[index_1]
      alpha3_star[index_2] <- alpha3[index_1]
      eta2_star[index_2] <-   eta2[index_1]
      gamma_T_star[index_2] <- gamma_T_MLE[index_1]
      gamma_E_star[index_2] <- gamma_E[index_1]
      lambda_1_star[index_2] <-lambda_1[index_1]
      lambda_2_star[index_2] <- lambda_2[index_1]
      lambda_3_star[index_2] <- lambda_3[index_1]
      gamma_d_star[,index_2]<-gamma_d[,index_1]
      
      i=index_2
      id=subgroup==i
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      data.ET=list(Y_T[id],Y_E[id],d[id])
      data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
      
      accept_pr<-log.likelihood.Y_ET(c(alpha0_star[i],alpha1_star[i],alpha2_star[i],alpha3_star[i],beta0_star[i],beta1_star[i],eta2_star[i],sigma12),data=data.ET,dose.level = dose.level)+
        log.likelihood.Y_S(c(gamma_T_star[i],gamma_E_star[i],lambda_1_star[i],lambda_2_star[i],lambda_3_star[i],gamma_d_star[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log.likelihood.Y_ET(c(alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],eta2[i],sigma12),data=data.ET,dose.level = dose.level)-
        log.likelihood.Y_S(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log(4/3*dnorm(beta0[i], mu.p.beta0, sqrt(sigma2.p.beta0))/dnorm(beta0[i], beta0_MLE[i], sqrt(sigma2_star)) *
              dtnorm(beta1[i], mu.p.beta1, sqrt(sigma2.p.beta1), lower = 0, upper = Inf)/dtnorm(beta1[i], beta1_MLE[i], sqrt(sigma2_star), lower = 0, upper = Inf)*
              dnorm(alpha0[i], mu.p.alpha0, sqrt(sigma2.p.alpha0))/dnorm(alpha0[i], alpha0_MLE[i], sqrt(sigma2_star))*
              dnorm(alpha1[i], mu.p.alpha1, sqrt(sigma2.p.alpha1))/dnorm(alpha1[i], alpha1_MLE[i], sqrt(sigma2_star))*
              dtnorm(alpha2[i], mu.p.alpha2, sqrt(sigma2.p.alpha2), lower = 0, upper = 1)/dtnorm(alpha2[i], alpha2_MLE[i], sqrt(0.1), lower = 0, upper = 1)*
              dtnorm(alpha3[i], mu.p.alpha3, sqrt(sigma2.p.alpha3), lower = 0, upper = 5)/dtnorm(alpha3[i], alpha3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 5)*
              dunif(eta2[i], min = 0, max = 10)/dtnorm(eta2[i], eta2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 10)*
              dtnorm(gamma_T[i], mu.p.gamma_T, sqrt(sigma2.p.gamma_T), lower = -5, upper = 5)/dtnorm(gamma_T[i], gamma_T_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dtnorm(gamma_E[i], mu.p.gamma_E, sqrt(sigma2.p.gamma_E), lower = -5, upper = 5)/dtnorm(gamma_E[i], gamma_E_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dgamma(lambda_1[i], alpha.p.lambda_1, beta.p.lambda_1)/(pgamma(3, alpha.p.lambda_1, beta.p.lambda_1)*dtnorm(lambda_1[i], lambda_1_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_2[i], alpha.p.lambda_2, beta.p.lambda_2)/(pgamma(3, alpha.p.lambda_2, beta.p.lambda_2)*dtnorm(lambda_2[i], lambda_2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_3[i], alpha.p.lambda_3, beta.p.lambda_3)/(pgamma(3, alpha.p.lambda_3, beta.p.lambda_3)*dtnorm(lambda_3[i], lambda_3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              ifelse(length(match(d.subgroup[-1],d.all[-1]))>0,prod(dtnorm(gamma_d[match(d.subgroup[-1],d.all[-1]), i], mu.p.gamma_d, sqrt(sigma2.p.gamma_d), lower = -5, upper = 5)/dtnorm(gamma_d[match(d.subgroup[-1],d.all[-1]), i],
                                                                                                                                                                                            gamma_d_MLE[match(d.subgroup[-1],d.all[-1]), i], sqrt(sigma2_star), lower = -5, upper = 5)),1))
      accept_pr<-min(exp(accept_pr),1)
      # print(accept_pr)
      dec<-rbinom(1, size = 1, prob =accept_pr)
      if(dec==1){
        beta0 <- beta0_star;beta1 <- beta1_star;alpha0 <- alpha0_star;alpha1 <- alpha1_star;alpha2 <- alpha2_star;alpha3 <- alpha3_star
        eta2 <- eta2_star;gamma_T <- gamma_T_star;gamma_E <- gamma_E_star;lambda_1 <- lambda_1_star;lambda_2 <- lambda_2_star;lambda_3 <- lambda_3_star
        gamma_d <- gamma_d_star;z_g<-z_g_star;omega_g<-omega_g_star
      }
    }
    if(sum(omega_g==1)==3){
      beta0_star<-beta0;beta1_star<-beta1;alpha0_star<-alpha0;alpha1_star<-alpha1
      alpha2_star<-alpha2;alpha3_star<-alpha3;eta2_star<-eta2
      gamma_T_star<-gamma_T;gamma_E_star<-gamma_E;lambda_1_star<-lambda_1
      lambda_2_star<-lambda_2;lambda_3_star<-lambda_3;gamma_d_star<-gamma_d
      index_cluster<-c(2,3)
      z_g_star<-z_g;omega_g_star<-omega_g
      index_1<-min(index_cluster);index_2<-max(index_cluster)
      z_g_star[index_2]<-index_1;omega_g_star[index_2]<-0
      
      beta0_star[index_2]<-beta0[index_1]
      beta1_star[index_2] <-  beta1[index_1]
      alpha0_star[index_2] <- alpha0[index_1]
      alpha1_star[index_2] <- alpha1[index_1]
      alpha2_star[index_2] <- alpha2[index_1]
      alpha3_star[index_2] <- alpha3[index_1]
      eta2_star[index_2] <-   eta2[index_1]
      gamma_T_star[index_2] <- gamma_T_MLE[index_1]
      gamma_E_star[index_2] <- gamma_E[index_1]
      lambda_1_star[index_2] <-lambda_1[index_1]
      lambda_2_star[index_2] <- lambda_2[index_1]
      lambda_3_star[index_2] <- lambda_3[index_1]
      gamma_d_star[,index_2]<-gamma_d[,index_1]
      
      i=index_2
      id=subgroup==i
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      data.ET=list(Y_T[id],Y_E[id],d[id])
      data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
      
      accept_pr<-log.likelihood.Y_ET(c(alpha0_star[i],alpha1_star[i],alpha2_star[i],alpha3_star[i],beta0_star[i],beta1_star[i],eta2_star[i],sigma12),data=data.ET,dose.level = dose.level)+
        log.likelihood.Y_S(c(gamma_T_star[i],gamma_E_star[i],lambda_1_star[i],lambda_2_star[i],lambda_3_star[i],gamma_d_star[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log.likelihood.Y_ET(c(alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],eta2[i],sigma12),data=data.ET,dose.level = dose.level)-
        log.likelihood.Y_S(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log(4/3*dnorm(beta0[i], mu.p.beta0, sqrt(sigma2.p.beta0))/dnorm(beta0[i], beta0_MLE[i], sqrt(sigma2_star)) *
              dtnorm(beta1[i], mu.p.beta1, sqrt(sigma2.p.beta1), lower = 0, upper = Inf)/dtnorm(beta1[i], beta1_MLE[i], sqrt(sigma2_star), lower = 0, upper = Inf)*
              dnorm(alpha0[i], mu.p.alpha0, sqrt(sigma2.p.alpha0))/dnorm(alpha0[i], alpha0_MLE[i], sqrt(sigma2_star))*
              dnorm(alpha1[i], mu.p.alpha1, sqrt(sigma2.p.alpha1))/dnorm(alpha1[i], alpha1_MLE[i], sqrt(sigma2_star))*
              dtnorm(alpha2[i], mu.p.alpha2, sqrt(sigma2.p.alpha2), lower = 0, upper = 1)/dtnorm(alpha2[i], alpha2_MLE[i], sqrt(0.1), lower = 0, upper = 1)*
              dtnorm(alpha3[i], mu.p.alpha3, sqrt(sigma2.p.alpha3), lower = 0, upper = 5)/dtnorm(alpha3[i], alpha3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 5)*
              dunif(eta2[i], min = 0, max = 10)/dtnorm(eta2[i], eta2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 10)*
              dtnorm(gamma_T[i], mu.p.gamma_T, sqrt(sigma2.p.gamma_T), lower = -5, upper = 5)/dtnorm(gamma_T[i], gamma_T_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dtnorm(gamma_E[i], mu.p.gamma_E, sqrt(sigma2.p.gamma_E), lower = -5, upper = 5)/dtnorm(gamma_E[i], gamma_E_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dgamma(lambda_1[i], alpha.p.lambda_1, beta.p.lambda_1)/(pgamma(3, alpha.p.lambda_1, beta.p.lambda_1)*dtnorm(lambda_1[i], lambda_1_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_2[i], alpha.p.lambda_2, beta.p.lambda_2)/(pgamma(3, alpha.p.lambda_2, beta.p.lambda_2)*dtnorm(lambda_2[i], lambda_2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_3[i], alpha.p.lambda_3, beta.p.lambda_3)/(pgamma(3, alpha.p.lambda_3, beta.p.lambda_3)*dtnorm(lambda_3[i], lambda_3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              ifelse(length(match(d.subgroup[-1],d.all[-1]))>0,prod(dtnorm(gamma_d[match(d.subgroup[-1],d.all[-1]), i], mu.p.gamma_d, sqrt(sigma2.p.gamma_d), lower = -5, upper = 5)/dtnorm(gamma_d[match(d.subgroup[-1],d.all[-1]), i],
                                                                                                                                                                                            gamma_d_MLE[match(d.subgroup[-1],d.all[-1]), i], sqrt(sigma2_star), lower = -5, upper = 5)),1))
      accept_pr<-min(exp(accept_pr),1)
      # print(accept_pr)
      dec<-rbinom(1, size = 1, prob =accept_pr)
      if(dec==1){
        beta0 <- beta0_star;beta1 <- beta1_star;alpha0 <- alpha0_star;alpha1 <- alpha1_star;alpha2 <- alpha2_star;alpha3 <- alpha3_star
        eta2 <- eta2_star;gamma_T <- gamma_T_star;gamma_E <- gamma_E_star;lambda_1 <- lambda_1_star;lambda_2 <- lambda_2_star;lambda_3 <- lambda_3_star
        gamma_d <- gamma_d_star;z_g<-z_g_star;omega_g<-omega_g_star
      }
    }
    
    
    if(sum(omega_g==1)==1){
      beta0_star<-beta0;beta1_star<-beta1;alpha0_star<-alpha0;alpha1_star<-alpha1
      alpha2_star<-alpha2;alpha3_star<-alpha3;eta2_star<-eta2
      gamma_T_star<-gamma_T;gamma_E_star<-gamma_E;lambda_1_star<-lambda_1
      lambda_2_star<-lambda_2;lambda_3_star<-lambda_3;gamma_d_star<-gamma_d
      i<-1
      id=(subgroup==i)
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      z_g_star<-z_g;omega_g_star<-omega_g
      z_g_star[i]<-i;omega_g_star[i]<-1
      
      beta0_star[i] <- rnorm(1,beta0_MLE[i],sqrt(sigma2_star))
      beta1_star[i] <- rtnorm(1,beta1_MLE[i],sqrt(sigma2_star),lower=0,upper=Inf)
      alpha0_star[i] <- rnorm(1,alpha0_MLE[i],sqrt(sigma2_star))
      alpha1_star[i] <- rnorm(1,alpha1_MLE[i],sqrt(sigma2_star))
      alpha2_star[i] <- rtnorm(1,alpha2_MLE[i],sqrt(0.1),lower=0,upper=1)
      alpha3_star[i] <- rtnorm(1,alpha3_MLE[i],sqrt(sigma2_star),lower=0,upper=5)
      eta2_star[i] <- rtnorm(1,eta2_MLE[i],sqrt(sigma2_star),lower=0,upper=10)
      gamma_T_star[i] <- rtnorm(1,gamma_T_MLE[i],sqrt(sigma2_star),lower=-5,upper=5)
      gamma_E_star[i] <- rtnorm(1,gamma_E_MLE[i],sqrt(sigma2_star),lower=-5,upper=5)
      lambda_1_star[i] <-rtnorm(1,lambda_1_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      lambda_2_star[i] <- rtnorm(1,lambda_2_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      lambda_3_star[i] <- rtnorm(1,lambda_3_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      gamma_d_star[,i]<-0
      for(j in match(d.subgroup[-1],d.all[-1])){
        gamma_d_star[j,i]<-rtnorm(1,gamma_d_MLE[j,i],sqrt(sigma2_star),lower=-5,upper=5)
      }
      
      data.ET=list(Y_T[id],Y_E[id],d[id])
      data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
      
      accept_pr<-log.likelihood.Y_ET(c(alpha0_star[i],alpha1_star[i],alpha2_star[i],alpha3_star[i],beta0_star[i],beta1_star[i],eta2_star[i],sigma12),data=data.ET,dose.level = dose.level)+
        log.likelihood.Y_S(c(gamma_T_star[i],gamma_E_star[i],lambda_1_star[i],lambda_2_star[i],lambda_3_star[i],gamma_d_star[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log.likelihood.Y_ET(c(alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],eta2[i],sigma12),data=data.ET,dose.level = dose.level)-
        log.likelihood.Y_S(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)+
        log(3/4*dnorm(beta0_star[i], mu.p.beta0, sqrt(sigma2.p.beta0))/dnorm(beta0_star[i], beta0_MLE[i], sqrt(sigma2_star)) *
              dtnorm(beta1_star[i], mu.p.beta1, sqrt(sigma2.p.beta1), lower = 0, upper = Inf)/dtnorm(beta1_star[i], beta1_MLE[i], sqrt(sigma2_star), lower = 0, upper = Inf)*
              dnorm(alpha0_star[i], mu.p.alpha0, sqrt(sigma2.p.alpha0))/dnorm(alpha0_star[i], alpha0_MLE[i], sqrt(sigma2_star))*
              dnorm(alpha1_star[i], mu.p.alpha1, sqrt(sigma2.p.alpha1))/dnorm(alpha1_star[i], alpha1_MLE[i], sqrt(sigma2_star))*
              dtnorm(alpha2_star[i], mu.p.alpha2, sqrt(sigma2.p.alpha2), lower = 0, upper = 1)/dtnorm(alpha2_star[i], alpha2_MLE[i], sqrt(0.1), lower = 0, upper = 1)*
              dtnorm(alpha3_star[i], mu.p.alpha3, sqrt(sigma2.p.alpha3), lower = 0, upper = 5)/dtnorm(alpha3_star[i], alpha3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 5)*
              dunif(eta2_star[i], min = 0, max = 10)/dtnorm(eta2_star[i], eta2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 10)*
              dtnorm(gamma_T_star[i], mu.p.gamma_T, sqrt(sigma2.p.gamma_T), lower = -5, upper = 5)/dtnorm(gamma_T_star[i], gamma_T_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dtnorm(gamma_E_star[i], mu.p.gamma_E, sqrt(sigma2.p.gamma_E), lower = -5, upper = 5)/dtnorm(gamma_E_star[i], gamma_E_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dgamma(lambda_1_star[i], alpha.p.lambda_1, beta.p.lambda_1)/(pgamma(3, alpha.p.lambda_1, beta.p.lambda_1)*dtnorm(lambda_1_star[i], lambda_1_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_2_star[i], alpha.p.lambda_2, beta.p.lambda_2)/(pgamma(3, alpha.p.lambda_2, beta.p.lambda_2)*dtnorm(lambda_2_star[i], lambda_2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_3_star[i], alpha.p.lambda_3, beta.p.lambda_3)/(pgamma(3, alpha.p.lambda_3, beta.p.lambda_3)*dtnorm(lambda_3_star[i], lambda_3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              ifelse(length(match(d.subgroup[-1],d.all[-1]))>0,prod(dtnorm(gamma_d_star[match(d.subgroup[-1],d.all[-1]), i], mu.p.gamma_d, sqrt(sigma2.p.gamma_d), lower = -5, upper = 5)/dtnorm(gamma_d_star[match(d.subgroup[-1],d.all[-1]), i],
                                                                                                                                                                                                 gamma_d_MLE[match(d.subgroup[-1],d.all[-1]), i], sqrt(sigma2_star), lower = -5, upper = 5)),1))
      
      accept_pr<-min(exp(accept_pr),1)
      # print(accept_pr)
      dec<-rbinom(1, size = 1, prob =accept_pr)
      if(dec==1){
        beta0 <- beta0_star;beta1 <- beta1_star;alpha0 <- alpha0_star;alpha1 <- alpha1_star;alpha2 <- alpha2_star;alpha3 <- alpha3_star
        eta2 <- eta2_star;gamma_T <- gamma_T_star;gamma_E <- gamma_E_star;lambda_1 <- lambda_1_star;lambda_2 <- lambda_2_star;lambda_3 <- lambda_3_star
        gamma_d <- gamma_d_star;z_g<-z_g_star;omega_g<-omega_g_star
      }
    }
    if(sum(omega_g==1)==1){
      beta0_star<-beta0;beta1_star<-beta1;alpha0_star<-alpha0;alpha1_star<-alpha1
      alpha2_star<-alpha2;alpha3_star<-alpha3;eta2_star<-eta2
      gamma_T_star<-gamma_T;gamma_E_star<-gamma_E;lambda_1_star<-lambda_1
      lambda_2_star<-lambda_2;lambda_3_star<-lambda_3;gamma_d_star<-gamma_d
      i<-2
      id=(subgroup==i)
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      z_g_star<-z_g;omega_g_star<-omega_g
      z_g_star[i]<-i;omega_g_star[i]<-1
      
      beta0_star[i] <- rnorm(1,beta0_MLE[i],sqrt(sigma2_star))
      beta1_star[i] <- rtnorm(1,beta1_MLE[i],sqrt(sigma2_star),lower=0,upper=Inf)
      alpha0_star[i] <- rnorm(1,alpha0_MLE[i],sqrt(sigma2_star))
      alpha1_star[i] <- rnorm(1,alpha1_MLE[i],sqrt(sigma2_star))
      alpha2_star[i] <- rtnorm(1,alpha2_MLE[i],sqrt(0.1),lower=0,upper=1)
      alpha3_star[i] <- rtnorm(1,alpha3_MLE[i],sqrt(sigma2_star),lower=0,upper=5)
      eta2_star[i] <- rtnorm(1,eta2_MLE[i],sqrt(sigma2_star),lower=0,upper=10)
      gamma_T_star[i] <- rtnorm(1,gamma_T_MLE[i],sqrt(sigma2_star),lower=-5,upper=5)
      gamma_E_star[i] <- rtnorm(1,gamma_E_MLE[i],sqrt(sigma2_star),lower=-5,upper=5)
      lambda_1_star[i] <-rtnorm(1,lambda_1_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      lambda_2_star[i] <- rtnorm(1,lambda_2_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      lambda_3_star[i] <- rtnorm(1,lambda_3_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      gamma_d_star[,i]<-0
      for(j in match(d.subgroup[-1],d.all[-1])){
        gamma_d_star[j,i]<-rtnorm(1,gamma_d_MLE[j,i],sqrt(sigma2_star),lower=-5,upper=5)
      }
      
      data.ET=list(Y_T[id],Y_E[id],d[id])
      data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
      
      accept_pr<-log.likelihood.Y_ET(c(alpha0_star[i],alpha1_star[i],alpha2_star[i],alpha3_star[i],beta0_star[i],beta1_star[i],eta2_star[i],sigma12),data=data.ET,dose.level = dose.level)+
        log.likelihood.Y_S(c(gamma_T_star[i],gamma_E_star[i],lambda_1_star[i],lambda_2_star[i],lambda_3_star[i],gamma_d_star[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log.likelihood.Y_ET(c(alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],eta2[i],sigma12),data=data.ET,dose.level = dose.level)-
        log.likelihood.Y_S(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)+
        log(3/4*dnorm(beta0_star[i], mu.p.beta0, sqrt(sigma2.p.beta0))/dnorm(beta0_star[i], beta0_MLE[i], sqrt(sigma2_star)) *
              dtnorm(beta1_star[i], mu.p.beta1, sqrt(sigma2.p.beta1), lower = 0, upper = Inf)/dtnorm(beta1_star[i], beta1_MLE[i], sqrt(sigma2_star), lower = 0, upper = Inf)*
              dnorm(alpha0_star[i], mu.p.alpha0, sqrt(sigma2.p.alpha0))/dnorm(alpha0_star[i], alpha0_MLE[i], sqrt(sigma2_star))*
              dnorm(alpha1_star[i], mu.p.alpha1, sqrt(sigma2.p.alpha1))/dnorm(alpha1_star[i], alpha1_MLE[i], sqrt(sigma2_star))*
              dtnorm(alpha2_star[i], mu.p.alpha2, sqrt(sigma2.p.alpha2), lower = 0, upper = 1)/dtnorm(alpha2_star[i], alpha2_MLE[i], sqrt(0.1), lower = 0, upper = 1)*
              dtnorm(alpha3_star[i], mu.p.alpha3, sqrt(sigma2.p.alpha3), lower = 0, upper = 5)/dtnorm(alpha3_star[i], alpha3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 5)*
              dunif(eta2_star[i], min = 0, max = 10)/dtnorm(eta2_star[i], eta2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 10)*
              dtnorm(gamma_T_star[i], mu.p.gamma_T, sqrt(sigma2.p.gamma_T), lower = -5, upper = 5)/dtnorm(gamma_T_star[i], gamma_T_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dtnorm(gamma_E_star[i], mu.p.gamma_E, sqrt(sigma2.p.gamma_E), lower = -5, upper = 5)/dtnorm(gamma_E_star[i], gamma_E_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dgamma(lambda_1_star[i], alpha.p.lambda_1, beta.p.lambda_1)/(pgamma(3, alpha.p.lambda_1, beta.p.lambda_1)*dtnorm(lambda_1_star[i], lambda_1_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_2_star[i], alpha.p.lambda_2, beta.p.lambda_2)/(pgamma(3, alpha.p.lambda_2, beta.p.lambda_2)*dtnorm(lambda_2_star[i], lambda_2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_3_star[i], alpha.p.lambda_3, beta.p.lambda_3)/(pgamma(3, alpha.p.lambda_3, beta.p.lambda_3)*dtnorm(lambda_3_star[i], lambda_3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              ifelse(length(match(d.subgroup[-1],d.all[-1]))>0,prod(dtnorm(gamma_d_star[match(d.subgroup[-1],d.all[-1]), i], mu.p.gamma_d, sqrt(sigma2.p.gamma_d), lower = -5, upper = 5)/dtnorm(gamma_d_star[match(d.subgroup[-1],d.all[-1]), i],
                                                                                                                                                                                                 gamma_d_MLE[match(d.subgroup[-1],d.all[-1]), i], sqrt(sigma2_star), lower = -5, upper = 5)),1))
      
      accept_pr<-min(exp(accept_pr),1)
      # print(accept_pr)
      dec<-rbinom(1, size = 1, prob =accept_pr)
      if(dec==1){
        beta0 <- beta0_star;beta1 <- beta1_star;alpha0 <- alpha0_star;alpha1 <- alpha1_star;alpha2 <- alpha2_star;alpha3 <- alpha3_star
        eta2 <- eta2_star;gamma_T <- gamma_T_star;gamma_E <- gamma_E_star;lambda_1 <- lambda_1_star;lambda_2 <- lambda_2_star;lambda_3 <- lambda_3_star
        gamma_d <- gamma_d_star;z_g<-z_g_star;omega_g<-omega_g_star
      }
    }
    if(sum(omega_g==1)==1){
      beta0_star<-beta0;beta1_star<-beta1;alpha0_star<-alpha0;alpha1_star<-alpha1
      alpha2_star<-alpha2;alpha3_star<-alpha3;eta2_star<-eta2
      gamma_T_star<-gamma_T;gamma_E_star<-gamma_E;lambda_1_star<-lambda_1
      lambda_2_star<-lambda_2;lambda_3_star<-lambda_3;gamma_d_star<-gamma_d
      i<-3
      id=(subgroup==i)
      d.subgroup<-sort(unique(d[id&Y_E!=0]))
      z_g_star<-z_g;omega_g_star<-omega_g
      z_g_star[i]<-i;omega_g_star[i]<-1
      
      beta0_star[i] <- rnorm(1,beta0_MLE[i],sqrt(sigma2_star))
      beta1_star[i] <- rtnorm(1,beta1_MLE[i],sqrt(sigma2_star),lower=0,upper=Inf)
      alpha0_star[i] <- rnorm(1,alpha0_MLE[i],sqrt(sigma2_star))
      alpha1_star[i] <- rnorm(1,alpha1_MLE[i],sqrt(sigma2_star))
      alpha2_star[i] <- rtnorm(1,alpha2_MLE[i],sqrt(0.1),lower=0,upper=1)
      alpha3_star[i] <- rtnorm(1,alpha3_MLE[i],sqrt(sigma2_star),lower=0,upper=5)
      eta2_star[i] <- rtnorm(1,eta2_MLE[i],sqrt(sigma2_star),lower=0,upper=10)
      gamma_T_star[i] <- rtnorm(1,gamma_T_MLE[i],sqrt(sigma2_star),lower=-5,upper=5)
      gamma_E_star[i] <- rtnorm(1,gamma_E_MLE[i],sqrt(sigma2_star),lower=-5,upper=5)
      lambda_1_star[i] <-rtnorm(1,lambda_1_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      lambda_2_star[i] <- rtnorm(1,lambda_2_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      lambda_3_star[i] <- rtnorm(1,lambda_3_MLE[i],sqrt(sigma2_star),lower=0,upper=3)
      gamma_d_star[,i]<-0
      for(j in match(d.subgroup[-1],d.all[-1])){
        gamma_d_star[j,i]<-rtnorm(1,gamma_d_MLE[j,i],sqrt(sigma2_star),lower=-5,upper=5)
      }
      
      data.ET=list(Y_T[id],Y_E[id],d[id])
      data.S=list(Y_T[id],Y_E[id],d[id],Y_S.obs[id],status[id])
      
      accept_pr<-log.likelihood.Y_ET(c(alpha0_star[i],alpha1_star[i],alpha2_star[i],alpha3_star[i],beta0_star[i],beta1_star[i],eta2_star[i],sigma12),data=data.ET,dose.level = dose.level)+
        log.likelihood.Y_S(c(gamma_T_star[i],gamma_E_star[i],lambda_1_star[i],lambda_2_star[i],lambda_3_star[i],gamma_d_star[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)-
        log.likelihood.Y_ET(c(alpha0[i],alpha1[i],alpha2[i],alpha3[i],beta0[i],beta1[i],eta2[i],sigma12),data=data.ET,dose.level = dose.level)-
        log.likelihood.Y_S(c(gamma_T[i],gamma_E[i],lambda_1[i],lambda_2[i],lambda_3[i],gamma_d[match(d.subgroup[-1],d.all[-1]),i]),data=data.S,dose.level = dose.level)+
        log(3/4*dnorm(beta0_star[i], mu.p.beta0, sqrt(sigma2.p.beta0))/dnorm(beta0_star[i], beta0_MLE[i], sqrt(sigma2_star)) *
              dtnorm(beta1_star[i], mu.p.beta1, sqrt(sigma2.p.beta1), lower = 0, upper = Inf)/dtnorm(beta1_star[i], beta1_MLE[i], sqrt(sigma2_star), lower = 0, upper = Inf)*
              dnorm(alpha0_star[i], mu.p.alpha0, sqrt(sigma2.p.alpha0))/dnorm(alpha0_star[i], alpha0_MLE[i], sqrt(sigma2_star))*
              dnorm(alpha1_star[i], mu.p.alpha1, sqrt(sigma2.p.alpha1))/dnorm(alpha1_star[i], alpha1_MLE[i], sqrt(sigma2_star))*
              dtnorm(alpha2_star[i], mu.p.alpha2, sqrt(sigma2.p.alpha2), lower = 0, upper = 1)/dtnorm(alpha2_star[i], alpha2_MLE[i], sqrt(0.1), lower = 0, upper = 1)*
              dtnorm(alpha3_star[i], mu.p.alpha3, sqrt(sigma2.p.alpha3), lower = 0, upper = 5)/dtnorm(alpha3_star[i], alpha3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 5)*
              dunif(eta2_star[i], min = 0, max = 10)/dtnorm(eta2_star[i], eta2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 10)*
              dtnorm(gamma_T_star[i], mu.p.gamma_T, sqrt(sigma2.p.gamma_T), lower = -5, upper = 5)/dtnorm(gamma_T_star[i], gamma_T_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dtnorm(gamma_E_star[i], mu.p.gamma_E, sqrt(sigma2.p.gamma_E), lower = -5, upper = 5)/dtnorm(gamma_E_star[i], gamma_E_MLE[i], sqrt(sigma2_star), lower = -5, upper = 5)*
              dgamma(lambda_1_star[i], alpha.p.lambda_1, beta.p.lambda_1)/(pgamma(3, alpha.p.lambda_1, beta.p.lambda_1)*dtnorm(lambda_1_star[i], lambda_1_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_2_star[i], alpha.p.lambda_2, beta.p.lambda_2)/(pgamma(3, alpha.p.lambda_2, beta.p.lambda_2)*dtnorm(lambda_2_star[i], lambda_2_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              dgamma(lambda_3_star[i], alpha.p.lambda_3, beta.p.lambda_3)/(pgamma(3, alpha.p.lambda_3, beta.p.lambda_3)*dtnorm(lambda_3_star[i], lambda_3_MLE[i], sqrt(sigma2_star), lower = 0, upper = 3))*
              ifelse(length(match(d.subgroup[-1],d.all[-1]))>0,prod(dtnorm(gamma_d_star[match(d.subgroup[-1],d.all[-1]), i], mu.p.gamma_d, sqrt(sigma2.p.gamma_d), lower = -5, upper = 5)/dtnorm(gamma_d_star[match(d.subgroup[-1],d.all[-1]), i],
                                                                                                                                                                                                 gamma_d_MLE[match(d.subgroup[-1],d.all[-1]), i], sqrt(sigma2_star), lower = -5, upper = 5)),1))
      
      accept_pr<-min(exp(accept_pr),1)
      # print(accept_pr)
      dec<-rbinom(1, size = 1, prob =accept_pr)
      if(dec==1){
        beta0 <- beta0_star;beta1 <- beta1_star;alpha0 <- alpha0_star;alpha1 <- alpha1_star;alpha2 <- alpha2_star;alpha3 <- alpha3_star
        eta2 <- eta2_star;gamma_T <- gamma_T_star;gamma_E <- gamma_E_star;lambda_1 <- lambda_1_star;lambda_2 <- lambda_2_star;lambda_3 <- lambda_3_star
        gamma_d <- gamma_d_star;z_g<-z_g_star;omega_g<-omega_g_star
      }
    }
    # print(z_g)
  }
  
  ind <- thin*seq((N.burnin+1),N.post)
  # str(pi_t[,,,ind,])
  
  Uti_mcmc=function(x){
    re=sum(x*score)
    return(re) 
  }  
  uti_mcmc<-apply(pi_t[,,,ind,], MARGIN = c(3,4,5),Uti_mcmc)
  Uti_est<-t(apply(uti_mcmc,MARGIN = c(1,3),mean))
  Sur.p_est<-t(apply(sur.p_t[,ind,],MARGIN=c(1,3),mean))
  Sur.p.lb_est<-t(apply(sur.p_t[,ind,]>phi_S.lb,MARGIN=c(1,3),mean))
  return(list(Uti_est=Uti_est,Sur.p_est=Sur.p_est,Sur.p.lb_est=Sur.p.lb_est))
}

#########################################################################################################
##### PGen12.I implements stage I of the PGen I-II design ###############################################
#########################################################################################################
PGen12.I=function(Y_T.1,Y_E.1,Y_S.1.obs,status.1,enter.1,d.1,g.1,dose.level,score=matrix(c(20,50,100,0,30,60),2,3, byrow=TRUE),A1=1/3,ncohort1=10,cohortsize1=3,targetT=0.3,targetE=0.5,cutoff.eli.T=0.9,cutoff.eli.E=0.9){
  ##### Y_T.1: Individual toxicity of patients in stage I
  ##### Y_E.1: Individual efficacy of patients in stage I
  ##### Y_S.1.obs: Individual observed Y_S^o of patients in stage I
  ##### status.1: Individual censor indicator of patients in stage I (1 for no censoring, 0 for censoring)
  ##### enter.1: Individual enter time in stage I
  ##### d.1: Individual dose allocation of patients in stage I
  ##### g.1: Group of patients in stage I
  ##### dose.level: dose levels
  ##### score: the utility function used
  ##### A1: accrual rate per patient, 1 patient per A1 month in stage I
  ##### ncohort1: number of cohorts in stage I
  ##### cohortsize1: cohortsize in stage I
  ##### targetT: highest acceptable toxicity rate for dose-finding
  ##### targetE: lowest acceptable efficacy rate for dose-finding
  ##### cutoff.eli.T: cut-off probability for toxicity admissible
  ##### cutoff.eli.E: cut-off probability for efficacy admissible
  
  ndose<-length(dose.level)
  yT<-yE<-rep(0, ndose)      ## number of DLT/efficacy at each dose level
  for(i in 1:ndose){
    yT[i]=sum(Y_T.1[d.1==i]==1)
    yE[i]=sum(Y_E.1[d.1==i]==2)
  }
  y00<-y01<-y02<-y10<-y11<-y12<-rep(0,ndose) ## number of different outcomes at each dose level
  for(i in 1:ndose){
    y00[i]=sum(Y_T.1[d.1==i]==0&Y_E.1[d.1==i]==0)
    y01[i]=sum(Y_T.1[d.1==i]==0&Y_E.1[d.1==i]==1)
    y02[i]=sum(Y_T.1[d.1==i]==0&Y_E.1[d.1==i]==2)
    y10[i]=sum(Y_T.1[d.1==i]==1&Y_E.1[d.1==i]==0)
    y11[i]=sum(Y_T.1[d.1==i]==1&Y_E.1[d.1==i]==1)
    y12[i]=sum(Y_T.1[d.1==i]==1&Y_E.1[d.1==i]==2)
  }
  n.1<-rep(0, ndose)        ## number of patients treated at each dose level
  for(i in 1:ndose){
    n.1[i]=sum(d.1==i)
  }
  y.1<-rbind(y00,y10,y01,y11,y02,y12)
  earlystop=0
  if(length(d.1)==ncohort1*cohortsize1) {
    acset.1=adm(y.1,targetT,targetE,cutoff.eli.T,cutoff.eli.E)
    print("This is the end of stage I. Please use PGen12.II() for second stage of dose-finding.")
    cat("The stage 1 acceptable dose set is:",formatC(acset.1, digits=0, format="f"), sep=" ", "\n")
  } else{
    uu=uti_benchmark_PGen12(score,targetT,targetE)
    n.earlystop=ncohort1*cohortsize1
    boundary.temp=get.boundary(targetT, targetE, ncohort1, cohortsize1,cutoff.eli=cutoff.eli.T, cutoff.eli.E = cutoff.eli.E)
    b.e=boundary.temp[4,]   # escalation boundary
    b.d=boundary.temp[3,]   # deescalation boundary
    b.elim=boundary.temp[2,]  # elimination boundary
    b.elimE=boundary.temp[5,]
    A1C=A1*cohortsize1 ## 1 cohort per A1C month in stage I
    acset.1<-NULL
    elimi = rep(0, ndose);    ## indicate whether doses are eliminated due to toxicity/efficacy
    safe = 0
    posH<-rep(1-uu/100,ndose)
    N1=6;N2=9
    
    d<-d.1[length(d.1)]##Current dose
    nc<-n.1[d]/cohortsize1
    # determine whether current dose level is overly toxic
    if(!is.na(b.elim[nc]))
    {
      if(yT[d]>=b.elim[nc]) 
      {      
        elimi[d:ndose]=1;
        if(d==1) {earlystop=1; break;} 
      }
    }
    
    if(!is.na(b.elimE[nc]))
    {
      if(yE[d]<=b.elimE[nc]) 
      {      
        elimi[d]=1;
      }
    }
    
    if(sum(elimi==1)==ndose) {earlystop=1; break;} 
    
    u_curr<-(score[1,1]*y00[d]+score[1,2]*y01[d]+score[1,3]*y02[d]+score[2,1]*y10[d]+score[2,2]*y11[d]+score[2,3]*y12[d])/100
    
    posH[d] = 1-pbeta(uu/100,1+u_curr,n.1[d]-u_curr+1)
    
    posH <- posH*(1-elimi);
    if(n.1[d]>=N1){safe=1} else{safe=0}
    if(n.1[d]>=n.earlystop){break}
    
    if (yT[d]>=b.d[nc] && d!=1) {
      if(sum(elimi[1:(d-1)]==0)>0){d_opt=max(which(elimi[1:(d-1)]==0))} else {
        if(elimi[d]==1){earlystop=1;break} else{d_opt=d}}
    } else if (yT[d]>=b.d[nc] && d==1) {if(elimi[d]==0){d_opt=d} else{earlystop=1;break}
    } else{
      admi_set=d;
      if(d>1){
        if(sum(elimi[1:(d-1)]==0)>0){admi_set<-c(admi_set,max(which(elimi[1:(d-1)]==0)))} 
      }
      if(d<ndose){
        if(safe==0){
          if(sum(elimi[(d+1):ndose]==0)>0){admi_set<-c(admi_set,d+min(which(elimi[(d+1):ndose]==0)))}
        } else {
          if(yT[d]<=b.e[nc] & sum(elimi[(d+1):ndose]==0)>0){admi_set<-c(admi_set,d+min(which(elimi[(d+1):ndose]==0)))}
        }
      }
      #if(length(admi_set)>1 & eff_cut>0.5 & n.1[d]>=9){admi_set<-admi_set[-1]}
      temp.posH<-posH[admi_set]+runif(length(admi_set))*(10^-15)
      d_opt=admi_set[which.max(temp.posH)]
    }	
    
    if (elimi[d_opt]==1) {earlystop=1; break} 
    if (sum(elimi)==ndose) {earlystop=1; break}
    
    if (d<ndose){
      if(sum(elimi[(d+1):ndose]==0)>0){
        d.temp=d+min(which(elimi[(d+1):ndose]==0))
        if(n.1[d]>=N2 & n.1[min(d.temp,ndose)]==0 & yT[d]<b.d[n.1[d]/cohortsize1]){d_opt<-d.temp} 
      }  
    }
    if(earlystop==1) print("No dose is acceptable, and the trial is stopped.") else cat("The dose for next cohort of patients: ", formatC(d_opt, digits=0, format="f"), sep="  ", "\n")
  }
  
}

#########################################################################################################
##### PGen12.II implements stage II of the PGen I-II design #############################################
#########################################################################################################
PGen12.II=function(acset.1,Y_T.12,Y_E.12,Y_S.12.obs,status.12,enter.12,d.12,g.12,dose.level,score=matrix(c(20,50,100,0,30,60),2,3, byrow=TRUE),A2=1/3,ncohort2=30,cohortsize2=3,targetT=0.3,targetE=0.5,
                   cutoff.eli.T=0.9,cutoff.eli.E=0.9,phi_S.lb=0.4,ps.lb=0.1,rho=0.7){
  ##### acset.1: the stage 1 acceptable dose set
  ##### Y_T.12: Individual toxicity of patients in stages I and II
  ##### Y_E.12: Individual efficacy of patients in stages I and II
  ##### Y_S.12.obs: Individual observed Y_S^o of patients in stages I and II
  ##### status.12: Individual censor indicator of patients in stages I and II (1 for no censoring, 0 for censoring)
  ##### enter.12: Individual enter time in stages I and II
  ##### d.12: Individual dose allocation of patients in stages I and II
  ##### g.12: Group of patients in stages I and II
  ##### dose.level: dose levels
  ##### score: the utility function used
  ##### A2: accrual rate per patient, 1 patient per A2 month in stageII
  ##### ncohort2: number of cohorts in stageII
  ##### cohortsize2: cohortsize in stage II
  ##### targetT: highest acceptable toxicity rate for dose-finding
  ##### targetE: lowest acceptable efficacy rate for dose-finding
  ##### cutoff.eli.T: cut-off probability for toxicity admissible
  ##### cutoff.eli.E: cut-off probability for efficacy admissible
  ##### phi_S.lb: parameter \underline{\phi}_S in the definition of stage 2 subgroup-specific acceptable dose set
  ##### ps.lb: parameter p_{S,L} in the definition of stage 2 subgroup-specific acceptable dose set
  ##### rho: parameter \rho in the definition of stage 2 candidate dose set
  
  datate.12<-matrix(0,6,length(Y_T.12)) ## Individual toxicity-efficacy distribution in stages I and II
  for(i in 1:length(Y_T.12)){
    if(Y_T.12[i]==0&Y_E.12[i]==0)  datate.12[1,i]=1
    if(Y_T.12[i]==1&Y_E.12[i]==0)  datate.12[2,i]=1
    if(Y_T.12[i]==0&Y_E.12[i]==1)  datate.12[3,i]=1
    if(Y_T.12[i]==1&Y_E.12[i]==1)  datate.12[4,i]=1
    if(Y_T.12[i]==0&Y_E.12[i]==2)  datate.12[5,i]=1
    if(Y_T.12[i]==1&Y_E.12[i]==2)  datate.12[6,i]=1
  }
  acset.g.12<-list(integer(0),integer(0),integer(0))## Acceptable set for each subgroup in stages I and II
  id.acset.g.12<-rep(0,3)
  for(i in 1:3){
    y.12<-matrix(0,6,length(acset.1))
    for(j in 1:length(acset.1)){
      index<-g.12==i&d.12==acset.1[j]
      if(sum(index)==0) y.12[,j]=rep(0,6) else
        if(sum(index)==1) y.12[,j]=datate.12[,index]
        else y.12[,j]<-apply(datate.12[,g.12==i&d.12==acset.1[j]],1,sum)
    }
    acset.g.12[[i]]<-acset.1[adm(y.12,targetT,targetE,cutoff.eli.T,cutoff.eli.E)]
    id.acset.g.12[i]<-(length(acset.g.12[[i]])!=0)
  }
  if(sum(id.acset.g.12)==0) {
    print("No dose is acceptable, and the trial is stopped.") 
  }else{
    data.12=list(d=d.12,dose=dose.level[d.12],Y_T=Y_T.12,Y_E=Y_E.12,Y_S.obs=Y_S.12.obs,status=status.12,g=g.12)
    res<-rjmcmc.ETS(data.12,dose.level,score,phi_S.lb)
    for(i in 1:3){
      if(length(acset.g.12[[i]])!=0) acset.g.12[[i]]<-acset.g.12[[i]][which(res[[1]][i,acset.g.12[[i]]]>rho*max(res[[1]][i,acset.g.12[[i]]])&res[[3]][i,acset.g.12[[i]]]>ps.lb)]
      if(length(acset.g.12[[i]])==0) earlystop[i]=1
    }
    print("This is the end of stage II. Please use PGen12.III() for third stage of dose-finding.")
    cat("The stage 2 candidate dose sets are:","\nsubgroup 1: ",formatC(acset.g.12[[1]], digits=0, format="f"),
        "\nsubgroup 2: ",formatC(acset.g.12[[2]], digits=0, format="f"),"\nsubgroup 3: ",formatC(acset.g.12[[3]], digits=0, format="f"), sep=" ", "\n")
  }
  
}







#########################################################################################################
##### PGen12.III implements stage III of the PGen I-II design ###########################################
#########################################################################################################
PGen12.III=function(acset.12,Y_T.123,Y_E.123,Y_S.123.obs,status.123,enter.123,d.123,g.123,dose.level,score=matrix(c(20,50,100,0,30,60),2,3, byrow=TRUE),A3=1/3,ncohort3=10,cohortsize3=3,targetT=0.3,targetE=0.5,
                    cutoff.eli.T=0.9,cutoff.eli.E=0.9,phi_S.lb=0.4,ps.lb=0.5){
  ##### acset.12: the stage 2 candidate dose sets
  ##### Y_T.123: Individual toxicity of patients in stages I, II and III
  ##### Y_E.123: Individual efficacy of patients in stages I, II and III
  ##### Y_S.123.obs: Individual observed Y_S^o of patients in stages I, II and III
  ##### status.123: Individual censor indicator of patients in stages I, II and III (1 for no censoring, 0 for censoring)
  ##### enter.123: Individual enter time in stages I, II and III
  ##### d.123: Individual dose allocation of patients in stages I, II and III
  ##### g.123: Group of patients in stages I, II and III
  ##### dose.level: dose levels
  ##### score: the utility function used
  ##### A3: accrual rate per patient, 1 patient per A2 month in stageIII
  ##### ncohort3: number of cohorts in stageIII
  ##### cohortsize3: cohortsize in stage III
  ##### targetT: highest acceptable toxicity rate for dose-finding
  ##### targetE: lowest acceptable efficacy rate for dose-finding
  ##### cutoff.eli.T: cut-off probability for toxicity admissible
  ##### cutoff.eli.E: cut-off probability for efficacy admissible
  ##### phi_S.lb: parameter \underline{\phi}_S in the definition of stage 2 subgroup-specific acceptable dose set
  ##### ps.lb: parameter p_{S,L} in the definition of stage 3 updated acceptable dose set
  ##### rho: parameter \rho in the definition of stage 3 updated candidate dose set
  
  datate.123<-matrix(0,6,length(Y_T.123)) ## Individual toxicity-efficacy distribution in stages I, II and III
  for(i in 1:length(Y_T.123)){
    if(Y_T.123[i]==0&Y_E.123[i]==0)  datate.123[1,i]=1
    if(Y_T.123[i]==1&Y_E.123[i]==0)  datate.123[2,i]=1
    if(Y_T.123[i]==0&Y_E.123[i]==1)  datate.123[3,i]=1
    if(Y_T.123[i]==1&Y_E.123[i]==1)  datate.123[4,i]=1
    if(Y_T.123[i]==0&Y_E.123[i]==2)  datate.123[5,i]=1
    if(Y_T.123[i]==1&Y_E.123[i]==2)  datate.123[6,i]=1
  }
  acset.g.123<-list(integer(0),integer(0),integer(0)) ## Acceptable set for each subgroup
  id.acset.g.123<-rep(0,3)
  for(i in 1:3){
    if(length(acset.12[[i]])>0){
      y.123<-matrix(0,6,length(acset.12[[i]]))
      for(j in 1:length(acset.12[[i]])){
        if(length(acset.12[[i]][j])==0) index<-NULL else index<-g.123==i&d.123==acset.12[[i]][j]
        if(sum(index)==0) y.123[,j]=rep(0,6) else 
          if(sum(index)==1) y.123[,j]=datate.123[,index] 
          else y.123[,j]<-apply(datate.123[,g.123==i&d.123==acset.12[[i]][j]],1,sum)
      }
      acset.g.123[[i]]=acset.12[[i]][adm_tox(y.123,targetT,cutoff.eli.T)]
    } else acset.g.123[[i]]=acset.12[[i]]
    id.acset.g.123[i]<-(length(acset.g.123[[i]])!=0)
  }
  if(sum(id.acset.g.123)==0) {
    print("No dose is acceptable, and the trial is stopped.") 
  }else{
    data.123=list(d=d.123,dose=dose.level[d.123],Y_T=Y_T.123,Y_E=Y_E.123,Y_S.ob=Y_S.123.obs,status=status.123,g=g.123)
    res<-rjmcmc.ETS(data.123,dose.level,score,phi_S.lb)
    for(i in 1:3){
      if(length(acset.g.123[[i]])!=0) {acset.g.123[[i]]<-acset.g.123[[i]][which(res[[3]][i,acset.g.123[[i]]]>ps.lb)];acset.g.123[[i]]<-acset.g.123[[i]][which.max(res[[2]][i,acset.g.123[[i]]])]}
      if(length(acset.g.123[[i]])==0) earlystop[i]=1
    }
    print("This is the end of PGen12 design.")
    cat("The subgroup-specific optimal dose is:","\nsubgroup 1: ",formatC(acset.g.123[[1]], digits=0, format="f"),
        "\nsubgroup 2: ",formatC(acset.g.123[[2]], digits=0, format="f"),"\nsubgroup 3: ",formatC(acset.g.123[[3]], digits=0, format="f"), sep=" ", "\n")
  }
  
}


##############*********************************************************************################################
##############*********************************************************************################################
##############*** Example: PGen I-II design implementation as shown in Table 3***##################################
##############*********************************************************************################################
##############*********************************************************************################################


#########################################################################################################
##### Implement stage I of the PGen I-II design ###############################################
#########################################################################################################
set.seed(1e3)
ndose<-4##number of dose
t1=1   ## t1: short-term outcome follow-up time
t2=6  ## t2: long-term outcome follow-up time
dose.level<-as.vector(scale(log(5*c(10^6,10^7,10^8,10^9)),center=FALSE,scale=TRUE))
dose.level<-dose.level/max(dose.level)## dose levels

Y_T.1 <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0) ## Y_T.1: Individual toxicity of patients in stage I
Y_E.1 <- c(1, 1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 1, 0, 2, 1, 2) ## Y_E.1: Individual efficacy of patients in stage I
Y_S.1.obs <- c(0.84, 0.49, 5.00, 2.74, 5.00, 1.33, 0.24, 4.06, 5.00, 3.19, 4.00, 0.58, 
               3.00, 0.69, 1.21, 0.14, 1.06, 2.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00) ## Y_S.1.obs: Individual observed Y_S^o of patients in stage I
status.1 <- c(1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0) ## status.1: Individual censor indicator of patients in stage I (1 for no censoring, 0 for censoring)
enter.1 <- rep(0:7, each = 3) ## enter.1: Individual enter time in stage I
d.1 <- c(rep(1, 3), rep(2, 3), rep(3, 9), rep(4, 6), rep(3, 3)) ## d.1: Individual dose allocation of patients in stage I
g.1 <- c(rep(3, 3), rep(2, 3), rep(1, 6), rep(3, 6), rep(2, 6)) ## g.1: Group of patients in stage I

PGen12.I(Y_T.1,Y_E.1,Y_S.1.obs,status.1,enter.1,d.1,g.1,dose.level)


#########################################################################################################
##### Implement stage II of the PGen I-II design ###############################################
#########################################################################################################


t_l=3## t_l: the assessment time between the phase II and III
acset.1<-c(1,2,3,4) ## acset.1: the acceptable dose set of Stage I
Y_T.12 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0, 0, 0, 0, 0, 1) ## Y_T.12: Individual toxicity of patients in stages I and II
Y_E.12 <- c(1, 0, 2, 2, 0, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 1, 2, 
            2, 2, 2, 1, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 
            2, 2, 2, 0, 1, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 0, 2, 2, 2, 
            1, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 1, 2, 1, 2, 1, 2, 0, 1, 2, 
            2, 1, 1, 1, 2, 2, 2, 2) ## Y_E.12: Individual efficacy of patients in stages I and II
Y_S.12.obs <- c(1.28, 0.00, 1.10, 5.00, 0.00, 5.00, 5.00, 0.37, 5.00, 0.06, 5.00, 5.00, 
                0.51, 5.00, 5.00, 5.00, 4.96, 5.00, 5.00, 5.00, 0.09, 0.81, 0.18, 5.00, 
                5.00, 5.00, 0.42, 5.00, 5.00, 3.70, 3.26, 0.17, 2.07, 1.90, 1.94, 0.73, 
                5.00, 2.45, 0.80, 5.00, 1.21, 5.00, 1.20, 0.74, 5.00, 5.00, 5.00, 1.88, 
                3.63, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 4.61, 5.00, 4.20, 1.32, 2.69, 
                5.00, 0.00, 0.70, 5.00, 2.36, 0.87, 3.30, 5.00, 0.34, 5.00, 1.58, 3.78, 
                0.99, 3.06, 0.05, 5.00, 0.43, 1.95, 5.00, 5.00, 5.00, 5.00, 0.00, 5.00, 
                0.74, 2.95, 2.21, 0.00, 2.05, 4.54, 5.00, 5.00, 1.86, 3.52, 5.00, 5.00, 
                4.12, 0.77, 5.00, 5.00, 0.00, 1.34, 5.00, 3.08, 0.21, 3.91, 5.00, 5.00, 
                5.00, 0.00, 0.04, 0.91, 4.00, 3.84, 0.93, 2.32, 3.00, 2.00, 0.53, 1.80) ## Y_S.12.obs: Individual observed Y_S^o of patients in stages I and II
status.12 <- c(1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 
               0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0,
               0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1,
               1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 
               1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1) ## status.12: Individual censor indicator of patients in stages I and II (1 for no censoring, 0 for censoring)
enter.12 <- rep(0:39, each = 3) ## enter.12: Individual enter time in stages I and II
d.12 <- c(1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
          2, 3, 1, 4, 1, 4, 3, 1, 2, 2, 4, 4, 4, 3, 4, 2, 4, 1, 2, 4, 2, 3, 4, 4, 2, 3, 1, 3, 2, 4, 
          3, 2, 1, 2, 1, 2, 1, 3, 1, 4, 3, 3, 4, 4, 1, 3, 2, 3, 2, 2, 3, 1, 1, 2, 4, 3, 2, 1, 2, 3, 
          2, 2, 1, 3, 4, 1, 3, 1, 4, 1, 4, 3, 4, 2, 1, 2, 4, 3, 2, 1, 3, 2, 4, 3, 4, 2, 3, 4, 1, 4) ## d.12: Individual dose allocation of patients in stages I and II
g.12 <- c(2, 2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 1, 1, 1, 3, 3, 3, 
          2, 3, 3, 1, 2, 1, 3, 2, 2, 3, 1, 2, 3, 3, 2, 3, 3, 2, 2, 3, 3, 2, 1, 2, 3, 3, 3, 1, 2, 3, 
          3, 3, 1, 3, 3, 1, 3, 1, 2, 2, 2, 2, 1, 3, 3, 3, 2, 2, 3, 3, 2, 3, 1, 3, 2, 2, 2, 3, 1, 2, 
          1, 1, 2, 1, 1, 1, 3, 3, 1, 2, 2, 1, 3, 3, 2, 2, 3, 3, 1, 1, 3, 2, 1, 2, 3, 2, 3, 1, 2, 3) ## g.12: Group of patients in stages I and II

PGen12.II(acset.1,Y_T.12,Y_E.12,Y_S.12.obs,status.12,enter.12,d.12,g.12,dose.level)

#########################################################################################################
##### Implement stage III of the PGen I-II design ###############################################
#########################################################################################################

acset.12<-list(c(1,2,3,4), c(2,3,4), c(1,2,3,4))## acset.12: the stage 2 candidate dose sets
Y_T.123 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
             0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
             0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1) ## Y_T.123: Individual toxicity of patients in stages I, II and III
Y_E.123 <- c(1,0,2,2,0,2,2,1,1,1,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,2,2,1,1,2,2,2,2,1,2,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,0,1,2,2,1,2,2,2,1,2,2,2,
             2,1,2,2,2,2,2,1,2,0,2,2,2,1,0,2,2,2,2,2,2,2,2,2,2,2,2,0,2,2,1,2,1,2,1,2,0,1,2,2,1,1,1,2,2,2,2,2,1,2,1,2,2,0,2,2,2,2,2,2,1,2,1,1,2,1,2,2,1,2,2,0,2,
             2,1,1,2) ## Y_E.123: Individual efficacy of patients in stages I, II and III
Y_S.123.obs <- c(1.28,0.00,1.10,5.00,0.00,5.00,5.00,0.37,5.00,0.06,5.00,5.00,0.51,5.00,5.00,5.00,4.96,5.00,5.00,5.00,0.09,0.81,0.18,5.00,5.00,5.00,0.42,5.00,5.00,
                 3.70,3.26,0.17,2.07,1.90,1.94,0.73,5.00,2.45,0.80,5.00,1.21,5.00,1.20,0.74,5.00,5.00,5.00,1.88,3.63,5.00,5.00,5.00,5.00,5.00,5.00,4.61,5.00,4.20,
                 1.32,2.69,5.00,0.00,0.70,5.00,2.36,0.87,3.30,5.00,0.34,5.00,1.58,3.78,0.99,3.06,0.05,5.00,0.43,1.95,5.00,5.00,5.00,5.00,0.00,5.00,0.74,2.95,2.21,
                 0.00,2.05,4.54,5.00,5.00,1.86,3.52,5.00,5.00,4.12,0.77,5.00,5.00,0.00,1.34,5.00,3.08,0.21,3.91,5.00,5.00,5.00,0.00,0.04,0.91,5.00,3.84,0.93,2.32,
                 5.00,5.00,0.53,1.80,3.85,0.64,2.63,4.98,5.00,5.00,0.00,5.00,5.00,5.00,5.00,5.00,5.00,4.04,5.00,0.01,5.00,1.90,0.52,5.00,5.00,5.00,2.71,1.39,0.00,
                 5.00,5.00,5.00,0.37,5.00) ## Y_S.123.obs: Individual observed Y_S^o of patients in stages I, II and III
status.123 <- c(1,1,1,0,1,0,0,1,0,1,0,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,0,1,1,1,1,1,1,1,0,1,1,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,0,1,0,1,1,1,0,1,1,0,1,1,1,0,1,0,1,1,1,
                1,1,0,1,1,0,0,0,0,1,0,1,1,1,1,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,1,1,1,0,0,0,1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,1,0,1,0,1,1,0,0,0,1,1,1,0,
                0,0,1,0) ## status.123: Individual censor indicator of patients in stages I, II and III (1 for no censoring, 0 for censoring)
enter.123 <- c( 0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15,
                16,16,16,17,17,17,18,18,18,19,19,19,20,20,20,21,21,21,22,22,22,23,23,23,24,24,24,25,25,25,26,26,26,27,27,27,28,28,28,29,29,29,30,30,30,31,31,31,
                32,32,32,33,33,33,34,34,34,35,35,35,36,36,36,37,37,37,38,38,38,39,39,39,42,42,42,43,43,43,44,44,44,45,45,45,46,46,46,47,47,47,48,48,48,49,49,49,
                50,50,50,51,51,51) ## enter.123: Individual enter time in stages I, II and III
d.123 <- c(1,1,1,2,2,2,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,2,3,1,4,1,4,3,1,2,2,4,4,4,3,4,2,4,1,2,4,2,3,4,4,2,3,1,3,2,4,3,2,1,2,1,2,1,3,1,4,3,3,4,
           4,1,3,2,3,2,2,3,1,1,2,4,3,2,1,2,3,2,2,1,3,4,1,3,1,4,1,4,3,4,2,1,2,4,3,2,1,3,2,4,3,4,2,3,4,1,4,3,2,4,2,1,4,2,3,4,3,2,3,4,3,4,3,4,2,1,4,1,4,2,2,1,1,
           4,2,2,4) ## d.123: Individual dose allocation of patients in stages I, II and III
g.123 <- c(2,2,2,1,1,1,3,3,3,3,3,3,2,2,2,1,1,1,3,3,3,3,3,3,1,1,1,3,3,3,2,3,3,1,2,1,3,2,2,3,1,2,3,3,2,3,3,2,2,3,3,2,1,2,3,3,3,1,2,3,3,3,1,3,3,1,3,1,2,2,2,2,1,
           3,3,3,2,2,3,3,2,3,1,3,2,2,2,3,1,2,1,1,2,1,1,1,3,3,1,2,2,1,3,3,2,2,3,3,1,1,3,2,1,2,3,2,3,1,2,3,2,2,1,1,1,3,2,1,2,1,3,3,1,3,2,1,2,2,3,3,1,2,1,2,1,3,
           2,2,2,1) ## g.123: Group of patients in stages I, II and III

PGen12.III(acset.12,Y_T.123,Y_E.123,Y_S.123.obs,status.123,enter.123,d.123,g.123,dose.level)
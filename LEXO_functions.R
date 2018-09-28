#### Lagged Exact Bayesian Online Changepoint Detection with Parameter Estimation
#### Code submission

######## Poisson setting ########
# Estimated posterior mean of parameters
# prior = vector (mu0, kappa0, alpha0, beta0)
calculatePosteriorPoisson = function(x, prior){
  alpha = prior[1]
  beta = prior[2]
  n = length(x)
  alpha_n = sum(x) + alpha
  beta_n = beta + n
  return(list(Mean=(n/(n+beta)* mean(x) + (beta)/(n+beta) * alpha/beta), SM=alpha_n*(alpha_n+1)/beta_n^2))
}

calculatePPDPoisson = function(x, alpha, beta){
  #nu and chi is a vector
  mapply(dnbinom,x=x,size=alpha, prob=beta/(beta+1))
}

#############################################
### Beginning algorithm
LEXO_poisson = function(x,alpha0=10^-4,beta0=10^-4,H,ell){
  prior = c(alpha0, beta0)
  totalTime=length(x)
  # Rlist is list of R matrix; each element is R matrix for one lag; R[[1]] = EXO
  # Initialize R matrix and posterior mean matrix
  R = list(ell+1)
  # Theta is list of mean posterior P(theta_t | r_t, x_1:t+\ell), each element corresponds to one lag
  # SM is list of second moment posterior P(theta_t | r_t, x_1:t+\ell), each element corresponds to one lag
  Theta = SM = list(ell+1)
  for (i in 1:(ell+1)){
    Theta[[i]] = SM[[i]] = matrix(0, nrow=totalTime, ncol=totalTime)
    R[[i]] = matrix(0, totalTime, totalTime)
    R[[i]][1,1] = 1
  }
  # Constant hazard rate
  # predprobs: Keeping track of predictive probabilities P(x_t | r_t, x_{1:t-1}) at each time point 
  # P(r_t|x_1:t+1) \propto P(x_t+1 | r_t, x_1:t) P(r_t, x_1:t)
  # reg1 = P(x_t+1 | r_t, x_1:t): regularization term for lexo-1
  predprobs = list()
  # We do not do anything for t=1
  predprobs[[1]] <- NULL
  # Initialize prior hyperparameter
  alphaT = alpha0;
  betaT = beta0;
  
  m=Sys.time()
  for (t in 1:totalTime){
    if (t==1) {
      # Update the sufficient statistics only
      alphaT0 = c(alpha0, alphaT+x[t])
      betaT0 = c(beta0, betaT+1)
      alphaT = alphaT0;
      betaT = betaT0
    }
    else{ #t>=2
      # evaluate predictive probabilities P(x_t | r_t, x_1:t-1)
      predprobs[[t]] = calculatePPDPoisson(x[t], alphaT, betaT)
      # Calculate growth probabilities
      R[[1]][2:t,t]= R[[1]][1:(t-1),t-1]*predprobs[[t]][-1]*(1-H)
      # Calculate changepoint probabilities
      # P(r_t=0, x_1:t)  \sum_{r_t-1} P(x_1:t-1, r_t-1) P(x_t|r_t=0, x_1:t-1) P(r_t|r_t-1)
      # Note that in the above sum, the sum is literally over the first term only
      R[[1]][1,t] = sum(R[[1]][1:(t-1),t-1])*predprobs[[t]][1]*H
      # Normalize
      R[[1]][,t] = R[[1]][,t]/sum(R[[1]][,t])
      # Update sufficient statistics
      alphaT0 = c(alpha0, alphaT+x[t])
      betaT0 = c(beta0, betaT+1)
      alphaT = alphaT0;
      betaT = betaT0
    }
    # Calculating EXO posterior mean parameter estimate
    Q=sapply(0:(t-1), function(j) calculatePosteriorPoisson(x[(t-j):t], prior))
    Theta[[1]][1:t,t]= unlist(Q["Mean",])
    SM[[1]][1:t,t]= unlist(Q["SM",])
    
    # LEXO algorithm
    for (j in 1:ell){
      # Calculate run length
      if (t>=j+1){ # LEXO-1 begins from t=2, LEXO-2 begins from t=3,...
        R[[j+1]][1:(t-j),t-j]=R[[j]][2:(t-j+1),t-j+1]+R[[j]][1,t-j+1]*R[[1]][1:(t-j),t-j]
        # Calculate posterior mean parameter estimate
        #### First calculate alpha 
        # Unnormalize version
        # \alpha = P(r_{t+1}=0 | r_t, x_{1:t\ell})
        alpha1 = 10^200 * R[[1]][1:(t-j),t-j]*R[[j]][1,t-j+1]#avoid underflow
        alpha2 = 10^200 * R[[j]][2:(t-j+1),t-j+1]
        alpha = alpha1/(alpha1+alpha2)
        alpha = replace(alpha, which(alpha1==0),1) # If alpha1 is too small, replace it with 0
        Theta[[j+1]][1:(t-j),t-j]=alpha*Theta[[1]][1:(t-j),t-j]+(1-alpha)*Theta[[j]][2:(t-j+1),t-j+1]
        SM[[j+1]][1:(t-j),t-j]=alpha*SM[[1]][1:(t-j),t-j]+(1-alpha)*SM[[j]][2:(t-j+1),t-j+1]
      }
    }
  }  
  MeanPosterior = mapply(function(x,y) apply(x*y,2,sum), x=R, y=Theta)
  SMPosterior = mapply(function(x,y) apply(x*y,2,sum), x=R, y=SM)
  VarPosterior = SMPosterior-MeanPosterior^2
  m1=Sys.time()-m
  return(list(RunLength=R, PosMean=list(Mean=MeanPosterior, Variance=VarPosterior,runTime=m1)))
}
######### End algorithm ####################################

##### Data Generation for Poisson model in the simulation
data_gen_poisson <- function(n.samp, n.cp,totalTime, amt_change){
  #amt_change = how much the true rate changes from one regime to the next
  the.samples = list()
  the.changeps <- matrix(NA,nrow=n.samp,ncol=n.cp)
  lambda <- matrix(NA,nrow=n.samp,ncol=n.cp+1)
  lambda[,1] = 1;
  for (j in 2:ncol(lambda)){
    lambda[,j] = lambda[,j-1]+amt_change
  }
  the.cp <- seq(floor(totalTime/(n.cp)),totalTime,by=floor(totalTime/(n.cp+1)))[-(n.cp+1)]
  true_seg=seg(the.cp, totalTime)
  start_point = c(1,the.cp)
  end_point = c(the.cp,totalTime+1)
  true_theta =rep(lambda[1,],end_point-start_point)
  # one sample
  for (i in 1:n.samp){
    #the.cp <- sort(sample(10:990,n.cp)) ; 
    the.cp <- seq(floor(totalTime/(n.cp)),totalTime,by=floor(totalTime/(n.cp+1)))[-(n.cp+1)]
    the.changeps[i,] <- the.cp # Choose where cp happens
    the.samples[[i]] = numeric(totalTime)
    the.samples[[i]][1:the.cp[1]] = rpois(the.cp[1],1) 
    for (k in 1:n.cp){
      if (k==n.cp) the.samples[[i]][(the.cp[k]+1):length(the.samples[[i]])] = rpois(length(the.samples[[i]])-the.cp[k], lambda[i,k+1])
      else the.samples[[i]][(the.cp[k]+1):the.cp[k+1]] = rpois(the.cp[k+1]-the.cp[k], lambda[i,k+1])
    }
  } 
  return(list(Samples=the.samples,CP=the.changeps, the.cp=the.cp, true_theta=true_theta, theta=lambda, true_seg=true_seg))  
}
#####################################################
########### Normal setting ##########################
# prior = (mu0, kappa0, alpha0, beta0) 
# For precision shift model when mean is assumed to be known, put kappa0 to be a very large value
calculatePosteriorNormal = function(x, prior){
  mu0 = prior[1]
  kappa0 = prior[2]
  alpha0 = prior[3]
  beta0 = prior[4]
  n = length(x)
  # Update posterior
  mu_n = mu0 + (kappa0*mu0+sum(x))/(kappa0+n)
  kappa_n = kappa0 + n
  alpha_n = alpha0 +n/2
  beta_n = beta0 + 1/2*sum((x-mean(x))^2)+(kappa0*n)*(mean(x)-mu0)^2/(2*(kappa0+n))
  # Calculate
  MeanMean = mu_n
  MeanPrecision = alpha_n/beta_n
  VarianceMean = 2*alpha_n/(2*alpha_n-2) * beta_n/(alpha_n*kappa_n)
  VariancePrecision=alpha_n/beta_n^2
  SMMean=MeanMean^2+VarianceMean
  SMPrecision=MeanPrecision^2 + VariancePrecision
  return(list(MeanMean=MeanMean, MeanPrecision=MeanPrecision, SMMean=SMMean, SMPrecision=SMPrecision))
}

# Calculate PPD
#### Functions that are used to calculate PPD
studentpdf <- function(x, df, mu0, var){
  c = exp(lgamma(df/2 + 0.5) - lgamma(df/2)) * (df*pi*var)^(-0.5);
  p = c * (1 + 1/(df*var)*(x-mu0)^2)^(-(df+1)/2);
  if (is.nan(p)) print("NaN")
  return(p)
}  
calculatePPDNormal = function(x, mu, kappa, alpha, beta){
  #nu and chi is a vector
  mapply(studentpdf,x=x,df=2*alpha, mu0=mu, var=beta*(kappa+1)/(alpha*kappa))
}

#############################################
### Beginning algorithm
# prior = (mu0, kappa0, alpha0, beta0)
LEXO_normal = function(x, mu0, kappa0,alpha0,beta0, H, ell){
  prior =c(mu0, kappa0,alpha0,beta0)
  totalTime=length(x)
  # Rlist is list of R matrix; each element is R matrix for one lag; R[[1]] = EXO
  # Initialize R matrix and posterior mean matrix
  R = list(ell+1)
  # Theta is list of mean posterior P(theta_t | r_t, x_1:t+\ell), each element corresponds to one lag
  # SM is list of second moment posterior P(theta_t | r_t, x_1:t+\ell), each element corresponds to one lag
  MeanMean = SMMean = MeanPrecision=SMPrecision=list(ell+1)
  
  for (i in 1:(ell+1)){
    MeanMean[[i]] = SMMean[[i]]=MeanPrecision[[i]] = SMPrecision[[i]] =
      matrix(0, nrow=totalTime, ncol=totalTime)
    R[[i]] = matrix(0, totalTime, totalTime)
    R[[i]][1,1] = 1
  }
  # Constant hazard rate
  # predprobs: Keeping track of predictive probabilities P(x_t | r_t, x_{1:t-1}) at each time point 
  # P(r_t|x_1:t+1) \propto P(x_t+1 | r_t, x_1:t) P(r_t, x_1:t)
  # reg1 = P(x_t+1 | r_t, x_1:t): regularization term for lexo-1
  predprobs = list()
  # We do not do anything for t=1
  predprobs[[1]] <- NULL
  # Initialize prior hyperparameter
  muT = mu0; kappaT = kappa0
  alphaT = alpha0; betaT=beta0;
  
  m=Sys.time()
  for (t in 1:totalTime){
    if (t==1) {
      # Update the sufficient statistics only
      muT0 = c(mu0, (kappaT*muT+x[t])/(kappaT+1))
      kappaT0 = c(kappa0, kappaT+1)
      alphaT0 = c(alpha0, alphaT+1/2)
      betaT0 = c(beta0, betaT + kappaT*(x[t]-muT)^2/(2*(kappaT+1)))
      muT = muT0;
      kappaT = kappaT0
      alphaT = alphaT0
      betaT = betaT0
    }
    else{ #t>=2
      # evaluate predictive probabilities P(x_t | r_t, x_1:t-1)
      predprobs[[t]] = calculatePPDNormal(x[t], muT, kappaT, alphaT, betaT)
      # Calculate growth probabilities
      R[[1]][2:t,t]= R[[1]][1:(t-1),t-1]*predprobs[[t]][-1]*(1-H)
      # Calculate changepoint probabilities
      # P(r_t=0, x_1:t)  \sum_{r_t-1} P(x_1:t-1, r_t-1) P(x_t|r_t=0, x_1:t-1) P(r_t|r_t-1)
      # Note that in the above sum, the sum is literally over the first term only
      R[[1]][1,t] = sum(R[[1]][1:(t-1),t-1])*predprobs[[t]][1]*H
      # Normalize
      R[[1]][,t] = R[[1]][,t]/sum(R[[1]][,t])
      # Update sufficient statistics
      muT0 = c(mu0, (kappaT*muT+x[t])/(kappaT+1))
      kappaT0 = c(kappa0, kappaT+1)
      alphaT0 = c(alpha0, alphaT+1/2)
      betaT0 = c(beta0, betaT + kappaT*(x[t]-muT)^2/(2*(kappaT+1)))
      muT = muT0;
      kappaT = kappaT0
      alphaT = alphaT0
      betaT = betaT0
    }
    # Calculating EXO posterior mean parameter estimate
    Q=sapply(0:(t-1), function(j) calculatePosteriorNormal(x[(t-j):t], prior))
    MeanMean[[1]][1:t,t]= unlist(Q["MeanMean",])
    SMMean[[1]][1:t,t]= unlist(Q["SMMean",])
    MeanPrecision[[1]][1:t,t]= unlist(Q["MeanPrecision",])
    SMPrecision[[1]][1:t,t]= unlist(Q["SMPrecision",])
    #SM[[1]][1:t,t]= sapply(0:(t-1), function(j) calculatePosteriorSecondMoment(x[(t-j):t], prior=c(10^-4, 10^-4)))
    
    # LEXO algorithm
    for (j in 1:ell){
      # Calculate run length
      if (t>=j+1){ # LEXO-1 begins from t=2, LEXO-2 begins from t=3,...
        R[[j+1]][1:(t-j),t-j]=R[[j]][2:(t-j+1),t-j+1]+R[[j]][1,t-j+1]*R[[1]][1:(t-j),t-j]
        # Calculate posterior mean parameter estimate
        #### First calculate alpha 
        # Unnormalize version
        #\alpha = P(r_{t+1}=0 | r_t, x_{1:t\ell})
        alpha1 = 10^200 * R[[1]][1:(t-j),t-j]*R[[j]][1,t-j+1]#avoid underflow
        alpha2 = 10^200 * R[[j]][2:(t-j+1),t-j+1]
        alpha = alpha1/(alpha1+alpha2)
        alpha = replace(alpha, which(alpha1==0),0)
        #alpha = replace(alpha,which(alpha<1e-300),0)
        MeanMean[[j+1]][1:(t-j),t-j]=alpha*MeanMean[[1]][1:(t-j),t-j]+(1-alpha)*MeanMean[[j]][2:(t-j+1),t-j+1]
        MeanPrecision[[j+1]][1:(t-j),t-j]=alpha*MeanPrecision[[1]][1:(t-j),t-j]+(1-alpha)*MeanPrecision[[j]][2:(t-j+1),t-j+1]
        SMMean[[j+1]][1:(t-j),t-j]=alpha*SMMean[[1]][1:(t-j),t-j]+(1-alpha)*SMMean[[j]][2:(t-j+1),t-j+1]
        SMPrecision[[j+1]][1:(t-j),t-j]=alpha*SMPrecision[[1]][1:(t-j),t-j]+(1-alpha)*SMPrecision[[j]][2:(t-j+1),t-j+1]
      }
    }
  }  
  MeanMeanPosterior = mapply(function(x,y) apply(x*y,2,sum), x=R, y=MeanMean)
  SMMeanPosterior = mapply(function(x,y) apply(x*y,2,sum), x=R, y=SMMean)
  
  MeanPrecisionPosterior = mapply(function(x,y) apply(x*y,2,sum), x=R, y=MeanPrecision)
  SMPrecisionPosterior = mapply(function(x,y) apply(x*y,2,sum), x=R, y=SMPrecision)
  
  VarMeanPosterior = SMMeanPosterior-MeanMeanPosterior^2
  VarPrecisionPosterior=SMPrecisionPosterior-MeanPrecisionPosterior^2
  m1=Sys.time()-m
  return(list(RunLength=R, PosPrecision=list(Mean=MeanPrecisionPosterior, Variance=VarPrecisionPosterior),
              PosMean=list(Mean=MeanMeanPosterior,Variance=VarMeanPosterior), RunTime=m1))
}

data_gen_normal <- function(n.samp, n.cp, totalTime, amt){
  # For mean shift model
  # amt_change = change in the mean from one to the other
  the.samples = list()
  the.changeps <- matrix(NA,nrow=n.samp,ncol=n.cp)
  mu <- matrix(NA,nrow=n.samp,ncol=n.cp+1)
  mu[,1] = 0;
  for (j in 2:ncol(mu)){
    mu[,j] = mu[,j-1]+amt
  }
  the.cp <- seq(floor(totalTime/(n.cp)),totalTime,by=floor(totalTime/(n.cp+1)))[-(n.cp+1)]
  true_seg=seg(the.cp, totalTime)
  start_point = c(1,the.cp)
  end_point = c(the.cp,totalTime+1)
  true_theta =rep(mu[1,],end_point-start_point)
  sigma <- matrix(NA,nrow=n.samp,ncol=n.cp+1)
  sigma[,1] = 1;
  for (j in 2:ncol(sigma)){
    sigma[,j]=rep(1,n.samp)
    #p = runif(n.samp,0,1)
    #sigma[,j] = (p<0.5)*sigma[,j-1]+ (p>0.5)*runif(n.samp,sigma[,j-1]*1.5, sigma[,j-1]*2)
  }
  # one sample
  for (i in 1:n.samp){
    the.changeps[i,] <- the.cp # Choose where cp happens
    the.samples[[i]] = numeric(1000)
    the.samples[[i]][1:(the.cp[1]-1)] = rnorm((the.cp[1]-1),mu[i,1],sigma[i,1]) 
    for (k in 1:n.cp){
      if (k==n.cp) the.samples[[i]][(the.cp[k]):length(the.samples[[i]])] = rnorm(length(the.samples[[i]])-the.cp[k]+1, mu[i,k+1],sqrt(sigma[i,k+1]))
      else the.samples[[i]][(the.cp[k]):(the.cp[k+1]-1)] = rnorm(the.cp[k+1]-the.cp[k], mu[i,k+1], sqrt(sigma[i,k+1]))
    }
  } 
  return(list(Samples=the.samples,CP=the.changeps, the.cp=the.cp, true_theta=true_theta, theta=mu, true_seg=true_seg))  
}

data_gen_normal_variance <- function(n.samp, n.cp, totalTime, amt){
  # amt_change = change_in_precision from one regime to the next
  the.samples = list()
  the.changeps <- matrix(NA,nrow=n.samp,ncol=n.cp)
  mu <- matrix(NA,nrow=n.samp,ncol=n.cp+1)
  mu[,1] = 0;
  for (j in 2:ncol(mu)){
    mu[,j] = mu[,j-1]
  }
  the.cp <- seq(floor(totalTime/(n.cp)),totalTime,by=floor(totalTime/(n.cp+1)))[-(n.cp+1)]
  true_seg=seg(the.cp, totalTime)
  start_point = c(1,the.cp)
  end_point = c(the.cp,totalTime+1)
  sigma <- matrix(NA,nrow=n.samp,ncol=n.cp+1)
  sigma[,1] = 1/16;
  for (j in 2:ncol(sigma)){
    sigma[,j]=sigma[,j-1]*amt # sigma: variance
  }
  true_theta =rep(1/sigma[1,],end_point-start_point)
  # one sample
  for (i in 1:n.samp){
    the.changeps[i,] <- the.cp # Choose where cp happens
    the.samples[[i]] = numeric(1000)
    the.samples[[i]][1:(the.cp[1]-1)] = rnorm((the.cp[1]-1),mu[i,1],sqrt(sigma[i,1])) 
    for (k in 1:n.cp){
      if (k==n.cp) the.samples[[i]][(the.cp[k]):length(the.samples[[i]])] = rnorm(length(the.samples[[i]])-the.cp[k]+1, mu[i,k+1],sqrt(sigma[i,k+1]))
      else the.samples[[i]][(the.cp[k]):(the.cp[k+1]-1)] = rnorm(the.cp[k+1]-the.cp[k], mu[i,k+1], sqrt(sigma[i,k+1]))
    }
  } 
  return(list(Samples=the.samples,CP=the.changeps, the.cp=the.cp, true_theta=true_theta, theta=mu, true_seg=true_seg))  
}
######### General Function
seg<- function(cp,bigT){
  if (length(cp)==0) seg=rep(1,bigT) else{
    cp=append(cp, bigT)
    seg = rep(1, cp[1]-1)
    for (i in 2:length(cp)){
      seg=c(seg,rep(i,cp[i]-cp[i-1]))
    }
    seg = c(seg, seg[bigT-1])
  }
  return(seg)  
}

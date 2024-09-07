# Y_1 = \beta_1 * Z + \omega_1 * A * Z_pred + \gamma_1 * A  
# Y_2 = \beta_2 * Z + \omega_2 * A * Z_pred + \gamma_2 * A
# Z is a vector of covariates
# \beta_1 and \beta_2: vector of coefficients for Z, length=10
# \omega_1 and \omega_2: vector of coefficients for Z*A, length=10 (non-zero \omega also indicates the predictive covariates)
# \gamma_1, \gamma_2 
# nobs: number of observation in the dataset
# expis the rate parameter for exponential distribution of censoring times
# p is the parameter for indicator of causes, p=Pr(T<=\infty,\theta=1|Y(A,Z)=0) 

FG.sim.10<- function(nobs, beta1, beta2, omega1, omega2, gamma, scaleC, exp, p,lambda0,censor_beta) {
  Z<-matrix(nrow=nobs, ncol=10)
  Z[,1:5]<-rbinom(nobs*5,1,0.5) # baseline binary covariates
  Z[,6:10]<-rnorm(nobs*5,0,1)
  
  A<-rbinom(nobs, 1, 1 / (1 + exp(-Z%*%gamma))) # treatments 1=treatment, -1=control logit(p(A=1))=Z
  Inter<-Z*A
  
  Y1.trt<-Z%*%beta1+0.5*(1-0.5)
  Y1.cnt<-Z%*%omega1+0.5*(0-0.5)
  Y2.trt<-Z%*%beta2+0.25*(1-0.5)
  Y2.cnt<-Z%*%omega2+0.25*(0-0.5)
  Y1<-ifelse(A==1, Y1.trt, Y1.cnt)
  Y2<-ifelse(A==1, Y2.trt, Y2.cnt)
  
  #- Generate indicator for cause
  theta <- 1 + rbinom(nobs, 1, prob = (1 - p)^exp(Y1))
  U <-runif(Y1[theta==1],0,1)
  T1.cnt<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.trt<- (-log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)) # event type 1
  summary(T1.cnt)
  summary(T1.trt)
  T1<-ifelse(A[theta==1]==1,T1.trt,T1.cnt)
  T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] )) # event type 1
  C <- (-log(runif(nobs,0,1))/(lambda0*exp(Z%*%matrix(censor_beta)))) #simulate censoring times

  event.time<-numeric(nobs)
  event.time[theta==1]<-T1
  event.time[theta==2]<-T2
  event.time<-pmin(event.time,C)
  event.type<-ifelse(event.time==C,0,1)
  event.type<-event.type*theta
  return(data.frame(A=A,event.time=event.time,event.type=event.type,Z=Z,Inter=Inter))
}
FG.sim.30<- function(nobs, beta1, beta2, omega1, omega2, gamma, scale1, shape,scaleC, exp, p,lambda0,censor_beta) {
  Z<-matrix(nrow=nobs, ncol=30)
  Z[,1:5]<-rbinom(nobs*5,1,0.5) # baseline binary covariates
  Z[,6:30]<-rnorm(nobs*25,0,1)
  
  A<-rbinom(nobs, 1, 1 / (1 + exp(-Z%*%gamma))) # treatments 1=treatment, -1=control logit(p(A=1))=Z
  Inter<-Z*A
  
  Y1.trt<-Z%*%beta1+0.5*(1-0.5)
  Y1.cnt<-Z%*%omega1+0.5*(0-0.5)
  Y2.trt<-Z%*%beta2+0.25*(1-0.5)
  Y2.cnt<-Z%*%omega2+0.25*(0-0.5)
  Y1<-ifelse(A==1, Y1.trt, Y1.cnt)
  Y2<-ifelse(A==1, Y2.trt, Y2.cnt)
  
  #- Generate indicator for cause
  theta <- 1 + rbinom(nobs, 1, prob = (1 - p)^exp(Y1))
  #Conditional on cause indicators, we simulate the model.
  U <-runif(Y1[theta==1],0,1)
  T1.cnt<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.trt<- (-(1/scale1)*log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p))^(1/shape) # event type 1
  summary(T1.cnt)
  summary(T1.trt)
  T1<-ifelse(A[theta==1]==1,T1.trt,T1.cnt)
  T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] )) # event type 1
  C <- (-log(runif(nobs,0,1))/(lambda0*exp(Z%*%matrix(censor_beta)))) #simulate censoring times

  event.time<-numeric(nobs)
  event.time[theta==1]<-T1
  event.time[theta==2]<-T2
  event.time<-pmin(event.time,C)
  event.type<-ifelse(event.time==C,0,1)
  event.type<-event.type*theta
  return(data.frame(A=A,event.time=event.time,event.type=event.type,Z=Z,Inter=Inter))
}

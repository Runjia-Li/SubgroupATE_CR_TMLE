# Y_1 = \beta_1 * Z + \omega_1 * A * Z_pred + \gamma_1 * A  
# Y_2 = \beta_2 * Z + \omega_2 * A * Z_pred + \gamma_2 * A
# Z is a d-dimention vector of covariates
# \beta_1 and \beta_2: vector of coefficients for Z, length=10
# \omega_1 and \omega_2: vector of coeffcients for Z*A, length=10 (non-zero \omega also inidcates the preditive covariates)
# \gamma_1, \gamma_2 

# nobs: number of observation in the dataset
# scaleC is the scale paramenter for exponential distribution of censoring times
# expis the rate paramenter for exponential distribution of censoring times
# p is the parameter for indicator of causes, p=Pr(T<=\infty,\theta=1|Y(A,Z)=0) 

# ----------------d: number of covariates
# ----------------binary: index of binary covaraites
FG.sim.6<- function(nobs, beta1, beta2, omega1, omega2, gamma, scale1, shape,scaleC, exp, p,lambda0,censor_beta) {
  Z<-matrix(nrow=nobs, ncol=6)
  Z[,1:3]<-rbinom(nobs*3,1,0.5) # baseline binary covariates
  # Z[,11:20]<-rnorm(nobs*10,0,0.5) # baseline continuous covariates
  Z[,4:6]<-rnorm(nobs*3,0,1)
  
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
  # table(theta)
  #Conditional on cause indicators, we simulate the model.
  U <-runif(Y1[theta==1],0,1)
  # T1<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.cnt<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.trt<- (-(1/scale1)*log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p))^(1/shape) # event type 1
  summary(T1.cnt)
  summary(T1.trt)
  T1<-ifelse(A[theta==1]==1,T1.trt,T1.cnt)
  T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] )) # event type 1
  C <- (-log(runif(nobs,0,1))/(lambda0*exp(Z%*%matrix(censor_beta)))) #simulate censoring times
  # C<-runif(nobs, min = min(T1,T2), max = quantile(c(T1,T2),0.95))  
  
  event.time<-numeric(nobs)
  event.time[theta==1]<-T1
  event.time[theta==2]<-T2
  event.time<-pmin(event.time,C)
  event.type<-ifelse(event.time==C,0,1)
  event.type<-event.type*theta
  return(data.frame(A=A,event.time=event.time,event.type=event.type,Z=Z,Inter=Inter))
}

FG.sim.8<- function(nobs, beta1, beta2, omega1, omega2, gamma, scale1, shape,scaleC, exp, p,lambda0,censor_beta) {
  Z<-matrix(nrow=nobs, ncol=8)
  Z[,1:4]<-rbinom(nobs*4,1,0.5) # baseline binary covariates
  # Z[,11:20]<-rnorm(nobs*10,0,0.5) # baseline continuous covariates
  Z[,5:8]<-rnorm(nobs*4,0,1)
  
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
  # table(theta)
  #Conditional on cause indicators, we simulate the model.
  U <-runif(Y1[theta==1],0,1)
  # T1<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.cnt<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.trt<- (-(1/scale1)*log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p))^(1/shape) # event type 1
  summary(T1.cnt)
  summary(T1.trt)
  T1<-ifelse(A[theta==1]==1,T1.trt,T1.cnt)
  T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] )) # event type 1
  C <- (-log(runif(nobs,0,1))/(lambda0*exp(Z%*%matrix(censor_beta)))) #simulate censoring times
  # C<-runif(nobs, min = min(T1,T2), max = quantile(c(T1,T2),0.95))  
  
  event.time<-numeric(nobs)
  event.time[theta==1]<-T1
  event.time[theta==2]<-T2
  event.time<-pmin(event.time,C)
  event.type<-ifelse(event.time==C,0,1)
  event.type<-event.type*theta
  return(data.frame(A=A,event.time=event.time,event.type=event.type,Z=Z,Inter=Inter))
}


FG.sim.10<- function(nobs, beta1, beta2, omega1, omega2, gamma, scaleC, exp, p,lambda0,censor_beta) {
  Z<-matrix(nrow=nobs, ncol=10)
  Z[,1:5]<-rbinom(nobs*5,1,0.5) # baseline binary covariates
  # Z[,11:20]<-rnorm(nobs*10,0,0.5) # baseline continuous covariates
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
  # table(theta)
  #Conditional on cause indicators, we simulate the model.
  U <-runif(Y1[theta==1],0,1)
  # T1<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.cnt<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.trt<- (-log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)) # event type 1
  summary(T1.cnt)
  summary(T1.trt)
  T1<-ifelse(A[theta==1]==1,T1.trt,T1.cnt)
  T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] )) # event type 1
  C <- (-log(runif(nobs,0,1))/(lambda0*exp(Z%*%matrix(censor_beta)))) #simulate censoring times
  # C<-runif(nobs, min = min(T1,T2), max = quantile(c(T1,T2),0.95))  
  
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
  # Z[,11:20]<-rnorm(nobs*10,0,0.5) # baseline continuous covariates
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
  # table(theta)
  #Conditional on cause indicators, we simulate the model.
  U <-runif(Y1[theta==1],0,1)
  # T1<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.cnt<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.trt<- (-(1/scale1)*log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p))^(1/shape) # event type 1
  summary(T1.cnt)
  summary(T1.trt)
  T1<-ifelse(A[theta==1]==1,T1.trt,T1.cnt)
  T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] )) # event type 1
  C <- (-log(runif(nobs,0,1))/(lambda0*exp(Z%*%matrix(censor_beta)))) #simulate censoring times
  # C<-runif(nobs, min = min(T1,T2), max = quantile(c(T1,T2),0.95))  
  
  event.time<-numeric(nobs)
  event.time[theta==1]<-T1
  event.time[theta==2]<-T2
  event.time<-pmin(event.time,C)
  event.type<-ifelse(event.time==C,0,1)
  event.type<-event.type*theta
  return(data.frame(A=A,event.time=event.time,event.type=event.type,Z=Z,Inter=Inter))
}

FG.sim.40<- function(nobs, beta1, beta2, omega1, omega2, gamma, scale1, shape,scaleC, exp, p,lambda0,censor_beta) {
  Z<-matrix(nrow=nobs, ncol=40)
  Z[,1:5]<-rbinom(nobs*5,1,0.5) # baseline binary covariates
  # Z[,11:20]<-rnorm(nobs*10,0,0.5) # baseline continuous covariates
  Z[,6:40]<-rnorm(nobs*35,0,1)
  
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
  # table(theta)
  #Conditional on cause indicators, we simulate the model.
  U <-runif(Y1[theta==1],0,1)
  # T1<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.cnt<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.trt<- (-(1/scale1)*log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p))^(1/shape) # event type 1
  summary(T1.cnt)
  summary(T1.trt)
  T1<-ifelse(A[theta==1]==1,T1.trt,T1.cnt)
  T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] )) # event type 1
  C <- (-log(runif(nobs,0,1))/(lambda0*exp(Z%*%matrix(censor_beta)))) #simulate censoring times
  # C<-runif(nobs, min = min(T1,T2), max = quantile(c(T1,T2),0.95))  
  
  event.time<-numeric(nobs)
  event.time[theta==1]<-T1
  event.time[theta==2]<-T2
  event.time<-pmin(event.time,C)
  event.type<-ifelse(event.time==C,0,1)
  event.type<-event.type*theta
  return(data.frame(A=A,event.time=event.time,event.type=event.type,Z=Z,Inter=Inter))
}
FG.sim.20<- function(nobs, beta1, beta2, omega1, omega2, gamma, scale1, shape,scaleC, exp, p,lambda0,censor_beta) {
  Z<-matrix(nrow=nobs, ncol=20)
  Z[,1:10]<-rbinom(nobs*10,1,0.5) # baseline binary covariates
  # Z[,11:20]<-rnorm(nobs*10,0,0.5) # baseline continuous covariates
  Z[,11:20]<-rnorm(nobs*10,0,1)
  
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
  # table(theta)
  #Conditional on cause indicators, we simulate the model.
  U <-runif(Y1[theta==1],0,1)
  # T1<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.cnt<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.trt<- (-(1/scale1)*log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p))^(1/shape) # event type 1
  summary(T1.cnt)
  summary(T1.trt)
  T1<-ifelse(A[theta==1]==1,T1.trt,T1.cnt)
  T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] )) # event type 1
  C <- (-log(runif(nobs,0,1))/(lambda0*exp(Z%*%matrix(censor_beta)))) #simulate censoring times
  # C<-runif(nobs, min = min(T1,T2), max = quantile(c(T1,T2),0.95))  
  
  event.time<-numeric(nobs)
  event.time[theta==1]<-T1
  event.time[theta==2]<-T2
  event.time<-pmin(event.time,C)
  event.type<-ifelse(event.time==C,0,1)
  event.type<-event.type*theta
  return(data.frame(A=A,event.time=event.time,event.type=event.type,Z=Z,Inter=Inter))
}
FG.sim.40<- function(nobs, beta1, beta2, omega1, omega2, gamma, scale1, shape,scaleC, exp, p,lambda0,censor_beta) {
  Z<-matrix(nrow=nobs, ncol=40)
  Z[,1:10]<-rbinom(nobs*10,1,0.5) # baseline binary covariates
  # Z[,11:20]<-rnorm(nobs*10,0,0.5) # baseline continuous covariates
  Z[,11:40]<-rnorm(nobs*30,0,1)
  
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
  # table(theta)
  #Conditional on cause indicators, we simulate the model.
  U <-runif(Y1[theta==1],0,1)
  # T1<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.cnt<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
  T1.trt<- (-(1/scale1)*log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p))^(1/shape) # event type 1
  summary(T1.cnt)
  summary(T1.trt)
  T1<-ifelse(A[theta==1]==1,T1.trt,T1.cnt)
  T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] )) # event type 1
  C <- (-log(runif(nobs,0,1))/(lambda0*exp(Z%*%matrix(censor_beta)))) #simulate censoring times
  # C<-runif(nobs, min = min(T1,T2), max = quantile(c(T1,T2),0.95))  
  
  event.time<-numeric(nobs)
  event.time[theta==1]<-T1
  event.time[theta==2]<-T2
  event.time<-pmin(event.time,C)
  event.type<-ifelse(event.time==C,0,1)
  event.type<-event.type*theta
  return(data.frame(A=A,event.time=event.time,event.type=event.type,Z=Z,Inter=Inter))
}
# ggplot(data.frame(T1,A[theta==1]),aes(T1, fill =as.factor( `A.theta....1.`))) + geom_density(alpha = 0.2)
# summary(C)
# FG.sim.20<- function(nobs, beta1, beta2, scaleC, exp, p) {
#   Z<-matrix(nrow=nobs, ncol=10)
#   Z[,1:5]<-rbinom(nobs*5,1,0.5) # baseline binary covariates
#   Z[,6:10]<-rnorm(nobs*5,0,1) # baseline continuous covariates
#   A<-rbinom(nobs, 1, 1 / (1 + exp(-0.1*Z[,1]-0.1*Z[,6]))) # treatments 1=treatment, -1=control logit(p(A=1))=Z
#   Inter<-Z*A
#   X.matrix<-cbind(A,Z,Inter)
#   colnames(X.matrix)<-c("A",paste0(rep("Z",10),c(1:10)),paste0(rep("Inter",10),c(1:10)))
#   
#   Y1<-X.matrix%*%beta1
#   Y2<-X.matrix%*%beta2
#   
#   #- Generate indicator for cause
#   theta <- 1 + rbinom(nobs, 1, prob = (1 - p)^exp(Y1))
#   # table(theta)
#   #Conditional on cause indicators, we simulate the model.
#   U <-runif(Y1[theta==1],0,1)
#   T1<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
#   T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] ))
#   C <- (-log(runif(nobs,0,1))/(scaleC*exp(exp))) #simulate censoring times
#   
#   event.time<-numeric(nobs)
#   event.time[theta==1]<-T1
#   event.time[theta==2]<-T2
#   event.time<-pmin(event.time,C)
#   event.type<-ifelse(event.time==C,0,1)
#   event.type<-event.type*theta
#   return(data.frame(A=A,event.time=event.time,event.type=event.type,Z=Z,Inter=Inter))
# }
# 
# ########################
# nobs=1000
# beta1=c(-0.1, 0.2, 0, rep(0,8), -0.1,0,rep(0,8)); beta2=beta1/2;
# # beta1=c(-0.01, 0.02,0.01,rep(0,8), -0.01,0,rep(0,8)); beta2=beta1/2;
# scaleC=1/10; exp=0; p=0.5
# Z<-matrix(nrow=nobs, ncol=10)
# Z[,1:5]<-rbinom(nobs*5,1,0.5) # baseline binary covariates
# Z[,6:10]<-rnorm(nobs*5,0,1) # baseline continuous covariates
# A<-rbinom(nobs, 1, 1 / (1 + exp(-0.1*Z[,1]-0.1*Z[,6]))) # treatments 1=treatment, -1=control logit(p(A=1))=Z
# Inter<-Z*A
# X.matrix<-cbind(A,Z,Inter)
# dat1 <- simulateTwoCauseFineGrayModel(nobs, beta1, beta2, X.matrix, u.min = 1, u.max = 10, p = 0.5)
# dat.Sim<-data.frame(A=A,event.time=dat1$ftime,event.type=dat1$fstatus,Z=Z,Inter=Inter)
# table1(~factor(event.type)+event.time+factor(A),data=dat.Sim)
# 
# data.train=dat.Sim
# lam.path <- 10^seq(log10(0.01), log10(0.000001), length = 100)
# (crr.Sim<-crr(data.train$event.time, data.train$event.type, cbind(A,Z,Inter))$coef)
# # crr.Sim2<-crr(data.train$event.time, data.train$event.type, cbind(A,Z[,c(1,2)],Inter[,c(1)]))$coef
# (crr.Sim2<-crr(data.train$event.time, data.train$event.type, cbind(A,Z[,c(1)],Inter[,c(1)]))$coef)
# 
# weight <- 1/abs(crr.Sim)
# names(weight)[1]<-"A"
# model<-as.formula(paste("Crisk(data.train$event.time, data.train$event.type, cencode = 0, failcode = 1)~A+", paste(paste0(rep("Z.",10),c(1:10)), collapse= "+"),"+",paste(paste0(rep("Inter.",10),c(1:10)), collapse= "+")))
# fit.al.Sim <- fastCrrp(model,data=data.train, penalty="LASSO",lambda = lam.path,penalty.factor=weight)
# coef.Sim<-fit.al.Sim$coef
# View(coef.Sim)                   
# 
# ########################
# n.train=1000;
# beta1=c(0.2,0,rep(0,8)); beta2=beta1/2;
# omega1=c(-0.1,0,rep(0,8));omega2=omega1/2;gamma1=-0.1; gamma2=gamma1/2; 
# scaleC=1/10; exp=0; p=0.5
# data.train<-FG.sim.10(n.train, beta1, beta2, omega1, omega2, gamma1, gamma2, scaleC, exp, p) # shape=1 is exponential
# table1(~factor(event.type)+event.time+factor(A),data=data.train)
# 
# lam.path <- 10^seq(log10(0.01), log10(0.000001), length = 100)
# (crr.man<-crr(data.train$event.time, data.train$event.type, cbind(data.train$A,data.train[,c(4:13)],data.train[,c(14:23)]))$coef)
# (crr.man2<-crr(data.train$event.time, data.train$event.type, cbind(data.train$A,data.train$Z.1,data.train$Inter.1))$coef)
# 
# weight <- 1/abs(crr.man)
# names(weight)[1]<-"A"
# colnames(X.matrix)<-c("A",paste0(rep("Z.",10),c(1:10)),paste0(rep("Inter.",10),c(1:10)))
# fit.al.man<- fastCrrp(model,data=data.train, penalty="LASSO",lambda = lam.path,penalty.factor=weight) 
# coef.man<-fit.al.man$coef
# View(coef.man)
# 
# Y1_1<-as.matrix(data.train[,c(4:13)])%*%beta1+as.matrix(data.train[,c(14:23)])%*%omega1+data.train$A*1
# Y1_0<-as.matrix(data.train[,c(4:13)])%*%beta1
# trueCIF_1<-sapply(tp,function(x) (1-(1-p*(1-exp(-x)))^exp(Y1_1)))
# trueCIF_0<-sapply(tp,function(x) (1-(1-p*(1-exp(-x)))^exp(Y1_0)))
# 
# predict(crr.man)



# nobs=1000; 
# beta1=c(0.02,0.01,rep(0,8)); beta2=c(0.01,0.005,rep(0,8));
# omega1=c(-0.01,0,rep(0,8));omega2=c(-0.005,0,rep(0,8));gamma1=-0.01; gamma2=-0.005; 
# scaleC=1/10; exp=0; p=0.5
# train<-FG.sim.10(nobs, beta1, beta2, omega1, omega2, gamma1, gamma2, scaleC, exp, p)
# table1(~factor(event.type)+event.time+factor(A),data=train)

# ######### check propensity score
# g_true <- 1 / (1 + exp(-Z[,1]))
# sl_libs <- c( 'SL.glmnet',  'SL.glm') #'SL.earth','SL.ranger',
# g <- SuperLearner(Y = A, # outcome is the A (treatment) vector
#                   X = as.data.frame(Z), # W is a matrix of predictors
#                   family=binomial(), # treatment is a binomial outcome
#                   SL.library=sl_libs) # using same candidate learners; could use different learners
# g_w <- as.vector(predict(g)$pred) # Pr(A=1|W)
# 
# prop.sc<-as.data.frame(A=A,g_w)
# ggplot(prop.sc, aes(x=g_w)) +                       # Draw overlaying histogram
#   geom_histogram(data=subset(prop.sc,A ==0),fill = "red", alpha = 0.2,bins = 50) +
#   geom_histogram(data=subset(prop.sc,A ==1), fill = "blue", alpha = 0.2,bins = 50)
# ggsave("hist_scenario2.png", width=4, height=6, dpi=600)#units="cm",

# FG.sim.30.10<- function(nobs, beta1, beta2, omega1, omega2, gamma, scale1, shape,scaleC, exp, p,lambda0,censor_beta) {
#   Z<-matrix(nrow=nobs, ncol=30)
#   Z[,1:4]<-rbinom(nobs*4,1,0.5) # baseline binary covariates
#   # Z[,11:20]<-rnorm(nobs*10,0,0.5) # baseline continuous covariates
#   Z[,5:30]<-rnorm(nobs*26,0,1)
#   
#   A<-rbinom(nobs, 1, 1 / (1 + exp(-Z%*%gamma))) # treatments 1=treatment, -1=control logit(p(A=1))=Z
#   Inter<-Z*A
#   
#   Y1.trt<-Z%*%beta1+0.5*(1-0.5)
#   Y1.cnt<-Z%*%omega1+0.5*(0-0.5)
#   Y2.trt<-Z%*%beta2+0.25*(1-0.5)
#   Y2.cnt<-Z%*%omega2+0.25*(0-0.5)
#   Y1<-ifelse(A==1, Y1.trt, Y1.cnt)
#   Y2<-ifelse(A==1, Y2.trt, Y2.cnt)
#   
#   #- Generate indicator for cause
#   theta <- 1 + rbinom(nobs, 1, prob = (1 - p)^exp(Y1))
#   # table(theta)
#   #Conditional on cause indicators, we simulate the model.
#   U <-runif(Y1[theta==1],0,1)
#   # T1<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
#   T1.cnt<- -log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p)
#   T1.trt<- (-(1/scale1)*log(1 - (1 - (1 - U * (1 - (1 - p)^exp(Y1[theta==1] )))^(1 / exp(Y1[theta==1] ))) / p))^(1/shape) # event type 1
#   summary(T1.cnt)
#   summary(T1.trt)
#   T1<-ifelse(A[theta==1]==1,T1.trt,T1.cnt)
#   T2 <- rexp(length(Y2[theta==2]), rate = exp(Y2[theta==2] )) # event type 1
#   C <- (-log(runif(nobs,0,1))/(lambda0*exp(Z%*%matrix(censor_beta)))) #simulate censoring times
#   # C<-runif(nobs, min = min(T1,T2), max = quantile(c(T1,T2),0.95))  
#   
#   event.time<-numeric(nobs)
#   event.time[theta==1]<-T1
#   event.time[theta==2]<-T2
#   event.time<-pmin(event.time,C)
#   event.type<-ifelse(event.time==C,0,1)
#   event.type<-event.type*theta
#   return(data.frame(A=A,event.time=event.time,event.type=event.type,Z=Z,Inter=Inter))
# }

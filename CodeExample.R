rm(list = ls())
source("functions/SimuData2model.R") # function for generating data
source("functions/SubgroupATE_CR_TMLE.R") # functions for estimating Subgroup ATE


########################################################################
###### Parameter used for generating data 
########################################################################
seed=2024
n.train=3000 # sample size of the simulated dataset
beta1=c(-0.4,-0.9,0.7,-0.1,0.1,0.1,-0.1,0,0,0); beta2=beta1/2; # Coefficients of baseline covaraites in the model for generating event time
omega1=c(0.6,-0.9,0,-0.1,0.1,0.1,-0.1,0,0,0);omega2=omega1/2; #  Coefficients of baseline covaraites in the model for generating event time
p=0.7;exp=0 # Parameters for determining event time
gamma=c(-0.2,1.5,-0.2,0.1,-0.1,0.1,0,-0.1,0,0); # Coefficients of baseline covaraites in the model for treatment
censor_beta=c(0.2,-0.2,-0.2,-0.1,-0.1,0.1,0,0,0.1,0);lambda0=0.3 #  Coefficients of baseline covaraites in the model for censoring 
########################################################################
###### Generate Data 
########################################################################
set.seed(seed);
data<-FG.sim.10(n.train, beta1, beta2, omega1, omega2, gamma, scaleC, exp, p,lambda0,censor_beta) 


########################################################################
###### Model for estimating subgroup ATE
########################################################################
# Algorithm used in SuperLearner for survival outcomes
treat.SL.library <-  c('SL.glm','SL.gam','SL.glmnet')
censor.SL.library <- lapply(c("survSL.coxph", "survSL.expreg",  "survSL.loglogreg"), function(alg) {
  c(alg, "survscreen.glmnet")
})
# Time percentile of interests
tq= c(0.25,0.5,0.75)
## Model input 
# Predictive variables used for dertermining subgroups
pred.var=c("Z.1","Z.3")
# Prognostic variables used in the model
prog.var=c("Z.2","Z.4","Z.5","Z.6","Z.7")
# Variables used for modeling treatment
treat.var=c("Z.2","Z.4","Z.5","Z.6","Z.8")
# Variables used for modeling censoring
censor.var=c("Z.2","Z.4","Z.5","Z.6","Z.9")
########################################################################
###### Execution
########################################################################
tp<-quantile(data$event.time,probs = tq) #Time corresponding to the percentile of interests
time<-sort(data$event.time[which(data$event.type==1)]) # Unique event time in the data

#true subgroup ATE and ITE of data
truth<-true_cate(tp=tp, dat=data,beta1, beta2, omega1, omega2, gamma, p)
trueITE<-truth$ite 
trueCATE<-truth$cate

#TMLE, with S-learner for the initial estimate, subdistribution hazard model for CIF
TMLE.S<-TMLE_FG(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                treat.var=treat.var,censor.var=censor.var,learner="S",SL=FALSE,LS=FALSE,ITE.out=FALSE)
#TMLE, with T-learner for the initial estimate, subdistribution hazard model for CIF
TMLE.T<-TMLE_FG(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                treat.var=treat.var,censor.var=censor.var,learner="T",SL=FALSE,LS=FALSE,ITE.out=FALSE)
#TMLE, with S-learner for the initial estimate, with LASSO-penalized subdistribution hazard model for CIF
TMLE.S.LS<-TMLE_FG(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                   treat.var=treat.var,censor.var=censor.var,learner="S",SL=FALSE,LS=TRUE,ITE.out=FALSE)
#TMLE, with T-learner for the initial estimate, with LASSO-penalized subdistribution hazard model for CIF
TMLE.T.LS<-TMLE_FG(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                   treat.var=treat.var,censor.var=censor.var,learner="T",SL=FALSE,LS=TRUE,ITE.out=FALSE)
#TMLE, with S-learner for the initial estimate, Super-learner for treatment and censoring
TMLE.S.SL<-TMLE_FG(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                   treat.var=treat.var,censor.var=censor.var,learner="S",SL=TRUE,LS=FALSE,ITE.out=FALSE)
#TMLE, with T-learner for the initial estimate, Super-learner for treatment and censoring
TMLE.T.SL<-TMLE_FG(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                   treat.var=treat.var,censor.var=censor.var,learner="T",SL=TRUE,LS=FALSE,ITE.out=FALSE)
#TMLE, with S-learner for the initial estimate, with LASSO-penalized subdistribution hazard model for CIF, Super-learner for treatment and censoring
TMLE.S.SL.LS<-TMLE_FG(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                      treat.var=treat.var,censor.var=censor.var,learner="S",SL=TRUE,LS=TRUE,ITE.out=FALSE)
#TMLE, with T-learner for the initial estimate, with LASSO-penalized subdistribution hazard model for CIF, Super-learner for treatment and censoring
TMLE.T.SL.LS<-TMLE_FG(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                      treat.var=treat.var,censor.var=censor.var,learner="T",SL=TRUE,LS=TRUE,ITE.out=FALSE)

#S-Learner, with subdistribution hazard model for CIF
SL<-FG_S(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
         treat.var=treat.var,censor.var=censor.var, LS=FALSE,ITE.out=FALSE)
#T-Learner, with subdistribution hazard model for CIF
TL<-FG_T(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
         treat.var=treat.var,censor.var=censor.var, LS=FALSE,ITE.out=FALSE)
#S-learner, with LASSO-penalized subdistribution hazard model for CIF
SL.LS<-FG_S(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
            treat.var=treat.var,censor.var=censor.var, LS=TRUE,ITE.out=FALSE)
#T-learner, with LASSO-penalized subdistribution hazard model for CIF
TL.LS<-FG_T(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
            treat.var=treat.var,censor.var=censor.var, LS=TRUE,ITE.out=FALSE)
#One-step Estimator, with S-learner for the initial estimate, subdistribution hazard model for CIF
OS.S<-one_step(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
               treat.var=treat.var,censor.var=censor.var,learner="S",SL=FALSE,LS=FALSE,ITE.out=FALSE)
#One-step Estimator, with T-learner for the initial estimate, with subdistribution hazard model for CIF
OS.T<-one_step(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
               treat.var=treat.var,censor.var=censor.var,learner="T",SL=FALSE,LS=FALSE,ITE.out=FALSE)
#One-step Estimator, with S-learner for the initial estimate, LASSO-penalized subdistribution hazard model for CIF
OS.S.LS<-one_step(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                  treat.var=treat.var,censor.var=censor.var,learner="S",SL=FALSE,LS=TRUE,ITE.out=FALSE)
#One-step Estimator, with T-learner for the initial estimate, LASSO-penalized subdistribution hazard model for CIF
OS.T.LS<-one_step(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                  treat.var=treat.var,censor.var=censor.var,learner="T",SL=FALSE,LS=TRUE,ITE.out=FALSE)
#One-step Estimator, with S-learner for the initial estimate, Super-learner for treatment and censoring
OS.S.SL<-one_step(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                  treat.var=treat.var,censor.var=censor.var,learner="S",SL=TRUE,LS=FALSE,ITE.out=FALSE)
#One-step Estimator, with T-learner for the initial estimate, Super-learner for treatment and censoring
OS.T.SL<-one_step(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                  treat.var=treat.var,censor.var=censor.var,learner="T",SL=TRUE,LS=FALSE,ITE.out=FALSE)
#One-step Estimator, with S-learner for the initial estimate, Super-learner for treatment and censoring, LASSO-penalized subdistribution hazard model for CIF
OS.S.SL.LS<-one_step(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                     treat.var=treat.var,censor.var=censor.var,learner="S",SL=TRUE,LS=TRUE,ITE.out=FALSE)
#One-step Estimator, with T-learner for the initial estimate, Super-learner for treatment and censoring, LASSO-penalized subdistribution hazard model for CIF
OS.T.SL.LS<-one_step(data=data, tp=tp, event.time="event.time",event.type="event.type", treat="A",pred.var=pred.var, prog.var=prog.var,
                     treat.var=treat.var,censor.var=censor.var,learner="T",SL=TRUE,LS=TRUE,ITE.out=FALSE)

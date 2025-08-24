### Estimating CATEs in Competing Risk Data--- 
#----------------------------------------------------------------------
## Author: Runjia Li 
## Created: 
## Version: 
## Last-Updated: 
#----------------------------------------------------------------------

#### Load required packages
library(BuyseTest) # perform GPC analysis with inference
library(cmprsk) # perform Wald test based on Fine-Gray model
library(FSA) # plot bias results (Figures 1-4)
library(ggplot2) # plot bias results (Figures 1-4)
library(ggpubr) # plot bias results (Figures 1-4)
library(gridExtra) # plot bias results (Figures 1-4)
# library(tcltk) # display progression bar during simulations
library(xtable) # put power results in Latex table (Table 3)
library(table1)
library(survival)
library(survminer)
library(dplyr)
library(fastcmprsk)
library(ggplot2)
library(gbm)
library(randomForest)
library(reshape)
library(latex2exp)
library(riskRegression)
library(earth)
library(data.table)
library(zoo)
library(nleqslv)
library(doParallel)
library(doRNG)
library(SuperLearner)
library(survSuperLearner)
library(reshape2)

# --------------------------------------------------------------------------
# Simulate data 
# --------------------------------------------------------------------------
# truth
true_cate<-function(tp, data.test, beta1, beta2, omega1, omega2, gamma, p){
  Z<-as.matrix(data.test[,which(substr(colnames(data.test), 1, 2)== "Z.")])
  A<-data.test$A
  Inter<-as.matrix(data.test[,which(substr(colnames(data.test), 1, 5)== "Inter")])
  Y1_1<-Z%*%beta1+0.5*(1-0.5)
  Y1_0<-Z%*%omega1+0.5*(0-0.5)
  
  trueCIF_1<-sapply(tp,function(x) (1-(1-p*(1-exp(-scale1*(x^shape))))^exp(Y1_1)))
  trueCIF_0<-sapply(tp,function(x) (1-(1-p*(1-exp(-x)))^exp(Y1_0)))
  
  trueITE<-trueCIF_1-trueCIF_0

  imp.Z.bin<-which(beta1!=omega1)
  CATE.lst<-as.data.frame(cbind(trueITE,Z[,imp.Z.bin]))
  colnames(CATE.lst)<-c(colnames(trueITE),colnames(Z)[imp.Z.bin])
  CATE.matrix<-aggregate(as.formula(paste(".~",paste(paste0(rep("Z.",length(imp.Z.bin)),c(imp.Z.bin)),collapse= "+"))),data=CATE.lst,mean)
  colnames(CATE.matrix)<-c(paste0(rep("Z.",length(imp.Z.bin)),imp.Z.bin),paste0(rep("DF_",length(tp)),c(1:length(tp))))

  return(list(ite=CATE.lst, cate=CATE.matrix,  trueCIF_1=trueCIF_1, trueCIF_0=trueCIF_0))
}

TMLE_FG<-function(data.train, data.test, tp=tp, Q=2,learner,Z_ind=NULL, Inter_ind=NULL, 
                  treat_Z_ind=NULL,censor_Z_ind=NULL,SL=FALSE,LS=FALSE,ITE=FALSE){
  data.train.all<-data.train
  data.train$id<-c(1:nrow(data.train))
  train.lst<-data.train.all%>%group_by(data.train.all[,which(colnames(data.train.all)%in%paste0("Z.",Inter_ind))])%>%group_split(.)
  ITE.train<-NULL
  for (l in 1: length(train.lst)){
    data.train<-train.lst[[l]]; data.train$id<-c(1:nrow(data.train))
    Z.train<-as.matrix(data.train[,which(substr(colnames(data.train), 1, 1)== "Z")])
    A.train<-data.train$A
    Inter.train<-as.matrix(data.train[,which(substr(colnames(data.train), 1, 5)== "Inter")])
    unique.times <- c(0,sort(unique(data.train$event.time)))
    bhaz1 <- data.table(event.time=unique.times)
    bhaz0 <- data.table(event.time=unique.times)
    
    ############################### Step 1.1: predict ##################################################
    # S-learner
    if (learner=="S"){
      if (length(Z_ind)==length(Inter_ind)){
        model.fg<-as.formula(paste0("Crisk(event.time, event.type, cencode = 0, failcode = 1)~A"))
      }else{
        model.fg<-as.formula(paste0("Crisk(event.time, event.type, cencode = 0, failcode = 1)~A+",
                                    paste(paste0(rep("Z.",length(Z_ind[which(!Z_ind%in%Inter_ind)])),Z_ind[which(!Z_ind%in%Inter_ind)]),collapse= "+")))
      }
      if (LS==FALSE){ 
        fit.fg <- fastCrr(model.fg,data=data.train,variance = FALSE)
        coef<-fit.fg$coef
        breslow<-fit.fg$breslowJump ## dH(t),CIF=(1-exp(-cumsum(exp(X\beta)*dH(t))))
      }else{
        lam.path <- 10^seq(log10(0.01), log10(0.000001), length = 100)
        weight <- 1/abs(fastCrr(model.fg,data=data.train,variance = FALSE)$coef)
        fit.fg <- fastCrrp(model.fg,data=data.train, penalty="LASSO",lambda = lam.path,penalty.factor=weight)
        coef.ind<-which.min(AIC(fit.fg,k=2))
        coef<-fit.fg$coef[,coef.ind]
        breslow<-fit.fg$breslowJump[,c(1,coef.ind+1)]
      }
      colnames(breslow)<-c("event.time","hazard")
      data.train$expxbeta1<-apply(cbind(1,Z.train[,Z_ind[which(!Z_ind%in%Inter_ind)]]),1,function(x) exp(sum(x * coef)) )
      data.train$expxbeta0<-apply(cbind(0,Z.train[,Z_ind[which(!Z_ind%in%Inter_ind)]]),1,function(x) exp(sum(x * coef)) )
      bhaz1 <- merge(bhaz1, rbind(data.table(event.time=0, hazard=0),breslow),by=c("event.time"), all.x=TRUE)
      bhaz1[, hazard:=ifelse(is.na(hazard),0,hazard)]
      bhaz0 <- merge(bhaz0, rbind(data.table(event.time=0, hazard=0),breslow),by=c("event.time"), all.x=TRUE)
      bhaz0[, hazard:=ifelse(is.na(hazard),0,hazard)]

      t.ind.fg<-vector()
      for (j in 1:length(tp)){t.ind.fg[j]<-ifelse(min(unique.times)>tp[j],NA,max(which(unique.times<=tp[j])))}
    } 
    if (learner=="T"){
      model.fg<-as.formula(paste0("Crisk(event.time, event.type, cencode = 0, failcode = 1)~",
                                  paste(paste0(rep("Z.",length(Z_ind[which(!Z_ind%in%Inter_ind)])),Z_ind[which(!Z_ind%in%Inter_ind)]),collapse= "+")
      ))        
      if (LS==FALSE){ 
        fit.fg.trt<- fastCrr(model.fg,data=data.train[which(data.train$A==1),],variance = FALSE)
        coef.trt<-fit.fg.trt$coef
        breslow.trt<-fit.fg.trt$breslowJump ## dH(t),CIF=(1-exp(-cumsum(exp(X\beta)*dH(t))))
        fit.fg.cnt<- fastCrr(model.fg,data=data.train[which(data.train$A==0),],variance = FALSE)
        coef.cnt<-fit.fg.cnt$coef
        breslow.cnt<-fit.fg.cnt$breslowJump ## dH(t),CIF=(1-exp(-cumsum(exp(X\beta)*dH(t))))
      }else{
        lam.path <- 10^seq(log10(0.01), log10(0.000001), length = 100)
        weight.trt <- 1/abs(fastCrr(model.fg,data=data.train[which(data.train$A==1),],variance = FALSE)$coef)
        weight.cnt <- 1/abs(fastCrr(model.fg,data=data.train[which(data.train$A==0),],variance = FALSE)$coef)
        fit.fg.trt <- fastCrrp(model.fg,data=data.train[which(data.train$A==1),], penalty="LASSO",lambda = lam.path,penalty.factor=weight.trt)
        coef.ind.trt<-which.min(AIC(fit.fg.trt,k=2))
        coef.trt<-fit.fg.trt$coef[,coef.ind.trt]
        breslow.trt<-fit.fg.trt$breslowJump[,c(1,coef.ind.trt+1)]
        fit.fg.cnt <- fastCrrp(model.fg,data=data.train[which(data.train$A==0),], penalty="LASSO",lambda = lam.path,penalty.factor=weight.cnt)
        coef.ind.cnt<-which.min(AIC(fit.fg.cnt,k=2))
        coef.cnt<-fit.fg.cnt$coef[,coef.ind.cnt]
        breslow.cnt<-fit.fg.cnt$breslowJump[,c(1,coef.ind.cnt+1)]
      }
      colnames(breslow.trt)<-colnames(breslow.cnt)<-c("event.time","hazard")
      if (length(which(!Z_ind%in%Inter_ind))==1){
        data.train$expxbeta1<-sapply(Z.train[,Z_ind[which(!Z_ind%in%Inter_ind)]],function(x) exp(sum(x * coef.trt)) )
        data.train$expxbeta0<-sapply(Z.train[,Z_ind[which(!Z_ind%in%Inter_ind)]],function(x) exp(sum(x * coef.cnt)) )
      }else{        
        data.train$expxbeta1<-apply(Z.train[,Z_ind[which(!Z_ind%in%Inter_ind)]],1,function(x) exp(sum(x * coef.trt)) )
        data.train$expxbeta0<-apply(Z.train[,Z_ind[which(!Z_ind%in%Inter_ind)]],1,function(x) exp(sum(x * coef.cnt)) )
      }
      
      bhaz1 <- merge(bhaz1, rbind(data.table(event.time=0, hazard=0),breslow.trt),by=c("event.time"), all.x=TRUE)
      bhaz1[, hazard:=ifelse(is.na(hazard),0,hazard)]
      bhaz1<-cbind(bhaz1,sapply(data.train$expxbeta1,function(x) 1-exp(-cumsum(bhaz1$hazard*x))))
      bhaz0 <- merge(bhaz0, rbind(data.table(event.time=0, hazard=0),breslow.cnt),by=c("event.time"), all.x=TRUE)
      bhaz0[, hazard:=ifelse(is.na(hazard),0,hazard)]
      bhaz0<-cbind(bhaz0,sapply(data.train$expxbeta0,function(x) 1-exp(-cumsum(bhaz0$hazard*x))))
      colnames(bhaz0)[3:ncol(bhaz0)]<-colnames(bhaz1)[3:ncol(bhaz1)]<-paste0("F_",c(1:nrow(data.train)))
      
      t.ind.fg<-vector()
      for (j in 1:length(tp)){t.ind.fg[j]<-ifelse(min(bhaz0$event.time)>tp[j],NA,max(which(bhaz0$event.time<=tp[j])))}
      
    }

    ######################## Step 1.2: Estimate the Probability of Treatment ###############################
    # # Logistic Regression
    if (length(Z_ind)==length(Inter_ind)){
      g_w<-ifelse(data.train$A==1,sum(data.train$A)/nrow(data.train),sum(1-data.train$A)/nrow(data.train))
    }else{
      if (SL==TRUE){
        #Super learner
        Z<-as.data.frame(cbind(Z.train[,treat_ind[!treat_ind%in% Inter_ind]]))
        colnames(Z)<-paste0("Z.",treat_ind[!treat_ind%in% Inter_ind])
        g.sl <- SuperLearner(Y = A.train, # outcome is the A (treatment) vector
                             X = Z, # Z is a dataframe of predictors
                             family=binomial(), # treatment is a binomial outcome
                             SL.library=treat.SL.library) # using same candidate learners; could use different learners
        g_w.sl <- as.vector(predict(g.sl)$pred) # Pr(A=1|W)
        g_w<-g_w.sl
      } else{
        g.glm<-glm(as.formula(paste("A~", paste(colnames(Z.train)[treat_ind[!treat_ind%in% Inter_ind]], collapse= "+"))),data=data.train, family = "binomial") #checked, same with the formulat for generating A
        Z<-as.data.frame(cbind(Z.train[,treat_ind[!treat_ind%in% Inter_ind]]))
        colnames(Z)<-paste0("Z.",treat_ind[!treat_ind%in% Inter_ind])
        g_w.glm<- predict(g.glm,Z,"response")# Pr(A=1|W)
        g_w<-g_w.glm
      }
    }
    H_1<--1/(g_w); H_0<-1/(1-g_w)
    H_A<-ifelse(A.train==1,H_1, H_0)
    ######################## Step 1.3: Censoring and hazard  ###############################
    # censoring at all time points for A=1 and A=0
    if (SL==FALSE){
      if (is.null(censor_Z_ind[!treat_ind%in% Inter_ind])){
       censor_model<-as.formula("Surv(event.time, event.type==0)~1")} else{                 
       censor_model<-as.formula(paste("Surv(event.time, event.type==0)~",paste(paste0(rep("Z.",length(censor_Z_ind[!treat_ind%in% Inter_ind])),censor_Z_ind[!treat_ind%in% Inter_ind]),collapse= "+")))
       }
      fit.cox0 <- coxph(censor_model, data=data.train)
      exp.Xbeta0<-predict(fit.cox0,type="risk")
      
      tmp0 <- suppressWarnings(setDT(basehaz(fit.cox0, centered=FALSE)))
      colnames(tmp0)<-c("hazard","event.time")
      bhaz.cox <- data.table(event.time=unique.times)
      bhaz.cox <- merge(bhaz.cox, rbind(data.table(event.time=0, hazard=0),tmp0),
                        by=c("event.time"), all.x=TRUE)
      bhaz.cox[, hazard0:=na.locf(hazard)]
      bhaz.cox[, dhaz0:=c(0, diff(hazard0))]
      setnames(bhaz.cox, "hazard0", "chaz0")
      bhaz.cox[, chaz0.1:=c(0, chaz0[-.N])]
      surv.C<-exp(-matrix(bhaz.cox$chaz0.1)%*%exp.Xbeta0)
    }else{
      fit.G<-survSuperLearner(time = data.train$event.time, event = data.train$event.type==0,
                              X = data.train%>%dplyr::select(paste0(rep("Z.",length(censor_Z_ind)),censor_Z_ind)),
                              newX = data.train%>%dplyr::select(paste0(rep("Z.",length(censor_Z_ind)),censor_Z_ind)),
                              new.times = unique.times, # unique.list
                              event.SL.library =censor.SL.library, cens.SL.library = censor.SL.library, verbose = FALSE)
      surv.C<-t(fit.G$event.SL.predict)
    }
    ########### Step 2.1: Update the initial estimate  ###############
    ccpart1<-ccpart1_1<-ccpart1_0<-list()
    for (j in 1:length(tp)){
      ccpart1[[j]]<-matrix(rep(H_A,length(unique.times)),nrow=length(unique.times),byrow=TRUE)/surv.C
      ccpart1_1[[j]]<-matrix(rep(H_1,length(unique.times)),nrow=length(unique.times),byrow=TRUE)/surv.C
      ccpart1_0[[j]]<-matrix(rep(H_0,length(unique.times)),nrow=length(unique.times),byrow=TRUE)/surv.C
    }

    ############################
    # initial haza and F
    # haza is the matrix of hazard (\lambda_t,id=\lambda_0,t,id * exp(X_id*\beta))
    w.Y.1<-do.call("rbind",lapply(unique.times,function(tt) (tt<=data.train$event.time)))
    w.Y.2<-do.call("rbind",lapply(unique.times,function(tt) (tt>data.train$event.time)*(data.train$event.type==2)))*
      surv.C/matrix(rep(diag(surv.C[-1,]),nrow(surv.C)),byrow = TRUE,nrow = nrow(surv.C))
    w.Y<-w.Y.1+w.Y.2

    eq.t.obs<-do.call("rbind",lapply(unique.times,function(tt) (tt==data.train$event.time)))
    bhaz.m1<-matrix(rep(bhaz1$hazard,nrow(data.train)),byrow=FALSE,ncol=nrow(data.train))
    bhaz.m0<-matrix(rep(bhaz0$hazard,nrow(data.train)),byrow=FALSE,ncol=nrow(data.train))
    
    haz1.0<-matrix(bhaz1$hazard)%*%data.train$expxbeta1
    haz0.0<-matrix(bhaz0$hazard)%*%data.train$expxbeta0
    haza.0<-haz1.0; haza.0[,A.train==0]=haz0.0[,A.train==0]
    F1.0<-apply(haz1.0, 2,function(x) 1-exp(-cumsum(x)))
    F0.0<-apply(haz0.0, 2,function(x) 1-exp(-cumsum(x)))
    Fa.0<-F1.0; Fa.0[,A.train==0]=F0.0[,A.train==0]
    ATE.0<-apply(F1.0[t.ind.fg,]-F0.0[t.ind.fg,],1,mean)
    F1_final<-F0_final<- list();ATE_final<-NULL
    F1_tp<-F0_tp<-EIF<-NULL
    
    start_time<- Sys.time()
    for (j in 1:length(tp)){
      tk=t.ind.fg[j]; tk.time=tp[j]
      haz1<-haz1.0;haz0<-haz0.0; haza<-haza.0;F1<-F1.0;F0<-F0.0;Fa<-Fa.0
      le.tk.time<-matrix(rep((unique.times<=tk.time) ,nrow(data.train)),byrow=FALSE,ncol=nrow(data.train))
      event1<-matrix(rep(data.train$event.type==1,length(unique.times)),byrow=TRUE,nrow=length(unique.times))
      exp.update1<-matrix(rep(data.train$expxbeta1,length(unique.times)),byrow = TRUE,nrow=length(unique.times))
      exp.update0<-matrix(rep(data.train$expxbeta0,length(unique.times)),byrow = TRUE,nrow=length(unique.times))
      eps.hat<-1
      while (abs(eps.hat)>=0.001){
        h1<-ccpart1_1[[j]]*(matrix(rep(1-F1[tk,],nrow(F1)),nrow=nrow(F1),byrow=TRUE)/(1-F1))
        h0<-ccpart1_0[[j]]*(matrix(rep(1-F0[tk,],nrow(F0)),nrow=nrow(F0),byrow=TRUE)/(1-F0))
        h1[F1>0.99]=ccpart1_1[[j]][F1>0.99];h0[F0>0.99]=ccpart1_0[[j]][F0>0.99]; #no big difference
        h1[h1>500]=500;h0[h0>500]=500;h1[h1<(-500)]=-500;h0[h0<(-500)]=-500;
        ha<-h1;ha[,A.train==0]<-h0[,A.train==0]
        
        # each id has one equation
        eps.hat<-nleqslv(0.01, function(eps){
          mean(colSums(le.tk.time*w.Y* ha*
                         (event1*eq.t.obs-haza*exp(eps*ha))))# score function of each subject
        })$x
        # update haza and F
        exp.update1<-exp.update1*exp(eps.hat*h1)
        exp.update0<-exp.update0*exp(eps.hat*h0)
        exp.update1[exp.update1>100]=100;exp.update0[exp.update0>100]=100 #no big difference
        
        haz1<-bhaz.m1*exp.update1
        haz0<-bhaz.m0*exp.update0
        haza<-haz1;haza[,A.train==0]=haz0[,A.train==0]
        F1<-apply(haz1, 2,function(x) 1-exp(-cumsum(x))); 
        F0<-apply(haz0, 2,function(x) 1-exp(-cumsum(x)))
        Fa<-apply(haza, 2,function(x) 1-exp(-cumsum(x)))
        ATE<-mean(F1[tk,]-F0[tk,])
      }
      F1_final[[j]]<-F1;   F0_final[[j]]<-F0; 
      ATE_final<-c(ATE_final,ATE)
      
      F1_tp<-cbind(F1_tp,F1_final[[j]][t.ind.fg[j],])
      F0_tp<-cbind(F0_tp,F0_final[[j]][t.ind.fg[j],])
      h1<-ccpart1_1[[j]]*(matrix(rep(1-F1_final[[j]][t.ind.fg[j],],nrow(F1)),nrow=nrow(F1),byrow=TRUE)/(1-F1_final[[j]]))
      h0<-ccpart1_0[[j]]*(matrix(rep(1-F0_final[[j]][t.ind.fg[j],],nrow(F0)),nrow=nrow(F0),byrow=TRUE)/(1-F0_final[[j]]))
      h1[F1_final[[j]]>0.99]=ccpart1_1[[j]][F1_final[[j]]>0.99];h0[F0_final[[j]]>0.99]=ccpart1_0[[j]][F0_final[[j]]>0.99]; #no big difference
      h1[h1>500]=500;h0[h0>500]=500;h1[h1<(-500)]=-500;h0[h0<(-500)]=-500;
      ha<-h1;ha[,A.train==0]<-h0[,A.train==0]
      EIF<-cbind(EIF,colSums(le.tk.time*w.Y* ha*(event1*eq.t.obs-haza))+
                   (F1_final[[j]][t.ind.fg[j],]-F0_final[[j]][t.ind.fg[j],])-mean((F1_final[[j]][t.ind.fg[j],]-F0_final[[j]][t.ind.fg[j],])))
    }
    end_time2 <- Sys.time()
    (arraytime2<- start_time-end_time2) # 4s
    
    ############################### Step 3.1: ITE ###############################
    data.train<-cbind(data.train,F1_tp,F0_tp,F1_tp-F0_tp,EIF)
    colnames(data.train)[(ncol(data.train)-11):ncol(data.train)]<-
      c(paste0("F1_",1:length(tp)),paste0("F0_",1:length(tp)),paste0("DF_",1:length(tp)),
        paste0("EIF_",1:length(tp)))
    ITE.train<-rbind(ITE.train,data.train)
  }
  ############################### Step 3.2: CATE
  imp.Z.bin<-which(beta1!=omega1)
  CATE.train.matrix<-merge(aggregate(as.formula(paste(".~",paste(paste0(rep("Z.",length(imp.Z.bin)),c(imp.Z.bin)),collapse= "+"))),
                               data=ITE.train[,which(colnames(ITE.train)%in%c(paste0("Z.",imp.Z.bin),paste0("DF_",c(1:length(tp)))))],mean),
                           aggregate(as.formula(paste(".~",paste(paste0(rep("Z.",length(imp.Z.bin)),c(imp.Z.bin)),collapse= "+"))),
                                     data=ITE.train[,which(colnames(ITE.train)%in%c(paste0("Z.",imp.Z.bin),paste0("EIF_",c(1:length(tp)))))],
                                     function(xx) sqrt(mean(xx^2)/length(xx))),
                           by=c(paste0(rep("Z.",length(imp.Z.bin)),c(imp.Z.bin))))
  colnames(CATE.train.matrix)<-c(paste0(rep("Z.",length(imp.Z.bin)),imp.Z.bin),
                                 paste0(rep("DF_",length(tp)),c(1:length(tp))),
                                 paste0(rep("SE_",length(tp)),c(1:length(tp))))
  if (ITE==TRUE){
    return(list(cate=CATE.train.matrix,ITE=ITE.train))
  }else{
    return(list(cate=CATE.train.matrix))
  }
  
}

FG_S<-function(data.train, data.test, tp=tp,Z_ind=NULL, Inter_ind=NULL,LS=FALSE,SL=FALSE,ITE=FALSE){ # ? time points
  # model
  Z.train<-as.matrix(data.train[,which(substr(colnames(data.train), 1, 1)== "Z")])
  A.train<-data.train$A
  Inter.train<-as.matrix(data.train[,which(substr(colnames(data.train), 1, 5)== "Inter")])
  model.fg<-as.formula(paste0("Crisk(event.time, event.type, cencode = 0, failcode = 1)~A+",
                              paste(paste0(rep("Z.",length(Z_ind)),c(Z_ind)),collapse= "+"),
                              "+",paste(paste0(rep("Inter.",length(Inter_ind)),c(Inter_ind)),collapse= "+")))
  if (LS==FALSE){ 
    fit.fg <- fastCrr(model.fg,data=data.train,variance = FALSE)
    coef<-fit.fg$coef
    breslow<-fit.fg$breslowJump[,2]
  }else{
    lam.path <- 10^seq(log10(0.01), log10(0.000001), length = 100)
    weight <- 1/abs(fastCrr(as.formula(paste0("Crisk(event.time, event.type, cencode = 0, failcode = 1)~A+",
                                              paste(paste0(rep("Z.",length(Z_ind)),c(Z_ind)),collapse= "+"),
                                              "+",paste(paste0(rep("Inter.",length(Inter_ind)),c(Inter_ind)),collapse= "+"))),data=data.train,variance=FALSE)$coef)
    fit.fg <- fastCrrp(model.fg,data=data.train, penalty="LASSO",lambda = lam.path,penalty.factor=weight)
    coef.ind<-which.min(AIC(fit.fg,k=2))
    coef<-fit.fg$coef[,coef.ind]
    breslow<-fit.fg$breslowJump[,coef.ind+1]
  }
  
  time.fg<-fit.fg$breslowJump[,1]
  t.ind<-vector()
  for (j in 1:length(tp)){t.ind[j]<-max(which(time.fg<=tp[j]))}
  
  Z.test<-as.matrix(data.test[,which(substr(colnames(data.test), 1, 1)== "Z")])
  A.test<-data.test$A
  Inter.test<-as.matrix(data.test[,which(substr(colnames(data.test), 1, 5)== "Inter")])
  
  S.CIF.hat.trt<-apply(cbind(1,Z.test[,Z_ind],Z.test[,Inter_ind]),1,function(x) (1 - exp(-cumsum(exp(sum(x * coef)) * breslow)))) #If every observation received the treatment.
  S.CIF.hat.cnt<-apply(cbind(0,Z.test[,Z_ind],0*Z.test[,Inter_ind]),1,function(x) (1 - exp(-cumsum(exp(sum(x * coef)) * breslow)))) # If every observation received the control.
  imp.Z.bin<-which(beta1!=omega1)
  
  S.CIF.trt<-t(S.CIF.hat.trt[t.ind,])
  S.CIF.cnt<-t(S.CIF.hat.cnt[t.ind,])
  
  ite.s<-S.CIF.trt-S.CIF.cnt
  CATE.s.lst<-as.data.frame(cbind(ite.s,Z.test[,imp.Z.bin]))
  colnames(CATE.s.lst)<-c(paste0(rep("DF_",length(tp)),c(1:length(tp))),colnames(Z.test)[imp.Z.bin])
  CATE.s.matrix<-merge(aggregate(as.formula(paste(".~",paste(paste0(rep("Z.",length(imp.Z.bin)),c(imp.Z.bin)),collapse= "+"))),
                           data=CATE.s.lst[,which(colnames(CATE.s.lst)%in%c(paste0("Z.",imp.Z.bin),paste0("DF_",c(1:length(tp)))))],mean),
              aggregate(as.formula(paste(".~",paste(paste0(rep("Z.",length(imp.Z.bin)),c(imp.Z.bin)),collapse= "+"))),
                        data=CATE.s.lst[,which(colnames(CATE.s.lst)%in%c(paste0("Z.",imp.Z.bin),paste0("DF_",c(1:length(tp)))))],
                        function(xx) sqrt(var(xx))), 
              by=c(paste0(rep("Z.",length(imp.Z.bin)),c(imp.Z.bin))))
  colnames(CATE.s.matrix)<-c(paste0(rep("Z.",length(imp.Z.bin)),imp.Z.bin),
                             paste0(rep("DF_",length(tp)),c(1:length(tp))),
                             paste0(rep("SE_",length(tp)),c(1:length(tp))))

  if (ITE==TRUE){
    return(list(cate=CATE.s.matrix,ITE=ite.s))
  }else{
    return(list(cate=CATE.s.matrix))
  }
} 


FG_T<-function(data.train, data.test, tp=tp,Z_ind=NULL,Inter_ind=NULL,LS=FALSE,SL=FALSE,ITE=FALSE){ # ? time points
  T.data.train<-data.train[which(data.train$A==1),]
  C.data.train<-data.train[which(data.train$A==0),]
  Z.train<-as.matrix(data.train[,which(substr(colnames(data.train), 1, 1)== "Z")])
  T.Z.train<-as.matrix(T.data.train[,which(substr(colnames(T.data.train), 1, 1)== "Z")])
  T.A.train<-T.data.train$A
  T.Inter.train<-as.matrix(T.data.train[,which(substr(colnames(T.data.train), 1, 5)== "Inter")])
  C.Z.train<-as.matrix(C.data.train[,which(substr(colnames(C.data.train), 1, 1)== "Z")])
  C.A.train<-C.data.train$A
  Inter.C.train<-as.matrix(C.data.train[,which(substr(colnames(C.data.train), 1, 5)== "Inter")])
  
  # model
  model.fg<-as.formula(paste0("Crisk(event.time, event.type, cencode = 0, failcode = 1)~",
                              paste(paste0(rep("Z.",length(Z_ind)),c(Z_ind)),collapse= "+")
  ))
  coef.trt<-coef.cnt<-rep(0,ncol(Z.train))
  
  if (LS==FALSE){ 
    fit.fg.trt<- fastCrr(model.fg,data=data.train[which(data.train$A==1),],variance = FALSE)
    breslow.trt<-fit.fg.trt$breslowJump[,2]
    fit.fg.cnt<- fastCrr(model.fg,data=data.train[which(data.train$A==0),],variance = FALSE)
    breslow.cnt<-fit.fg.cnt$breslowJump[,2]
    coef.trt[Z_ind]<-fit.fg.trt$coef
    coef.cnt[Z_ind]<-fit.fg.cnt$coef
  }else{
    lam.path <- 10^seq(log10(0.01), log10(0.000001), length = 100)
    weight.trt <- 1/abs(fastCrr(as.formula(paste0("Crisk(event.time, event.type, cencode = 0, failcode = 1)~",
                                                  paste(paste0(rep("Z.",length(Z_ind)),c(Z_ind)),collapse= "+"))),data=data.train[data.train$A==1,],variance=FALSE)$coef)
    weight.cnt <- 1/abs(fastCrr(as.formula(paste0("Crisk(event.time, event.type, cencode = 0, failcode = 1)~",
                                                  paste(paste0(rep("Z.",length(Z_ind)),c(Z_ind)),collapse= "+"))),data=data.train[data.train$A==0,],variance=FALSE)$coef)
    
    fit.fg.trt<- fastCrrp(model.fg,data=data.train[data.train$A==1,], penalty="LASSO",lambda = lam.path,penalty.factor=weight.trt)
    coef.ind.trt<-which.min(AIC(fit.fg.trt,k=2))
    coef.trt[Z_ind]<-fit.fg.trt$coef[,coef.ind.trt]
    breslow.trt<-fit.fg.trt$breslowJump[,coef.ind.trt+1]
    fit.fg.cnt<- fastCrrp(model.fg,data=data.train[data.train$A==0,], penalty="LASSO",lambda = lam.path,penalty.factor=weight.cnt)
    coef.ind.cnt<-which.min(AIC(fit.fg.cnt,k=2))
    coef.cnt[Z_ind]<-fit.fg.cnt$coef[,coef.ind.cnt]
    breslow.cnt<-fit.fg.cnt$breslowJump[,coef.ind.cnt+1]
  }
  time.fg.trt<-fit.fg.trt$breslowJump[,1]
  time.fg.cnt<-fit.fg.cnt$breslowJump[,1]
  t.ind.trt<-vector()
  t.ind.cnt<-vector()
  
  for (j in 1:length(tp)){
    t.ind.trt[j]<-max(which(time.fg.trt<=tp[j]))
    t.ind.cnt[j]<-max(which(time.fg.cnt<=tp[j]))
  }
  
  #### predict test data
  Z.test<-as.matrix(data.test[,which(substr(colnames(data.test), 1, 1)== "Z")])
  A.test<-data.test$A
  Inter.test<-as.matrix(data.test[,which(substr(colnames(data.test), 1, 5)== "Inter")])
  
  T.CIF.hat.trt<-apply(cbind(Z.test),1,function(x) (1 - exp(-cumsum(exp(sum(x * coef.trt)) * breslow.trt)) )) #If every observation received the treatment.
  T.CIF.hat.cnt<-apply(cbind(Z.test),1,function(x) (1 - exp(-cumsum(exp(sum(x * coef.cnt)) * breslow.cnt)))) # If every observation received the control.
  imp.Z.bin<-which(beta1!=omega1)
  
  T.CIF.trt<-t(T.CIF.hat.trt[t.ind.trt,])
  T.CIF.cnt<-t(T.CIF.hat.cnt[t.ind.cnt,])
  ite.t<-(T.CIF.trt-T.CIF.cnt)
  imp.Inter<-which(omega1!=0)
  CATE.t.lst<-as.data.frame(cbind(ite.t,Z.test[,imp.Z.bin]))
  colnames(CATE.t.lst)<-c(paste0(rep("DF_",length(tp)),c(1:length(tp))),colnames(Z.test)[imp.Z.bin])
  
  CATE.t.matrix<-merge(aggregate(as.formula(paste(".~",paste(paste0(rep("Z.",length(imp.Z.bin)),c(imp.Z.bin)),collapse= "+"))),
                                 data=CATE.t.lst[,which(colnames(CATE.t.lst)%in%c(paste0("Z.",imp.Z.bin),paste0("DF_",c(1:length(tp)))))],mean),
                       aggregate(as.formula(paste(".~",paste(paste0(rep("Z.",length(imp.Z.bin)),c(imp.Z.bin)),collapse= "+"))),
                                 data=CATE.t.lst[,which(colnames(CATE.t.lst)%in%c(paste0("Z.",imp.Z.bin),paste0("DF_",c(1:length(tp)))))],
                                 function(xx) sqrt(var(xx))),
                       by=c(paste0(rep("Z.",length(imp.Z.bin)),c(imp.Z.bin))))
  colnames(CATE.t.matrix)<-c(paste0(rep("Z.",length(imp.Z.bin)),imp.Z.bin),
                             paste0(rep("DF_",length(tp)),c(1:length(tp))),
                             paste0(rep("SE_",length(tp)),c(1:length(tp))))

  if (ITE==TRUE){
    return(list(cate=CATE.t.matrix,ITE=ite.t))
  }else{
    return(list(cate=CATE.t.matrix))
  }
} 

#############################################################
#############################################################
#' Performs the following analyses: 
#' Full cohort; Standard case-cohort (Prentice, Lin-Ying);
#' Case-cohort plus full cohort, using MI approx;  Case-cohort plus full cohort, using MI SMC
#' Case-cohort plus superset1, using MI approx (Prentice, Lin-Ying); 
#' Case-cohort plus superset1, using MI SMC (Prentice, Lin-Ying); Case-cohort plus superset1;
#' Case-cohort plus superset2, using MI approx (Prentice, Lin-Ying); 
#' Case-cohort plus superset2, using MI SMC (Prentice, Lin-Ying); Case-cohort plus superset2.
#' The Lin-Ying option refers to the IPW estimator.
#' 
#' This file differs from analysis_casecohort_interactionF.R only in that the MI-SMC analyses assume that log(X) is conditionally normally distributed instead of X being conditionally normally distributed.
#' This alternative MI-SMC analysis was considered when data are generated using misspec=T.
#############################################################
#############################################################

#===============================
# setup
#===============================

#predictors used in the two MI approaches
predictors_approx = c("z1","z2","d","chaz")
predictors_smcfcs = c("z1","z2")

#number of imputations used in MI analyses
nimp=10

#number of iterations used in MI-approx analyses
nit = 100

#number of iterations used in MI-SMC analyses
nit.smc = 500

#functions for implementing MI-SMC: a simplified version of smcfcs
source("smcfcs_reduced.R")

#===============================
# generate data
# Can choose n=5000 or n=25000
# Must have misspec=T.
# Must have interaction=F.
#===============================

n=5000
misspec=T
interaction=F

source("generate_data_casecohort.R")

#===============================
# Compute Nelson-Aalen estimate of the cumulative hazard for full cohort
# This is used in the imputation analyses below
#===============================

temp=survfit(coxph(Surv(t,d)~1,data=cohort),type="fleming-harrington",censor = T)
temp.step=stepfun(temp$time,c(0,temp$cumhaz))
cohort$chaz=temp.step(cohort$t)

#add cumulative hazard into case-cohort data
cohort.merge<-cohort[,c("id","chaz")] 

cohort.caco<-merge(cohort.caco,cohort.merge,by.x="id")
cohort.super1.caco<-merge(cohort.super1.caco,cohort.merge,by.x="id")
cohort.super2.caco<-merge(cohort.super2.caco,cohort.merge,by.x="id")

#===============================
# cox analysis using full cohort data
#===============================

model=coxph(Surv(t,d)~x+z1+z2,data=cohort)
cohort.coef=model$coefficients
cohort.se=sqrt(diag(model$var))

#===============================
# traditional case-cohort analysis 
#===============================

#Prentice
model=cch(Surv(t,d)~x+z1+z2,data=caco, subcoh=~subco, id=~id, method="Prentice", cohort.size=n)

cacoP.coef = model$coefficients
cacoP.se = sqrt(diag(model$var))

#LinYing
model=cch(Surv(t,d)~x+z1+z2,data=caco, subcoh=~subco, id=~id, method="LinYing", cohort.size=n,robust=T)

cacoLY.coef= model$coefficients
cacoLY.se= sqrt(diag(model$var))

#===============================
#case-cohort
#MI-approx:  full-cohort approach
#===============================

#predictor matrix which determines the imputation model for x
pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
colnames(pred.mat)=names(cohort.caco)
rownames(pred.mat)=names(cohort.caco)
pred.mat["x",predictors_approx]=1

#method of imputation for x
method.vec=rep("",dim(cohort.caco)[2])
method.vec[which(colnames(cohort.caco)=="x")]="norm"

#perform the imputation 
imp<-mice(cohort.caco, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = nit, diagnostics = FALSE, printFlag = F)

# Fit the analysis model in each imputed data set
models<-with(imp,coxph(Surv(t,d)~x+z1+z2))

# Combine estimates across the imputed data sets using Rubin's Rules
summary_aprx = summary(pool(models))

caco.MIapprox.full.coef= summary_aprx[,"estimate"]
caco.MIapprox.full.se= summary_aprx[,"std.error"]

#===============================
#case-cohort
#MI-SMC: full-cohort approach
#===============================

minx<-min(cohort.caco$x,na.rm=T)
cohort.caco$logx<-ifelse(rep(minx,n)<0,log(cohort.caco$x-minx+0.1),log(cohort.caco$x))

#predictor matrix which determines the imputation models for x
pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
colnames(pred.mat)=names(cohort.caco)
rownames(pred.mat)=names(cohort.caco)
pred.mat["logx",predictors_smcfcs]=1

#method of imputation for x
method.vec=rep("",dim(cohort.caco)[2])
method.vec[which(colnames(cohort.caco)=="x")]="exp(logx)"
method.vec[which(colnames(cohort.caco)=="logx")]="norm"

#initial values for x
inits<-complete(imp,1)$x

#perform the imputation
time1<-proc.time()[3]
imp <- smcfcs.red(cohort.caco, smtype="coxph", smformula="Surv(t,d)~x+z1+z2",method=method.vec,
                  predictorMatrix=pred.mat,m = nimp, numit = nit.smc, rjlimit = 10000,noisy=F)

# obtain estimates across the imputed data sets and combine using Rubin's Rules
impobj <- imputationList(imp$impDatasets)
models <- with(impobj, coxph(Surv(t,d)~x+z1+z2))

caco.MIsmc.full.coef = MIcombine(models)$coefficients
caco.MIsmc.full.se= sqrt(diag(MIcombine(models)$variance))

#===============================
#case-cohort with superset (super1)
#MI-approx
#===============================

#predictor matrix which determines the imputation models for x
pred.mat=matrix(0,nrow=dim(cohort.super1.caco)[2],ncol=dim(cohort.super1.caco)[2])
colnames(pred.mat)=names(cohort.super1.caco)
rownames(pred.mat)=names(cohort.super1.caco)
pred.mat["x",predictors_approx]=1

#method of imputation for x
method.vec=rep("",dim(cohort.super1.caco)[2])
method.vec[which(colnames(cohort.super1.caco)=="x")]="norm"

#perform the imputation 
imp<-mice(cohort.super1.caco, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = nit, diagnostics = FALSE, printFlag = F)

#-----
#Prentice

# Fit the analysis model in each imputed data set
models <- vector("list", nimp)
for (k in 1:nimp){
  model=cch(Surv(t,d)~x+z1+z2,data=complete(imp,k),subcoh=~subco.super1, id=~id, method="Prentice", cohort.size=n)
  models[[k]] = model
}

# Combine estimates across the imputed data sets using Rubin's Rules
cacoP.MIapprox.super1.coef = MIcombine(models)$coef
cacoP.MIapprox.super1.se = sqrt(diag(MIcombine(models)$variance))

#-----
#LinYing

# Fit the analysis model in each imputed data set
models <- vector("list", nimp)
for (k in 1:nimp){
  model=cch(Surv(t,d)~x+z1+z2,data=complete(imp,k),subcoh=~subco.super1, id=~id, method="LinYing", cohort.size=n)
  models[[k]] = model
}

# Combine estimates across the imputed data sets using Rubin's Rules
cacoLY.MIapprox.super1.coef = MIcombine(models)$coef
cacoLY.MIapprox.super1.se= sqrt(diag(MIcombine(models)$variance))

#===============================
#case-cohort with superset (super1)
#MI-SMC
#===============================

minx<-min(cohort.super1.caco$x,na.rm=T)
cohort.super1.caco$logx<-ifelse(rep(minx,dim(cohort.super1.caco)[1])<0,log(cohort.super1.caco$x-minx+0.1),log(cohort.super1.caco$x))

#predictor matrix which determines the imputation models for x
pred.mat=matrix(0,nrow=dim(cohort.super1.caco)[2],ncol=dim(cohort.super1.caco)[2])
colnames(pred.mat)=names(cohort.super1.caco)
rownames(pred.mat)=names(cohort.super1.caco)
pred.mat["logx",predictors_smcfcs]=1

#method of imputation for x
method.vec=rep("",dim(cohort.super1.caco)[2])
method.vec[which(colnames(cohort.super1.caco)=="x")]="exp(logx)"
method.vec[which(colnames(cohort.super1.caco)=="logx")]="norm"

my.sampfrac = sum(cohort.super1.caco$subco.super1==1)/n 

#initial values for x
inits<-complete(imp,1)$x

#-----
#Prentice

#perform the imputation
imp <- smcfcs.casecohort.red.Prentice(cohort.super1.caco,smformula="Surv(t,d)~x+z1+z2",sampfrac=my.sampfrac,in.subco="subco.super1",
                                      method=method.vec,cohort.size=n,predictorMatrix=pred.mat,m=nimp,numit=nit.smc,rjlimit=10000,noisy=FALSE)

# Fit the analysis model in each imputed data set
impobj <- imputationList(imp$impDatasets)
models <- vector("list", nimp)
for (k in 1:nimp){
  model <- cch(Surv(t,d)~x+z1+z2,subcoh=~subco.super1, id=~id, method="Prentice", cohort.size=n,data=impobj$imputations[[k]])
  models[[k]] = model
}

cacoP.MIsmc.super1.coef = MIcombine(models)$coefficients
cacoP.MIsmc.super1.se = sqrt(diag(MIcombine(models)$variance))

#-----
#LinYing

#perform the imputation
imp <- smcfcs.casecohort.red.LinYing(cohort.super1.caco,smformula="Surv(t,d)~x+z1+z2",sampfrac=my.sampfrac,in.subco="subco.super1",
                                     method=method.vec,cohort.size=n,predictorMatrix=pred.mat,m=nimp,numit=nit.smc,rjlimit=10000,noisy=FALSE)

# Fit the analysis model in each imputed data set
impobj <- imputationList(imp$impDatasets)
models <- vector("list", nimp)
for (k in 1:nimp){
  model <- cch(Surv(t,d)~x+z1+z2,subcoh=~subco.super1, id=~id, method="LinYing", cohort.size=n,data=impobj$imputations[[k]])
  models[[k]] = model
}

cacoLY.MIsmc.super1.coef = MIcombine(models)$coefficients
cacoLY.MIsmc.super1.se = sqrt(diag(MIcombine(models)$variance))

#===============================
#case-cohort with superset (super2)
#MI-approx
#===============================

#predictor matrix which determines the imputation models for x
pred.mat=matrix(0,nrow=dim(cohort.super2.caco)[2],ncol=dim(cohort.super2.caco)[2])
colnames(pred.mat)=names(cohort.super2.caco)
rownames(pred.mat)=names(cohort.super2.caco)
pred.mat["x",predictors_approx]=1

#method of imputation for x
method.vec=rep("",dim(cohort.super2.caco)[2])
method.vec[which(colnames(cohort.super2.caco)=="x")]="norm"

#perform the imputation 
imp<-mice(cohort.super2.caco, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = nit, diagnostics = FALSE, printFlag = F)

#-----
#Prentice

# Fit the analysis model in each imputed data set
models <- vector("list", nimp)
for (k in 1:nimp){
  model=cch(Surv(t,d)~x+z1+z2,data=complete(imp,k),subcoh=~subco.super2, id=~id, method="Prentice", cohort.size=n)
  models[[k]] = model
}

# Combine estimates across the imputed data sets using Rubin's Rules
cacoP.MIapprox.super2.coef = MIcombine(models)$coef
cacoP.MIapprox.super2.se = sqrt(diag(MIcombine(models)$variance))

#-----
#LinYing

# Fit the analysis model in each imputed data set
models <- vector("list", nimp)
for (k in 1:nimp){
  model=cch(Surv(t,d)~x+z1+z2,data=complete(imp,k),subcoh=~subco.super2, id=~id, method="LinYing", cohort.size=n)
  models[[k]] = model
}

# Combine estimates across the imputed data sets using Rubin's Rules
cacoLY.MIapprox.super2.coef = MIcombine(models)$coef
cacoLY.MIapprox.super2.se = sqrt(diag(MIcombine(models)$variance))


#===============================
#case-cohort with superset (super2)
#MI-SMC
#===============================

minx<-min(cohort.super2.caco$x,na.rm=T)
cohort.super2.caco$logx<-ifelse(rep(minx,dim(cohort.super2.caco)[1])<0,log(cohort.super2.caco$x-minx+0.1),log(cohort.super2.caco$x))

#predictor matrix which determines the imputation models for x
pred.mat=matrix(0,nrow=dim(cohort.super2.caco)[2],ncol=dim(cohort.super2.caco)[2])
colnames(pred.mat)=names(cohort.super2.caco)
rownames(pred.mat)=names(cohort.super2.caco)
pred.mat["logx",predictors_smcfcs]=1

#method of imputation for x
method.vec=rep("",dim(cohort.super2.caco)[2])
method.vec[which(colnames(cohort.super2.caco)=="x")]="exp(logx)"
method.vec[which(colnames(cohort.super2.caco)=="logx")]="norm"

my.sampfrac = sum(cohort.super2.caco$subco.super2==1)/n 

#initial values for x
inits<-complete(imp,1)$x

#-----
#Prentice

#perform the imputation
imp <- smcfcs.casecohort.red.Prentice(cohort.super2.caco,smformula="Surv(t,d)~x+z1+z2",sampfrac=my.sampfrac,in.subco="subco.super2",
                                      method=method.vec,cohort.size=n,predictorMatrix=pred.mat,m=nimp,numit=nit.smc,rjlimit=10000,noisy=FALSE)

# Fit the analysis model in each imputed data set
impobj <- imputationList(imp$impDatasets)
models <- vector("list", nimp)
for (k in 1:nimp){
  model <- cch(Surv(t,d)~x+z1+z2,subcoh=~subco.super2, id=~id, method="Prentice", cohort.size=n,data=impobj$imputations[[k]])
  models[[k]] = model
}

cacoP.MIsmc.super2.coef = MIcombine(models)$coefficients
cacoP.MIsmc.super2.se = sqrt(diag(MIcombine(models)$variance))

#-----
#LinYing

#perform the imputation
imp <- smcfcs.casecohort.red.LinYing(cohort.super2.caco,smformula="Surv(t,d)~x+z1+z2",sampfrac=my.sampfrac,in.subco="subco.super2",
                                     method=method.vec,cohort.size=n,predictorMatrix=pred.mat,m=nimp,numit=nit.smc,rjlimit=10000,noisy=FALSE)

# Fit the analysis model in each imputed data set
impobj <- imputationList(imp$impDatasets)
models <- vector("list", nimp)
for (k in 1:nimp){
  model <- cch(Surv(t,d)~x+z1+z2,subcoh=~subco.super2, id=~id, method="LinYing", cohort.size=n,data=impobj$imputations[[k]])
  models[[k]] = model
}

cacoLY.MIsmc.super2.coef = MIcombine(models)$coefficients
cacoLY.MIsmc.super2.se = sqrt(diag(MIcombine(models)$variance))

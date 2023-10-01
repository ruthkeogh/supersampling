#############################################################
#############################################################
#' Performs the following analyses: 
#' 
#' Full cohort; Standard nested case-control (NCC);
#' NCC plus full cohort, using MI approx;  NCC plus full cohort, using MI SMC
#' 
#' NCC plus superset1, using MI approx; 
#' NCC plus superset1, using MI SMC; 
#' 
#' NCC plus superset1, using MI approx; 
#' NCC plus superset2, using MI SMC. 
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
# Can choose misspec=F or T
# Must have interaction=F. For interaction=T use analysis_ncc_interactionT.R
#===============================

n=5000
misspec=F
interaction=F

source("generate_data_ncc.R")

#===============================
# Compute Nelson-Aalen estimate of the cumulative hazard for full cohort
# This is used in the imputation analyses below
#===============================

temp=survfit(coxph(Surv(t,d)~1,data=cohort),type="fleming-harrington",censor = T)
temp.step=stepfun(temp$time,c(0,temp$cumhaz))
cohort$chaz=temp.step(cohort$t)

#add cumulative hazard into ncc data
cohort.merge<-cohort[,c("id","chaz")] 
cohort.ncc<-merge(cohort.ncc,cohort.merge,by.x="id")

#===============================
# define variable for time of case in ncc.super1, ncc.super2 
#===============================

#order by set then case/controls
ncc.super1<-ncc.super1[order(ncc.super1$setno,-ncc.super1$case),]
ncc.super2<-ncc.super2[order(ncc.super2$setno,-ncc.super2$case),]

ncc$t.case=ifelse(ncc$case==1,ncc$t,NA)
ncc=ncc%>%group_by(setno)%>%fill(t.case,.direction="down")
ncc<-as.data.frame(ncc)

ncc.super1$t.case=ifelse(ncc.super1$case==1,ncc.super1$t,NA)
ncc.super1=ncc.super1%>%group_by(setno)%>%fill(t.case,.direction="down")
ncc.super1<-as.data.frame(ncc.super1)

ncc.super2$t.case=ifelse(ncc.super2$case==1,ncc.super2$t,NA)
ncc.super2=ncc.super2%>%group_by(setno)%>%fill(t.case,.direction="down")
ncc.super2<-as.data.frame(ncc.super2)

#for ncc.super1 and ncc.super2 we use the cumulative hazard at the time of the case
ncc.super1$chaz<-temp.step(ncc.super1$t.case)
ncc.super2$chaz<-temp.step(ncc.super2$t.case)

#===============================
# cox analysis using full cohort data
#===============================

model=coxph(Surv(t,d)~x+z1+z2,data=cohort)
cohort.coef=model$coefficients
cohort.se=sqrt(diag(model$var))

#===============================
# traditional NCC analysis 
#===============================

model=coxph(Surv(t.case,case)~x+z1+z2+strata(setno),data=ncc)
ncc.coef = model$coefficients
ncc.se = sqrt(diag(model$var))

#===============================
#NCC
#MI-approx:  full-cohort approach
#===============================

#predictor matrix which determines the imputation model for x
pred.mat=matrix(0,nrow=dim(cohort.ncc)[2],ncol=dim(cohort.ncc)[2])
colnames(pred.mat)=names(cohort.ncc)
rownames(pred.mat)=names(cohort.ncc)
pred.mat["x",predictors_approx]=1

#method of imputation for x
method.vec=rep("",dim(cohort.ncc)[2])
method.vec[which(colnames(cohort.ncc)=="x")]="norm"

#perform the imputation 
imp<-mice(cohort.ncc, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = nit, diagnostics = FALSE, printFlag = F)

# Fit the analysis model in each imputed data set
models<-with(imp,coxph(Surv(t,d)~x+z1+z2))

# Combine estimates across the imputed data sets using Rubin's Rules
summary_aprx = summary(pool(models))

ncc.MIapprox.full.coef = summary_aprx[,"estimate"]
ncc.MIapprox.full.se = summary_aprx[,"std.error"]

#===============================
#NCC
#MI-SMC: full-cohort approach
#===============================

#predictor matrix which determines the imputation models for x
pred.mat=matrix(0,nrow=dim(cohort.ncc)[2],ncol=dim(cohort.ncc)[2])
colnames(pred.mat)=names(cohort.ncc)
rownames(pred.mat)=names(cohort.ncc)
pred.mat["x",predictors_smcfcs]=1

#method of imputation for x
method.vec=rep("",dim(cohort.ncc)[2])
method.vec[which(colnames(cohort.ncc)=="x")]="norm"

#perform the imputation
imp <- smcfcs.red(cohort.ncc, smtype="coxph", smformula="Surv(t,d)~x+z1+z2",method=method.vec,
              predictorMatrix=pred.mat,m = nimp, numit = nit.smc, rjlimit = 10000,noisy=F)

# obtain estimates across the imputed data sets and combine using Rubin's Rules
impobj <- imputationList(imp$impDatasets)
models <- with(impobj, coxph(Surv(t,d)~x+z1+z2))

ncc.MIsmc.full.coef = MIcombine(models)$coefficients
ncc.MIsmc.full.se = sqrt(diag(MIcombine(models)$variance))

#===============================
#NCC with superset (super1)
#MI-approx
#===============================

#predictor matrix which determines the imputation models for x
pred.mat=matrix(0,nrow=dim(ncc.super1)[2],ncol=dim(ncc.super1)[2])
colnames(pred.mat)=names(ncc.super1)
rownames(pred.mat)=names(ncc.super1)
pred.mat["x",predictors_approx]=1

#method of imputation for x
method.vec=rep("",dim(ncc.super1)[2])
method.vec[which(colnames(ncc.super1)=="x")]="norm"

#perform the imputation 
imp<-mice(ncc.super1, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = nit, diagnostics = FALSE, printFlag = F)

# Fit the analysis model in each imputed data set
models <- vector("list", nimp)
for (k in 1:nimp){
  model=coxph(Surv(t.case,case)~x+z1+z2+strata(setno),data=complete(imp,k))
  models[[k]] = model
}

# Combine estimates across the imputed data sets using Rubin's Rules
ncc.MIapprox.super1.coef = MIcombine(models)$coef
ncc.MIapprox.super1.se = sqrt(diag(MIcombine(models)$variance))

#===============================
#case-cohort with superset (super1)
#MI-SMC
#===============================

# Compute number at risk at each event time using the full cohort data
nrisk.fit<-survfit(Surv(t,d)~1,data=cohort)
ord.t.d1<-order(cohort$t[cohort$d==1])
numrisk<-summary(nrisk.fit,censored=F)$n.risk #this gives the number at risk at each unique event time

#add numbers at risk at each event time into the nested case-control data set
#note that the way the data are generated means that the ordering of cases is the same in cohort and in ncc.super1
ncc.super1$numrisk<-NA
ncc.super1$numrisk[ncc.super1$case==1][ord.t.d1]<-numrisk

#In each matched set: assign number at risk at the case's event time to every individual in the set
ncc.super1$numrisk<-ave(ncc.super1$numrisk, ncc.super1$setno, FUN = function(x) sum(x, na.rm=T))

#predictor matrix which determines the imputation models for x
pred.mat=matrix(0,nrow=dim(ncc.super1)[2],ncol=dim(ncc.super1)[2])
colnames(pred.mat)=names(ncc.super1)
rownames(pred.mat)=names(ncc.super1)
pred.mat["x",predictors_smcfcs]=1

#method of imputation for x
method.vec=rep("",dim(ncc.super1)[2])
method.vec[which(colnames(ncc.super1)=="x")]="norm"

#perform the imputation
imp<-smcfcs.nestedcc.red(originaldata=ncc.super1,smformula="Surv(t.case,case)~x+z1+z2+strata(setno)",
                          set="setno",event="case",nrisk="numrisk",
                          method=method.vec, predictorMatrix=pred.mat,
                          m=nimp,numit=nit.smc,rjlimit=10000,noisy=FALSE) 

# Fit the analysis model in each imputed data set
impobj <- imputationList(imp$impDatasets)
models <- vector("list", nimp)
for (k in 1:nimp){
  model <- coxph(Surv(t.case,case)~x+z1+z2+strata(setno),data=impobj$imputations[[k]])
  models[[k]] = model
}

ncc.MIsmc.super1.coef = MIcombine(models)$coefficients
ncc.MIsmc.super1.se = sqrt(diag(MIcombine(models)$variance))

#===============================
#NCC with superset (super2)
#MI-approx
#===============================

#predictor matrix which determines the imputation models for x
pred.mat=matrix(0,nrow=dim(ncc.super2)[2],ncol=dim(ncc.super2)[2])
colnames(pred.mat)=names(ncc.super2)
rownames(pred.mat)=names(ncc.super2)
pred.mat["x",predictors_approx]=1

#method of imputation for x
method.vec=rep("",dim(ncc.super2)[2])
method.vec[which(colnames(ncc.super2)=="x")]="norm"

#perform the imputation 
imp<-mice(ncc.super2, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = nit, diagnostics = FALSE, printFlag = F)

# Fit the analysis model in each imputed data set
models <- vector("list", nimp)
for (k in 1:nimp){
  model=coxph(Surv(t.case,case)~x+z1+z2+strata(setno),data=complete(imp,k))
  models[[k]] = model
}

# Combine estimates across the imputed data sets using Rubin's Rules
ncc.MIapprox.super2.coef = MIcombine(models)$coef
ncc.MIapprox.super2.se = sqrt(diag(MIcombine(models)$variance))

#===============================
#case-cohort with superset (super2)
#MI-SMC
#===============================

# Compute number at risk at each event time using the full cohort data
nrisk.fit<-survfit(Surv(t,d)~1,data=cohort)
ord.t.d1<-order(cohort$t[cohort$d==1])
numrisk<-summary(nrisk.fit,censored=F)$n.risk #this gives the number at risk at each unique event time

#add numbers at risk at each event time into the nested case-control data set
#note that the way the data are generated means that the ordering of cases is the same in cohort and in ncc.super2
ncc.super2$numrisk<-NA
ncc.super2$numrisk[ncc.super2$case==1][ord.t.d1]<-numrisk

#In each matched set: assign number at risk at the case's event time to every individual in the set
ncc.super2$numrisk<-ave(ncc.super2$numrisk, ncc.super2$setno, FUN = function(x) sum(x, na.rm=T))

#predictor matrix which determines the imputation models for x
pred.mat=matrix(0,nrow=dim(ncc.super2)[2],ncol=dim(ncc.super2)[2])
colnames(pred.mat)=names(ncc.super2)
rownames(pred.mat)=names(ncc.super2)
pred.mat["x",predictors_smcfcs]=1

#method of imputation for x
method.vec=rep("",dim(ncc.super2)[2])
method.vec[which(colnames(ncc.super2)=="x")]="norm"

#perform the imputation
imp<-smcfcs.nestedcc.red(originaldata=ncc.super2,smformula="Surv(t.case,case)~x+z1+z2+strata(setno)",
                           set="setno",event="case",nrisk="numrisk",
                           method=method.vec, predictorMatrix=pred.mat,
                           m=nimp,numit=nit.smc,rjlimit=10000,noisy=FALSE) 

# Fit the analysis model in each imputed data set
impobj <- imputationList(imp$impDatasets)
models <- vector("list", nimp)
for (k in 1:nimp){
  model <- coxph(Surv(t.case,case)~x+z1+z2+strata(setno),data=impobj$imputations[[k]])
  models[[k]] = model
}

ncc.MIsmc.super2.coef = MIcombine(models)$coefficients
ncc.MIsmc.super2.se = sqrt(diag(MIcombine(models)$variance))

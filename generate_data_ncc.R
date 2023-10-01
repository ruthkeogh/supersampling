#############################################################
#############################################################
#' Simulates data for a full cohort.
#' Creates a standard nested case-control (NCC) sample and two superset samples. 
#' Need to specify n=5000 or n=25000.
#' Need to specify interaction=F or T (F means no interaction between X and Z1 in the hazard model. T means there is an interaction between X and Z1 in the hazard model).
#' Need to specify misspec=F or T (misspec=F means that X|Z1,Z2 is normal. misspec=T means that X|Z1,Z2 is non-normal).
#############################################################
#############################################################

#===============================
#packages
#===============================

library(survival)
library(mice)
library(mitools)
library(splitstackshape)
library(tidyverse)
library(smcfcs)

#===============================
# Set parameter values, etc
#===============================

close.time=15  # Maximum follow-up time

p_z1 = 0.5   # Bernoulli probability for Z1
mu_z2 = 0    # mean for normal Z2

a_x = 0      # constant for X
b_x=0.25   # influence of Z1 on X
c_x=0.5   # influence of Z2 on X

beta.x=1   # log hazard ratio for X
beta.z1=0.5   # log hazard ratio for continuous Z1
beta.z2=1     # log hazard ratio for binary Z2
kappa = 4  # Weibull baseline shape 
kappa_drop=4     # Weibull baseline shape for dropout time

sigma.sq<-0.694^2 #used in generating X from a non-normal distribution - for mis-specified imputation model scenario

#scenarios with misspec=F

if(misspec==F){
  if(n==5000 & interaction==F){
    lambda=3.73*10^(-7) # Weibull baseline scale for event of interest=
    lambda_drop = 1.92*10^(-5)    # Weibull baseline scale dropout time
    beta.x.z1=0     # log hazard ratio for interaction X*Z1
  }
  
  if(n==25000 & interaction==F){
    lambda=0.672*10^(-7) # Weibull baseline scale for event of interest=
    lambda_drop = 2.04*10^(-5)    # Weibull baseline scale dropout time
    beta.x.z1=0     # log hazard ratio for interaction X*Z1
  }
  
  if(n==5000 & interaction==T){
    lambda=2.35*10^(-7) # Weibull baseline scale for event of interest=
    lambda_drop = 1.93*10^(-5)    # Weibull baseline scale dropout time
    beta.x.z1=0.5     # log hazard ratio for interaction X*Z1
  }
  
  if(n==25000 & interaction==T){
    lambda=0.292*10^(-7) # Weibull baseline scale for event of interest=
    lambda_drop = 2.04*10^(-5)    # Weibull baseline scale dropout time
    beta.x.z1=0.5     # log hazard ratio for interaction X*Z1
  }
}

#scenarios with misspec=T
  
if(misspec==T){
  if(n==5000 & interaction==F){
    lambda=3.5*10^(-7) # Weibull baseline scale for event of interest
    lambda_drop = 1.94*10^(-5)    # Weibull baseline scale dropout time
    beta.x.z1=0
  }
  
  if(n==25000 & interaction==F){
    lambda=0.44*10^(-7) # Weibull baseline scale for event of interest
    lambda_drop = 2.1*10^(-5)    # Weibull baseline scale dropout time
    beta.x.z1=0
  }
}

# Generate id numbers
id=seq(1,n)

#===============================
# Generate covariates 
#===============================

z1=rnorm(n,mu_z2,1)
z2=rbinom(n,1,p_z1)

if(misspec==F){
  x=rnorm(n,a_x +b_x*z1+c_x*z2,1)   
}

if(misspec==T){
  y=rnorm(n,0,sqrt(sigma.sq))
  eps<-exp(y)-exp(sigma.sq/2)
  x=a_x+b_x*z1+c_x*z2+eps
}

#===============================
# Generate potential event times
#===============================

u=runif(n,0,1)
t.event=(-log(u)*(1/lambda)*
           exp(-(beta.x*x+beta.z1*z1+beta.z2*z2+beta.x.z1*x*z1)))^(1/kappa)

#===============================
# Generate potential drop-out time
#===============================

u=runif(n,0,1)
t.drop=(-log(u)*(1/lambda_drop))^(1/kappa_drop)

#===============================
# Generate time for event or drop out
#===============================

t=pmin(t.event,t.drop,close.time)
cause=1*(t==t.event)+2*(t==t.drop)+3*(t==close.time) 
# 1: event, 2: drop out, 3: administrative censoring
d=ifelse(cause==1,1,0) 

prop.table(table(d))

#===============================
# the full-cohort data with no missingness
#===============================

cohort=data.frame(id,t,d,x,z1,z2,cause)

if(interaction==T){
  cohort$xz1int<-cohort$x*cohort$z1 #used when applying smcfcs
}

#===============================
# GENERATE NCC DATA WITH 12 CONTROLS PER CASE 
#LATER WE DROP SOME CONTROLS. 12 IS THE MAXIMUM NUMBER OF CONTROLS USED
#===============================

n.controls=12

ncc.full=NULL
n.case=sum(cohort$d==1)
case=NULL
setid=NULL
setno=NULL    
setno.count=0
for (i in which(cohort$d==1))
{
  # Select control(s) for nested case-control
  possible.controls=which(cohort$t>=cohort$t[i])
  if (length(possible.controls)>=n.controls){
    controls=sample(possible.controls,n.controls)
    ncc.full=rbind(ncc.full,cohort[i,])
    ncc.full=rbind(ncc.full,cohort[controls,])
    
    case=c(case,c(1,rep(0,n.controls))) #case-control indicator within the set (different from d)
    setid=c(setid,rep(1:(n.controls+1))) #id number of individual within set. Case always has setid=1 and controls have 2,3,....
    setno.count=setno.count+1 
    setno=c(setno,rep(setno.count,n.controls+1)) #c-c set identifier
  } else if(length(possible.controls)<n.controls){
    #we do this in case there are not enough controls available
    controls=possible.controls
    ncc.full=rbind(ncc.full,cohort[i,])
    ncc.full=rbind(ncc.full,cohort[controls,])
    
    case=c(case,c(1,rep(0,length(controls))))
    setid=c(setid,rep(1:(length(controls)+1)))
    setno.count=setno.count+1
    setno=c(setno,rep(setno.count,length(controls)+1))
  }
}

ncc.full$case=case
ncc.full$setid=setid
ncc.full$setno=setno

#===============================
# Generate data-set which is just the ncc substudy with 1 control per case
#===============================

ncc=ncc.full[ncc.full$setid<=2,]

#===============================
# Generate data-set which is the full cohort but with x missing for those not in the ncc study
#===============================

cohort.ncc=cohort
cohort.ncc$in.ncc <- cohort.ncc$id%in%ncc$id
cohort.ncc$x<-ifelse(cohort.ncc$in.ncc==1,cohort.ncc$x,NA)

if(interaction==T){
  cohort.ncc$xz1int<-ifelse(cohort.ncc$in.ncc==1,cohort.ncc$xz1int,NA)
}

#===============================
# Generate ncc superset data
#===============================

# ------------------------------
# Superset 1: extra 3 controls per case (i.e. 4 controls in total)

ncc.super1=ncc.full[ncc.full$setid<=5,]

#set x to missing for controls 2:4
ncc.super1$x<-ifelse(ncc.super1$setid>2 & ncc.super1$d==0,NA,ncc.super1$x)

if(interaction==T){
  ncc.super1$xz1int<-ifelse(ncc.super1$setid>2 & ncc.super1$d==0,NA,ncc.super1$xz1int)
}

# ------------------------------
# Superset 2: extra 11 controls per case (i.e. 12 controls in total)

ncc.super2=ncc.full

#set x to missing for controls 2:12
ncc.super2$x<-ifelse(ncc.super2$setid>2 & ncc.super2$d==0,NA,ncc.super2$x)

if(interaction==T){
  ncc.super2$xz1int<-ifelse(ncc.super2$setid>2 & ncc.super2$d==0,NA,ncc.super2$xz1int)
}

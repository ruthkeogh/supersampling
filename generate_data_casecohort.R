#############################################################
#############################################################
#' Simulates data for a full cohort.
#' Creates a standard case-cohort sample and two superset samples. 
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
    lambda=3.73*10^(-7) # Weibull baseline scale for event of interest
    lambda_drop = 1.92*10^(-5)    # Weibull baseline scale dropout time
    beta.x.z1=0     # log hazard ratio for interaction X*Z1
  }
  
  if(n==25000 & interaction==F){
    lambda=0.672*10^(-7) # Weibull baseline scale for event of interest
    lambda_drop = 2.04*10^(-5)    # Weibull baseline scale dropout time
    beta.x.z1=0     # log hazard ratio for interaction X*Z1
  }
  
  if(n==5000 & interaction==T){
    lambda=2.35*10^(-7) # Weibull baseline scale for event of interest
    lambda_drop = 1.93*10^(-5)    # Weibull baseline scale dropout time
    beta.x.z1=0.5     # log hazard ratio for interaction X*Z1
  }
  
  if(n==25000 & interaction==T){
    lambda=0.292*10^(-7) # Weibull baseline scale for event of interest
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
# GENERATE CASE-COHORT DATA
#===============================

cohort.caco=cohort

# ------------------------------
# Generate subcohort

if(n==5000){
  n.subco=0.05*n
}

if(n==25000){
  n.subco=0.01*n
}

cohort.caco$subco<-c(rep(1,n.subco),rep(0,n-n.subco))

# ------------------------------
#make x missing in those outside the case-cohort sample

cohort.caco$x<-ifelse(cohort.caco$subco==1|cohort.caco$d==1,cohort.caco$x,NA)

if(interaction==T){
  cohort.caco$xz1int<-ifelse(cohort.caco$subco==1|cohort.caco$d==1,cohort.caco$xz1int,NA)
}

# ------------------------------
# Generate data-set which is just the case-cohort substudy

caco=cohort.caco[cohort.caco$subco==1|cohort.caco$d==1,]

caco$entertime=ifelse(caco$d==1 & caco$subco==0,caco$t-0.001,0)

#===============================
# Generate case-cohort superset data
#===============================

# ------------------------------
# Superset 1

if(n==5000){
  n.subco.super1 = n.subco+0.15*n
}

if(n==25000){
  n.subco.super1 = n.subco+0.03*n
}

cohort.caco$subco.super1 <- c(rep(1,n.subco.super1),rep(0,n-n.subco.super1))

cohort.super1.caco=cohort.caco[cohort.caco$subco.super1==1|cohort.caco$d==1,]
cohort.super1.caco$subco=NULL
cohort.super1.caco$subco.super2=NULL

cohort.super1.caco$entertime=ifelse(cohort.super1.caco$d==1 &
                                      cohort.super1.caco$subco.super1==0,
                                    cohort.super1.caco$t-0.001,0)

# ------------------------------
# Superset 2

if(n==5000){
  n.subco.super2 = n.subco+0.55*n
}

if(n==10000){
  n.subco.super2 = n.subco+0.275*n
}

if(n==25000){
  n.subco.super2 = n.subco+0.11*n
}

cohort.caco$subco.super2 <- c(rep(1,n.subco.super2),rep(0,n-n.subco.super2))

cohort.super2.caco=cohort.caco[cohort.caco$subco.super2==1|cohort.caco$d==1,]
cohort.super2.caco$subco=NULL
cohort.super2.caco$subco.super1=NULL

cohort.super2.caco$entertime=ifelse(cohort.super2.caco$d==1 &
                                         cohort.super2.caco$subco.super2==0,
                                       cohort.super2.caco$t-0.001,0)

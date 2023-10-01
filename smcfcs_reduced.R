#' smcfcs imputes missing values of covariates using the Substantive Model Compatible
#' Fully Conditional Specification multiple imputation approach proposed by
#' Bartlett \emph{et al} 2015 (see references).
#'
#' The main function smcfcs in the package smcfcs enables imputation under different analysis models, including 
#' linear regression, logistic regression, Poisson regression, and Cox regression.
#' In this reduced version of the smcfcs function, we reduce it to the elements required for Cox regression, case-cohort studies and nested case-control studies. 
#' 
#' The main function smcfcs in the package smcfcs enables imputation under different imputation model types, including 
#' normal linear regression, logistic, proportional odds, multinomial logistic. 
#' In this reduced version of the smcfcs function, we reduce it to the elements required for normal linear regression. method="norm"
#' 
#' @export
smcfcs.red <- function(originaldata,smtype,smformula,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,errorProneMatrix=NULL) {
  #call core  smcfcs function, passing through arguments
  smcfcs.core(originaldata,smtype,smformula,method,predictorMatrix,m,numit,rjlimit,noisy,errorProneMatrix=errorProneMatrix)
}
#'
#' @export
smcfcs.casecohort.red.Prentice <- function(originaldata,smformula,sampfrac,in.subco,method,cohort.size,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,
                                           errorProneMatrix=NULL) {
  smcfcs.core(originaldata,smtype="casecohort.Prentice",smformula,method,predictorMatrix,m,numit,rjlimit,noisy,sampfrac=sampfrac,cohort.size=cohort.size,in.subco=in.subco,
              errorProneMatrix=errorProneMatrix)
}
smcfcs.casecohort.red.LinYing <- function(originaldata,smformula,sampfrac,in.subco,method,cohort.size,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,
                                          errorProneMatrix=NULL) {
  smcfcs.core(originaldata,smtype="casecohort.LinYing",smformula,method,predictorMatrix,m,numit,rjlimit,noisy,sampfrac=sampfrac,cohort.size=cohort.size,in.subco=in.subco,
              errorProneMatrix=errorProneMatrix)
}
smcfcs.nestedcc.red <- function(originaldata,smformula,set,event,nrisk,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,errorProneMatrix=NULL) {
  smcfcs.core(originaldata,smtype="nestedcc",smformula,method,predictorMatrix,m,numit,rjlimit,noisy,set=set,event=event,nrisk=nrisk,
              errorProneMatrix=errorProneMatrix)
}


#this is the core of the smcfcs function, called by wrapper functions for certain different substantive models
smcfcs.core <- function(originaldata,smtype,smformula,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,errorProneMatrix=NULL,
                        ...) {
  
  #get extra arguments passed in ...
  extraArgs <- list(...)
  
  stopifnot(is.data.frame(originaldata))
  if (ncol(originaldata)!=length(method)) stop("Method argument must have the same length as the number of columns in the data frame.")
  
  n <- dim(originaldata)[1]
  
  #create matrix of response indicators
  r <- 1*(is.na(originaldata)==0)
  
  #find column numbers of partially observed, fully observed variables, and outcome
  if (smtype=="coxph") {
    
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[,dCol]
    
    nullMod <- survival::coxph(Surv(originaldata[,timeCol],originaldata[,dCol])~1,
                               control = survival::coxph.control(timefix = FALSE))
    basehaz <- survival::basehaz(nullMod)
    H0indices <- match(originaldata[,timeCol], basehaz[,2])
    rm(nullMod)
  } else if (smtype=="casecohort.Prentice"|smtype=="casecohort.LinYing") {
    
    subcoCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% extraArgs$in.subco]
    #subcoMembers is a vector of row numbers of those in the subcohort
    subcoMembers <- which(originaldata[,subcoCol]==1)
    
    #generate weights for use in later analysis which we use to obtain baseline cumulative hazard
    #assign a weight of /samp.frac to individuals in the subcohort and 0 to those outside the subcohort
    subco.weight<-ifelse(originaldata[,subcoCol]==1,1/extraArgs$sampfrac,0)
    
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[,dCol]
    
    #list of unique event times - used in calculation of baseline cumulative hazard
    list.times=sort(unique(originaldata[,timeCol][originaldata[,dCol]==1]))  #RUTH 21/03/17: ADDED THIS LINE
    
    smcfcsid <- 1:n
  } else if (smtype=="nestedcc") {
    
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    setCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$set)]
    nriskCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$nrisk)]
    eventCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$event)] #this may or may not be distinct from the dCol!
    
    #the below command creates "Surv(t,case)~x" from "Surv(t,case)~x+strata(setno)" (for example) (i.e. it removes the strata part of the formula)
    #this is used when obtaining outmodxb
    #note this is done slightly oddly, but this is because as.formula does not work well for long formulas as it splits across lines
    exp1=as.formula(paste(smformula))[[2]]
    exp2=as.formula(smformula)[[3]][[2]]
    smformula2<-paste(deparse(exp1),"~",deparse(exp2,width.cutoff = 500L))
    
    #This is the indicator of whether an individual ever has the event (regardless of whether they are sometimes used as a control and sometimes (one) as a case)
    d <- originaldata[,eventCol]
    
    #number of individuals in each sampled risk set (matched set)
    num.sampriskset<-ave(rep(1,dim(originaldata)[1]), originaldata[,setCol], FUN = function(x) sum(x))
    
    controls<- which(originaldata[,eventCol]==0)
  }
  
  if (smtype=="casecohort.Prentice") {
    #generate weights for use in later analysis which we use to obtain baseline cumulative hazard
    #assign a weight of /samp.frac to individuals in the subcohort and 0 to those outside the subcohort
    subco.weight<-ifelse(originaldata[,subcoCol]==1,1/extraArgs$sampfrac,0)
  }  else if (smtype=="casecohort.LinYing") {
    subco.weight<-ifelse(originaldata[,subcoCol]==1,1/extraArgs$sampfrac,1)
  } 
  
  if (smtype=="coxph"|smtype=="casecohort.Prentice"|smtype=="casecohort.LinYing") {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")
  } else if (smtype=="nestedcc") {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")[-length(attr(terms(as.formula(smformula)), "term.labels"))]
  }
  smcovcols <- (1:ncol(originaldata))[colnames(originaldata) %in% smcovnames]
  
  #partial vars are those variables for which an imputation method has been specified among the available regression types
  partialVars <- which((method=="norm"))
  
  if (length(partialVars)==0) stop("You have not specified any valid imputation methods in the method argument.")
  
  #fully observed vars are those that are fully observed and are covariates in the substantive model
  fullObsVars <- which((colSums(r)==n) & (colnames(originaldata) %in% smcovnames))
  
  #passive variables
  passiveVars <- which((method!="") & (method!="norm"))
  
  print(paste("Outcome variable(s):", paste(colnames(originaldata)[outcomeCol],collapse=',')))
  print(paste("Passive variables:", paste(colnames(originaldata)[passiveVars],collapse=',')))
  print(paste("Partially obs. variables:", paste(colnames(originaldata)[partialVars],collapse=',')))
  print(paste("Fully obs. substantive model variables:", paste(colnames(originaldata)[fullObsVars],collapse=',')))
  
  imputations <- list()
  for (imp in 1:m) {
    imputations[[imp]] <- originaldata
  }
  
  rjFailCount <- 0
  
  for (imp in 1:m) {
    
    print(paste("Imputation ",imp))
    
    #initial imputation of each partially observed variable based on observed values
    for (var in 1:length(partialVars)) {
      targetCol <- partialVars[var]
      imputations[[imp]][r[,targetCol]==0,targetCol] <- sample(imputations[[imp]][r[,targetCol]==1,targetCol], size=sum(r[,targetCol]==0), replace=TRUE)
    }
    
    for (cyclenum in 1:numit) {
      
      if (noisy==TRUE) {
        print(paste("Iteration ",cyclenum))
      }
      #update passive variable(s)
      imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
      
      for (var in 1:length(partialVars)) {
        targetCol <- partialVars[var]
        if (is.null(predictorMatrix)) {
          predictorCols <- c(partialVars[! partialVars %in% targetCol], fullObsVars)
        } else {
          predictorCols <- which(predictorMatrix[targetCol,]==1)
          #ensure that user has not included outcome variable(s) here
          predictorCols <- predictorCols[! predictorCols %in% outcomeCol]
        }
        if ((imp==1) & (cyclenum==1)) {
          print(paste("Imputing: ",colnames(imputations[[imp]])[targetCol]," using ",paste(colnames(imputations[[imp]])[predictorCols],collapse=',')," plus outcome",collapse=','))
        }
        if (length(predictorCols)>0) {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~", paste(colnames(imputations[[imp]])[predictorCols], collapse="+"),sep=""))
        } else {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~1",sep=""))
        }
        if (smtype=="coxph") {
          xmoddata <- imputations[[imp]]
        } else if(smtype=="casecohort.Prentice"|smtype=="casecohort.LinYing") {
          xmoddata <- imputations[[imp]][subcoMembers,]
        } else if (smtype=="nestedcc"){
          xmoddata <- imputations[[imp]][controls,]
        } 
        if (method[targetCol]=="norm") {
          #estimate parameters of covariate model
          xmod <- lm(xmodformula, data=xmoddata)
          #take draw from posterior of covariate model parameters
          beta <- xmod$coef
          sigmasq <- summary(xmod)$sigma^2
          newsigmasq <- (sigmasq*xmod$df) / rchisq(1,xmod$df)
          covariance <- (newsigmasq/sigmasq)*vcov(xmod)
          newbeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
          #calculate fitted values
          if (smtype=="casecohort.Prentice"|smtype=="casecohort.LinYing"|smtype=="nestedcc") {
            xfitted <- model.matrix(xmodformula, data=imputations[[imp]]) %*% newbeta
          } else {
            xfitted <- model.matrix(xmod) %*% newbeta
          }
        }
        if (noisy==TRUE) {
          print(summary(xmod))
        }
        
        #estimate parameters of substantive model
        if (smtype=="coxph") {
          ymod <- survival::coxph(as.formula(smformula), imputations[[imp]],
                                  control = survival::coxph.control(timefix = FALSE))
          outcomeModBeta <- modPostDraw(ymod)
          ymod$coefficients <- outcomeModBeta
          basehaz <- survival::basehaz(ymod, centered=FALSE)[,1]
          H0 <- basehaz[H0indices]
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        } else if (smtype=="casecohort.Prentice") {
          ymod <- survival::cch(as.formula(smformula),id=~id,subcoh=imputations[[imp]][,subcoCol],cohort.size=extraArgs$cohort.size, data=imputations[[imp]],method="Prentice")
          outcomeModBeta <- modPostDraw(ymod)
          
          cumhaz.denom.elements=exp(model.matrix(as.formula(smformula),imputations[[imp]])[,-1] %*% outcomeModBeta)
          cumhaz.denom=sapply(list.times,function(x){sum(cumhaz.denom.elements[which(originaldata[,timeCol]>=x)]*subco.weight[which(originaldata[,timeCol]>=x)])})
          exp.func.denom=cumsum(1/cumhaz.denom)
          H0.fun=stepfun(list.times,c(0,exp.func.denom))
          H0=H0.fun(originaldata[,timeCol])
          
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        } else if (smtype=="casecohort.LinYing") {
          ymod <- survival::cch(as.formula(smformula),id=~id,subcoh=imputations[[imp]][,subcoCol],cohort.size=extraArgs$cohort.size, data=imputations[[imp]],method="LinYing")
          outcomeModBeta <- modPostDraw(ymod)
          
          cumhaz.denom.elements=exp(model.matrix(as.formula(smformula),imputations[[imp]])[,-1] %*% outcomeModBeta)
          cumhaz.denom=sapply(list.times,function(x){sum(cumhaz.denom.elements[which(originaldata[,timeCol]>=x)]*subco.weight[which(originaldata[,timeCol]>=x)])})
          exp.func.denom=cumsum(1/cumhaz.denom)
          H0.fun=stepfun(list.times,c(0,exp.func.denom))
          H0=H0.fun(originaldata[,timeCol])
          
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        } else if (smtype=="nestedcc") {
          ymod <- survival::coxph(as.formula(smformula), imputations[[imp]],
                                  control = survival::coxph.control(timefix = FALSE))
          outcomeModBeta <- modPostDraw(ymod)
          
          explan.matrix<-model.matrix(ymod)
          cumbasehaz.denom<-exp(matrix(outcomeModBeta,nrow=1)%*%t(explan.matrix))*originaldata[,nriskCol]/num.sampriskset
          cumbasehaz.denom<-ave(cumbasehaz.denom, originaldata[,setCol], FUN = sum)[originaldata[,dCol]==1] #this is the denominator of the contribution to the cumulative baseline hazard at each event time
          cumbasehaz.t<-originaldata[,timeCol][originaldata[,dCol]==1] #times to which the baseline cumulative hazards refer
          H0<-unlist(lapply(originaldata[,timeCol],function(x) {sum((1/cumbasehaz.denom)[cumbasehaz.t<=x])}))
          
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        
        if ((imp==1) & (cyclenum==1) & (var==1)) {
          smCoefIter <- array(0, dim=c(m, length(outcomeModBeta), numit))
        }
        
        if (var==length(partialVars)) {
          #then we have reached end of a cycle
          smCoefIter[imp,,cyclenum] <- outcomeModBeta
        }
        
        #impute x, either directly where possibly, or using rejection sampling otherwise
        imputationNeeded <- (1:n)[r[,targetCol]==0]
        
        
        #use rejection sampling
        #first draw for all subjects who need imputing, using a small number of attempts
        firstTryLimit <- 25
        j <- 1
        
        while ((length(imputationNeeded)>0) & (j<firstTryLimit)) {
          #sample from covariate model
          if ((method[targetCol]=="norm")) {
            imputations[[imp]][imputationNeeded,targetCol] <- rnorm(length(imputationNeeded),xfitted[imputationNeeded],newsigmasq^0.5)
          }
          
          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
          
          #accept/reject
          uDraw <- runif(length(imputationNeeded))
          if ((smtype=="coxph") | (smtype=="casecohort.Prentice")| (smtype=="casecohort.LinYing")) {
            outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]])
            outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
            s_t = exp(-H0[imputationNeeded]* exp(outmodxb[imputationNeeded]))
            prob = exp(1 + outmodxb[imputationNeeded] - (H0[imputationNeeded]* exp(outmodxb[imputationNeeded])) ) * H0[imputationNeeded]
            prob = d[imputationNeeded]*prob + (1-d[imputationNeeded])*s_t
            reject = 1*(uDraw > prob )
          } else if (smtype=="nestedcc") {
            outmodxb <-  model.matrix(as.formula(smformula2),imputations[[imp]])
            outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
            s_t = exp(-H0[imputationNeeded]* exp(outmodxb[imputationNeeded]))
            prob = exp(1 + outmodxb[imputationNeeded] - (H0[imputationNeeded]* exp(outmodxb[imputationNeeded])) ) * H0[imputationNeeded]
            prob = d[imputationNeeded]*prob + (1-d[imputationNeeded])*s_t
            reject = 1*(uDraw > prob )
          }
          imputationNeeded <- imputationNeeded[reject==1]
          
          j <- j+1
        }
        
        #now, for those remaining, who must have low acceptance probabilities, sample by subject
        for (i in imputationNeeded) {
          
          tempData <- imputations[[imp]][i,]
          tempData <- tempData[rep(1,rjlimit),]
          if (method[targetCol]=="norm") {
            tempData[,targetCol] <- rnorm(rjlimit,xfitted[i],newsigmasq^0.5)
          }
          #passively impute
          tempData <- updatePassiveVars(tempData, method, passiveVars)
          
          #accept reject
          uDraw <- runif(rjlimit)
          
          if ((smtype=="coxph") | (smtype=="casecohort.Prentice")| (smtype=="casecohort.LinYing")) {
            outmodxb <-  model.matrix(as.formula(smformula),tempData)
            outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
            s_t = exp(-H0[i]* exp(outmodxb))
            prob = exp(1 + outmodxb - (H0[i]* exp(outmodxb)) ) * H0[i]
            prob = d[i]*prob + (1-d[i])*s_t
            reject = 1*(uDraw > prob )
          }  else if (smtype=="nestedcc") {
            outmodxb <-  model.matrix(as.formula(smformula2),tempData)
            outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
            s_t = exp(-H0[i]* exp(outmodxb))
            prob = exp(1 + outmodxb - (H0[i]* exp(outmodxb)) ) * H0[i]
            prob = d[i]*prob + (1-d[i])*s_t
            reject = 1*(uDraw > prob )
          }
          if (sum(reject)<rjlimit) {
            imputations[[imp]][i,targetCol] <- tempData[reject==0,targetCol][1]
          } else {
            rjFailCount <- rjFailCount + 1
          }
        }
        #update passive variables
        imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
      }
      
    }
    
  }
  
  if (rjFailCount>0) {
    warning(paste("Rejection sampling failed ",rjFailCount," times (across all variables, iterations, and imputations). You may want to increase the rejection sampling limit.",sep=""))
  }
  
  # Added smformula and smtype to metadata, and make "smcfcs class"
  res <- list(
    impDatasets = imputations,
    smCoefIter = smCoefIter,
    smInfo = list("smtype" = smtype, "smformula" = smformula)
  )
  class(res) <- "smcfcs"
  
  return(res)
}

updatePassiveVars <- function(data, method, passivecols) {
  for (i in passivecols) {
    data[,i] <- with(data, eval(parse(text=method[i])))
  }
  data
}

modPostDraw <- function(modobj) {
  beta <- modobj$coef
  varcov <- vcov(modobj)
  beta + MASS::mvrnorm(1, mu=rep(0,ncol(varcov)), Sigma=varcov)
}


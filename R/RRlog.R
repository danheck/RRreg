# Randomized response as dependent variable in logistic regression

# generic function to allow for a formula interface:
#' Logistic randomized response regression
#' 
#' A dichotomous variable, measured by a randomized response method, serves as dependent variable using one or more continuous and/or categorical predictors
#' @param formula specifying the regression model, see \code{\link{formula}}
#' @param data \code{data.frame}, in which variables can be found (optional)
#' @param model Available RR models: \code{"Warner","UQTknown","UQTunknown","Mangat","Kuk","FR","Crosswise","CDM","CDMsym","SLD"}. See \code{vignette("RRreg")} for details.
#' @param p randomization probability/probabilities (depending on model)
#' @param group vector specifying group membership. Can be omitted for single-group RR designs (e.g., Warner). For two-group RR designs (e.g., \code{CDM} or \code{SLD}), use 1 and 2 to indicate the group membership, matching the respective randomization probabilities \code{p[1], p[2]}. If an RR design and a direct question (DQ) were both used in the study, the group indices are set to 0 (DQ) and 1 (RR; 1 or 2 for two-group RR designs). This can be used to test, whether the RR design leads to a different prevalence estimate by including a dummy variable for the question format (RR vs. DQ) as predictor. If the corresponding regression coefficient is significant, the prevalence estimates differ between RR and DQ. Similarly, interaction hypotheses can be tested (e.g., the correlation between a sensitive attribute and a predictor is only found using the RR but not the DQ design). Hypotheses like this can be tested by including the interaction of the DQ-RR-dummy variable and the predictor in \code{formula} (e.g., \code{RR ~ dummy*predictor}).
#' @param LR.test test regression coefficients by a likelihood ratio test, i.e., fitting the model repeatedly while excluding one parameter at a time
# @param intercept should  the model contain an intercept?
#' @param fit.n Minimum and maximum number of fitting replications using random starting values to avoid local maxima (only if \code{start=NULL})
#' @param fit.bound The model is fitted repeatedly either until the absolute parameter estimates are below \code{fit.bound} or the maximum number of fitting replication is reached. Thereby, stability of the estimates is increased. \code{fit.bound} should be increased if extreme parameter estimates are to be expected.
#' @param start starting values for optimization. Might be useful if model does not converge with default starting values.
#' @param maxit Maximum number of iterations within each run of \code{optim}
#' @param ... ignored
#' @author Daniel W. Heck
#' @seealso \code{vignette('RRreg')} or \url{https://dl.dropboxusercontent.com/u/21456540/RRreg/index.html} for a detailed description of the RR models and the appropriate definition of \code{p} 
#' @return Returns an object \code{RRlog} which can be analysed by the generic method \code{\link{summary}}
#' @references van den Hout, A., van der Heijden, P. G., & Gilchrist, R. (2007). The logistic regression model with response variables subject to randomized response. \emph{Computational Statistics & Data Analysis, 51}, 6060-6069. 
#' @examples
#' # generate data set without biases
#' dat <- RRgen(1000,pi=.3,"Warner",p=.8)
#' dat$covariate <- rnorm(1000)
#' dat$covariate[dat$true==1] <- rnorm(sum(dat$true==1),.4,1)
#' # analyse
#' ana <- RRlog(response~covariate,dat,"Warner",.8, fit.n = c(1,5))
#' summary(ana)
# @rdname RRlog
#' @export
RRlog <- function(formula, data,model, p,group, LR.test=TRUE, 
                  fit.n=c(10,100),fit.bound=10, maxit=1000, 
                  start=NULL, ...) UseMethod("RRlog")

# choose model and construct S3 method 'RRlog' 

#' @export
RRlog.default <-function(formula,data,model,p,group, LR.test=TRUE, fit.n=c(10,100), fit.bound=10, maxit=1000, start=NULL, ...){
  # not very nice: avoiding CMD CHECK errors by using the same parameter names as in the generic function
  x <- formula;
  y <- data;
  # construct design matrix on predictor side (Intercept first)
  #   if (!intercept){
  #     temp <- colnames(x)[colnames(x) != "(Intercept)"]
  #     x <- x[, colnames(x) != "(Intercept)"]
  #     x <- as.matrix(x)
  #     colnames(x) <- temp
  #   }
  x <- as.matrix(x)
  y <- as.numeric(y)
  model <- match.arg(model,c("Warner","UQTknown","UQTunknown","Mangat","Kuk","FR","Crosswise","CDM","CDMsym","SLD"))
  
  # no DQ format included
  if (missing(group) || is.null(group)){
    group=rep(1,length(y))
  }
  
  RRcheck.log(model,y,p,group,names(y)[1])   ## for FR and Kuk: n model check for dichot. response
  
  # get estimates repeatedly (because of local minima and starting points)
  est <- list(logLik = -Inf, coef=Inf)
  # given starting values
  if (! ( missing(start) || is.null(start) || is.na(start))){
    start <- RRcheck.start(model,x,start) 
    cnt <- fit.n[2]-1
  }else{  
    cnt <- 0
  }
  
  while (cnt <fit.n[1] | 
           ( (est$logLik==-Inf |max(abs(est$coef))>fit.bound) & cnt <fit.n[2])){
    cnt <- cnt +1
    if (cnt == 1){
      glmcoef <- coef(glm.fit(x,y, family=binomial(link = "logit")))
      start <- glmcoef
      if(is2group(model)){
        uni <- RRuni(response=y,model=model,p=p,group=group)
        start <- c(start, summary(uni)$coef[2,1])
      }
    }else{
      start <- runif (ncol(x),-1e-4, 1e-4)*(3+25*cnt/fit.n[2])^3
      if(is2group(model)){
        start <- c(start, runif(1, .35,.95))
      }
    }
    
      switch(model,
             "Warner" = est2 <- RRlog.Warner(x,y,p,start,group, maxit=maxit) ,
             "UQTknown" = est2 <- RRlog.UQTknown(x,y,p,start,group, maxit=maxit),
             "UQTunknown" = est2 <- RRlog.UQTunknown(x,y,p,start,group, maxit=maxit),
             "Mangat" = est2 <- RRlog.Mangat(x,y,p,start,group, maxit=maxit),
             "Kuk" = est2 <- RRlog.Kuk(x,y,p,start,max(y),group, maxit=maxit),
             "FR" = est2 <- RRlog.FR(x,y,p,start,group, maxit=maxit),
             "Crosswise" = {
               est2 <- RRlog.Warner(x,y,p,start,group, maxit=maxit)
               est2$model="Crosswise"},
             "CDM" = est2 <- RRlog.CDM(x,y,p,start,group, maxit=maxit),
             "CDMsym" = est2 <- RRlog.CDMsym(x,y,p,start,group, maxit=maxit),
             "SLD" = est2 <- RRlog.SLD(x,y,p,start,group, maxit=maxit)
             )
      
    if (!is.na(est2$logLik) && est2$logLik > est$logLik) 
        est <- est2
  }
#   if (cnt == fit.n[2]) warning(paste0("Maximum number of fitting replications reached (fit.n=",fit.n[2],"). This could indicate extreme and/or unstable parameter estimates. Consider re-fitting the model (e.g., using fit.n=c(5,1000) and/or fit.bound=25)"))
  
  est$n <- length(y)
  est$n.dq <- sum(group==0)
  est$npar <- length(est$param)
  names(est$coefficients) <- est$param
  try({
    est$vcov <- solve(-est$hessian)
   names(est$gradient) <- est$param
    colnames(est$hessian) <- est$param
   rownames(est$hessian) <- est$param
    colnames(est$vcov) <- est$param
    rownames(est$vcov) <- est$param
  }, silent=T)
  est$start <- start
  names(est$start) <- est$param
  
  # LR test for each parameter
  if (LR.test){
    ncoef <- ncol(x)
    deltaLogLik <- rep(NA,est$npar)
    start <- est$coefficients
    for (i in 1:ncoef){
      xx <- x
      xx[,i] <- rep(0,length(y))
      switch(model,
             "Warner" = est2 <- RRlog.Warner(xx,y,p,start,group) ,
             "UQTknown" = est2 <- RRlog.UQTknown(xx,y,p,start,group),
             "UQTunknown" = est2 <- RRlog.UQTunknown(xx,y,p,start,group),
             "Mangat" = est2 <- RRlog.Mangat(xx,y,p,start,group),
             "Kuk" = est2 <- RRlog.Kuk(xx,y,p,start,max(y),group),
             "FR" = est2 <- RRlog.FR(xx,y,p,start,group),
             "Crosswise" = est2 <- RRlog.Warner(xx,y,p,start,group),
             "CDM" = est2 <- RRlog.CDM(xx,y,p,start,group),
             "CDMsym" = est2 <- RRlog.CDMsym(xx,y,p,start,group),
             "SLD" = est2 <- RRlog.SLD(xx,y,p,start,group)
      )
      deltaLogLik[i] <- est2$logLik - est$logLik
    }
    # multi group models: additional parameter
    if (model %in% c("SLD","CDM","CDMsym","UQTunknown")){
      switch(model,
             "SLD" = est2 <- RRlog.SLD(x,y,p,start,group,setT=T),
             "CDM" = est2 <- RRlog.CDM(x,y,p,start,group,setGamma=T),
             "CDMsym" = est2 <- RRlog.CDMsym(x,y,p,start,group,setGamma=T),
             "UQTunknown" = est2 <- RRlog.UQTunknown(x,y,p,start,group,setPiUQ=T) )
      deltaLogLik[est$npar] <- est2$logLik - est$logLik
    }
    
    prob <- pchisq( -2*deltaLogLik,1,lower.tail =FALSE)
    est$prob <- prob
    est$deltaLogLik <- deltaLogLik
    names(est$prob) <- est$param
    names(est$deltaLogLik) <- est$param
  }
  
  
  coef <- est$coefficients
  if (model %in% c("SLD","CDM","CDMsym","UQTunknown")){
    coef <- coef[-est$npar]
  }
  try({
    est$fitted.values <- as.vector( x %*% coef)
    ## pi schätzen
    e <- exp(est$fitted.values)
    est$pi <- mean(e/(1+e))
    est$fit.n <- cnt
  }, silent=T)
  
  # SE: Scheers & Dayton 1988: propagation of error  page 970 
  # ableitung von pi = e/(1+e) nach den coeffizienten (gradient)
  #   ncoef <- length(coef)
  #   pi.grad <- rep(NA,ncoef)
  #   for (i in 1:length(coef)){
  #     pi.grad[i] <- mean( e*x[,i] / (1+ e)^2)
  #   }
  #   est$piSE <- sqrt(  sum( pi.grad %*% est$vcov[1:ncoef,1:ncoef] %*% pi.grad))
  
  # Scheers: df
  # multiple group models: andere df
  #   if (model %in% c("SLD","CDM","CDMsym","UQTunknown")){
  #     est$df <- 2*length(table(x))- (est$npar-1) -3 
  #   }else{
  #     est$df <- length(table(x)) - est$npar -1 
  #   }
  
  #   print("pi: first calc odds: e/(1+e) ; then mean() . BETTER ESTIMATE!")
  #   print(est$pi)
  #   print("pi: first mean of fitted.values, then odds e/(1+e)")
  #   e <- mean(est$fitted.values)
  #   print(exp(e)/(1+exp(e)))
  est$call <- match.call()  
  class(est) <- "RRlog"
  return(est)  
}




# formula interface: from formula to design matrix

#' @export
RRlog.formula <- function(formula,data=list(),model,p,group, ...){
  model <- match.arg(model,c("Warner","UQTknown","UQTunknown","Mangat","Kuk","FR","Crosswise","CDM","CDMsym","SLD"))
  
  if ( model %in% c("UQTunknown","SLD","CDM","CDMsym") &&  !missing(data) ){
    try({data <- as.data.frame(data)
         group <-  eval(substitute(group),data, parent.frame())
    },silent=T)
  }
  mf <- model.frame(formula=formula, data=data,na.action=na.omit)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  
  names(y) <- all.vars(formula)[1]
  
  # send to default function:
  est <- RRlog.default(x,y,model,p,group, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}


# predict values in response variable (different for logisitc regression ?!)
# predict.RRlog <- function(object, newdata=NULL, ...)
# {
#   if(is.null(newdata))
#     y <- fitted(object)
#   else{
#     if(!is.null(object$formula)){
#       ## model has been fitted using formula interface
#       x <- model.matrix(object$formula, newdata)
#     }
#     else{
#       x <- newdata
#     }
#     y <- as.vector(x %*% coef(object))
#     y <- exp(y)/(1+exp(y))
#   }
#   y
# }

#' @aliases RRlog
#' @method print RRlog
#' @export
print.RRlog <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients,5))
}


#' @aliases RRlog
#' @method summary RRlog
#' @export
summary.RRlog <- function(object, ...)
{
  se <- sqrt(abs(diag(object$vcov)))
  #   tval <- coef(object) / se
  wald_chi <- (object$coef/se)^2
  if (object$model  == "SLD"){
    wald_chi[object$npar] <- ((1-object$coef[object$npar])/se[object$npar])^2 
  }
  TAB <- cbind(Estimate = object$coef,
               StdErr = se,
               "Wald test"=wald_chi,
               "Pr(>Chi2,df=1)"=1-pchisq(wald_chi, 1)
               #                oddsRatio = exp(coef(object)),
               #                StdErr = exp(se)
  )
  if (!is.null(object$prob)){
    TAB <- cbind(TAB,
                 "deltaG2"=-2*object$deltaLogLik,
                 "Pr(>deltaG2)" = object$prob)
  }
  #   index <- pmatch("(Intercept)",object$param)
  #   if (!is.na(index) && !is.null(object$prob)){
  #     TAB[index,5] <- NA
  #     TAB[index,6] <- NA
  #   }
  # next lines: only if odds-ratios are printed
  #   if (object$model %in% c("UQTunknown","SLD","CDM","CDMsym")){
  #     TAB[length(se),3]=NA
  #     TAB[length(se),4]=NA
  #   }
  fitInfo <- cbind(n=object$n, 
                   logLik= object$logLik) 
  #                    df=object$df, 
  #                    AIC=-2*object$logLik+2*object$npar,
  #                    BIC=-2*object$logLik+log(object$n)*object$npar)
  colnames(fitInfo)=c("n","logLik") #,"df","AIC","BIC")
  rownames(fitInfo)=c("")
  
  modelInfo <- paste(object$model," with ",object$pString,sep="")
  if (object$n.dq>0){
    modelInfo <- paste0(modelInfo," (n=",object$n-object$n.dq,") combined with DQ (n=",object$n.dq,")")
  }
  if (object$model=="Kuk"){
    modelInfo <- paste(modelInfo," (",object$rep," repetition/s)",sep="")
  }
  res <- list(call=object$call,modelInfo=modelInfo,
              coefficients=TAB,fitInfo=fitInfo,
              model=object$model,pi=object$pi,piSE=object$piSE)
  class(res) <- "summary.RRlog"
  res
}


#' @aliases RRlog
#' @method print summary.RRlog
#' @export
print.summary.RRlog <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nModel:\n")
  write(x$modelInfo,"")
  cat("\nModel fit:\n")
  print(x$fitInfo)
  cat("\n")
  printCoefmat( round(x$coefficients,5))
  cat("\n")
  if(x$model == "SLD") 
    cat("Note that the parameter t is tested against the null hypothesis that all participants answered truthfully (i.e., H0: t=1).")
  #   cat("\nEstimate of pi:\n")
  #   write(paste0("pi = ",round(x$pi,6))) #" (use RRuni to get standard errors"))
  #                piSE = ",round(x$piSE,6),") 
}

#' @aliases RRlog
#' @method logLik RRlog
#' @export
logLik.RRlog <- function(object, ...){
  return(object$logLik)
}

#' @aliases RRlog
#' @method vcov RRlog
#' @export
vcov.RRlog <- function(object, ...){
  return(object$vcov)
}
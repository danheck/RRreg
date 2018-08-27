#' Logistic randomized response regression
#' 
#' A dichotomous variable, measured once or more per person by a randomized 
#' response method, serves as dependent variable using one or more continuous 
#' and/or categorical predictors.
#' 
#' @param formula specifying the regression model, see \code{\link{formula}}
#' @param data \code{data.frame}, in which variables can be found (optional)
#' @param model Available RR models: \code{"Warner"}, \code{"UQTknown"}, 
#'   \code{"UQTunknown"}, \code{"Mangat"}, \code{"Kuk"}, \code{"FR"}, 
#'   \code{"Crosswise"}, \code{"Triangular"}, \code{"CDM"}, \code{"CDMsym"}, 
#'   \code{"SLD"}, \code{"custom"}. See \code{vignette("RRreg")} for details.
#\code{"custom"} (custom: the argument \code{p} defines the sensitivity and specificity of the binary RR response, i.e., the probabilities P(1 | 1) and P(0 | 0))
#'
#' @param p randomization probability/probabilities (depending on model, see 
#'   \code{\link{RRuni}} for details)
#' @param group vector specifying group membership. Can be omitted for single-group RR designs (e.g., Warner). For two-group RR designs (e.g., \code{CDM} or \code{SLD}), use 1 and 2 to indicate the group membership, matching the respective randomization probabilities \code{p[1], p[2]}. If an RR design and a direct question (DQ) were both used in the study, the group indices are set to 0 (DQ) and 1 (RR; 1 or 2 for two-group RR designs). This can be used to test, whether the RR design leads to a different prevalence estimate by including a dummy variable for the question format (RR vs. DQ) as predictor. If the corresponding regression coefficient is significant, the prevalence estimates differ between RR and DQ. Similarly, interaction hypotheses can be tested (e.g., the correlation between a sensitive attribute and a predictor is only found using the RR but not the DQ design). Hypotheses like this can be tested by including the interaction of the DQ-RR-dummy variable and the predictor in \code{formula} (e.g., \code{RR ~ dummy*predictor}).
#' @param n.response number of responses per participant, e.g., if a participant 
#'   responds to 5 RR questions with the same randomization probability \code{p} 
#'   (either a single number if all participants give the same number of responses or a vector)
#' @param LR.test test regression coefficients by a likelihood ratio test, i.e., 
#'   fitting the model repeatedly while excluding one parameter at a time (each 
#'   nested model is fitted only once, which can result in local maxima).
#'   The likelihood-ratio test statistic \eqn{G^2(df=1)} is reported in the table 
#'   of coefficiencts as \code{deltaG2}.
#' @param fit.n Number of fitting replications using random starting values to avoid local maxima
#' @param EM.max maximum number of iterations of the EM algorithm. 
#'   If \code{EM.max=0}, the EM algorithm is skipped.
#' @param optim.max Maximum number of iterations within each run of \code{optim}
# @param start starting values for optimization. Might be useful if model does not converge with default starting values.
#' @param ... ignored
#' 
#' @details The logistic regression model is fitted first by an EM algorithm, 
#'   in which the dependend RR variable is treated as a misclassified binary variable 
#'   (Magder & Hughes, 1997). The results are used as starting values for a 
#'   Newton-Raphson based optimization by \code{\link{optim}}.
#' @author Daniel W. Heck
#' @seealso \code{vignette('RRreg')} or 
#'   \url{http://www.dwheck.de/separate_content/RRregManual/index.html} for a 
#'   detailed description of the RR models and the appropriate definition of \code{p} 
#' @return Returns an object \code{RRlog} which can be analysed by the generic 
#'   method \code{\link{summary}}. In the table of coefficients, the column \code{Wald} 
#'   refers to the Chi^2 test statistic which is computed as Chi^2 = z^2 = Estimate^2/StdErr^2.
#'   If \code{LR.test = TRUE}, the test statistic \code{deltaG2} is
#'   the likelihood-ratio-test statistic, which is computed by fitting a nested logistic
#'   model without the corresponding predictor.
#'   
#' 
#' @references van den Hout, A., van der Heijden, P. G., & Gilchrist, R. (2007). 
#'   The logistic regression model with response variables subject to randomized response. 
#'   \emph{Computational Statistics & Data Analysis, 51}, 6060-6069. 
#' 
#' @examples
#' # generate data set without biases
#' dat <- RRgen(1000,pi=.3,"Warner",p=.9)
#' dat$covariate <- rnorm(1000)
#' dat$covariate[dat$true==1] <- rnorm(sum(dat$true==1),.4,1)
#' # analyse
#' ana <- RRlog(response~covariate,dat,"Warner", p=.9, fit.n = 1)
#' summary(ana)
#' # check with true, latent states:
#' glm(true~covariate, dat, family=binomial(link="logit"))
#' @export
RRlog <- function(formula, data, model, p, group, n.response=1, LR.test=TRUE, 
                  fit.n=3, EM.max=1000, optim.max=500,  ...) UseMethod("RRlog")

# choose model and construct S3 method 'RRlog' 

#' @export
RRlog.default <-function(formula, data, model, p, group, n.response=1, LR.test=TRUE, fit.n=1, 
                         EM.max=1000, optim.max=500, ...){
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
  n <- length(y)
  model <- match.arg(model,c("Warner","UQTknown","UQTunknown","Mangat",
                             "Kuk","FR","Crosswise","Triangular",
                             "CDM","CDMsym","SLD","custom"))
  
  checkPerfectSeparation(y,x,n.response)
  
  # no DQ format included
  if (missing(group) || is.null(group)){
    group=rep(1,n)
  }
  if(length(n.response == 1)){
    n.response <- rep(n.response, n)
  }
  
  
  RRcheck.log(model,y,p,group, n.response, names(y)[1])   ## for FR and Kuk: n model check for dichot. response
  
  # get estimates repeatedly (because of local minima and starting points)
  est <- list(logLik = -Inf, coef=Inf)
  
  ## EM: get sensitivity P(1 | 1) and specificity P(0 | 0) from missclassification matrices (separately for groups)
  # not working for n.resposne > 1 currently
  if(all(n.response <= 1)){
    uni <- RRuni(response=y, model=model, p=p, group=group)
    
    sens <- rep(1, n)
    spec <- rep(1, n)
    #   if(model == "custom"){
    #     sens[group == 1] <- p[1]
    #     sens[group == 1] <- p[2]
    if (!(is2group(model))){
      P <- getPW(model, p, group=1)
      sens[group == 1] <- P[2,2]
      spec[group == 1] <- P[1,1]
    }else{
      P1 <- getPW(model, p, group=1, summary(uni)$coef[2,1])
      sens[group == 2] <- P1[2,2]
      spec[group == 2] <- P1[1,1]
      P2 <- getPW(model, p, group=2, summary(uni)$coef[2,1])
      sens[group == 2] <- P1[2,2]
      spec[group == 2] <- P1[1,1]
    }
  }
  
  for (cnt in 1:fit.n){
    # | ( (est$logLik==-Inf |max(abs(est$coef))>fit.bound) & cnt <fit.n[2]) # old scheme: fit.n=c(5,20)
    
    
    ####################### EM Algorithm: Magder & Hughes (1997)
    if(all(n.response <= 1)){
      # random dependent variable to get different starting values for EM:
      yy <- ifelse(rbinom(n, 1, max(.5,cnt/fit.n)) == 1, 
                   rbinom(n, 1, .5), y)
      #     print(sum(yy == y) /n)   # proporiton of identical RR outcomes for starting values
      glm.mod <- glm(cbind(yy, n.response-yy) ~ x+0, family = binomial(link = "logit"))
      
      EM.cnt <- 0
      repeat{
        EM.cnt <- EM.cnt+1
        eb <- glm.mod$fitted.values
        beta <- glm.mod$coef
        
        # Expectation: not corret for n.response > 1
        EY <- ifelse(y==0,
                     (1-sens)*eb/((1-sens)*eb+    spec*(1-eb)),
                     sens*eb    /(sens*eb    +(1-spec)*(1-eb)))

        # Maximization:
        suppressWarnings(glm.mod <- #glm(cbind(EY, n.response-EY)~x+0, family = binomial(link = "logit"),
                           #   control=list(maxit=5000), start=beta))
                           glm.fit(x, EY, family = binomial(link="logit"), 
                                   control=list(maxit=5000), start=beta))
        # stopping criterium
        if(EM.cnt >= EM.max | sqrt(sum((glm.mod$coef - beta)^2)) < 1e-5) break
      }
       
      # starting values for optim
      start <- coef(glm.mod)
      if(is2group(model)){
        start <- c(start, summary(uni)$coef[2,1])   # par2: e.g. t for SLD
      }
    }else{
      start <- rnorm(ncol(x), 0, .5)
      if(is2group(model))
        start <- c(start, runif(1, .4,.8))
    }
    
    ####################### optim estimation of full likelihood
    try(est2 <- RRlog.fit(model, x, y, n.response, p, start, group, maxit=optim.max))
    
    if (!is.na(est2$logLik) && est2$logLik > est$logLik){ 
      est <- est2
    }
  }
  if(is.null(est$convergence) || est$convergence != 0)
    warning("Model did not convergence in ML estimation using optim.")
  
  try(if(est$convergence != 0)
    warning(paste0("optim$convergence=",est$convergence,
                   ". Check convergence of model (e.g. by refitting using fit.n=c(20,20).")), silent=T)
  
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
    
    if(ncoef > 1){
      start <- est$coefficients
      for (i in 1:ncoef){
        xx <- x[,-i,drop=FALSE]
        # xx[,i] <- rep(0,length(y))
        try({
          est.res <- RRlog.fit(model, xx, y, n.response, p, start[-i], 
                               group, setPar2=-1, maxit=optim.max)
         
          # [bugfix]: sanity check for local maxima -> crosscheck LR test with Wald test
          deltaLogLik[i] <- est.res$logLik - est$logLik
          chi2_wald <- est$coef[i]^2 / est$vcov[i,i]
          if(abs(deltaLogLik[i] - chi2_wald) >  1) {
            est.res2 <- RRlog.fit(model, xx, y, n.response, p, start[-i]*.2, 
                                  group, setPar2=-1, maxit=optim.max)
            if(est.res2$logLik > est.res$logLik)
              deltaLogLik[i] <- est.res2$logLik - est$logLik
          }
        })
      }
      # multi group models: additional parameter
      if (is2group(model)){
        par2fix <- fix.par2(model)
        try(est.res <- RRlog.fit(model, x, y, n.response, p, start[1:ncol(x)], 
                                 group, setPar2=par2fix, maxit=optim.max))
        cnt <- 0
        while(is.na(est.res$logLik) && cnt<5){
          cnt <- cnt + 1
          try(est.res <- RRlog.fit(model, x, y, n.response, p, rnorm(ncol(x),0,.2), 
                                   group, setPar2=par2fix, maxit=optim.max))}
        deltaLogLik[est$npar] <- est.res$logLik - est$logLik
      }
    }
    prob <- pchisq( -2*deltaLogLik,1,lower.tail =FALSE)
    est$prob <- prob
    est$deltaLogLik <- deltaLogLik
    names(est$prob) <- est$param
    names(est$deltaLogLik) <- est$param
  }
  
  
  coef <- est$coefficients
  if (is2group(model)){
    coef <- coef[-est$npar]
  }
  try({
    latent.values <- as.vector( x %*% coef)
    if (!is2group(model)){
      P <- getPW(model, p)
      est$fitted <- P[2,1]+(P[2,2] - P[2,1])/(1+exp(-latent.values))
    } else{
      est$fitted <- rep(NA, length(latent.values))
      for (i in 1:2){
        P <- getPW(model, p, par2fix, group = i)
        est$fitted[group == i] <- P[2,1]+(P[2,2] - P[2,1])/(1+exp(-latent.values[group == i]))     
      }
    }
    ## pi estimate
    e <- exp(latent.values)
    est$pi <- mean(e/(1+e))  #mean cheating rate
    est$fit.n <- fit.n
  }, silent = TRUE)
  
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
RRlog.formula <- function(formula, data=list(), model, p, group, n.response=1,...){
  model <- match.arg(model,c("Warner","UQTknown","UQTunknown","Mangat",
                             "Kuk","FR","Crosswise","Triangular",
                             "CDM","CDMsym","SLD", "custom"))
  
  if (!missing(data) ){
    try({data <- as.data.frame(data)
    group <-  eval(substitute(group),data, parent.frame())
    },silent=T)
    try({data <- as.data.frame(data)
    n.response <-  eval(substitute(n.response),data, parent.frame())
    },silent=T)
  }
  mf <- model.frame(formula=formula, data=data,na.action=na.omit)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  
  names(y) <- all.vars(formula)[1]
  
  # send to default function:
  est <- RRlog.default(x, y, model, p, group, n.response, ...)
  est$call <- match.call()
  est$formula <- formula
  est$model.frame <- mf
  est
}


#' Predict Individual Prevalences of the RR Attribute
#' 
#' Predictions of the RR logistic regression model for the individual probabilities (or logits) 
#' of having the sensitive RR attribute, or of the probability of the RR responses.
#' 
#' @param object A fitted \code{\link{RRlog}} model
#' @param newdata An optional vector, matrix, or data.frame with values on the 
#'     predictor variables. Note that for matrices, the order of predictors 
#'     should match the order of predictors in the formula. Uses the fitted 
#'     values of the model if omitted.
#' @param type \code{"response"} returns predicted probabilities for the (observable) RR responses,
#'     \code{"link"} returns predicted logit-values for the (latent) sensitive attribute, 
#'     and \code{"attribute"} returns predicted probabilities of having the (latent) sensitive attribute.
#' @param se.fit Return standard errors for the predicted values in addition to confidence intervals.
#'     SEs on the logit scale are computed using the observed Fisher information from the fitted model. 
#'     Standard errors for the probability scale are computed using the delta method.
#' @param ci Confidence level for confidence interval. If \code{ci=FALSE}, no 
#'     confidence interval is returned. Confidence intervals  on the probability 
#'     scale (if \code{type="response"} or \code{type="attribute"}) are computed 
#'     based on the CI on the logit-scale using the inverse link function
#'     (hence, the CI will in general not be symmetric). 
#' @param ... ignored
#' 
#' @return either a vector of predicted values or a matrix with columns for the 
#'     point estimates, confidence interval, and standard errors (if \code{se.fit=TRUE} and \code{ci=.95}).
#' @export
predict.RRlog <- function (object, newdata = NULL, type = c("link", "response", "attribute"), 
                           se.fit = FALSE, ci = 0.95, ...)
{
  if (is.null(newdata)){
    x <- model.matrix(object$formula, object$model.frame)
  } else {
    if (!is.null(object$formula)){
      newdata <- data.frame(newdata, predict=1)
      ## model has been fitted using formula interface
      pred.formula <- update(object$formula, predict~.)
      x <- model.matrix(pred.formula, newdata, na.action=na.exclude)
    }
    else{
      x <- newdata
    }
  }
  # logit values:
  ypred <- predict <- as.vector(x %*% coef(object))
  vcov <- vcov(object)
  
  if(is2group(object$model)){
    k <- length(object$coef)
    vcov <- vcov[1:k,1:k]
  }
  
  predict.vcov <- x %*% vcov %*% t(x)
  predict.se <- sqrt(diag(predict.vcov))
  zcrit <- qnorm((1-ci)/2, lower.tail=F)
  ci.lower <- predict - predict.se*zcrit
  ci.upper <- predict + predict.se*zcrit
  
  type  <- match.arg(type, c("link", "response", "attribute"))
  
  # predict probability for RR response ("yes" or "no" etc.)
  if (type == "response"){
    # get complete link function (logit + RR response)
    p <- getPW(object$model, object$p)[2,]
    linkfun <- RRloglink(p)
    predict <- linkfun$linkinv(ypred)
    ci.lower <- linkfun$linkinv(ci.lower)
    ci.upper <- linkfun$linkinv(ci.upper)  
    
    # standard errors for probability scale (RR responses): delta method
    # https://cran.r-project.org/web/packages/modmarg/vignettes/delta-method.html
    # c <- p[1]
    d <- p[2] - p[1]
    deriv <- d * exp(-ypred)/(1+exp(-ypred))^2   # derivation of:    g(y) = c + d/(1 + exp(-y))
    j <- deriv * x                       # Jacobian for each prediction
    vcov_prob <- j %*% vcov %*% t(j)
    predict.se <- sqrt(diag(vcov_prob)) # computed SE refer to response probability scale
  }
  
  # predict sensitive attribute  (one step between "link" / "response")
  if (type == "attribute"){
    predict <- plogis(ypred)
    ci.lower <- plogis(ci.lower)
    ci.upper <- plogis(ci.upper)
    
    # standard errors for probability scale (latent attribute): delta method
    # https://cran.r-project.org/web/packages/modmarg/vignettes/delta-method.html
    deriv <- exp(-ypred)/(1+exp(-ypred))^2  # derivation of:    g(y) = 1/(1+exp(-y))
    j <- deriv * x                # Jacobian for each prediction
    vcov_prob <- j %*% vcov %*% t(j)
    predict.se <- sqrt(diag(vcov_prob)) # computed SE refer to latent-attribute probability scale
  } 
  
  # output
  if (!se.fit && ci == 0)
    return(predict)
  else if (se.fit && ci != 0)
    return (cbind(predict=predict, se=predict.se, ci.lower=ci.lower, ci.upper=ci.upper))
  else if (se.fit)
    return (cbind(predict=predict, se=predict.se))
  else
    return (cbind(predict=predict, ci.lower=ci.lower, ci.upper=ci.upper))
}

#' @export
print.RRlog <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients,5))
}


#' @export
summary.RRlog <- function(object, ...)
{
  se <- sqrt(abs(diag(object$vcov)))
  if(length(se) == 0)
    se <- rep(NA, length(object$coef))
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


#' @export
print.summary.RRlog <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nRR Model:\n")
  write(x$modelInfo,"")
  cat("\nModel fit:\n")
  print(x$fitInfo)
  cat("\n")
  printCoefmat( round(x$coefficients,5))
  cat("\n")
  if(x$model == "SLD") 
    cat("Note that the parameter t is tested against the null hypothesis 
        that all participants answered truthfully (i.e., H0: t=1).\n")
}

#' @export
fitted.RRlog <- function(object, ...){
  return(object$fitted)
}

#' @export
logLik.RRlog <- function(object, ...){
  return(object$logLik)
}

#' @export
vcov.RRlog <- function(object, ...){
  return(object$vcov)
}

#' Plot Logistic RR Regression
#' 
#' Plot predicted logit values/probabilities of a randomized response logistic regression model.
#' 
#' @param x a fitted \link{RRlog} object
#' @param predictor character name of a predictor of the model to be fitted
#' @param center.preds whether to compute predictions by assuming that all other 
#'     predictors are at their respective mean values (if \code{FALSE}: all other 
#'     predictors are set to zero)
#' @param plot.mean whether to plot the mean of the predictor as a vertical line
#' @param ci level for confidence intervals. Use \code{ci=0} to omit.
#' @param xlim if provided, these boundaries are used for the predictor on the x-axis
#' @param steps number of steps for plotting
#' @param ... other arguments passed to the function \link{plot} (e.g., \code{ylim=c(0,1)}).
#' @inheritParams predict.RRlog
#' 
#' @seealso \code{\link{predict.RRlog}}
#'  
#' @examples
#'  # generate data
#'  n <- 500
#'  x <- data.frame(x1=rnorm(n))
#'  pi.true <- 1/(1+exp(.3+1.5*x$x1))
#'  true <- rbinom(n, 1, plogis(pi.true))
#'  dat <- RRgen(n, trueState=true, model="Warner", p=.1)
#'  x$response <- dat$response
#'  
#'  # fit and plot RR logistic regression
#'  mod <- RRlog(response ~ x1, data=x, model="Warner", p=.1)
#'  plot(mod, "x1" ,ci=.95, type = "attribute", ylim = 0:1)
#'  
#' @export
plot.RRlog <- function(x, predictor = NULL, type = c("link", "response", "attribute"),
                       center.preds = TRUE, plot.mean = TRUE, ci = .95, 
                       xlim = NULL, steps = 50, ...){
  type <- match.arg(type, c("link", "response", "attribute"))
  
  # single predictor: choose automatically
  beta <- x$coef
  if (missing(predictor) | is.null(predictor)){
    if(length(beta) == 2 & names(beta)[1] == "(Intercept)"){
      predictor <- names(beta)[2]
    }else{
      stop("Please specify a predictor for the x-axis.")
    }
  }
  
  # predict values
  if (missing(xlim) | is.null(xlim)){
    predvals <- x$model.frame[,predictor]
    xlim <- c(min(predvals), max(predvals))
  }
  xx <- seq(xlim[1], xlim[2], length.out=steps)
  
  # construct new design matrix
  X.old <-  model.matrix(x$formula, x$model.frame)
  X.new <- colMeans(X.old[,colnames(X.old) != predictor, drop=F])
  X <- cbind(matrix(X.new, length(xx),length(beta)-1,  byrow=T), predictor=xx)
  colnames(X) <- c(names(X.new), predictor)
  
  pred <- predict(x, newdata = X, ci = ci, type = type) 
  if(is.null(dim(pred))){
    yy <- pred
  } else {
    yy <- pred[,"predict"]
  }
  
  plot(xx, yy, main=paste0("Predictions for RRlog: ",substitute(x), " (type=", type[1], ")"), type="l", 
       xlab=substitute(predictor), ylab=colnames(x$model.frame)[1], ...)
  if(!is.null(dim(pred))){
    polygon(c(xx, rev(xx)), c(pred[,2], rev(pred[,3])),border=NA,
            col=adjustcolor("gray", alpha.f=.5))
  }
  if(plot.mean){
    abline(v=mean(x$model.frame[,predictor]) , lty=3)
  }
  lines(x=xx, y=yy)
  
}
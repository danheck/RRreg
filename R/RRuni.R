#' Univariate analysis of randomized response data
#' 
#' Analyse a data vector \code{response} with a specified RR model (e.g., \code{Warner}) with known randomization probability \code{p}
#' 
#'  @param response either vector of responses containing 0 (No) and 1 (Yes) or name of response variable in \code{data}. In Kuk's card playing method (\code{Kuk}), the observed response variable gives the number of red cards. For the Forced Response (\code{FR}) model, response values are integers from 0 to (m-1), where 'm' is the number of response categories. 
#'  @param data optional \code{data.frame} containing the response variable
#'  @param model choose RR model. Available models: \code{"Warner","UQTknown","UQTunknown","Mangat","Kuk","FR","Crosswise","CDM","CDMsym","SLD", "mix.norm", "mix.exp","mix.unknown"}. See \code{vignette('RRreg')} for detailed specifications.
#'  @param p randomization probability. For the Cheating Detection Model (\code{CDM}) or the Stochastic Lie Detector (\code{SLD}): a vector with two values. For the Forced Response model (\code{FR}): a vector of the length of the number of categories
#'  @param group a group vector of the same length as \code{response} containing values 1 or 2, only required for two-group models, which specify different randomization probabilities for two groups, e.g., \code{CDM} or \code{SLD}. If a data.frame \code{data} is provided, the variable \code{group} is searched within it.
#' @param MLest if \code{TRUE}, least-squares estimates of pi outside of [0,1] are corrected to obtain maximum likelihood estimates
#' @return an \code{RRuni} object, can by analyzed by using \code{\link{summary}}
#' @seealso \code{vignette('RRreg')} or \url{https://dl.dropboxusercontent.com/u/21456540/RRreg/index.html} for a detailed description of the RR models and the appropriate definition of \code{p} 
#' @examples 
#' # Generate responses of 1000 people according to Warner's model
#' # with an underlying true proportion of .3
#' genData <- RRgen(n=1000, pi=.3, model="Warner", p=.7)
#' # Analyse univariate data to estimate 'pi'
#' analyse <- RRuni(response=genData$response, model="Warner", p=.7)
#' summary(analyse)
#' 
#' # Generate data in line with the Stochastic Lie Detector 
#' # assuming that 90% of the respondents answer truthfully
#' genData2 <- RRgen(n=1000, pi=.3, model="SLD", p=c(.2,.8), complyRates=c(.8,1),groupRatio=0.4)
#' analyse2 <- RRuni(response=genData2$response, model="SLD", p=c(.2,.8), group=genData2$group)
#' summary(analyse2)
#' @export
RRuni<-function(response, data, model, p,group = NULL, MLest=TRUE){
  # extract column 'response' from data.frame 'data'
  if ( !missing(data)){
    try( {data <- as.data.frame(data)
          response <-  eval(substitute(response),data, parent.frame())
          })
  }
  model <- match.arg(model, c("Warner","UQTknown","UQTunknown","Mangat","Kuk",
                       "FR","Crosswise","CDM","CDMsym","SLD", "mix.norm", 
                       "mix.exp","mix.unknown"))
  if ( is2group(model) &&  !missing(data) ){
      try({data <- as.data.frame(data)
         group <-  eval(substitute(group),data, parent.frame())
           },silent=T)
  }
  RRcheck.xpgroup(model,response,p,group,"response")
  
  switch(model,
         "Warner" = res <- RRuni.Warner(response,p),
         "UQTknown" = res <- RRuni.UQTknown(response,p),
         "UQTunknown" = res <- RRuni.UQTunknown(response,p,group),
         "Mangat" = res <- RRuni.Mangat(response,p),
         "Kuk" = res <- RRuni.Kuk(response,p),
         "FR" = res <- RRuni.FR(response,p),
         "Crosswise" = res <- RRuni.Crosswise(response,p),
         "SLD" = res <- RRuni.SLD(response,p,group),
         "CDM" = res <- RRuni.CDM(response,p,group),
         "CDMsym" = res <- RRuni.CDMsym(response,p,group),
         "mix.norm" = res <- RRuni.mix.norm(response,p),
         "mix.exp" = res <- RRuni.mix.exp(response,p),
         "mix.unknown" = res <- RRuni.mix.unknown(response,p, group)
         )
  if (MLest && ! (model %in% c("mix.norm","mix.exp", "mix.unknown"))){
    res$pi <- RRcheck.param(res$pi)
    if (model=="UQTunknown"){
      res$piUQ <- RRcheck.param(res$piUQ)
    }else if (model %in% c("CDM","CDMsym")){
      res$gamma <- RRcheck.param(res$gamma)
      res$beta <- RRcheck.param(res$beta)
    }else if (model == "SLD"){
      res$t <- RRcheck.param(res$t)
    }
  }
  class(res)="RRuni"
  return(res)
}

#' @aliases RRuni
#' @method print RRuni
#' @export
print.RRuni<-function(x,...){
  cat("RR model: \n")
  write(x$call,"")
  cat("\nEstimate of pi:\n")
  write( paste0(round(x$pi,6)," (",round(x$piSE,6),") "),"")
  if (x$model=="SLD"){
    cat("\nEstimate of t:\n")
    write(paste0(round(x$t,6), " (",round(x$tSE,6),")"),"")
  }else if (x$model %in% c("CDM","CDMsym")){
    cat("\nEstimate of gamma:\n")
    write(paste0(round(x$gamma,6), " (",round(x$gammaSE,6),")"),"")
  }else if (x$model == "UQTunknown"){
    cat("\nEstimate for prevalence of unrelated question:\n")
    write(paste(round(x$piUQ,6), " (",round(x$piUQSE,6),")",sep=""),"")
  }else if (x$model == "mix.unknown"){
    cat("\nEstimate for unrelated question:\n")
    write(paste(round(x$piUQ,6), " (",round(x$piUQSE,6),")",sep=""),"")
  }
}

#' @aliases RRuni
#' @method summary RRuni
#' @export
summary.RRuni<-function(object,...){
  zval <- object$pi/object$piSE
  TAB <- cbind(Estimate = object$pi,
               StdErr = object$piSE,
               z=zval,
               "Pr(>|z|)"=pnorm(zval,lower.tail=F))
  if (object$model != "FR"){
    rownames(TAB) <- "pi"
  }else{
    rownames(TAB) <- paste("pi",0:(length(object$pi)-1),sep="")
  }
  if (object$model=="UQTunknown"){
    zval_piUQ=object$piUQ/object$piUQSE
    TAB=rbind(TAB,
              cbind(object$piUQ, object$piUQSE,zval_piUQ,pnorm(zval_piUQ,lower.tail=F)))
    rownames(TAB)=c("pi","piUQ")
  }
  if (object$model=="SLD"){
    ifelse
    zval_t=(object$t-1)/object$tSE  ##teste  t gegen 1 !
    TAB=rbind(TAB,
              cbind(object$t, object$tSE,zval_t,pnorm(zval_t)))
    rownames(TAB)=c("pi","t")
  }
  if (object$model %in% c("CDM")){
    zval_b=object$beta/object$betaSE
    zval_g=object$gamma/object$gammaSE
    TAB=rbind(TAB,
              cbind(object$beta,object$betaSE,zval_b,pnorm(zval_b,lower.tail=F)),
              cbind(object$gamma,object$gammaSE,zval_g,pnorm(zval_g,lower.tail=F)))
    rownames(TAB)=c("pi","beta","gamma")
  }
  if (object$model %in% c("CDMsym")){
    zval_g=object$gamma/object$gammaSE
    TAB=rbind(TAB,
#               cbind(object$beta,object$betaSE,zval_b,pnorm(-abs(zval_b))),
              cbind(object$gamma,object$gammaSE,zval_g,pnorm(zval_g,lower.tail=F)))
    rownames(TAB)=c("pi","gamma")
  }
  if (object$model %in% c("mix.unknown")){
    zval_g=object$piUQ/object$piUQSE
    TAB=rbind(TAB,
              #               cbind(object$beta,object$betaSE,zval_b,pnorm(-abs(zval_b))),
              cbind(object$piUQ,object$piUQSE,zval_g,pnorm(zval_g,lower.tail=F)))
    rownames(TAB)=c("pi","piUQ")
  }
  res <- list(call=object$call,n=object$n,
              coefficients=TAB, model=object$model)
  class(res) <- "summary.RRuni"
  return(res)
}


#' @aliases RRuni
#' @method print summary.RRuni
#' @export
print.summary.RRuni<-function(x,...){
  cat("Call:\n")
  write(x$call,"")
  cat("Sample size: ")
  write(x$n,"")
  cat("\n")
  printCoefmat( round(x$coefficients,6))
  if (x$model=="SLD"){
    cat("\n(for the parameter t, i.e. probability of true responding of carriers, the test is H0: t=1; H1: t<1 and the one-sided probability value is given)")
  }
}

# plot.RRuni necessary?

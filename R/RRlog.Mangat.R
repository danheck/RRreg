
# Mangat's Model
RRlog.Mangat <- function(x,y,p,start,group, maxit=1000){
  grad <- rep(NA,ncol(x));
  logLik<- NA;
  coef<- rep(NA,ncol(x));
  iter<- NA;
  hessian<- matrix(NA,ncol=ncol(x),nrow=ncol(x))
  tryCatch(
{est <- optim(par=start,fn=RRlog.Mangat.ll,
              gr=RRlog.Mangat.llgrad, 
#               method="BFGS",
              control=list(fnscale=-1, maxit=maxit),hessian=T,  cov=x,y=y,prand=p,group=group)
 grad <- RRlog.Mangat.llgrad(est$par,x,y,p,group)
 logLik=est$value;
 coef=est$par;
 iter=est$counts;
 hessian=est$hessian
}
,error = function(e) {})

  
#   print(grad)
#   print(grad(RRlog.Mangat.ll,est$par,cov=x,y=y,prand=p,group=group))
  res <- list(model="Mangat",pString=paste0("p = ",round(p,3)),
              coefficients=coef, logLik=logLik,param=colnames(x),
              gradient=grad,hessian=hessian,iter=iter, convergence=est$convergence)
  return(res)
}


# loglikelihoodfunction for Mangat model
RRlog.Mangat.ll = function (param,cov,y,prand,group){
  p <- prand
  rr <- group==1
  e <- exp(cov%*%param)
  vec <- y*log( e+1-p )+(1-y)*log(p)-log(1+e)  ## RR likelihood
  ll <- sum(vec[rr])
  if (min(group)==0){
    vec2 <- y*log(e)-log(1+e)    ## DQ likelihood
    ll <- sum(vec[rr],vec2[!rr])
  }
  return(ll)
}

# gradient
RRlog.Mangat.llgrad=function (param,cov,y,prand,group){
  p <- prand
  rr <- group==1
  gradient <- rep(0,length(param))
  e <- exp(cov%*%param)
  g <- y/(e+1-p) - 1/(1+e)
  for (j in 1: length(gradient)){
    vec <- e[rr]*cov[rr,j]*g[rr]  
    gradient[j] <- sum(vec)
    if (min(group)==0){
      vec2 <- (y[!rr]-e[!rr]/(1+e[!rr])) *cov[!rr,j] 
      gradient[j] <- sum(vec,vec2)
    }
  } 
  return(gradient)  
}

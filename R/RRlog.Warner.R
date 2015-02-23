
# Warner's Model
RRlog.Warner <- function(x,y,p,start,group, maxit=1000){
  grad <- rep(NA,ncol(x));
  logLik<- NA;
  coef<- rep(NA,ncol(x));
  iter<- NA;
  hessian<- matrix(NA,ncol=ncol(x),nrow=ncol(x))
  tryCatch(
    {est <- optim(par=start,fn=RRlog.Warner.ll,
            gr=RRlog.Warner.llgrad, 
#             method="BFGS",
            control=list(fnscale=-1, maxit=maxit),hessian=T,  cov=x,y=y,prand=p,group=group);
     grad <- RRlog.Warner.llgrad(est$par,x,y,p,group);
     logLik=est$value;
     coef=est$par;
     iter=est$counts;
     hessian=est$hessian
    }
    ,error = function(e) {})
#   print(RRlog.Warner.llgrad(est$par,x,y,p,group))
#   print(grad(RRlog.Warner.ll,est$par,cov=x,y=y,prand=p,group=group))
  res <- list(model="Warner",pString=paste0("p = ",round(p,3)),
           coefficients=coef, logLik=logLik,param=colnames(x),
           gradient=grad,hessian=hessian,iter=iter, convergence=est$convergence)
  return(res)
}


# loglikelihoodfunction for Warner model
RRlog.Warner.ll = function (param,cov,y,prand,group){
  p <- prand
  e <- exp(cov%*%param)
  rr <- group==1
  vec <- y*log( p*e+(1-p) )+(1-y)*log((1-p)*e+p)-log(1+e)  ## RR likelihood per person
  ll <- sum(vec[rr])
  if (min(group)==0){     # DQ format
    vec2 <- y*log(e)-log(1+e)
    ll <- sum(vec[rr],vec2[!rr])
  }
  return(ll)
}

# gradient
RRlog.Warner.llgrad=function (param,cov,y,prand,group){
  p <- prand
  gradient <- rep(0,length(param))
  rr <- group==1
  e <- exp(cov%*%param)
  g <- (y*p)/(p*e+(1-p))+(1-y)*(1-p)/((1-p)*e+p)-1/(1+e)
  for (j in 1: length(gradient)){
    vec <- e[rr]*cov[rr,j]*g[rr]         ## RR likelihood
    gradient[j] <- sum(vec)
    if (min(group)==0){
      vec2 <- (y[!rr]-e[!rr]/(1+e[!rr])) *cov[!rr,j]  ##  DQ likelihood
      gradient[j] <- sum(vec,vec2)
    }
  } 
  return(gradient)  
}
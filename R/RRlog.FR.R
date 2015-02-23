
# Forced Response (FR) Model (only dichotomous)
RRlog.FR <- function(x,y,p,start,group, maxit=1000){
  grad <- rep(NA,ncol(x));
  logLik<- NA;
  coef<- rep(NA,ncol(x));
  iter<- NA;
  hessian<- matrix(NA,ncol=ncol(x),nrow=ncol(x))
  tryCatch(
{est=optim(par=start,fn=RRlog.FR.ll,
           gr=RRlog.FR.llgrad, 
#            method="L-BFGS-B",
           control=list(fnscale=-1, maxit=maxit),hessian=T,  cov=x,y=y,prand=p,group=group)
 grad <- RRlog.FR.llgrad(est$par,x,y,p,group)  ## falsch?!
 logLik=est$value;
 coef=est$par;
 iter=est$counts;
 hessian=est$hessian
}
,error = function(e) {})

  
#   print(grad)
#   print(grad(RRlog.FR.ll,est$par,cov=x,y=y,prand=p,group=group))
  res=list(model="FR",pString=paste("p1 = ",round(p[1],3),"; p2 = ",round(p[2],3),sep=""),
           coefficients=coef,logLik=logLik,param=colnames(x),
           gradient=grad, hessian=hessian,iter=iter, convergence=est$convergence)
  return(res);
}

# loglikelihoodfunction for Forced Response model
RRlog.FR.ll = function (param,cov,y,prand,group){
  p1 <- prand[1]  # forced no
  p2 <- prand[2]  # forced yes
  rr <- group==1
  e <- exp(cov%*%param)
  vec <- y*log( (1-p1)*e+p2) + (1-y)*log( p1*e+ 1-p2) - log(1+e)
  ll <- sum(vec[rr])
  if (min(group)==0){
    vec2 <- y*log(e)-log(1+e)
    ll <- ll+sum(vec2[!rr])
  }
  return(ll)
}

# gradient   ### FALSCH
RRlog.FR.llgrad=function (param,cov,y,prand,group){
  p1 <- prand[1];
  p2 <- prand[2];
  rr <- group==1
  gradient <- rep(0,length(param))
  e <- exp(cov%*%param)
  g <- y*(1-p1)/((1-p1)*e+p2) + (1-y)*p1/(p1*e+1-p2) - 1/(1+e)
  for (j in 1: length(gradient)){
    vec <- e[rr]*cov[rr,j]*g[rr]      
    gradient[j] <- sum(vec)
    if (min(group)==0){
      vec2 <- (y[!rr]-e[!rr]/(1+e[!rr])) *cov[!rr,j]  ##  DQ likelihood
      gradient[j] <- sum(vec,vec2)
    }
  }
  return(gradient)  
}
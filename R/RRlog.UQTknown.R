# UQTknown's Model 
RRlog.UQTknown <- function(x,y,p,start,group, maxit=1000){
  grad <- rep(NA,ncol(x));
  logLik<- NA;
  coef<- rep(NA,ncol(x));
  iter<- NA;
  hessian<- matrix(NA,ncol=ncol(x),nrow=ncol(x))
  tryCatch(
{est <- optim(par=start,fn=RRlog.UQTknown.ll,
              gr=RRlog.UQTknown.llgrad, 
#               method="BFGS",
              control=list(fnscale=-1, maxit=maxit),hessian=T,  cov=x,y=y,prand=p,group=group)
 grad <- RRlog.UQTknown.llgrad(est$par,x,y,p,group)
 logLik=est$value;
 coef=est$par;
 iter=est$counts;
 hessian=est$hessian
}
,error = function(e) {})
  
#   print(grad)
#   print(grad(RRlog.UQTknown.ll,est$par,cov=x,y=y,prand=p,group=group))
  res <- list(model="UQTknown", pString=paste0("probability of answering sensitive question = ",p[1],"; prevalence of irrelevant question = ",p[2]),param=colnames(x), coefficients=coef, logLik=logLik,
           gradient=grad,hessian=hessian,
           iter=iter)
  return(res)
}

# loglikelihoodfunction for UQTknown model
RRlog.UQTknown.ll = function (param,cov,y,prand,group){
  p1 <- prand[1]
  p2 <- prand[2]
  rr <- group==1
  e <- exp(cov%*%param)
  vec <- y*log( (p1+p2-p1*p2)*e +p2-p1*p2 ) +
     (1-y)*log( (1-p1-p2+p1*p2)*e+ 1+p1*p2-p2 ) - log(1+e)  ## RR likelihood
  ll <- sum(vec[rr])
  if (min(group)==0){
    vec2 <- y*log(e)-log(1+e)  ## DQ Likelihood
    ll <- sum(vec[rr],vec2[!rr])
  }
  return(ll)
}

# gradient
RRlog.UQTknown.llgrad=function (param,cov,y,prand,group){
  p1 <- prand[1]
  p2 <- prand[2]
  rr <- group==1
  gradient <- rep(0,length(param))
  e <- exp(cov%*%param)
  g <- ((p1+p2-p1*p2)*y)/((p1+p2-p1*p2)*e+p2-p1*p2)+
    (1-y)*(1-p1-p2+p1*p2)/((1-p1-p2+p1*p2)*e+1+p1*p2-p2)-1/(1+e)
  for (j in 1: length(gradient)){
    vec=e[rr]*cov[rr,j]*g[rr]  
    gradient[j] <- sum(vec)
    if (min(group)==0){
      vec2 <- (y[!rr]-e[!rr]/(1+e[!rr])) *cov[!rr,j] 
      gradient[j] <- sum(vec,vec2)
    }
  }
  return(gradient)  
}
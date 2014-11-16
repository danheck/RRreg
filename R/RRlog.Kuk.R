# Kuk's Model 
RRlog.Kuk <- function(x,y,p,start,rep,group, maxit=1000){
  grad <- rep(NA,ncol(x));
  logLik<- NA;
  coef<- rep(NA,ncol(x));
  iter<- NA;
  hessian<- matrix(NA,ncol=ncol(x),nrow=ncol(x))
  tryCatch(
{est <- optim(par=start,fn=RRlog.Kuk.ll,
              gr=RRlog.Kuk.llgrad, 
#               method="BFGS",
              control=list(fnscale=-1, maxit=maxit),hessian=T,  cov=x,y=y,prand=p,rep=rep,group=group)
 grad <- RRlog.Kuk.llgrad(est$par,x,y,p,rep,group)
 logLik=est$value;
 coef=est$par;
 iter=est$counts;
 hessian=est$hessian
}
,error = function(e) {})

  
#   print(grad)
  
  res <- list(model="Kuk",pString=paste("p1 = ",p[1],"; p2 = ",p[2],sep=""), coefficients=coef, logLik=logLik,param=colnames(x),
           gradient=grad,hessian=hessian,rep=rep,
           iter=iter)
  return(res);
}

# loglikelihoodfunction for Kuk's model
RRlog.Kuk.ll = function (param,cov,y,prand,rep,group){
  p1 <- prand[1]
  p2 <- prand[2]
  rr <- group==1
  e <- exp(cov%*%param)  
  vec <- -log(1+e)+log(p1^y*(1-p1)^(rep-y)* e+
               p2^y*(1-p2)^(rep-y)) +log(choose(rep,y))
  ll <- sum(vec[rr])
  if (min(group)==0){
    vec2 <- y*log(e)-log(1+e)
    ll <- sum(vec[rr],vec2[!rr])
  }
  return(ll)
}

# gradient
RRlog.Kuk.llgrad=function (param,cov,y,prand,rep,group){
  p1 <- prand[1]
  p2 <- prand[2]
  rr <- group==1
  gradient <- rep(0,length(param))
  e <- exp(cov%*%param)
  g <- (p1^y*(1-p1)^(rep-y)) / 
    (p1^y*(1-p1)^(rep-y)*e+p2^y*(1-p2)^(rep-y)) - 1/(1+e)
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
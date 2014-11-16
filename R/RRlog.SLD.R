# SLD
RRlog.SLD <- function(x,y,p,start,group,setT=FALSE, maxit=1000){
  grad <- rep(NA,ncol(x)+1);
  logLik<- NA;
  coef<- rep(NA,ncol(x));
  iter<- NA;
  hessian<- matrix(NA,ncol=ncol(x),nrow=ncol(x))
  try({est <- optim(par=start,fn=RRlog.SLD.ll,
                    gr=RRlog.SLD.llgrad, 
                    method="L-BFGS-B",
                    lower = c(rep(-Inf,ncol(x)),0), 
                    upper = c(rep( Inf,ncol(x)),1),
                    control=list(fnscale=-1, maxit=maxit),hessian=T,  
                    cov=x,y=y,prand=p,group=group,setT=setT)
       grad <- RRlog.SLD.llgrad(est$par,x,y,p,group,setT)
       logLik=est$value;
       coef=est$par;
       iter=est$counts;
       hessian=est$hessian
  },silent=T)
  
#   print(grad)
#   print(grad(func=RRlog.SLD.ll,x=est$par,cov=x,y=y,prand=p,group=group,setT=setT))
  res <- list(model="SLD",
              pString=paste("p1 = ",p[1],"; p2 = ",p[2],sep=""),
              coefficients=coef,
              logLik=logLik,param=c(colnames(x),"t"),
              gradient=grad,hessian=hessian,iter=iter)
  return(res)
}

# loglikelihoodfunction for SLD: (param[1] is t)
RRlog.SLD.ll = function (param,cov,y,prand,group,setT){
  p1 <- prand[1]
  p2 <- prand[2]
  m <- length(param)
  t <- param[m]
  if (setT) {t <- 1}
  e <- exp(cov%*%param[-m])
  s1 <- group==1
  s2 <- group==2
  vec1 <- y*log( t*e + (1-p1) ) + (1-y)*log( (1-t)*e + p1 )-log(1+e)  ## RR likelihood
  vec2 <- y*log( t*e + (1-p2) ) + (1-y)*log( (1-t)*e + p2 )-log(1+e)
  ll <- sum(vec1[s1],vec2[s2])
  if (min(group)==0){ ## additional group with direct questioning: group==0
    dq <- group==0
    vec3 <- y*log(e)-log(1+e)
    ll <- ll + sum(vec3[dq])
  }
  return(ll)
}

# gradient
RRlog.SLD.llgrad <- function(param,cov,y,prand,group,setT){
  p1 <- prand[1]
  p2 <- prand[2]
  gradient <- rep(0,length(param))
  m <- length(param)
  t <- param[m]
  s1 <- group==1
  s2 <- group==2
  e <- exp(cov%*%param[-m])
  
  if (setT) {
    t <- 1
  }else{    # d/dt  
    vec1 <- (e*y) /(t*e-p1+1)-(e*(1-y))/(e*(1-t)+p1)
    vec2 <- (e*y) /(t*e-p2+1)-(e*(1-y))/(e*(1-t)+p2)
    gradient[m] <- sum(vec1[s1],vec2[s2])
  }
  
  g1 <- ( (t-1)*(1-y))/(e*(t-1)-p1) + (t*y)/(e*t-p1+1) - 1/(e+1)
  g2 <- ( (t-1)*(1-y))/(e*(t-1)-p2) + (t*y)/(e*t-p2+1) - 1/(e+1)
  # d/d betai
  for (j in 1:(m-1)){
    vec1 <- e[s1] *cov[s1,j]*g1[s1]
    vec2 <- e[s2] *cov[s2,j]*g2[s2]
    gradient[j] <- sum(vec1,vec2)
    if (min(group)==0){ ## additional group with direct questioning: group==0
      dq <- group==0
      vec3 <- (y[dq]-e[dq]/(1+e[dq])) *cov[dq,j] 
      gradient[j] <- gradient[j] + sum(vec3)
    }
  }
  return(gradient)
}
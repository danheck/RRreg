# UQTunknown
RRlog.UQTunknown <- function(x,y,p,start,group,setPiUQ=FALSE, maxit=1000){
  grad <- rep(NA,ncol(x)+1);
  logLik<- NA;
  coef<- rep(NA,ncol(x));
  iter<- NA;
  hessian<- matrix(NA,ncol=ncol(x),nrow=ncol(x))
  tryCatch(
{est <- optim(par=start,fn=RRlog.UQTunknown.ll,
              gr=RRlog.UQTunknown.llgrad, 
#                method="L-BFGS-B",
              lower = c(rep(-Inf,ncol(x)),0), 
              upper = c(rep( Inf,ncol(x)),1),
              control=list(fnscale=-1, maxit=maxit),hessian=T,  
              cov=x,y=y,prand=p,group=group,setPiUQ=setPiUQ) 
 grad <- RRlog.UQTunknown.llgrad(est$par,cov=x,y,p,group,setPiUQ)
 logLik=est$value;
 coef=est$par;
 iter=est$counts;
 hessian=est$hessian
}
,error = function(e) {})

#   print(grad)
#   print(grad(RRlog.UQTunknown.ll,est$par,cov=x,y=y,prand=p,group=group,setPiUQ=setPiUQ))
  res <- list(model="UQTunknown",
              pString=paste0("probability of answering sensitive question = ",p[1],"/",p[2]," (group 1/2)"),
              coefficients=coef,
              logLik=logLik, param=c(colnames(x),"piUQ"),
              gradient=grad,hessian=hessian,iter=iter)
  return(res)
}

# loglikelihoodfunction for UQTunknown: (param[m] is piUQ)
RRlog.UQTunknown.ll = function (param,cov,y,prand,group,setPiUQ){
  p1 <- prand[1]
  p2 <- prand[2]
  m <- length(param)
  piUQ <- param[m]  
  s1 <- group==1
  s2 <- group==2
  if (setPiUQ) piUQ <- 0
  
  e <- exp(cov%*%param[-m])
  vec1 <- y*log( (p1+piUQ-piUQ*p1)*e+piUQ-piUQ*p1 )+
    (1-y)*log( (1-piUQ-p1+p1*piUQ)*e+1-piUQ+p1*piUQ ) -log(1+e)
  vec2 <- y*log( (p2+piUQ-piUQ*p2)*e+piUQ-piUQ*p2 )+
    (1-y)*log( (1-piUQ-p2+p2*piUQ)*e+1-piUQ+p2*piUQ ) -log(1+e)
  ll <- sum(vec1[s1],vec2[s2])
  if (min(group)==0){
    dq <- group==0
    vec3 <- y*log(e)-log(1+e)
    ll <- ll + sum(vec3[dq])
  }
  return(ll)
}

# gradient
RRlog.UQTunknown.llgrad <- function(param,cov,y,prand,group,setPiUQ){
  p1 <- prand[1]
  p2 <- prand[2]
  m <- length(param)
  gradient <- rep(0,m)
  piUQ <- param[m]  # greek gamma =: gg  
  s1 <- group==1
  s2 <- group==2
  if (setPiUQ) piUQ <- 0
  
  e <- exp(cov%*%param[-m])
  g1 <- y*(p1+piUQ-piUQ*p1)/((p1+piUQ-piUQ*p1)*e+piUQ-piUQ*p1)+
    (1-y)*(1-piUQ-p1+p1*piUQ)/((1-piUQ-p1+p1*piUQ)*e+1-piUQ+p1*piUQ)-1/(1+e)
  g2 <- y*(p2+piUQ-piUQ*p2)/((p2+piUQ-piUQ*p2)*e+piUQ-piUQ*p2)+
    (1-y)*(1-piUQ-p2+p2*piUQ)/((1-piUQ-p2+p2*piUQ)*e+1-piUQ+p2*piUQ)-1/(1+e)  
  # d/d betai
  for (j in 1:(m-1)){
    vec1 <- e[s1]*cov[s1,j]*g1[s1]
    vec2 <- e[s2]*cov[s2,j]*g2[s2]
    gradient[j] <- sum(vec1,vec2)
    if (min(group)==0){ ## additional group with direct questioning: group==0
      dq <- group==0
      vec3 <- (y[dq]-e[dq]/(1+e[dq])) *cov[dq,j] 
      gradient[j] <- gradient[j] + sum(vec3)
    }
  }
  # d/d piUQ
  vec1 <-  ( (1-p1)*e+1-p1)*(y/((p1+piUQ-piUQ*p1)*e+piUQ-piUQ*p1)-
                        (1-y)/((1-piUQ-p1+p1*piUQ)*e+1-piUQ+p1*piUQ))
  vec2 <- ( (1-p2)*e+1-p2)*(y/((p2+piUQ-piUQ*p2)*e+piUQ-piUQ*p2)-
                        (1-y)/((1-piUQ-p2+p2*piUQ)*e+1-piUQ+p2*piUQ))
  gradient[m] <- sum(vec1[s1],vec2[s2])
  if (setPiUQ) { gradient[m] <- 0}
  return(gradient)
}


# CDM
RRlog.CDM <- function(x,y,p,start,group,setGamma=FALSE, maxit=1000){
  grad <- rep(NA,ncol(x)+1);
  logLik<- NA;
  coef<- rep(NA,ncol(x));
  iter<- NA;
  hessian<- matrix(NA,ncol=ncol(x),nrow=ncol(x))
  tryCatch(
  {est <- optim(par=start,fn=RRlog.CDM.ll,
              gr=RRlog.CDM.llgrad, 
#               method="L-BFGS-B",
              lower = c(rep(-Inf,ncol(x)),0), 
              upper = c(rep( Inf,ncol(x)),1),
              control=list(fnscale=-1, maxit=maxit),hessian=T,  
              cov=x,y=y,prand=p,group=group,setGamma=setGamma)
 grad <- RRlog.CDM.llgrad(est$par,x,y,p,group,setGamma)
 logLik=est$value;
 coef=est$par;
 iter=est$counts;
 hessian=est$hessian
}
,error = function(e) {})

  # print(grad)
  #print(grad(RRlog.CDM.ll,est$par,cov=x,y=y,prand=p,group=group,setGamma=setGamma))
  res <- list(model="CDM",
              pString=paste0("forced Yes = ",p[1],"/",p[2], " (group 1/2)"),
              coefficients=coef,
              logLik=logLik, param=c(colnames(x),"gamma"),
              gradient=grad,hessian=hessian,iter=iter)
  return(res)
}

# loglikelihoodfunction for CDM: (param[1] is gamma=gg)
RRlog.CDM.ll = function (param,cov,y,prand,group,setGamma){
  p <- prand
  m <- length(param)
  gg <- param[m]  # greek gamma =: gg
  s1 <- group==1
  s2 <- group==2
  if (setGamma) gg <- 0
  
  e <- exp(cov%*%param[-m])
  vec1 <- y*log((1-p[1]*gg)*e+p[1]-p[1]*gg)+
    (1-y)*log( p[1]*gg*e+1-p[1]+p[1]*gg)-log(1+e)
  vec2 <- y*log((1-p[2]*gg)*e+p[2]-p[2]*gg)+
    (1-y)*log( p[2]*gg*e+1-p[2]+p[2]*gg)-log(1+e)
  ll <- sum(vec1[s1],vec2[s2])
  if (min(group)==0){       ## additional group with DQ
    dq <- group==0
    vec3 <- y*log(e)-log(1+e)
    ll <- ll + sum(vec3[dq])
  }
  return(ll)
}

# gradient
RRlog.CDM.llgrad <- function(param,cov,y,prand,group,setGamma){
  p1 <- prand[1]
  p2 <- prand[2]
  gradient <- rep(0,length(param))
  m <- length(param)
  gg <- param[m]  # greek gamma =: gg
  s1 <- group==1
  s2 <- group==2
  e <- exp(cov%*%param[-m])
  if (setGamma) {
    gg <- 0
  } else{# d/d gb
    vec1 <- p1*(e+1)*(-y/((1-p1*gg)*e+p1-p1*gg)+(1-y)/(p1*gg*e+1-p1+p1*gg))
    vec2 <- p2*(e+1)*(-y/((1-p2*gg)*e+p2-p2*gg)+(1-y)/(p2*gg*e+1-p2+p2*gg))
    gradient[m] <- sum(vec1[s1],vec2[s2])
  }
  
  g1 <- y*(1-p1*gg)/((1-p1*gg)*e+p1-p1*gg)+
    (1-y)*(p1*gg)/(p1*gg*e+1-p1+p1*gg)-1/(1+e)
  g2 <- y*(1-p2*gg)/((1-p2*gg)*e+p2-p2*gg)+
    (1-y)*(p2*gg)/(p2*gg*e+1-p2+p2*gg)-1/(1+e)  
  # d/d betai
  for (j in 1:(m-1)){
    vec1 <- e[s1] *cov[s1,j]*g1[s1]
    vec2 <- e[s2]*cov[s2,j]*g2[s2]
    gradient[j] <- sum(vec1,vec2)
    if (min(group)==0){     ## additional group with DQ
      dq <- group==0
      vec3 <- (y[dq]-e[dq]/(1+e[dq])) *cov[dq,j] 
      gradient[j] <- gradient[j] + sum(vec3)
    }
  }
  return(gradient)
}


###############################
###############################
# CDMsym
RRlog.CDMsym <- function(x,y,p,start,group,setGamma=FALSE, maxit=1000){
  grad <- rep(NA,ncol(x)+1);
  logLik<- NA;
  coef<- rep(NA,ncol(x));
  iter<- NA;
  hessian<- matrix(NA,ncol=ncol(x),nrow=ncol(x))
  tryCatch(
{est <- optim(par=start,fn=RRlog.CDMsym.ll,
              #                gr=RRlog.CDMsym.llgrad, 
              method="L-BFGS-B",
              lower = c(rep(-Inf,ncol(x)),0), 
              upper = c(rep( Inf,ncol(x)),1),
              control=list(fnscale=-1, maxit=maxit),hessian=T,  
              cov=x,y=y,prand=p,group=group,setGamma=setGamma)
#  grad <- RRlog.CDM.llgrad(est$par,x,y,p,group,setGamma)
 logLik=est$value;
 coef=est$par;
 iter=est$counts;
 hessian=est$hessian
}
,error = function(e) {})

#   grad <- RRlog.CDMsym.llgrad(est$par,x,y,p,group,setGamma)
#   print(grad)
#   print(grad(RRlog.CDMsym.ll,est$par,cov=x,y=y,prand=p,group=group,setGamma=setGamma))
res <- list(model="CDMsym",
            pString=paste0("forced Yes/No=",p[1],"/",p[2]," (group 1) ; forced Yes/No=",p[3],"/",p[4]," (group 2)"),
            coefficients=coef,
            logLik=logLik, param=c(colnames(x),"gamma"),
            gradient=grad, hessian=hessian,iter=iter)
return(res)
}

# loglikelihoodfunction for CDMsym: (param[1] is gamma=gg)
RRlog.CDMsym.ll = function (param,cov,y,prand,group,setGamma){
  p <- prand
  m <- length(param)
  gg <- param[m]  # greek gamma =: gg
  e <- exp(cov%*%param[-m])
  s1 <- group==1
  s2 <- group==2
  if (setGamma) gg <- 0
  pi <- e/(1+e)
  
  #   ll <- RRuni.CDMsym.ll(c(pi,gg),x,prand,group)
  vec1 <- y* log(pi*(1-p[2])+(1-pi-gg)*p[1]) +
          (1-y)* log(pi*p[2]    +(1-pi-gg)*(1-p[1])+gg)
  vec2 <- y* log(pi*(1-p[4])+(1-pi-gg)*p[3]) +
          (1-y)* log(pi*p[4]    +(1-pi-gg)*(1-p[3])+gg)
  ll <- sum(vec1[s1],vec2[s2])
  if (min(group)==0){ ## additional group with direct questioning: group==0
    dq <- group==0
    vec3 <- y*log(e)-log(1+e)
    ll <- ll + sum(vec3[dq])
  }
  return(ll)
}

# gradient
RRlog.CDMsym.llgrad <- function(param,cov,y,prand,group,setGamma){
  p1 <- prand[1]
  p2 <- prand[2]
  p3 <- prand[3]
  p4 <- prand[4]
  gradient <- rep(0,length(param))
  m <- length(param)
  gg <- param[m]  # greek gamma =: gg
  s1 <- group==1
  s2 <- group==2
  e <- exp(cov%*%param[-m])
  
  if (setGamma){
    gg <- 0
  }else{ # d/d gamma
    vec1 <- -p1*(e+1)*y/( (1-p2-p1*gg)*e-p1*gg)     +(1-y)*(p1*e-1+p1)/( (p1*gg+1+p2)*e-gg+gg*p1)
    vec2 <- -p3*(e+1)*y/( (1-p4-p3*gg)*e-p3*gg)     +(1-y)*(p3*e-1+p3)/( (p3*gg+1+p4)*e-gg+gg*p3)
    gradient[m] <- sum(vec1[s1],vec2[s2])
  }
  
  g1 <- y*(1-p1*gg)/((1-p1*gg)*e+p1-p1*gg)+
    (1-y)*(p1*gg)/(p1*gg*e+1-p1+p1*gg)-1/(1+e)
  g2 <- y*(1-p2*gg)/((1-p2*gg)*e+p2-p2*gg)+
    (1-y)*(p2*gg)/(p2*gg*e+1-p2+p2*gg)-1/(1+e)  
  # d/d betai
  for (j in 1:(m-1)){
    vec1 <- e[s1] *cov[s1,j]*g1[s1]
    vec2 <- e[s2]*cov[s2,j]*g2[s2]
    gradient[j] <- sum(vec1,vec2)
    if (min(group)==0){ ## additional group with direct questioning: group==0
      dq <- group==0
      vec3 <- (y[dq]-e[dq]/(1+e[dq])) *cov[dq,j] 
      gradient[j] <- gradient[j] + sum(vec3)
    }
  }
  
  return(gradient)
}
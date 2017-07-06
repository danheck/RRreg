
# loglikelihood for univariate RR models
RRuni.ll <- function(par, model, response, pp, group, ncat, Kukrep=1){
  y <- by(factor(response, levels = 1:ncat-1), group, table)
  if (Kukrep > 1) y <- by(factor(response, levels = 0:Kukrep), group, table)
  n <- sapply(y, sum)
  
  if(is2group(model)){
    P <- list(getPW(model, pp, group= 1, par2=par[2]),
               getPW(model, pp, group= 2, par2=par[2]))
  }else{
    P <- list(getPW(model, pp, Kukrep = Kukrep))
  }
  
  pi <- c(1-sum(par[1:(ncat-1)]), par[1:(ncat-1)])
  if(model %in% c("CDM", "CDMsym"))
    pi <- c(par, 1-sum(par))
  
  #pi <- pi/sum(pi)
  prob <- lapply(P, function(PPP) PPP %*% pi)
  
  # problems for categorical RR: sum(pi) must be 1 !
  prob.negative <- any(pi<0) | any(sapply(prob, function(ppp) ppp<0))
  if(is.na(prob.negative) || prob.negative)
    return(1e6)
  
  ll <- -sum(mapply(dmultinom, y, n, prob, log=TRUE))

  return(ll)
}
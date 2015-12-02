# 
# #' Beta-Binomial ANOVA for RR Data
# #'
# #' The model assumes that true prevalence rates \eqn{\pi_ij} of individuals within groups follow beta-distributions with separate parameters \eqn{\alpha_j} and \eqn{\beta_j} for each group. 
# #' @param response vector providing the number of affirmative responses to the RR questions
# #' @param group vector of same length indicating group membership
# #' @param n.per.person either a single number (>1) if each respondent gave the same number of responses or a vector of the same length as \code{response} indicating the number of responses per participants
# #' @param data an optional data.frame from which the variables are taken
# #' @param model RR design (only one-group designs allowed)
# #' @param p randomiation probability
# #' @param ci credibility region for the individual prevalence estimates
# #' @param LR.test whether to run likelihood ratio tests
# # @param fit.n number of fitting runs
# #' @details Data must be in the long format (one observation per group/person combination)
# #' 
# #' The group prevalence means and variances are estimated using maximum likelihood and tested for equality using likelihood ratio tests. The following models are included:
# #' \itemize{
# #'  \item{free: }{separate prevalence and variance in each group}
# #'  \item{restrict_pi: }{identical prevalence estimate across groups}
# #'  \item{restrict_var: }{identical variance across groups}
# #'  \item{restrict_both: }{identical prevalence AND variance across groups}
# #  \item{no_var: }{variance=0 in all groups}
# #' }
# #' 
# #' Individual prevalence rates are estimated by empirical Bayes: First, the hyperparameter \eqn{\alpha} and \eqn{\beta} are estimated in each group via maximum likelihood. Next, individual prevalence rates are obtained by using these estimates as hyperpriors. The standard error of these individual estimates are adjusted for the additional uncertainty in estimating \eqn{\alpha} and \eqn{\beta} by a Taylor series approximation.
# #' 
# #' @return a list of the class \code{RRanova} containing:
# #' \itemize{
# #'  \item{\code{estimates}: }{separate \code{prevalence} estimates for all models; \code{parameter} (\eqn{\alpha} and \eqn{\beta})}
# #'  \item{\code{individual}: estimates for unrestricted model}
# #'  \item{\code{LR.test}: }{results of likelihood-ratio test}
# #  \item{\code{data}: }{response, group, and id vector}
# #'  \item{\code{model}: }{information which RR design was used}
# #' }
# #' @examples 
# #' # generate data with 3 groups, 200 participants per group, 10 responses each:
# #' pi_mean <- c(.2, .7, .5)
# #' pi_var <- c(.15,.15,.15)
# #' # transform into beta-parameters:
# #' alpha <- -(pi_mean*(pi_mean^2-pi_mean+pi_var))/pi_var
# #' beta <- ((pi_mean-1)*(pi_mean^2-pi_mean+pi_var))/pi_var
# #' group <- rep(LETTERS[1:3], each=200*10)
# #' id <-  rep(1:(3*200), each=10)
# #' # separate true prevalance for each person:
# #' pi <- c(rbeta(200, alpha[1], beta[1]),
# #'         rbeta(200, alpha[2], beta[2]),
# #'         rbeta(200, alpha[3], beta[3]))
# #' true <- rbinom(3*200*10, 1, rep(pi, each=10))
# #' 
# #' # RR procedure
# #' model <- "FR"
# #' p <- c(.1,.15)
# #' response <- RRgen(trueState=true, model=model, p=p)$response
# #' 
# #' # fit model
# #' mod1 <- RRanova(response, group=rep("all", 3*200*10), id=id, model=model, p=p)
# #' mod1
# #' mod2 <- RRanova(response, group, id=id, model=model, p=p)
# #' mod2
# #' @references Fox, J.-P. (2008). Beta-Binomial ANOVA for Multivariate Randomized Response Data. British Journal of Mathematical and Statistical Psychology, 61(2), 453–470.
# 
# #' @export
# RRanova <- function(response, group, n.per.person, data, model, p, ci=.9,  LR.test = TRUE){
#   
#   ##################### check routines
#   if(!missing(data)){
#     try({
#       data <- as.data.frame(data)
#       response <-  eval(substitute(response),data, parent.frame())
#       group <-  eval(substitute(group),data, parent.frame())
#       n.per.person <-  eval(substitute(n.per.person),data, parent.frame())
#     },silent=T)
#     if(is.null(response))
#       stop(paste0("Variable ",substitute(response)," missing in 'data'."))
#     if(is.null(group))
#       stop(paste0("Variable ",substitute(group)," missing in 'data'."))
#     if(is.null(n.per.person))
#       stop(paste0("Variable ",substitute(n.per.person)," missing in 'data'."))
#     
#   }
#   J <- length(unique(group))
#   N <- length(response)
#   if(length(n.per.person) == 1)
#     n.per.person <- rep(n.per.person, N)
#   
#   if(N != length(group) | length(n.per.person)!=N)
#     stop("Input vectors 'response', 'group', and/or 'n.per.person' do not have the same length!")
#   if(min(response) < 0 | any(response>n.per.person) | any(round(response) != response))
#     stop("The dependent variable 'response' should only contain the values between 0 and 'n.per.person'.")
#   
#   # TODO if()
#   # stop("Participants not nested within groups!")
#   
#   
#   ###################### check feasible RR designs and get linear transformation funktion
#   
#   model <- match.arg(model, modelnames())
#   if(is2group(model) | isContinuous(model) )
#     stop("Only one-group, dichotomous RR models allowed at the moment (see ?RRanova)")
#   
#   # misclassification matrix:
#   pw <- getPW(model, p)[2,]
#   # p(1) = phi1*pi + (1-phi1)*phi2    # in FR design
#   # p(1) = p(1|0)*(1-pi) + p(1|1)*pi  # as misclassification probabilities
#   # p(1) = [p(1|1)-p(1|0)]*pi + p(1|0)
#   phi1 <- pw[2] - pw[1]
#   phi2 <- pw[1] / (1-phi1)
#   # linear transformation of true prevalence to observed frequency: (and inverse)
#   deltaRR <- function(pi) phi1*pi + pw[1]
#   deltaRRinv <- function(p) (p - pw[1])/phi1
#   
#   ###################### 1. estimate hyperpriors alpha[1...J] and beta[1...J] 
#   
#   aggr <- aggregate(response, by=list(group), length)
#   groupnames <- aggr$Group
#   n <-  aggr$x
#   
#   # (b) ML estimate for alpha[j], beta[j]
#   # par = c(alpha_j, beta_j)
#   # y_j = number of 1s per participant  (in group j)
#   # n_j per participant (in group j)
#   loglik <- function(par, y_j, n_j){
#     a <- par[1]
#     b <- par[2]
#     # implement beta-binomial directly:
#     ll <- sum(lgamma(a+y_j) + lgamma(n_j+b-y_j) - lgamma(n_j+a+b)) + 
#       length(n_j)*(lgamma(a+b) - lgamma(a)-lgamma(b))
# #     loglik.pp <- function(yy,nn)  
# #       sum(log(a+0:yy))+sum(log(b+0:(nn-yy-1)))-sum(log(a+b+0:(nn-1)))
# #     ll <- mapply(loglik.pp, y_j, n_j)
#     return(-sum(ll))
#   }
#   
#   individual <- data.frame(group=group, response=response, 
#                            n.per.person=n.per.person, pi=NA, piSE=NA)
#   y.list <- n.list <- vector("list", J)
#   alpha <- runif(J, .7,1.3)
#   beta <- runif(J, .7,1.3)
#   BF_var <- rep(NA, J)
#   parameter <- cbind(alpha=alpha, alpha.SE=NA, beta=beta, beta.SE=NA)
#   rownames(parameter) <- groupnames
#   for(j in 1:J){
#     groupname <- groupnames[j]
#     selG <- group == groupname
#     y.list[[j]] <- y_j <- response[selG]  # aggregate(response[selG], list(id[selG]), sum)$x
#     n.list[[j]] <- n_j <- n.per.person[selG] #aggregate(response[selG], list(id[selG]), length)$x
#     
# #     ym <- mean(y_j)
# #     ys <- var(y_j)
# #     astart <- ym*(y_j*(n_j-ym)-ys)/(n_j*ys-ym*(n_j-ym))
# #     bstart <- (n_j-ym)*(ym*(n_j-ym)-ys)/(n_j*ys-ym*(n_j-ym))
#     # estimate alpha and beta using ML:
# #     oo <- list(value=Inf)
# #     for(i in 1:fit.n){
#       try({
#         oo <- optim(par=c(parameter[j,"alpha"], parameter[j,"beta"]), 
#                         fn=loglik, lower=1e-10, method="L-BFGS-B", hessian=TRUE,
#                         y_j=y.list[[j]], n_j=n.list[[j]])
# #         if(oo.new$value<oo$value)
# #           oo <- oo.new
#       })
#     # }
#     parameter[j, "alpha"] <- alpha[j] <- oo$par[1]
#     parameter[j, "beta"] <- beta[j] <- oo$par[2]
#     try({
#       parSE <- sqrt(diag(solve(oo$hessian)))
#       parameter[j, c("alpha.SE", "beta.SE")] <- parSE
#     })
#     
#     ###################### 2. use estimates to compute posterior for p_ij
#     
#     # selGG <- individual$group==groupname
#     # individual$id[selGG] <- aggregate(response[selG], list(id[selG]), sum)$Group
#     # Eq. 4:
#     individual$pi[selG] <- deltaRRinv((y_j+alpha[j])/(n_j+alpha[j]+beta[j]))
#     # Eq. 5:
#     piVar1 <- (y_j+alpha[j])*(n_j-y_j+beta[j])/(n_j+alpha[j]+beta[j]+1) / (phi1*(n_j+alpha[j]+beta[j]))^2
#     # Eq. 7:
#     # corrected variance via tayler expansion
#     # expectation: E[p | y, a, b] = ( (y+a)/(n+a+b) - pw[1])/phi1
#     # gradient (alpha):   (beta[j]+n_j-y_j)/(phi1*(alpha[j]+beta[j]+n_j)^2)  ### FOR RR EXPECTATION
#     # gradient (beta):    -(alpha[j]+y_j)/(phi1*(alpha[j]+beta[j]+n_j)^2)  ### FOR RR EXPECTATION
#     grad_alpha <- (beta[j]+n_j-y_j)/(alpha[j]+beta[j]+n_j)^2# for beta-binomial w/o RR
#     grad_beta <-  -(alpha[j]+y_j)/(alpha[j]+beta[j]+n_j)^2  # for beta-binomial w/o RR
#     gr <- cbind(grad_alpha, grad_beta)
#     
#     piVar2 <- 0
#     for(c in 1:2){
#       for(d in 1:2){
#         try(piVar2 <- piVar2 + solve(oo$hessian[c,d])*gr[,c]*gr[,d] / phi1^2)
#       }
#     }
#     # print(cbind(piVar1, piVar2))
#     individual$piSE[selG] <- sqrt(piVar1 + piVar2)
#     
#     # Eq. 9  : df2=alpha-1 <0 in F-distribution!! 
#     nu <- (1-ci)/2
#     # posterior estimates for each person:
#     as <- alpha[j] + y_j
#     bs <- beta[j]-y_j+n_j
#     # lower and upper probabilities on RR-scrambled scale
#     
#     for(i in 1:n[j]){
#       if(as[i]<=1){
#         individual$pi.lower[selG][i] <-0
#       }else{
#         fL <- function(low) pbeta(deltaRR(low),as[i], bs[i])- nu
#         try(individual$pi.lower[selG][i] <- uniroot(fL, c(-1,2))$root, silent=T)
#       }
#       if(bs[i]<=1){
#         individual$pi.upper[selG][i] <-1
#       }else{
#         fU <- function(up) pbeta(deltaRR(up),as[i], bs[i])- (1-nu)
#         try(individual$pi.upper[selG][i] <- uniroot(fU, c(-1,2))$root, silent=T)
#       }
#     }
#     
#     # approximation by fox (does not work)
# #     suppressWarnings({
# #       pL <- ifelse(as<=1, 0, 1/(1+(bs+1)*(as-1)* pf(1-nu/2, 2*(bs+1), 2*(as-1))) )
# #       pU <- ifelse(bs<=1, 1, (as/bs*pf(1-nu/2, 2*as, 2*bs))/
# #                      (1+as/bs*pf(1-nu/2, 2*as, 2*bs)))
# #     })
# #     individual$pi.lower[selG] <-deltaRRinv(pL)
# #     individual$pi.upper[selG] <-deltaRRinv(pU)
#     ff <- Vectorize(function(ksi) exp(sum(dbinom(y_j,n_j,ksi, log=T))),"ksi")
#     p.var.num <- integrate(ff,0,1)$value
#     ff <- Vectorize(function(ksi) exp(sum(dbinom(y_j,n_j,ksi, log=T))),"ksi")
#     dhsig <- 1/(alpha[j]+beta[j])*sqrt(pi/2)
#     dhalfnormal <- function(x, sigma) 
#       exp(-x^2/(2*sigma^2))*sqrt(2/pi)/sigma
#     dbetabin <- Vectorize(function(ksi, w){
#       a <- ksi/w
#       b <- (1-w*a)/w
#       exp(sum(lchoose(n_j,y_j)+lgamma(a+y_j) + lgamma(n_j+b-y_j) - lgamma(n_j+a+b)) + 
#             length(n_j)*(lgamma(a+b) - lgamma(a)-lgamma(b)))
#     }, c("ksi", "w"))
#     ff.full <- Vectorize(function(ksi){
#       ftmp <- Vectorize(function(w) dhalfnormal(w, dhsig)*dbetabin(ksi, w),"w")
#       integrate(ftmp,0,10)$value
#     }, "ksi")
#     p.var.den <- integrate(ff.full,0,1)$value
#     BF_var[j] <- p.var.num / p.var.den
#     # print(p.var.num)
#   }
#   individual$pi.lower <- ifelse( is.na(individual$pi.lower) |  individual$pi.lower<0 |
#                                   response==0, 0, individual$pi.lower)
#   individual$pi.upper <- ifelse(is.na(individual$pi.upper) |individual$pi.upper>1 |
#                                 response == n.per.person, 1, individual$pi.upper)
#   
#   ###################### 3.  (bayes factor: no detailed procedure in paper)
#   # compute likelihood ratio test instead
#   
#   # par = c(mean, variance) for a single group  --- no RR!!
#   loglik2 <- function(par, y_j, n_j){
#     a <-  -(par[1]*(par[1]^2-par[1]+par[2]))/par[2]
#     b <-  ((par[1]-1)*(par[1]^2-par[1]+par[2]))/par[2]
#     # implement beta-binomial directly:
#     ll <- sum(lgamma(a+y_j) + lgamma(n_j+b-y_j) - lgamma(n_j+a+b)) + 
#       length(n_j)*(lgamma(a+b) - lgamma(a)-lgamma(b))
#     return(-ll)
#   }
#   
#   # par: vector of parameter values (means/variances depending on model)
#   # y.list, n.list = lists of length J with vectors containing response sums and number of responses per participant
#   totalLL <- function(par, y.list, n.list, restriction){
#     J <- length(y.list)
#     par <- switch(restriction,
#                   free = par,
#                   restrict_prev = c(rep(par[1], J), par[1:J+1]),
#                   restrict_var = c(par[1:J], rep(par[J+1], J)),
#                   restrict_both= rep(par, each=J),
#                   no_var=c(par, rep(0, J))  # problem: teilen durch null wenn beta/alpha berechnet werden!
#     )
#     ll <- 0
#     for(j in 1:J){
#       ll <- ll +loglik2(par[c(j, J+j)], y.list[[j]], n.list[[j]])
#     }
#     return(ll)  # already negative
#   }
#   
#   # get standard errors for mean and variance estimates:
#   mm <- alpha/(alpha+beta)
#   vv <- alpha*beta/((alpha+beta)^2*(alpha+beta+1))
#   
#   models <- c("free", "restrict_prev", "restrict_var", "restrict_both") # , "no_var")
#   n.models <- ifelse(J == 1, 1, 4)
#   est.pi <- vector("list", n.models)
#   names(est.pi) <- models[1:n.models]
#   fit <- data.frame(model=models[1:n.models], G2=NA, npar=NA, delta.df=NA, delta.G2=NA, p=NA)
#   
#   for(mod in 1:n.models){
#     # oo2 <- list(value=Inf)
#     restriction <- models[mod]
#     # get starting values
#     start <- switch(restriction,
#                     free = c(mm, vv),
#                     restrict_prev = c(mean(mm), vv),
#                     restrict_var = c(mm,mean(vv)),
#                     restrict_both= c(mean(mm), mean(vv)),
#                     no_var=mm)
#     n.pi <- switch(restriction,
#                    free = J,
#                    restrict_prev = 1,
#                    restrict_var = J,
#                    restrict_both= 1,
#                    no_var=J)
#     # fit model
#     try({
#       #       for(i in 1:fit.n){
#       #         if(i>1)
#       #           start <- runif(length(start), .3,.7)
#       oo2 <- optim(par=start, fn=totalLL, method="L-BFGS-B", hessian=TRUE,
#                    lower=c(rep(1e-6, n.pi),rep(1e-6, length(start)-n.pi)), 
#                    upper=c(rep(1-1e-6, n.pi),rep(.25-1e-6, length(start)-n.pi)), 
#                    y.list=y.list, n.list=n.list, restriction=restriction)
#       #       if(oo.new$value<oo2$value)
#       #         oo2 <- oo.new
#       
#       # transfer estiamtes back into table
#       par <- switch(restriction,
#                     free = oo2$par,
#                     restrict_prev = c(rep(oo2$par[1], J), oo2$par[1:J+1]),
#                     restrict_var = c(oo2$par[1:J], rep(oo2$par[J+1], J)),
#                     restrict_both = rep(oo2$par, each=J),
#                     no_var = c(oo2$par, rep(0, J)))
#       try({
#         tmp <- sqrt(diag(solve(oo2$hessian)))
#       se <- switch(restriction,
#                    free = tmp,
#                    restrict_prev = c(rep(tmp[1], J), tmp[1:J+1]),
#                    restrict_var = c(tmp[1:J], rep(tmp[J+1], J)),
#                    restrict_both= rep(tmp, each=J),
#                    no_var=c(tmp, rep(0, J)))
#       est.pi[[mod]] <- cbind(mean=deltaRRinv(par[1:J]), mean.SE=se[1:J], 
#                              variance=par[J+1:J], variance.SE=se[1:J+J
#                                                                                                  ])
#       rownames(est.pi[[mod]]) <- groupnames
#       })
#       fit[mod,"G2"] <- 2*oo2$value
#       fit[mod,"npar"] <- length(start)
#     })
#   }
#   fit$delta.G2 <- fit$G2 - fit$G2[1]
#   fit$delta.df <- fit$npar[1] - fit$npar
#   fit$p <- with(fit, pchisq(delta.G2, npar, lower=FALSE))
#   
#   res <- list(estimates=list(prevalence=est.pi, parameter=parameter),
#               BF_var=BF_var,
#               LR.test = fit, individual=individual,
#               model=paste0(model, " with p=", paste0(round(p, 3), collapse=", ")))
#   class(res) <- "RRanova"
#   return(res)
# }
# 
# 
# #' @export
# print.RRanova <- function(x,...){
#   cat("RR design:",x$model, ":\n")
#   
#   cat("\nEstimated means and variance of prevalence rates (unrestricted model):\n")
#   print(round(x$estimates$prevalence[[1]],3))
#   
#   if(nrow(x$LR.test) > 1){
#     cat("\nLikelihood ratio tests:\n")
#     LR.test <- x$LR.test
#     LR.test$p <- round(LR.test$p, 4)
#     LR.test$p[LR.test$p == 0] <- "<.0001"
#     print(LR.test)
#   }
# }

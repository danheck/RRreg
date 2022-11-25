#' Mixed Effects Logistic Regression for RR Data
#'
#' Uses the package \code{\link{lme4}} to fit a generalized linear mixed model
#' (GLMM) with an adjusted link funciton.
#'
#' @param formula two-sided formula including random and fixed effects (see
#'   below or \code{\link{glmer}} for details)
#' @param data an optional data frame with variables named in formula
#' @param model type of RR design. Only 1-group RR designs are supported at the
#'   moment (i.e., \code{"Warner"}, \code{"FR"}, \code{"UQTknown"},
#'   \code{"Crosswise"}, \code{"Triangular"}, \code{"Kuk"}, \code{"Mangat"},
#'   \code{"custom"}). See \code{\link{RRuni}} or \code{vignette(RRreg)} for
#'   details.
#' @param p randomization probability
#' @param const the RR link function is not defined for small and/or large
#'   probabilities (the boundaries depend on \code{model} and \code{p}). To
#'   increase robustness of the estimation, these probabilities smaller and
#'   larger than the respective boundaries (plus/minus a constant defined via
#'   \code{const}) are smoothed by separate logit-link functions.
#' @param adjust_control whether to adjust the control arguments for
#'   \code{glmer}, which might help in case of convergence issues for some
#'   models. \code{\link[lme4]{lmerControl}} settings are changed to
#'   \code{nAGQ0initStep = FALSE} and \code{optimizer = "bobyqa"}.
#' @param ... further arguments passed to \code{\link{glmer}}
#'
#' @details Some examples for formula: 
#' \itemize{ 
#'   \item{random intercept: }
#'        {\code{response ~ covariate + (1 | group)}} 
#'   \item{random slope: }
#'        {\code{response ~ covariate + (0 + covariate | group)}} 
#'   \item{both random slope and intercept: }
#'        {\code{response ~ covariate +(covariate | group)}}
#'   \item{level-2 predictor (must have constant values within groups!):}
#'        {\code{response ~ lev2 + (1|group)}} 
#' }
#'
#' Note that parameter estimation will be unstable and might fail if the
#' observed responses are not in line with the model. For instance, a
#' Forced-Response model (\code{model="FR"}) with \code{p=c(0,1/4)} requires
#' that expected probabilities for responses are in the interval [.25,1.00]. If
#' the observed proportion of responses is very low (<<.25), intercepts will be
#' estimated to be very small (<<0) and/or parameter estimation might fail. See
#' \code{\link[lme4]{glmer}} for setting better starting values and
#' \code{\link[lme4]{lmerControl}} for further options to increase stability.
#'
#' @return an object of class \code{glmerMod}
#'
#' @references  van den Hout, A., van der Heijden, P. G., & Gilchrist, R.
#'   (2007). The Logistic Regression Model with Response Variables Subject to
#'   Randomized Response. Computational Statistics & Data Analysis, 51,
#'   6060–6069.
#'
#' @examples
#' # generate data with a level-1 predictor
#' set.seed(1234)
#' d <- data.frame(
#'   group = factor(rep(LETTERS[1:20], each = 50)),
#'   cov = rnorm(20 * 50)
#' )
#' # generate dependent data based on logistic model (random intercept):
#' d$true <- simulate(~ cov + (1 | group),
#'   newdata = d,
#'   family = binomial(link = "logit"),
#'   newparams = list(
#'     beta = c("(Intercept)" = -.5, cov = 1),
#'     theta = c("group.(Intercept)" = .8)
#'   )
#' )[[1]]
#' # scramble responses using RR:
#' model <- "FR"
#' p <- c(true0 = .1, true1 = .2)
#' d$resp <- RRgen(model = "FR", p = p, trueState = d$true)$response
#' # fit model:
#' mod <- RRmixed(resp ~ cov + (1 | group), data = d, model = "FR", p = p)
#' summary(mod)
#' @export
RRmixed <- function(formula, data, model, p, const = .0001,
                    adjust_control = FALSE, ...) {
  ######### check if model is allowed
  model <- match.arg(model, modelnames())
  if (is2group(model) | isContinuous(model)) {
    stop("Only one-group, dichotomous RR models allowed at the moment (see ?RRmixed)")
  }

  ######### get link function
  p <- getPW(model, p)[2, ]

  ######### fit model using lme4
  args <- list(...)
  if (is.null(args$mustart)) {
    args$mustart <- runif(nrow(data), min(p) + const, max(p) - const)
  }
  if (adjust_control) {
    args$control <- lme4::glmerControl(
      nAGQ0initStep = FALSE,
      optimizer = "bobyqa"
    )
  }
  if (is.null(args$start)) {
    glmer_formula <- lme4::glFormula(formula, data)
    args$start <- list(fixef = runif(ncol(glmer_formula$X), -.3, .3))
  }

  input <- c(args, list(
    formula = formula, data = data,
    family = binomial(link = RRloglink(p = p, const = const))
  ))
  mod <- do.call("glmer", args = input)
  mod@call$data <- substitute(data)
  mod@call$control <- list(...)$control
  mod
}




################ LINK FUNCTION #############################

# p = c(0,1)  = probability to respond 1 given a true state of 0 or 1 respectively
# gegeben: p = c(p(1|0), p(1|1))
# gesucht: p(resp==1) = c+d*pi
# p(1) = p(1|0)*(1-pi)     + p(1|1)*pi
# p(1) = p(1|0)- p(1|0)*pi + p(1|1)*pi
# p(1) = p(1|0) + (p(1|1)-p(1|0))*pi
#' @importFrom stats qlogis
RRloglink <- function(p = c(0, 1), const = .0001) {
  c <- p[1]
  d <- p[2] - p[1]
  linkfun <- function(mu) {
    suppressWarnings(eta <- log((mu - c) / (c + d - mu))) # stability if mu<p[1] or mu>p[2]

    below <- mu < min(p) + const
    above <- mu > max(p) - const
    slope <- ifelse(p[1] < p[2], 1, -1)
    eta[below] <- slope * (2 * log(const) + qlogis(mu[below] / (min(p) + 2 * const)))
    eta[above] <- slope * (-2 * log(const) + qlogis((mu[above] - max(p) + 2 * const) / (1 - max(p) + 2 * const)))
    eta
  }
  linkinv <- function(eta) c + d / (1 + exp(-eta))
  mu.eta <- function(eta) (exp(eta) * d) / (exp(eta) + 1)^2
  valideta <- function(eta) TRUE
  link <- paste0("logRR(", paste0(p, collapse = ","), ")")
  simulate <- structure(
    list(
      linkfun = linkfun, linkinv = linkinv,
      mu.eta = mu.eta, valideta = valideta, name = link
    ),
    class = "link-glm"
  )
}

# # Basic checks for link function:
# vv <- RRloglink(p=c(.1,.9))
# ## check invertibility
# xx <- seq(-30,30, length=100)
# plot(xx, sapply(xx, function(x) vv$linkfun(vv$linkinv(x))), type="l")
# abline(0,1, col="red")
#
# library(numDeriv)
# all.equal(grad(vv$linkinv,1),vv$mu.eta(1))  ## check derivative
# plot(xx, sapply(xx, vv$mu.eta), type="l")
# lines(xx, sapply(xx, function(x) grad(vv$linkinv,x)), col="red")
#
# # stability checks:
# curve(RRloglink(p = c(.2, .9))$linkfun(x))
# p <- getPW("Warner", p = .9)
# p <- getPW("UQTknown", p = c(.1,.8), par2 = .1)
# # c(p[1], p[2]-p[1])
# par(mfrow=c(1,2))
# curve(RRloglink(p = c(.1,.6))$linkfun(x))
# abline(v=c(.1,.6))
# curve(RRloglink(p = c(.6,.11))$linkfun(x))
# abline(v=c(.1,.6))

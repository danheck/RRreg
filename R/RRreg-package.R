#' Correlation and Regression Analyses for Randomized Response Designs

#' @description{
#' \if{html}{\figure{RRreg.png}{options: width=150 alt ="RRreg logo" style='float: right'}}
#' \if{latex}{\figure{RRreg.png}{options: width=0.5in}}
#' Univariate and multivariate methods for randomized response (RR) survey
#' designs (e.g., Warner, 1965). Univariate estimates of true proportions can be
#' obtained using \code{\link{RRuni}}. RR variables can be used in multivariate
#' analyses for correlations (\code{\link{RRcor}}), as dependent variable in a
#' logistic regression (\code{\link{RRlog}}), or as predictors in a linear
#' regression (\code{\link{RRlin}}). The function \code{\link{RRgen}} generates
#' single RR data sets, whereas \code{\link{RRsimu}} generates and analyzes RR
#' data repeatedly for simulation and bootstrap purposes. An overview of the
#' available RR designs and examples can be found in the package vignette by
#' \code{vignette('RRreg')}.
#' }
#' 
#' @details
#' In case of issues or questions, please refer to the GitHub repository: 
#' \url{https://github.com/danheck/RRreg}
#' 
#' An introduction with examples is available via \code{vignette("RRreg")} or 
#' at the website: \url{https://www.dwheck.de/vignettes/RRreg.html}
#' 
#' @author Daniel W. Heck \email{daniel.heck@@uni-marburg.de}
#'
#' @section Citation:
#' If you use \code{RRreg} in publications, please cite the package as follows:
#'
#' Heck, D. W., & Moshagen, M. (2018).
#' RRreg: An R package for correlation and regression analyses of randomized response data.
#' \emph{Journal of Statistical Software. 85 (2)}, 1-29. \doi{10.18637/jss.v085.i02}
#'
#' @references
#' Warner, S. L. (1965). Randomized response: A survey technique for eliminating
#' evasive answer bias. \emph{Journal of the American Statistical Association, 60}, 63–69.
#' 
#' @name RRreg-package
#' @docType package
#' @keywords package
#' 
#' @import stats
#' @import graphics
#' @importFrom grDevices adjustcolor
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster makeCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom lme4 glmer
"_PACKAGE"


#' Minaret Data
#'
#' Data by Radukic and Musch (2014)
#'
#'
#' The following variables are included:
#' \itemize{
#'   \item \code{age} in years
#'   \item \code{leftRight} political left-right orientation on a scale from -5 to 5
#'   \item \code{rrt} response to RR question (SLD with randomization probabilities 
#'      \code{p=c(2/12,10/12)})
#'   \item \code{condition} group membership in SLD (either randomization probability 
#'      \code{p[1]} or \code{p[2]})
#'   \item \code{RRdesign} whether the respondent answered to the RR question 
#'      (RRdesign=1) or to the direct question (RRdesign=-1)
#'   \item \code{leftRight.c} zero-centered political left-right orientation
#'   \item \code{age.c} zero-centered age
#' }
#'
#' @docType data
#' @keywords datasets
#' @name minarets
#' @usage data(minarets)
# @format A data frame with 1621 rows and 6 variables
"minarets"

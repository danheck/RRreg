[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/RRreg)](http://cran.r-project.org/package=RRreg)
[![Build Status](https://travis-ci.org/danheck/RRreg.svg?branch=master)](https://travis-ci.org/danheck/RRreg)
[![Licence](https://img.shields.io/badge/licence-GPL--2-green.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
[![monthly downloads](http://cranlogs.r-pkg.org/badges/RRreg)](http://cranlogs.r-pkg.org/badges/RRreg)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/RRreg)](http://cranlogs.r-pkg.org/badges/grand-total/RRreg)
<!--[![Research software impact](http://depsy.org/api/package/cran/RRreg/badge.svg)](http://depsy.org/package/r/RRreg)-->

RRreg
=====

The package RRreg provides univariate and multivariate methods for randomized response (RR) survey designs (e.g., Warner, 1965). Univariate estimates of true proportions can be obtained using RRuni. RR variables can be used in multivariate analyses for correlations (RRcor), as dependent variable in a logistic regression (RRlog) and as predictors in a linear regression (RRlin). The function RRgen generates single RR data sets, whereas RRsimu generates and analyzes RR data repeatedly for simulation and bootstrap purposes.

### Links

The package can be downloaded from CRAN by typing `install.packages("RRreg")` in an active R session.

To get the most recent version, `RRreg` can also directly be installed from GitHub by:
```
# install.packages("devtools")
devtools::install_github("danheck/RRreg")
```

The manual is available within R by typing `vignette('RRreg')`.

### Citation

If you use `RRreg` in publications, please cite the package as follows:

Heck, D. W., & Moshagen, M. (in press). 
RRreg: An R package for correlation and regression analyses of randomized response data. 
*Journal of Statistical Software.*

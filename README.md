[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/RRreg)](http://cran.r-project.org/package=RRreg)
[![Workflow](https://github.com/danheck/RRreg/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/danheck/RRreg/actions/workflows/check-standard.yaml)
[![Licence](https://img.shields.io/badge/licence-GPL--2-green.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
[![monthly downloads](http://cranlogs.r-pkg.org/badges/RRreg)](http://cranlogs.r-pkg.org/badges/RRreg)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/RRreg)](http://cranlogs.r-pkg.org/badges/grand-total/RRreg)
<!--[![Research software impact](http://depsy.org/api/package/cran/RRreg/badge.svg)](http://depsy.org/package/r/RRreg)-->

RRreg
=====

<img src="man/figures/RRreg.png" width="200" style="float: right">
The package RRreg provides univariate and multivariate methods for randomized response (RR) survey designs (e.g., Warner, 1965). Univariate estimates of true proportions can be obtained using RRuni. RR variables can be used in multivariate analyses for correlations (RRcor), as dependent variable in a logistic regression (RRlog) and as predictors in a linear regression (RRlin). The function RRgen generates single RR data sets, whereas RRsimu generates and analyzes RR data repeatedly for simulation and bootstrap purposes.

### Links

The package can be downloaded from CRAN by typing `install.packages("RRreg")` in an active R session.

To get the most recent version, `RRreg` can also directly be installed from GitHub by:
```
# install.packages("devtools")
devtools::install_github("danheck/RRreg")
```

The vignette is available within R by typing `vignette('RRreg')` or at https://dwheck.de/vignettes/RRreg.html.

### Citation

If you use `RRreg` in publications, please cite the package as follows:

- Heck, D. W., & Moshagen, M. (2018). RRreg: An R package for correlation and regression analyses of randomized response data. *Journal of Statistical Software, 85* (2), 1–29. https://doi.org/10.18637/jss.v085.i02


<!-- README.md is generated from README.Rmd. Please edit that file -->

# statgenMPP

[![R-CMD-check](https://github.com/Biometris/statgenMPP/workflows/R-CMD-check/badge.svg)](https://github.com/Biometris/statgenMPP/actions?workflow=R-CMD-check)
[![codecov](https://codecov.io/gh/Biometris/statgenMPP/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Biometris/statgenMPP)

The statgenMPP package is developed as an easy-to-use package for QTL
mapping in multi-parent populations. The package has many ways of
visualizing inputs and results. First Identity By Descent (IBD)
probabilities are computed using Hidden Markov Models. These
probabilities are then used in a mixed model approach for QTL Mapping as
described in Li et al. .

## Installation

-   Install from CRAN:

``` r
install.packages("statgenMPP")
```

-   Install latest development version from GitHub (requires
    [remotes](https://github.com/r-lib/remotes) package):

``` r
remotes::install_github("Biometris/statgenMPP", ref = "develop", dependencies = TRUE)
```

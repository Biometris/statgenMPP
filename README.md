
<!-- README.md is generated from README.Rmd. Please edit that file -->

# statgenMPP <img src="man/figures/logo.png" align="right" height="139" alt="" />

[![](https://www.r-pkg.org/badges/version/statgenMPP)](https://www.r-pkg.org/pkg/statgenMPP)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/statgenMPP)](https://www.r-pkg.org/pkg/statgenMPP)
[![R-CMD-check](https://github.com/Biometris/statgenMPP/workflows/R-CMD-check/badge.svg)](https://github.com/Biometris/statgenMPP/actions?workflow=R-CMD-check)
[![codecov](https://codecov.io/gh/Biometris/statgenMPP/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Biometris/statgenMPP)

The statgenMPP package is developed as an easy-to-use package for QTL
mapping in biparental and multi-parent populations. The package has many
ways of visualizing inputs and results. First Identity By Descent (IBD)
probabilities are computed using Hidden Markov Models. These
probabilities are then used in a mixed model approach for QTL Mapping as
described in [Li et
al. 2021](https://link.springer.com/article/10.1007/s00122-021-03919-7).

- Install from CRAN:

``` r
install.packages("statgenMPP")
```

- Install latest development version from GitHub (requires
  [remotes](https://github.com/r-lib/remotes) package):

``` r
remotes::install_github("Biometris/statgenMPP", ref = "develop", dependencies = TRUE)
```

## Example

Here we give a simple example, using a single biparental population. The
example contains simulated data for one F4DH population. The population
type F4DH is a cross between two parents, A and C, followed by 3
generations of selfings, followed by a DH generation, see
[statgenIBD](https://biometris.github.io/statgenIBD/) for details.

First we load the marker data and phenotypic data and calculate the IBDs
using the `calcIBDMPP` function:

``` r
library(statgenMPP)
#> Loading required package: statgenGWAS

markerFiles <- system.file("extdata/multipop", "AxC.txt", 
                             package = "statgenMPP")
mapFile <- system.file("extdata/multipop", "mapfile.txt",
                       package = "statgenMPP")
phenoDat <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                                   package = "statgenMPP"))

ACMPP <- calcIBDMPP(crossNames = c("AxC"),  
                     markerFiles = markerFiles,
                     pheno = phenoDat,
                     popType = "F4DH",
                     mapFile = mapFile,
                     evalDist = 5)
```

The population has the following simple structure, for more complicated
examples see the vignette.

``` r
plot(ACMPP, plotType = "pedigree")
```

<img src="man/figures/README-plotPACMPP-1.png" width="100%" />

The next step is to select QTLs using `selQTLMPP`:

``` r
ACMQM <- selQTLMPP(MPPobj = ACMPP, trait = "yield", threshold = 3)
```

The QTL-mapping profile and parental effects:

``` r
plot(ACMQM, plotType = "QTLProfileExt")
```

<img src="man/figures/README-plotQPEACMQM-1.png" width="100%" />

A summary of the QTL-analyis gives a short overview containing the total
number of markers and the number of QTLs found. Also for all QTL their
position on the chromosome is shown as well as the nearest marker on the
original map, the explained variance and the effects and the standard
errors of all parents:

``` r
## Print summary
summary(ACMQM)
#> Trait analysed: yield 
#> 
#> Data are available for 95 markers.
#> Threshold: 3 
#> 
#> Number of QTLs: 2 
#> 
#>   evalPos chr pos mrkNear minlog10p varExpl  eff_A  eff_C se_eff_A se_eff_C
#>  EXT_1_20   1  20    M1_3      9.25   0.280 -0.921  0.921    0.125    0.125
#>  EXT_3_65   3  65    M3_7     12.50   0.513  1.251 -1.251    0.138    0.138
```

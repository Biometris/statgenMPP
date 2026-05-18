# Multi round genome scans for QTL detection

Multi round genome scans for QTL detection.  
  
Several rounds of QTL detection are performed. First a model is fitted
without cofactors. If for at least one marker the \\-log10(p)\\ value is
above the threshold the marker with the lowest p-Value is added as
cofactor in the next round of QTL detection. This process continues
until there are no new markers with a \\-log10(p)\\ value above the
threshold or until the maximum number of cofactors is reached.

## Usage

``` r
selQTLMPP(
  MPPobj,
  trait = NULL,
  QTLwindow = 10,
  threshold = 3,
  maxCofactors = NULL,
  K = NULL,
  computeKin = FALSE,
  parallel = FALSE,
  verbose = FALSE
)
```

## Arguments

- MPPobj:

  An object of class gDataMPP, typically the output of either
  [`calcIBDMPP`](calcIBDmpp.md) or [`readRABBITMPP`](readRABBITMPP.md).

- trait:

  A character string indicating the trait QTL mapping is done for.

- QTLwindow:

  A numerical value indicating the window around a QTL that is
  considered as part of that QTL.

- threshold:

  A numerical value indicating the threshold for the \\-log10p\\ value
  of a marker to be considered a QTL.

- maxCofactors:

  A numerical value, the maximum number of cofactors to include in the
  model. If `NULL` cofactors are added until no new cofactors are found.

- K:

  A list of chromosome specific kinship matrices. If `NULL` and
  `computeKin = FALSE` no kinship matrix is included in the models.

- computeKin:

  Should chromosome specific kinship matrices be computed?

- parallel:

  Should the computation of variance components be done in parallel?
  This requires a parallel back-end to be registered. See examples.

- verbose:

  Should progress and intermediate plots be output?

## Value

An object of class `QTLMPP`

## Details

By default only family specific effects and residual variances and no
kinship relations are included in the model. It is possible to include
kinship relations by either specifying `computeKin = TRUE`. When doing
so the kinship matrix is computed by averaging \\Z Z^t\\ over all
markers, where \\Z\\ is the genotype x parents matrix for the marker. It
is also possible to specify a list of precomputed chromosome specific
kinship matrices in `K`. Note that adding a kinship matrix to the model
increases the computation time a lot, especially for large populations.

## See also

[`kinshipIBD`](kinshipIBD.md)

## Examples

``` r
## Read phenotypic data.
pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                               package = "statgenMPP"))
## Rename first column to genotype.
colnames(pheno)[1] <- "genotype"


## Compute IBD probabilities for simulated population - AxB, AxC.
ABC <- calcIBDMPP(crossNames = c("AxB", "AxC"),
                  markerFiles = c(system.file("extdata/multipop", "AxB.txt",
                                              package = "statgenMPP"),
                                  system.file("extdata/multipop", "AxC.txt",
                                              package = "statgenMPP")),
                  pheno = pheno,
                  popType = "F4DH",
                  mapFile = system.file("extdata/multipop", "mapfile.txt",
                                        package = "statgenMPP"),
                  evalDist = 5)

## Single-QTL Mapping.
ABC_SQM <- selQTLMPP(ABC, trait = "yield", maxCofactors = 0)

## Multi-QTL Mapping.
if (FALSE) { # \dontrun{
## Register parallel back-end with 2 cores.
doParallel::registerDoParallel(cores = 2)

## Run multi-QTL mapping.
ABC_MQM <- selQTLMPP(ABC, trait = "yield", parallel = TRUE)

## Run multi-QTL mapping - include kinship matrix.
ABC_MQM_kin <- selQTLMPP(ABC, trait = "yield", parallel = TRUE,
                        computeKin = TRUE)
} # }
```

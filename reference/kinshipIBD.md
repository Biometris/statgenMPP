# Compute kinship matrix for IBD probabilities

Compute a kinship matrix or a list of chromosome specific kinship
matrices for a 3D array of IBD probabilities. The kinship matrix is
computed by averaging \\Z Z^t\\ over all markers, where \\Z\\ is the
genotype x parents matrix for the marker. If `chrSpecific = TRUE`
chromosome specific kinship matrices are computed for each chromosome
based only on the markers on all other chromosomes.

## Usage

``` r
kinshipIBD(markers, map = NULL, chrSpecific = TRUE)
```

## Arguments

- markers:

  An n x m x p array with IBD probabilities with genotypes in the rows
  (n), markers in the columns (m), and parents in the 3rd dimension (p).

- map:

  A data.frame with columns `chr` for chromosome and `pos` for position.
  Positions should be in centimorgan (cM). They should not be cumulative
  over the chromosomes. Other columns are ignored. Marker names should
  be in the row names. These should match the marker names in the input
  file. Only required if `chrSpecific = TRUE`.

- chrSpecific:

  Should chromosome specific kinship matrices be computed?

## Value

A kinship matrix or a list of chromosome specific kinship matrices.

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

## Compute chromosome specific kinship matrices.
KChrSpec <- kinshipIBD(markers = ABC$markers, map = ABC$map)
```

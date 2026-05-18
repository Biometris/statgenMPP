# Create an object of class gDataMPP

Function for creating an object of class gDataMPP from an object of
class IBDprob (computed or imported using statgenIBD) and phenotypic
data.

## Usage

``` r
createGDataMPP(IBDprob, pheno)
```

## Arguments

- IBDprob:

  An object of class `IBDprob`.

- pheno:

  A data frame with at least columns "genotype" for the "genotype" and
  one or more numerical columns containing phenotypic information. A
  column "cross" can be used for indicating the cross the genotype comes
  from. This column is ignored if the `IBDprob` has an attribute
  genoCross containing this information (the default behaviour).

## Value

An object of class `gDataMPP`

## Examples

``` r
## Read phenotypic data.
pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                               package = "statgenMPP"))
## Rename first column to genotype.
colnames(pheno)[1] <- "genotype"

## Compute IBD probabilities for simulated population using statgenIBD - AxB
AB <- statgenIBD::calcIBD(popType = "F4DH",
                         markerFile = system.file("extdata/multipop", "AxB.txt",
                                                 package = "statgenMPP"),
                         mapFile = system.file("extdata/multipop", "mapfile.txt",
                                              package = "statgenMPP"),
                         evalDist = 5)

## Create an object of class gDataMPP from IBD probabilities and
## phenotypic data.
ABMPP <- createGDataMPP(IBDprob = AB, pheno = pheno)
```

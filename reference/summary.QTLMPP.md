# Summary function for the class `QTLMPP`

Gives a summary for an object of S3 class `QTLMPP`.

## Usage

``` r
# S3 method for class 'QTLMPP'
summary(object, ...)
```

## Arguments

- object:

  An object of class `QTLMPP`.

- ...:

  Not used.

## Examples

``` r
if (FALSE) { # \dontrun{
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

## Multi-QTL Mapping.
ABC_MQM <- selQTLMPP(ABC, trait = "yield")

## Print summary.
summary(ABC_MQM)
} # }
```

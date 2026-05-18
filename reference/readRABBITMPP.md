# Read IBD probabilities

Read a file with IBD probabilities computed by the RABBIT software
package. It is possible to additionally read the pedigree file that is
also used by RABBIT. Reading this file allows for plotting the pedigree.
Phenotypic data can be added from a data.frame.

## Usage

``` r
readRABBITMPP(infile, pedFile = NULL, pheno = NULL)
```

## Arguments

- infile:

  A character string, a link to a .csv file with IBD probabilities.

- pedFile:

  A character string, a link to a .csv file with pedigree information as
  used by RABBIT as input.

- pheno:

  A data frame with at least columns "genotype" for the "genotype" and
  one or more numerical columns containing phenotypic information. A
  column "cross" can be used for indicating the cross the genotype comes
  from.

## Value

An object of class `gDataMPP` with map and markers corresponding to the
imported information in the imported .csv file.

## References

Zheng, Chaozhi, Martin P Boer, and Fred A Van Eeuwijk. “Recursive
Algorithms for Modeling Genomic Ancestral Origins in a Fixed Pedigree.”
G3 Genes\|Genomes\|Genetics 8 (10): 3231–45.
https://doi.org/10.1534/G3.118.200340.

## Examples

``` r
if (FALSE) { # \dontrun{
## Read RABBIT data for barley.
genoFile <- system.file("extdata/barley", "barley_magicReconstruct.zip",
                       package = "statgenMPP")
barleyMPMPP <- readRABBITMPP(unzip(genoFile, exdir = tempdir()))
} # }
```

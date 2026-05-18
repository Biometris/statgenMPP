# IBD calculation for multi parental populations

IBD calculation for multi parental populations. Per cross IBD
probabilities are calculated using `calcIBD` in the statgenIBD package.
These probabilities are combined with optional phenotypic data and
stored in a single object of class `gDataMPP`.

## Usage

``` r
calcIBDMPP(
  crossNames,
  markerFiles,
  pheno,
  popType,
  mapFile,
  evalDist,
  grid = TRUE,
  verbose = FALSE
)
```

## Arguments

- crossNames:

  A character vector, the names of the crosses.

- markerFiles:

  A character vector indicating the locations of the files with
  genotypic information for the populations. The files should be in
  tab-delimited format with a header containing marker names.

- pheno:

  A data.frame or a list of data.frames with phenotypic data, with
  genotypes in the first column `genotype` and traits in the following
  columns. The trait columns should be numerical columns only. A list of
  data.frames can be used for replications, i.e. different trials.

- popType:

  A character string indicating the type of population. One of DH, Fx,
  FxDH, BCx, BCxDH, BC1Sx, BC1SxDH, C3, C3DH, C3Sx, C3SxDH, C4, C4DH,
  C4Sx, C4SxDH (see Details).

- mapFile:

  A character string indicating the location of the map file for the
  population. The file should be in tab-delimited format. It should
  consist of exactly three columns, marker, chromosome and position.
  There should be no header. The positions in the file should be in
  centimorgan.

- evalDist:

  A numeric value, the maximum distance in cM between evaluation points.

- grid:

  Should the extra markers that are added to assure the a maximum
  distince of `evalDist` be on a grid (`TRUE`) or in between marker
  existing marker positions (`FALSE`).

- verbose:

  Should progress be printed?

## Value

An object of class `gDataMPP` with the following components:

- `map`:

  a data.frame containing map data. Map is sorted by chromosome and
  position.

- `markers`:

  a 3D matrix containing IBD probabilities.

- `pheno`:

  data.frame or list of data.frames containing phenotypic data.

- `kinship`:

  a kinship matrix.

- `covar`:

  a data.frame with extra covariates (including the name of the cross).

## Details

IBD probabilities can be calculated for many different types of
populations. In the following table all supported populations are
listed. Note that the value of x in the population types is variable,
with its maximum value depicted in the last column.

|  |  |  |  |
|----|----|----|----|
| **Population type** | **Cross** | **Description** | **max. x** |
| DH | biparental | doubled haploid population |  |
| Fx | biparental | Fx population (F1, followed by x-1 generations of selfing) | 8 |
| FxDH | biparental | Fx, followed by DH generation | 8 |
| BCx | biparental | backcross, second parent is recurrent parent | 9 |
| BCxDH | biparental | BCx, followed by DH generation | 9 |
| BC1Sx | biparental | BC1, followed by x generations of selfing | 7 |
| BC1SxDH | biparental | BC1, followed by x generations of selfing and DH | 6 |
| C3 | three-way | three way cross: (AxB) x C |  |
| C3DH | three-way | C3, followed by DH generation |  |
| C3Sx | three-way | C3, followed by x generations of selfing | 7 |
| C3SxDH | three-way | C3, followed by x generations of selfing and DH generation | 6 |
| C4 | four-way | four-way cross: (AxB) x (CxD) |  |
| C4DH | four-way | C4, followed by DH generation |  |
| C4Sx | four-way | C4, followed by x generations of selfing | 6 |
| C4SxDH | four-way | C4, followed by x generations of selfing and DH generation | 6 |

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

summary(ABC)
#> map
#>  Number of markers: 95 
#>  Number of chromosomes: 5 
#> 
#> markers
#>  Number of markers: 95 
#>  Number of genotypes: 180 
#>  Parents: A, B, C 
#> pheno
#>  Number of traits: 1 
#>  Traitnames: yield 
#>  Number of genotypes: 180 
#> 
#> crosses          
#>  AxB:100  
#>  AxC: 80  
```

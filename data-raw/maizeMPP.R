## Use saved data to speed up vignette building.

## Define names of crosses.
crosses <- paste0("F353x", c("B73", "D06", "D09", "EC169", "F252", "F618",
                             "Mo17", "UH250", "UH304", "W117"))
head(crosses)

## Specify files containing crosses.
## Extract them in a temporary directory.
tempDir <- tempdir()
crossFiles <- unzip(system.file("extdata/maize/maize.zip", package = "statgenMPP"),
                    files = paste0(crosses, ".txt"), exdir = tempDir)

## Specify file containing map.
mapFile <- unzip(system.file("extdata/maize/maize.zip", package = "statgenMPP"),
                 files = "map.txt", exdir = tempDir)

## Read phenotypic data.
phenoFile <- unzip(system.file("extdata/maize/maize.zip", package = "statgenMPP"),
                   files = "EUmaizePheno.txt", exdir = tempDir)
phenoDat <- read.delim(phenoFile)
head(phenoDat[, 1:5])

## Perform IBD calculations.
maizeMPP <- calcIBDMPP(crossNames = crosses,
                       markerFiles = crossFiles,
                       pheno = phenoDat,
                       popType = "DH",
                       mapFile = mapFile,
                       evalDist = 5)

## Perform simple interval mapping.
maizeSIM <- selQTLMPP(MPPobj = maizeMPP,
                      trait = "mean_DtSILK",
                      maxCofactors = 0)

maizeCIM <- selQTLMPP(MPPobj = maizeMPP,
                      trait = "mean_DtSILK",
                      threshold = 5)

usethis::use_data(maizeSIM, maizeCIM, internal = FALSE, overwrite = TRUE)

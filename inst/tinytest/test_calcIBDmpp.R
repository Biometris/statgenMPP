### Test calcIBDmpp

## Define input files.

markerFiles <- c(system.file("extdata/multipop", "AxB.txt",
                             package = "statgenMPP"),
                 system.file("extdata/multipop", "AxC.txt",
                             package = "statgenMPP"))
mapFile <- system.file("extdata/multipop", "mapfile.txt",
                       package = "statgenMPP")

pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                                package = "statgenMPP"))
colnames(pheno)[1] <- "genotype"

## Checks for correct input.
expect_error(calcIBDmpp(crossNames = 1),
             "crossNames should be a character vector")
expect_error(calcIBDmpp(crossNames = "AxB", markerFiles = 1),
             "markerFiles should be a character vector")
expect_error(calcIBDmpp(crossNames = "AxB", markerFiles = markerFiles),
             "crossNames and markerFiles should have the same length")
expect_error(calcIBDmpp(crossNames = "AxB", markerFiles = "tst"),
             "The following files don't exist")
expect_error(calcIBDmpp(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = 1),
             "pheno should be a data.frame")
expect_error(calcIBDmpp(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno[, 2:4]),
             "pheno should have a column genotype")
expect_error(calcIBDmpp(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno, popType = 1),
             "popType should be a character string of length 1")
expect_error(calcIBDmpp(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno, popType = "F4DH", mapFile = 1),
             "mapFile should be a character string of length 1")
expect_error(calcIBDmpp(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno, popType = "F4DH", mapFile = mapFile,
                        evalDist = "a"),
             "evalDist should be a positive numerical value")
expect_error(calcIBDmpp(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno, popType = "F4DH", mapFile = mapFile,
                        evalDist = c(1, 2)),
             "evalDist should be a positive numerical value")

expect_silent(ABC <- calcIBDmpp(crossNames = c("AxB", "AxC"),
                                markerFiles = markerFiles,
                                pheno = pheno, popType = "F4DH",
                                mapFile = mapFile, evalDist = 5))


## General structure.
expect_inherits(ABC, "gDataMpp")

expect_inherits(ABC$map, "data.frame")
expect_equal(dim(ABC$map), c(95, 2))

expect_inherits(ABC$markers, "array")
expect_equal(dim(ABC$markers), c(180, 95, 3))

expect_inherits(ABC$pheno$pheno, "data.frame")

expect_inherits(ABC$covar, "data.frame")

expect_inherits(attr(ABC, "genoCross"), "data.frame")
expect_inherits(attr(ABC, "pedigree"), "data.frame")


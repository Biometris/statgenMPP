### Test calcIBDMPP

## Define input files.

markerFiles <- c(system.file("extdata/multipop", "AxB.txt",
                             package = "statgenMPP"),
                 system.file("extdata/multipop", "AxC.txt",
                             package = "statgenMPP"))
mapFile <- system.file("extdata/multipop", "mapfile.txt",
                       package = "statgenMPP")

pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                                package = "statgenMPP"))

phenoLst <- list(pheno[1:100, ], pheno[101:180, ])

## Checks for correct input.
# crossNames.
expect_error(calcIBDMPP(crossNames = 1),
             "crossNames should be a character vector")
expect_error(calcIBDMPP(crossNames = "AxB", markerFiles = markerFiles),
             "crossNames and markerFiles should have the same length")

# markerFiles.
expect_error(calcIBDMPP(crossNames = "AxB", markerFiles = 1),
             "markerFiles should be a character vector")
expect_error(calcIBDMPP(crossNames = "AxB", markerFiles = "tst"),
             "The following files don't exist")

# pheno
expect_error(calcIBDMPP(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = 1),
             "pheno should be a data.frame or a list of data.frames")
expect_error(calcIBDMPP(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno[, 2, drop = FALSE]),
             "pheno should have a column genotype")

# pheno list.
expect_error(calcIBDMPP(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = list(1, 2)),
             "pheno should be a data.frame or a list of data.frames")
expect_error(calcIBDMPP(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = list(pheno[1:100, 2, drop = FALSE],
                                     pheno[101:180, ])),
             "All data.frames in pheno should have a column genotype")

# popType
expect_error(calcIBDMPP(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno, popType = 1),
             "popType should be a character string of length 1")

# mapFile
expect_error(calcIBDMPP(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno, popType = "F4DH", mapFile = 1),
             "mapFile should be a character string of length 1")
expect_error(calcIBDMPP(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno, popType = "F4DH", mapFile = "tst.txt"),
             "mapFile doesn't exist")

# evalDist
expect_error(calcIBDMPP(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno, popType = "F4DH", mapFile = mapFile,
                        evalDist = "a"),
             "evalDist should be a positive numerical value")
expect_error(calcIBDMPP(crossNames = c("AxB", "AxC"), markerFiles = markerFiles,
                        pheno = pheno, popType = "F4DH", mapFile = mapFile,
                        evalDist = c(1, 2)),
             "evalDist should be a positive numerical value")

## Large evalDist for faster computations.
expect_silent(ABC <- calcIBDMPP(crossNames = c("AxB", "AxC"),
                                markerFiles = markerFiles,
                                pheno = pheno, popType = "F4DH",
                                mapFile = mapFile, evalDist = 25))
expect_silent(ABCPhenoLst <- calcIBDMPP(crossNames = c("AxB", "AxC"),
                                        markerFiles = markerFiles,
                                        pheno = phenoLst, popType = "F4DH",
                                        mapFile = mapFile, evalDist = 25))

## General structure.
expect_inherits(ABC, "gDataMPP")

expect_inherits(ABC$map, "data.frame")
expect_equal(dim(ABC$map), c(15, 2))

expect_inherits(ABC$markers, "array")
expect_equal(dim(ABC$markers), c(180, 15, 3))

expect_inherits(ABC$pheno$pheno, "data.frame")

expect_inherits(ABC$covar, "data.frame")

expect_inherits(attr(ABC, "genoCross"), "data.frame")
expect_inherits(attr(ABC, "pedigree"), "data.frame")

expect_equal(ABC, ABCPhenoLst)

## Verbose = TRUE
printOut <- capture.output(
  ABCVerb <- calcIBDMPP(crossNames = c("AxB", "AxC"),
                        markerFiles = markerFiles,
                        pheno = pheno, popType = "F4DH",
                        mapFile = mapFile, evalDist = 25,
                        verbose = TRUE))
expect_true("calculating IBD in cross: AxB." %in% printOut)
expect_true("reading data .............." %in% printOut)
expect_true("analysis of family ........" %in% printOut)

## Single cross.
ABCOne <- calcIBDMPP(crossNames = "AxB",
                     markerFiles = markerFiles[1],
                     pheno = pheno, popType = "F4DH",
                     mapFile = mapFile, evalDist = 25)
expect_true(all(ABCOne$covar$cross == "AxB"))

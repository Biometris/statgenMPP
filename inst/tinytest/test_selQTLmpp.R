### Test selQTLMPP

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

## Compute IBD probabilities.

# High evaldist for faster computations.
ABC <- calcIBDMPP(crossNames = c("AxB", "AxC"),
                  markerFiles = markerFiles,
                  pheno = pheno, popType = "F4DH",
                  mapFile = mapFile, evalDist = 25)

## Check inputs.

expect_error(selQTLMPP(MPPobj = 1),
             "MPPobj should be an object of class gDataMPP")
expect_error(selQTLMPP(MPPobj = ABC, trait = 1),
             "trait should be a character string of length one present in")
expect_error(selQTLMPP(MPPobj = ABC, trait = "tst"),
             "trait should be a character string of length one present in")
expect_error(selQTLMPP(MPPobj = ABC, trait = "yield", QTLwindow = "1"),
             "QTLwindow should be a positive numerical value")
expect_error(selQTLMPP(MPPobj = ABC, trait = "yield", QTLwindow = c(1, 2)),
             "QTLwindow should be a positive numerical value")
expect_error(selQTLMPP(MPPobj = ABC, trait = "yield", threshold = "1"),
             "threshold should be a positive numerical value")
expect_error(selQTLMPP(MPPobj = ABC, trait = "yield", threshold = c(1, 2)),
             "threshold should be a positive numerical value")
expect_error(selQTLMPP(MPPobj = ABC, trait = "yield", maxCofactors = "1"),
             "maxCofactors should be a positive numerical value")
expect_error(selQTLMPP(MPPobj = ABC, trait = "yield", maxCofactors = c(1, 2)),
             "maxCofactors should be a positive numerical value")

ABC_SQM <- selQTLMPP(MPPobj = ABC, trait = "yield", maxCofactors = 0)
ABC_MQM_max <- selQTLMPP(MPPobj = ABC, trait = "yield")

## General structure.
expect_inherits(ABC_MQM_max, "QTLMPP")
expect_inherits(ABC_MQM_max, "GWAS")

expect_equal(names(ABC_MQM_max),
             c("GWAResult", "signSnp", "kinship", "thr", "GWASInfo"))
expect_equal_to_reference(ABC_MQM_max, "ABC_MQM_max", tolerance = 1e-6)

## Option verbose.
expect_stdout(selQTLMPP(MPPobj = ABC, trait = "yield", maxCofactors = 1,
                        verbose = TRUE),
              "QTL scan for trait yield, 0 cofactors")

### Test selQTLmpp

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
ABC <- calcIBDmpp(crossNames = c("AxB", "AxC"),
                  markerFiles = markerFiles,
                  pheno = pheno, popType = "F4DH",
                  mapFile = mapFile, evalDist = 25)

## Check inputs.

expect_error(selQTLmpp(MPPobj = 1),
             "MPPobj should be an object of class gData")
expect_error(selQTLmpp(MPPobj = ABC, trait = 1),
             "trait should be a character string of length one present in")
expect_error(selQTLmpp(MPPobj = ABC, trait = "tst"),
             "trait should be a character string of length one present in")
expect_error(selQTLmpp(MPPobj = ABC, trait = "geno", QTLwindow = "1"),
             "QTLwindow should be a positive numerical value")
expect_error(selQTLmpp(MPPobj = ABC, trait = "geno", QTLwindow = c(1, 2)),
             "QTLwindow should be a positive numerical value")
expect_error(selQTLmpp(MPPobj = ABC, trait = "geno", threshold = "1"),
             "threshold should be a positive numerical value")
expect_error(selQTLmpp(MPPobj = ABC, trait = "geno", threshold = c(1, 2)),
             "threshold should be a positive numerical value")
expect_error(selQTLmpp(MPPobj = ABC, trait = "geno", maxCofactors = "1"),
             "maxCofactors should be a positive numerical value")
expect_error(selQTLmpp(MPPobj = ABC, trait = "geno", maxCofactors = c(1, 2)),
             "maxCofactors should be a positive numerical value")

ABC_SIM <- selQTLmpp(MPPobj = ABC, trait = "geno", CIM = FALSE)
ABC_CIM1 <- selQTLmpp(MPPobj = ABC, trait = "geno", maxCofactors = 0)
ABC_CIM_max <- selQTLmpp(MPPobj = ABC, trait = "geno")

## General structure.
expect_inherits(ABC_CIM_max, "QTLmpp")
expect_inherits(ABC_CIM_max, "GWAS")

expect_equal(names(ABC_CIM_max),
             c("GWAResult", "signSnp", "kinship", "thr", "GWASInfo"))
expect_equal_to_reference(ABC_CIM_max, "ABC_CIM_max", tolerance = 1e-6)

# CIM and SIM with no cofactors should be the same.
expect_equal(ABC_SIM, ABC_CIM1)

## Option verbose.
expect_stdout(selQTLmpp(MPPobj = ABC, trait = "geno", maxCofactors = 1,
                        verbose = TRUE),
              "QTL scan for trait geno, 1 cofactors")


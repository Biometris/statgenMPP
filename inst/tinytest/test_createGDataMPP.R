### Test createGDataMPP

markerFile <- system.file("extdata/multipop", "AxB.txt",
                          package = "statgenMPP")
mapFile <- system.file("extdata/multipop", "mapfile.txt",
                       package = "statgenMPP")

pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                                package = "statgenMPP"))
pheno2 <- pheno[, -1, drop = FALSE]
pheno3 <- pheno
pheno3[["pheno"]] <- as.character(pheno3[["pheno"]])


## Compute IBD probabilities using statgenIBD.
tstIBD <- statgenIBD::calcIBD(popType = "F4DH", markerFile = markerFile,
                              mapFile = mapFile, evalDist = 5)

expect_error(createGDataMPP(IBDprob = 1),
             "IBDprob should be an object of class IBDprob")
expect_error(createGDataMPP(IBDprob = tstIBD, pheno = 1),
             "pheno should be an data.frame")
expect_error(createGDataMPP(IBDprob = tstIBD, pheno = pheno2),
             "The following columns are missing in pheno")
expect_error(createGDataMPP(IBDprob = tstIBD, pheno = pheno3),
             "The following columns in pheno are not numeric")

tstMPP <- createGDataMPP(tstIBD, pheno)

expect_inherits(tstMPP, "gDataMPP")
expect_equal_to_reference(tstMPP, "gDataMPP")

## Check that cross is copied from pheno.
pheno4 <- pheno
pheno4[["cross"]] <- factor(1)

tstMPP_cross <- createGDataMPP(tstIBD, pheno4)
expect_inherits(tstMPP_cross$covar, "data.frame")
expect_null(tstMPP_cross$pheno$pheno[["cross"]])





### Test readRABBIT

## Define input files.

genoFile <- system.file("extdata/barley", "barley_magicReconstruct.zip",
                        package = "statgenMPP")
pedFile <- system.file("extdata/barley", "barley_pedInfo.csv",
                       package = "statgenMPP")

## Checks for correct input.
expect_error(readRABBIT(infile = 1),
             "infile should be a character string indicating a readable")
expect_error(readRABBIT(infile = "tst"),
             "infile should be a character string indicating a readable")
expect_error(readRABBIT(infile = "tst.csv"),
             "infile should be a character string indicating a readable")

expect_error(readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                        pedFile = 1),
             "pedFile should be a character string indicating a readable")
expect_error(readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                        pedFile = "tst"),
             "pedFile should be a character string indicating a readable")
expect_error(readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                        pedFile = "tst.csv"),
             "pedFile should be a character string indicating a readable")

expect_error(readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                        pheno  = 1),
             "pheno should be a data.frame")
expect_error(readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                        pheno  = barleyPheno["cross"]),
             "The following columns are missing in pheno")
barleyPhenoChr <- barleyPheno
barleyPhenoChr[["Awn_length"]] <- as.character(barleyPhenoChr[["Awn_length"]])
expect_error(readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                        pheno  = barleyPhenoChr),
             "The following columns in pheno are not numeric")

## Different combinations of inputs should give similar output.
expect_silent(barleyMPP <-
                readRABBIT(infile = unzip(genoFile, exdir = tempdir())))
expect_silent(barleyMPP2 <-
                readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                           pedFile = pedFile))
expect_silent(barleyMPP3 <-
                readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                           pheno = barleyPheno))
expect_silent(barleyMPP4 <-
                readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                           pedFile = pedFile, pheno = barleyPheno))

## General structure.
expect_inherits(barleyMPP4, "gDataMpp")

expect_inherits(barleyMPP4$map, "data.frame")
expect_equal(dim(barleyMPP4$map), c(355, 2))

expect_inherits(barleyMPP4$markers, "array")
expect_equal(dim(barleyMPP4$markers), c(916, 355, 5))

expect_inherits(barleyMPP4$pheno$pheno, "data.frame")

expect_inherits(barleyMPP4$covar, "data.frame")

expect_inherits(attr(barleyMPP2, "genoCross"), "data.frame")
expect_inherits(attr(barleyMPP2, "pedigree"), "data.frame")

## map and markers should be the same for all.
expect_equal(barleyMPP$map, barleyMPP4$map)
expect_equal(barleyMPP$markers, barleyMPP4$markers)

## covar read from pedFile.
expect_equal(barleyMPP2$covar, barleyMPP4$covar)
expect_equal(attr(barleyMPP2, "pedigree"), attr(barleyMPP4, "pedigree"))
expect_equal(attr(barleyMPP2, "genoCross"), attr(barleyMPP4, "genoCross"))

## Covar read from pheno should be identical to that read from pedFile.
expect_equal(barleyMPP3$pheno, barleyMPP4$pheno)
expect_equal(barleyMPP3$covar, barleyMPP4$covar)





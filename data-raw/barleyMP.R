## Read csv with phenoypic data.
barleyPheno <- read.csv("data-raw/BarleyMP_pheno.csv")
## Rename colums to match input format for readRABBIT.
colnames(barleyPheno)[1] <- "genotype"
colnames(barleyPheno)[3] <- "cross"
## Cross should be a factor for use in residual part of model.
barleyPheno[["cross"]] <- factor(barleyPheno[["cross"]])
## Restrict to relevant columns.
barleyPheno <- barleyPheno[, 1:3]

## Use saved data to speed up vignette building.
tempDir <- tempdir()
inFile <- unzip(system.file("extdata/barley/barley_magicReconstruct.zip",
                            package = "statgenMPP"), exdir = tempDir)

## Specify pedigree file.
pedFile <- system.file("extdata/barley/barley_pedInfo.csv",
                       package = "statgenMPP")

## Read phenotypic data.
data("barleyPheno")

## read RABBIT output.
barleyMPP <- readRABBIT(infile = inFile,
                        pedFile = pedFile,
                        pheno = barleyPheno)

## Perform Multi-QTL Mapping.
barleyMQM <- selQTLMPP(MPPobj = barleyMPP,
                       trait = "Awn_length",
                       threshold = 4)

usethis::use_data(barleyPheno, barleyMQM, overwrite = TRUE)

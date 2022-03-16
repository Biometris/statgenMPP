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

## Compute IBD probabilities.

# High evaldist for faster computations.
ABC <- calcIBDMPP(crossNames = c("AxB", "AxC"),
                  markerFiles = markerFiles,
                  pheno = pheno, popType = "F4DH",
                  mapFile = mapFile, evalDist = 25)

## QTL Detection.
ABC_MQM <- selQTLMPP(MPPobj = ABC, trait = "pheno")

## Summary.
sumABC <- capture.output(summary(ABC_MQM))
expect_true("Number of QTLs: 3 " %in% sumABC)
expect_true("Threshold: 3 " %in% sumABC)
expect_true(" EXT_1_25   1  25    M1_3     15.89   0.225 -1.402  0.809  0.593" %in% sumABC)
expect_true(" EXT_3_75   3  75    M3_8     11.00   0.318  1.056  0.607 -1.662" %in% sumABC)

## QTLProfile is a call to manhattan plot in statgenGWAS. Not much to test.
p1 <- plot(ABC_MQM, plotType = "QTLProfile")
expect_inherits(p1, "ggplot")

## QTLRegion is a call to geneticMap plot in statgenGWAS. Not much to test.
p2 <- plot(ABC_MQM, plotType = "QTLRegion")
expect_inherits(p2, "ggplot")

## Parental effects.
p3 <- plot(ABC_MQM, plotType = "parEffs")
expect_inherits(p3, "ggplot")

## Extended QTL profile.
p4 <- plot(ABC_MQM, plotType = "QTLProfileExt")
expect_inherits(p4, "gtable")


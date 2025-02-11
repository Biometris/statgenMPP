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
ABC_MQM <- selQTLMPP(MPPobj = ABC, trait = "yield")

## Summary.
sumABC <- capture.output(summary(ABC_MQM))
expect_true("Number of QTLs: 3 " %in% sumABC)
expect_true("Threshold: 3 " %in% sumABC)

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

## Confidence intervals around parental effects.
p5 <- plot(ABC_MQM, plotType = "parCIs")
expect_inherits(p5, "ggplot")





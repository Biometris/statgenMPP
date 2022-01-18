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

## QTL Detection.
ABC_CIM <- selQTLmpp(MPPobj = ABC, trait = "geno")

## QTLProfile is a call to manhattan plot in statgenGWAS. Not much to test.
p1 <- plot(ABC_CIM, plotType = "QTLProfile")
expect_inherits(p1, "ggplot")

## QTLRegion is a call to geneticMap plot in statgenGWAS. Not much to test.
p2 <- plot(ABC_CIM, plotType = "QTLRegion")
expect_inherits(p2, "ggplot")

## Parental effects.
p3 <- plot(ABC_CIM, plotType = "parEffs")
expect_inherits(p3, "ggplot")

## Extended QTL profile.
p4 <- plot(ABC_CIM, plotType = "QTLProfileExt")
expect_inherits(p4, "gtable")


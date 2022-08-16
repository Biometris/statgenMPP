## Specify files containing markers.
# One file for each of the two crosses.
markerFiles <- c(system.file("extdata/multipop", "AxB.txt",
                             package = "statgenMPP"),
                 system.file("extdata/multipop", "AxC.txt",
                             package = "statgenMPP"))

## Specify file containing map.
# Both crosses use the same map file.
mapFile <- system.file("extdata/multipop", "mapfile.txt",
                       package = "statgenMPP")

## Read phenotypic data
phenoDat <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                                   package = "statgenMPP"))
# Check contents.
head(phenoDat)

## Perform IBD calculations.
ABCMPP <- calcIBDMPP(crossNames = c("AxB", "AxC"),
                     markerFiles = markerFiles,
                     pheno = phenoDat,
                     popType = "F4DH",
                     mapFile = mapFile,
                     evalDist = 5) #1

## Register parallel back-end with 4 cores.
doParallel::registerDoParallel(cores = 4)

## Run multi-QTL mapping.
ABC_MQM <- selQTLMPP(ABCMPP, trait = "pheno", parallel = TRUE)

## Run multi-QTL mapping with kinship.
ABC_MQM2 <- selQTLMPP(ABCMPP, trait = "pheno", parallel = TRUE, computeKin = TRUE)

ABC_MQM$signSnp
ABC_MQM2$signSnp

plot(ABC_MQM, plotType = "QTLProfileExt")
plot(ABC_MQM2, plotType = "QTLProfileExt")



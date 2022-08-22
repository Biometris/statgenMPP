library(statgenMPP)

# 1. specify files in lists
## 1.1. popC4S3

par.names <- list(par.cross1 = c("BLACK", "BLUE", "GREEN", "RED")) # parent names
cross.names <- list(cross1.name = "popC4S3") # cross name
loc.names <- list(loc.cross1 = "testScripts/data/popC4S3/cross.loc") # loc file
qua.names <- list(loc.cross1 = "testScripts/data/popC4S3/cross.qua") # qua file
pop.types <- list(poptype.cross1 = "C4S3") # cross type
mapfile <- "testScripts/data/popC4S3/mapfile.map" # map file
evaldist <- 1  # distance for IBD calc

# 2. calculate IBDs and create the data frame
MPPobj <- calcIBDmppOrig(par.names, cross.names, loc.names, qua.names,
                         pop.types, mapfile, evaldist)
head(MPPobj$calcIBDres$IBDdata)[1:6,1:10] # pheno + design matrix for LMMsolve latter
head(MPPobj$calcIBDres$IBDmatrix)[1:6,1:6]

# 3. genome scan for multi-QTLs
MPPobj <- selQTLmppOrig(MPPobj, QTLwindow = 10, threshold = 3,
                        trait.name = "pheno", CIM = TRUE)
MPPobj$Result$QTLcandidates


##### Identical analysis using new functions.

crossNames <- "popC4S3" # cross name(s)
locFiles <- "testScripts/data/popC4S3/cross.loc" # loc files
quaFiles <- "testScripts/data/popC4S3/cross.qua" # qua files
popType <- "C4S3" # cross type(s)
mapFile <- "testScripts/data/popC4S3/mapfile.map" # map file names
evalDist <- 5 # distance for IBD calc

## Read phenotypic data.
quaDat <- lapply(quaFiles, read.table, header = TRUE)

# 2. calculate IBDs and create the data frame
MPPobj2 <- calcIBDmpp(crossNames = crossNames,
                      markerFiles = locFiles,
                      pheno = quaDat,
                      popType = popType,
                      mapFile = mapFile,
                      evalDist = evalDist,
                      grid = TRUE,
                      verbose = TRUE)

# 3. genome scan for multi-QTLs
MPPobj2a <- selQTLmpp(MPPobj2, QTLwindow = 10, threshold = 1.5,
                      trait = "yield", CIM = TRUE)
MPPobj2a$signSnp


plot(MPPobj2a)
plot(MPPobj2a, plotType = "parEffs")
plot(MPPobj2a, plotType = "QTLRegion")

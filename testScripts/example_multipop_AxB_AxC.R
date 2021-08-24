library(statgenMPP)

# load original functions Wenhao:
source("R/orig/calcIBDmppOrig.R")
source("R/orig/selQTLmppOrig.R")
source("R/orig/scanQTLOrig.R")
source("R/orig/utilsOrig.R")
source("R/orig/randomQTLmodelOrig.R")
source("R/orig/plotQTLscanOrig.R")
##

# 1. specify files in lists
## 1.1. AxBxC
par.names <- list(par.cross1 = c("A","B"),
                  par.cross2 = c("A","C")) # parents names per cross
cross.names <- list(cross1.name = "AxB",
                    cross2.name = "AxC") # cross name(s)
loc.names <- list(loc.cross1 = "testScripts/data/multipop/AxB.loc",
                  loc.cross2 = "testScripts/data/multipop/AxC.loc") # loc files
qua.names <- list(loc.cross1 = "testScripts/data/multipop/AxB.qua",
                  loc.cross2 = "testScripts/data/multipop/AxC.qua") # qua files
pop.types <- list(poptype.cross1 = "F4DH",
                  poptype.cross2 = "F4DH") # cross type(s)
mapfile <- "testScripts/data/multipop/mapfile.map" # map file names
evaldist <- 5 # distance for IBD calc

# 2. calculate IBDs and create the data frame
MPPobj1 <- calcIBDmppOrig(par.names, cross.names, loc.names, qua.names,
                         pop.types, mapfile, evaldist)
head(MPPobj1$calcIBDres$IBDdata)[1:6,1:10] # pheno + design matrix for LMMsolve latter
head(MPPobj1$calcIBDres$IBDmatrix)[1:6,1:6]

# 3. genome scan for multi-QTLs
MPPobj1a <- selQTLmppOrig(MPPobj1, QTLwindow = 10, threshold = 3,
                          trait.name = "pheno", CIM = TRUE)
MPPobj1a$Result$QTLcandidates


##### Identical analysis using new functions.

crossNames <- c("AxB", "AxC") # cross name(s)
locFiles <- c("testScripts/data/multipop/AxB.loc",
              "testScripts/data/multipop/AxC.loc") # loc files
quaFiles <- c("testScripts/data/multipop/AxB.qua",
              "testScripts/data/multipop/AxC.qua") # qua files
popType <- "F4DH" # cross type(s)
mapFile <- "testScripts/data/multipop/mapfile.map" # map file names
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
MPPobj2a <- selQTLmpp(MPPobj2, QTLwindow = 10, threshold = 3,
                      trait = "pheno", CIM = TRUE)
MPPobj2a$Result$QTLcandidates

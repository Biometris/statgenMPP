library(statgenMPP)

# 1. specify files in lists
## 1.1. popC3S4DH

par.names <- list(par.cross1 = c("BLACK", "BLUE", "GREEN")) # parent names
cross.names <- list(cross1.name = "popC3S4DH") # cross name
loc.names <- list(loc.cross1 = "testScripts/data/popC3S4DH/cross.loc") # loc file
qua.names <- list(loc.cross1 = "testScripts/data/popC3S4DH/cross.qua") # qua file
pop.types <- list(poptype.cross1 = "C3S4DH") # cross type
mapfile <- "testScripts/data/popC3S4DH/mapfile.map" # map file
evaldist <- 1  # distance for IBD calc

# 2. calculate IBDs and create the data frame
MPPobj <- calcIBDmppOrig(par.names, cross.names, loc.names, qua.names,
                         pop.types, mapfile, evaldist)
head(MPPobj$calcIBDres$IBDdata)[1:6,1:10] # pheno + design matrix for LMMsolve latter
head(MPPobj$calcIBDres$IBDmatrix)[1:6,1:6]

# 3. genome scan for multi-QTLs
MPPobj <- selQTLMPP(MPPobj, QTLwindow = 10, threshold = 3,
                    trait.name = "pheno", CIM = TRUE)
MPPobj$Result$QTLcandidates


crossNames <- "popC3S4DH"
locFiles <- "testScripts/data/popC3S4DH/cross.loc" # loc file
quaFiles <- "testScripts/data/popC3S4DH/cross.qua" # qua file
poptypes <- "C3S4DH" # cross type
mapFile <- "testScripts/data/popC3S4DH/mapfile.map" # map file
evaldist <- 1

quaDat <- lapply(quaFiles, read.table, header = TRUE)

MPPobj2 <- calcIBDmpp(crossNames, locFiles, quaDat,
                      poptypes, mapFile, evaldist)

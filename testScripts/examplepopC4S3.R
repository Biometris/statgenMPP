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
MPPobj <- calcIBDmpp(par.names, cross.names, loc.names, qua.names,
                     pop.types, mapfile, evaldist)
head(MPPobj$calcIBDres$IBDdata)[1:6,1:10] # pheno + design matrix for LMMsolve latter
head(MPPobj$calcIBDres$IBDmatrix)[1:6,1:6]

# 3. genome scan for multi-QTLs
MPPobj <- selQTLmpp(MPPobj, QTLwindow = 10, threshold = 3,
                    trait.name = "pheno", CIM = TRUE)
MPPobj$Result$QTLcandidates

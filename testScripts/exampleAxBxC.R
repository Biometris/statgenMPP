library(statgenMPP)

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
MPPobj <- calcIBDmpp(par.names, cross.names, loc.names, qua.names,
                     pop.types, mapfile, evaldist)
head(MPPobj$calcIBDres$IBDdata)[1:6,1:10] # pheno + design matrix for LMMsolve latter
head(MPPobj$calcIBDres$IBDmatrix)[1:6,1:6]

# 3. genome scan for multi-QTLs
MPPobj <- selQTLmpp(MPPobj, QTLwindow = 10, threshold = 3,
                    trait.name = "pheno", CIM = TRUE)
MPPobj$Result$QTLcandidates

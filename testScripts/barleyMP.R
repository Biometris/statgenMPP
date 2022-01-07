library(statgenMPP)

genoFile <- system.file("extdata", "BarleyMP_magicReconstruct_Summary.zip",
                        package = "statgenMPP")

barleyMPMPP <- readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                          pheno = barleyPheno)

barleyRes <- selQTLmpp(MPPobj = barleyMPMPP, trait = "Awn_length",
                       maxCofactors = 1)

plot(barleyRes, plotType = "manhattan")
plot(barleyRes, plotType = "parEffs")
plot(barleyRes, plotType = "QTLRegion")

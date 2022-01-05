genoFile <- system.file("extdata", "BarleyMP_magicReconstruct_Summary.zip",
                        package = "statgenMPP")

pheno <- read.csv(system.file("extdata", "BarleyMP_pheno.csv",
                              package = "statgenMPP"))
colnames(pheno)[1] <- "genotype"
colnames(pheno)[3] <- "cross"
pheno$cross <- factor(pheno$cross)
pheno <- pheno[, 1:3]


barleyMPMPP <- readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                          pheno = pheno)

barleyRes <- selQTLmpp(MPPobj = barleyMPMPP, trait = "Awn_length")

plot(barleyRes, plotType = "manhattan")
plot(barleyRes, plotType = "parEffs")
plot(barleyRes, plotType = "QTLRegion")

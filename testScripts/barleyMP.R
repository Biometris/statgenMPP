genoFile <- system.file("extdata", "BarleyMP_magicReconstruct_Summary.zip",
                        package = "statgenMPP")

pheno <- read.csv(system.file("extdata", "BarleyMP_pheno.csv",
                              package = "statgenMPP"))
phenoDat <- pheno[1:2]
colnames(phenoDat)[1] <- "genotype"

covar <- pheno[, 3, drop = FALSE]
covar$cross <- factor(paste0("cross", covar$pop))
rownames(covar) <- phenoDat$genotype

barleyMPMPP <- readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                          pheno = phenoDat, covar = covar)

barleyRes <- selQTLmpp(MPPobj = barleyMPMPP, trait = "Awn_length")

plot(barleyRes, plotType = "manhattan")
plot(barleyRes, plotType = "parEffs")
plot(barleyRes, plotType = "QTLRegion")

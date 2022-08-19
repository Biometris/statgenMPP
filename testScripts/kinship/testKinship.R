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
ABC_MQM <- selQTLMPP(ABCMPP, trait = "pheno", parallel = TRUE, maxCofactors = 0)

## Run multi-QTL mapping with kinship.
ABC_MQM2 <- selQTLMPP(ABCMPP, trait = "pheno", parallel = TRUE, computeKin = TRUE,
                      maxCofactors = 0)

ABC_MQM$signSnp
ABC_MQM2$signSnp

plot(ABC_MQM, plotType = "QTLProfileExt")
plot(ABC_MQM2, plotType = "QTLProfileExt")

# Example how to implement kinship using LMMsolver

library(LMMsolver)
eps <- 1.0e-10
chr <- 5
map <- ABCMPP$map
markers <- ABCMPP$markers

mrkNamesNonChr <- rownames(map)[map$chr != chr]
mrkNonChr <- apply(markers[, mrkNamesNonChr, ], MARGIN = 2,
                   FUN = tcrossprod, simplify = FALSE)
KChr <- Reduce(`+`, mrkNonChr) / length(mrkNonChr)
eigKChr <- eigen(KChr)

sel <- eigKChr$values > eps
U <- eigKChr$vectors[, sel]
d <- eigKChr$values[sel]
Usc <- U %*% diag(sqrt(d))
# Should be close to original kinship
range(KChr - tcrossprod(Usc))

colnames(Usc) <- paste0("eig", 1:ncol(Usc))
y <- ABCMPP$pheno$pheno$pheno
dat <- data.frame(y, Usc)

lGrp <- list(Eigen=2:ncol(dat))
obj <- LMMsolve(y~1, random=~grp(Eigen), group=lGrp, data=dat)
summary(obj)

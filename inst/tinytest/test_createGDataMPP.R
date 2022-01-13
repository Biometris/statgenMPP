### Test calcIBDmpp

## Define input files.

markerFiles <- c(system.file("extdata/multipop", "AxB.txt",
                             package = "statgenMPP"),
                 system.file("extdata/multipop", "AxC.txt",
                             package = "statgenMPP"))
mapFile <- system.file("extdata/multipop", "mapfile.txt",
                       package = "statgenMPP")

pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                                package = "statgenMPP"))
colnames(pheno)[1] <- "genotype"

## Compute IBD probabilities.
expect_silent(ABC <- calcIBDmpp(crossNames = c("AxB", "AxC"),
                                markerFiles = markerFiles,
                                pheno = pheno, popType = "F4DH",
                                mapFile = mapFile, evalDist = 5))

## Check summary.

sumABC <- summary(ABC)

expect_inherits(sumABC, "summary.gData")
expect_equal(names(sumABC), c("mapSum", "markerSum", "phenoSum", "covarSum"))

# Function is a copy of function from statgenGWAS except for markerSum bit.
# Only check that.

expect_inherits(sumABC$markerSum, "list")
expect_equal(names(sumABC$markerSum), c("nMarkers", "nGeno", "markerContent"))
expect_equal(sumABC$markerSum$markerContent, c(`parents:` = "A, B, C"))

## Check plot.

# All plots are called directly from statgenGWAS or statgenIBD.
# Not much can be tested here.

# genMap
expect_silent(p1 <- plot(ABC, plotType = "genMap"))
p1a <- plot(ABC, highlight = "EXT_3_30")

expect_equal(length(p1a$layers), length(p1$layers) + 2)


expect_silent(p2 <- plot(ABC, plotType = "allGeno"))
expect_silent(p3 <- plot(ABC, plotType = "pedigree"))


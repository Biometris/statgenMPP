## Read csv with phenoypic data.
barleyPheno <- read.csv("data-raw/BarleyMP_pheno.csv")
## Rename colums to match input format for readRABBIT.
colnames(barleyPheno)[1] <- "genotype"
colnames(barleyPheno)[3] <- "cross"
## Cross should be a factor for use in residual part of model.
barleyPheno[["cross"]] <- factor(barleyPheno[["cross"]])
## Restrict to relevant columns.
barleyPheno <- barleyPheno[, 1:3]

usethis::use_data(barleyPheno, overwrite = TRUE)

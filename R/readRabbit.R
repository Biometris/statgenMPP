#' Read IBD probabilities
#'
#' Read a file with IBD probabilities computed by the RABBIT software package.
#'
#' @param infile A character string, a link to a .csv file with IBD
#' probabilities.
#' @param pheno A data frame with at least columns "genotype" for the
#' "genotype", "cross" indicating the cross the genotype comes from and one
#' or more numerical columns containing phenotypic information.
#'
#' @return A \code{gData} object with map and markers corresponding to the
#' imported information in the imported .csv file.
#'
#' @examples
#' genoFile <- system.file("extdata", "BarleyMP_magicReconstruct_Summary.zip",
#'                        package = "statgenMPP")
#' barleyMPMPP <- readRABBIT(unzip(genoFile, exdir = tempdir()))
#'
#' @references Fine mapping of a major QTL for awn length in barley using a
#' multiparent mapping population. Liller CB, Walla A, Boer MP, Hedley P,
#' Macaulay M, Effgen S, von Korff M, van Esse GW, Koornneef M.
#' Theor Appl Genet. 2017 Feb;130(2):269-281. \doi{10.1007/s00122-016-2807-y}.
#'
#' @importFrom utils hasName
#' @export
readRABBIT <- function(infile,
                       pheno = NULL) {
  if (missing(infile) || !is.character(infile) || length(infile) > 1 ||
      file.access(infile, mode = 4) == -1 || tools::file_ext(infile) != "csv") {
    stop("infile should be a character string indicating a readable .csv file")
  }
  if (!is.null(pheno)) {
    if (!inherits(pheno, "data.frame")) {
      stop("pheno should be a data.frame.\n")
    }
    minCols <- c("genotype", "cross")
    missCols <- minCols[!sapply(X = minCols, FUN = function(minCol) {
      hasName(x = pheno, name = minCol)
    })]
    if (length(missCols) > 0) {
      stop("The following columns are missing in pheno:\n",
           paste(missCols, collapse = ", "))
    }
    trtCols <- colnames(pheno)[!colnames(pheno) %in% minCols]
    nonNumCols <- trtCols[!sapply(X = pheno[trtCols], FUN = is.numeric)]
    if (length(nonNumCols) > 0) {
      stop("The following columns in pheno are not numeric:\n",
           paste(nonNumCols, collapse = ", "))
    }
  }
  ## Read map and marker probabilities.
  markMap <- data.table::fread(infile, skip = "haploprob", fill = TRUE,
                               data.table = FALSE)
  ## Extract map.
  map <- data.frame(chr = as.numeric(markMap[3, -1]),
                    pos = as.numeric(markMap[4, -1]),
                    row.names = as.character(markMap[2, -1]))
  markMap <- markMap[5:(which(markMap[[2]] == "ibdprob") - 1), ]
  ## Get names of genotypes and compute number of founder alleles per genotype.
  genoNames <- unique(sapply(X = strsplit(x = markMap[, 1],
                                          split = "_haplotype"), FUN = "[[", 1))
  nAlleles = nrow(markMap) / length(genoNames)
  ## Convert markers to 3D array.
  markArr <- array(dim = c(length(genoNames), nrow(map), nAlleles))
  for (i in 1:nrow(map)) {
    markArr[, i, ] <- matrix(as.numeric(markMap[, i + 1]),
                             ncol = nAlleles, byrow = TRUE)
  }
  ## Read founder names from file.
  foundNames <- data.table::fread(infile, header = FALSE, nrows = nAlleles + 2,
                                  skip = "haplotypes in order",
                                  data.table = FALSE, select = 3)
  foundNames <- as.character(foundNames[3:nrow(foundNames), ])
  ## Add dimnames to markers: genotypes x markers x founders.
  dimnames(markArr) <- list(genoNames, rownames(map), foundNames)
  ## Split pheno in pheno and covar.
  covar <- pheno["cross"]
  rownames(covar) <- pheno[["genotype"]]
  pheno <- pheno[-which(colnames(pheno) == "cross")]
  ## Create gData object.
  res <- createGData(geno = markArr, map = map, pheno = pheno, covar = covar)
  attr(x = res, which = "popType") <- "RABBIT"

  # attr(x = res, which = "pedigree") <- crossIBD$pedigree
  # attr(x = res, which = "genoCross") <- attr(x = crossIBD, which = "genoCross")

  return(res)
}

#' Read IBD probabilities
#'
#' Read a file with IBD probabilities computed by the RABBIT software package.
#' It is possible to additionally read the pedigree file that is also used by
#' RABBIT. Reading this file allows for plotting the pedigree. Phenotypic data
#' can be added from a data.frame.
#'
#' @param infile A character string, a link to a .csv file with IBD
#' probabilities.
#' @param pedFile A character string, a link to a .csv file with pedigree
#' information as used by RABBIT as input.
#' @param pheno A data frame with at least columns "genotype" for the
#' "genotype" and one or more numerical columns containing phenotypic
#' information. A column "cross" can be used for indicating the cross the
#' genotype comes from.
#'
#' @returns An object of class \code{gDataMPP} with map and markers
#' corresponding to the imported information in the imported .csv file.
#'
#' @examples
#' ## Read RABBIT data for barley.
#' genoFile <- system.file("extdata/barley", "barley_magicReconstruct.zip",
#'                        package = "statgenMPP")
#' barleyMPMPP <- readRABBITMPP(unzip(genoFile, exdir = tempdir()))
#'
#' @references Zheng, Chaozhi, Martin P Boer, and Fred A Van Eeuwijk.
#' “Recursive Algorithms for Modeling Genomic Ancestral Origins in a Fixed
#' Pedigree.” G3 Genes|Genomes|Genetics 8 (10): 3231–45.
#' https://doi.org/10.1534/G3.118.200340.
#'
#' @importFrom utils hasName unzip
#' @export
readRABBITMPP <- function(infile,
                          pedFile = NULL,
                          pheno = NULL) {
  if (!is.null(pheno)) {
    if (!inherits(pheno, "data.frame")) {
      stop("pheno should be a data.frame.\n")
    }
    ## If pedFile is included cross will be read from there.
    minCols <- "genotype"
    missCols <- minCols[!sapply(X = minCols, FUN = function(minCol) {
      hasName(x = pheno, name = minCol)
    })]
    if (!is.null(pedFile) && hasName(x = pheno, name = "cross")) {
      pheno <- pheno[-which(colnames(pheno) == "cross")]
    }
    if (length(missCols) > 0) {
      stop("The following columns are missing in pheno:\n",
           paste(missCols, collapse = ", "))
    }
    trtCols <- colnames(pheno)[!colnames(pheno) %in% c(minCols, "cross")]
    nonNumCols <- trtCols[!sapply(X = pheno[trtCols], FUN = is.numeric)]
    if (length(nonNumCols) > 0) {
      stop("The following columns in pheno are not numeric:\n",
           paste(nonNumCols, collapse = ", "))
    }
  }
  rabbitRes <- statgenIBD::readRABBIT(infile, pedFile)
  covar <- NULL
  genoCross <- attr(x = rabbitRes, which = "genoCross")
  if (!is.null(pheno) && hasName(x = pheno, name = "cross")) {
    ## Split pheno in pheno and covar.
    covar <- pheno["cross"]
    rownames(covar) <- pheno[["genotype"]]
    pheno <- pheno[-which(colnames(pheno) == "cross")]
    genoCross <- covar
    genoCross[["geno"]] <- rownames(covar)
  } else if (!is.null(genoCross)) {
    covar <- genoCross
    rownames(covar) <- covar[["geno"]]
    covar[["geno"]] <- NULL
  }
  ## Create gDataMPP object.
  res <- createGDataMPPInternal(geno = rabbitRes$markers,
                                map = rabbitRes$map,
                                pheno = pheno,
                                covar = covar)
  attr(x = res, which = "popType") <- rabbitRes$popType
  attr(x = res, which = "genoCross") <- genoCross
  attr(x = res, which = "mapOrig") <- rabbitRes$map
  return(res)
}

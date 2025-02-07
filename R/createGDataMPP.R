#' Create an object of class gDataMPP
#'
#' Function for creating an object of class gDataMPP from an object of class
#' IBDprob (computed or imported using statgenIBD) and phenotypic data.
#'
#' @param IBDprob An object of class \code{IBDprob}.
#' @param pheno A data frame with at least columns "genotype" for the
#' "genotype" and one or more numerical columns containing phenotypic
#' information. A column "cross" can be used for indicating the cross the
#' genotype comes from. This column is ignored if the `IBDprob` has an
#' attribute genoCross containing this information (the default behaviour).
#'
#' @returns An object of class \code{gDataMPP}
#'
#' @examples
#' ## Read phenotypic data.
#' pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
#'                                package = "statgenMPP"))
#' ## Rename first column to genotype.
#' colnames(pheno)[1] <- "genotype"
#'
#' ## Compute IBD probabilities for simulated population using statgenIBD - AxB
#' AB <- statgenIBD::calcIBD(popType = "F4DH",
#'                          markerFile = system.file("extdata/multipop", "AxB.txt",
#'                                                  package = "statgenMPP"),
#'                          mapFile = system.file("extdata/multipop", "mapfile.txt",
#'                                               package = "statgenMPP"),
#'                          evalDist = 5)
#'
#' ## Create an object of class gDataMPP from IBD probabilities and
#' ## phenotypic data.
#' ABMPP <- createGDataMPP(IBDprob = AB, pheno = pheno)
#'
#' @importFrom utils packageVersion
#' @export
createGDataMPP <- function(IBDprob,
                           pheno) {

  if (!inherits(IBDprob, "IBDprob")) {
    stop("IBDprob should be an object of class IBDprob.\n")
  }
  if (!inherits(pheno, "data.frame")) {
    stop("pheno should be a data.frame.\n")
  }
  genoCross <- attr(IBDprob, "genoCross")
  markers <- IBDprob$markers
  ## If pedFile is included cross will be read from there.
  minCols <- "genotype"
  missCols <- minCols[!sapply(X = minCols, FUN = function(minCol) {
    hasName(x = pheno, name = minCol)
  })]
  if (length(missCols) > 0) {
    stop("The following columns are missing in pheno:\n",
         paste(missCols, collapse = ", "))
  }
  if (!is.null(genoCross) && hasName(x = pheno, name = "cross")) {
    pheno <- pheno[-which(colnames(pheno) == "cross")]
  }
  trtCols <- colnames(pheno)[!colnames(pheno) %in% c(minCols, "cross")]
  nonNumCols <- trtCols[!sapply(X = pheno[trtCols], FUN = is.numeric)]
  if (length(nonNumCols) > 0) {
    stop("The following columns in pheno are not numeric:\n",
         paste(nonNumCols, collapse = ", "))
  }
  if (hasName(x = pheno, name = "cross")) {
    ## Split pheno in pheno and covar.
    covar <- pheno["cross"]
    rownames(covar) <- pheno[["genotype"]]
    pheno <- pheno[-which(colnames(pheno) == "cross")]
  } else {
    ## Construct covar from genoCross.
    covar <- genoCross["cross"]
    rownames(covar) <- genoCross[["geno"]]
  }
  ## Create gDataMPP object.
  res <- createGDataMPPInternal(geno = markers, map = IBDprob$map,
                                pheno = pheno, covar = covar)
  attr(x = res, which = "popType") <- IBDprob$popType
  attr(x = res, which = "genoCross") <- genoCross
  attr(x = res, which = "mapOrig") <- IBDprob$map
  attr(x = res, which = "pedigree") <- IBDprob$pedigree
  return(res)
}

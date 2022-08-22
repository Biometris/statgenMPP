#' Compute kinship matrix for IBD probabilities
#'
#' Compute a kinship matrix or a list of chromosome specific kinship matrices
#' for a 3D array of IBD probabilities. The kinship matrix is computed by
#' averaging \eqn{Z Z^t} over all markers, where \eqn{Z} is the genotype x
#' parents matrix for the marker. If \code{chrSpecific = TRUE} chromosome
#' specific kinship matrices are computed for each chromosome based only on the
#' markers on all other chromosomes.
#'
#' @param markers An n x m x p array with IBD probabilities with genotypes in
#' the rows (n), markers in the columns (m), and parents in the 3rd dimension
#' (p).
#' @param map A data.frame with columns \code{chr} for chromosome and
#' \code{pos} for position. Positions should be in centimorgan (cM). They
#' should not be cumulative over the chromosomes. Other columns are ignored.
#' Marker names should be in the row names. These should match the marker names
#' in the input file. Only required if \code{chrSpecific = TRUE}.
#' @param chrSpecific Should chromosome specific kinship matrices be
#' computed?
#'
#' @return A kinship matrix or a list of chromosome specific kinship matrices.
#'
#' @examples
#' ## Read phenotypic data.
#' pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
#'                                package = "statgenMPP"))
#' ## Rename first column to genotype.
#' colnames(pheno)[1] <- "genotype"
#'
#'
#' ## Compute IBD probabilities for simulated population - AxB, AxC.
#' ABC <- calcIBDMPP(crossNames = c("AxB", "AxC"),
#'                   markerFiles = c(system.file("extdata/multipop", "AxB.txt",
#'                                               package = "statgenMPP"),
#'                                   system.file("extdata/multipop", "AxC.txt",
#'                                               package = "statgenMPP")),
#'                   pheno = pheno,
#'                   popType = "F4DH",
#'                   mapFile = system.file("extdata/multipop", "mapfile.txt",
#'                                         package = "statgenMPP"),
#'                   evalDist = 5)
#'
#' ## Compute chromosome specific kinship matrices.
#' KChrSpec <- kinshipIBD(markers = ABC$markers, map = ABC$map)
#'
#' @export
kinshipIBD <- function(markers,
                       map = NULL,
                       chrSpecific = TRUE) {
  if (!is.array(markers) && !length(dim(markers)) == 3) {
    stop("markers should be a 3 dimensional array.\n")
  }
  if (chrSpecific) {
    if (is.null(map) || !is.data.frame(map)) {
      stop("map should be a data.frame.\n")
    }
    if (!all(hasName(x = map, name = c("chr", "pos")))) {
      ## chr and pos are obligatory cols.
      stop("chr and pos should be columns in map.\n")
    }
    K <- sapply(X = as.character(unique(map$chr)), FUN = function(chr) {
      mrkNamesNonChr <- rownames(map)[map$chr != chr]
      mrkNonChr <- apply(markers[, mrkNamesNonChr, ], MARGIN = 2,
                         FUN = tcrossprod, simplify = FALSE)
      KChr <- Reduce(`+`, mrkNonChr) / length(mrkNonChr)
      rownames(KChr) <- colnames(KChr) <- rownames(markers)
      return(KChr)
    }, simplify = FALSE)
  } else {
    mrk <- apply(markers, MARGIN = 2, FUN = tcrossprod, simplify = FALSE)
    K <- Reduce(`+`, mrk) / ncol(markers)
    rownames(K) <- colnames(K) <- rownames(markers)
  }
  return(K)
}

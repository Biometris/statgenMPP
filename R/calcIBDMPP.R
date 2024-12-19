#' IBD calculation for multi parental populations
#'
#' IBD calculation for multi parental populations. Per cross IBD probabilities
#' are calculated using \code{calcIBD} in the statgenIBD package. These
#' probabilities are combined with optional phenotypic data and stored in a
#' single object of class \code{gDataMPP}.
#'
#' IBD probabilities can be calculated for many different types of populations.
#' In the following table all supported populations are listed. Note that the
#' value of x in the population types is variable, with its maximum value
#' depicted in the last column.
#'
#' | __Population type__ | __Cross__ | __Description__ | __max. x__ |
#' | ------ | ----------------- | -------------------------------------- | --- |
#' | DH | biparental | doubled haploid population | |
#' | Fx | biparental | Fx population (F1, followed by x-1 generations of selfing) | 8 |
#' | FxDH | biparental | Fx, followed by DH generation | 8 |
#' | BCx | biparental | backcross, second parent is recurrent parent | 9 |
#' | BCxDH | biparental | BCx, followed by DH generation | 9 |
#' | BC1Sx | biparental | BC1, followed by x generations of selfing | 7 |
#' | BC1SxDH | biparental | BC1, followed by x generations of selfing and DH | 6 |
#' | C3 | three-way | three way cross: (AxB) x C |  |
#' | C3DH | three-way | C3, followed by DH generation |  |
#' | C3Sx | three-way | C3, followed by x generations of selfing | 7 |
#' | C3SxDH | three-way | C3, followed by x generations of selfing and DH generation | 6 |
#' | C4 | four-way | four-way cross: (AxB) x (CxD)	| |
#' | C4DH | four-way | C4, followed by DH generation |  |
#' | C4Sx | four-way | C4, followed by x generations of selfing | 6 |
#' | C4SxDH | four-way | C4, followed by x generations of selfing and DH generation | 6 |
#'
#' @param crossNames A character vector, the names of the crosses.
#' @param markerFiles A character vector indicating the locations of the files
#' with genotypic information for the populations. The files should be in
#' tab-delimited format with a header containing marker names.
#' @param pheno A data.frame or a list of data.frames with phenotypic data,
#' with genotypes in the first column \code{genotype} and traits in the
#' following columns. The trait columns should be numerical columns only.
#' A list of data.frames can be used for replications, i.e. different
#' trials.
#' @param popType A character string indicating the type of population. One of
#' DH, Fx, FxDH, BCx, BCxDH, BC1Sx, BC1SxDH, C3, C3DH, C3Sx, C3SxDH, C4, C4DH,
#' C4Sx, C4SxDH (see Details).
#' @param mapFile A character string indicating the location of the map file
#' for the population. The file should be in tab-delimited format. It should
#' consist of exactly three columns, marker, chromosome and position. There
#' should be no header. The positions in the file should be in centimorgan.
#' @param evalDist A numeric value, the maximum distance in cM between
#' evaluation points.
#' @param grid Should the extra markers that are added to assure the a
#' maximum distince of \code{evalDist} be on a grid (\code{TRUE}) or in between
#' marker existing marker positions (\code{FALSE}).
#' @param verbose Should progress be printed?
#'
#' @returns An object of class \code{gDataMPP} with the following components:
#' \item{\code{map}}{a data.frame containing map data. Map is sorted by
#' chromosome and position.}
#' \item{\code{markers}}{a 3D matrix containing IBD probabilities.}
#' \item{\code{pheno}}{data.frame or list of data.frames containing phenotypic
#' data.}
#' \item{\code{kinship}}{a kinship matrix.}
#' \item{\code{covar}}{a data.frame with extra covariates (including the
#' name of the cross).}
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
#' summary(ABC)
#'
#' @importFrom utils read.delim packageVersion
#' @importFrom stats setNames
#' @export
calcIBDMPP <- function(crossNames,
                       markerFiles,
                       pheno,
                       popType,
                       mapFile,
                       evalDist,
                       grid = TRUE,
                       verbose = FALSE) {
  ## Checks.
  if (!is.character(crossNames)) {
    stop("crossNames should be a character vector.\n")
  }
  if (!is.character(markerFiles)) {
    stop("markerFiles should be a character vector.\n")
  }
  if (length(crossNames) != length(markerFiles)) {
    stop("crossNames and markerFiles should have the same length.\n")
  }
  missFiles <- markerFiles[!file.exists(markerFiles)]
  if (length(missFiles) > 0) {
    stop("The following files don't exist: \n",
         paste(missFiles, collapse = ", "))
  }
  if (!is.null(pheno) && !inherits(pheno, "data.frame") &&
      !inherits(pheno, "list")) {
    stop("pheno should be a data.frame or a list of data.frames.\n")
  }
  if (inherits(pheno, "data.frame") && !hasName(x = pheno, name = "genotype")) {
    stop("pheno should have a column genotype.\n")
  }
  if (inherits(pheno, "list")) {
    if (!all(sapply(X = pheno, FUN = inherits, "data.frame"))) {
      stop("pheno should be a data.frame or a list of data.frames.\n")
    }
    if (!all(sapply(X = pheno, FUN = function(x) {
      hasName(x = x, name = "genotype")
    }))) {
      stop("All data.frames in pheno should have a column genotype.\n")
    }
  }
  if (!is.character(popType) || length(popType) > 1) {
    stop("popType should be a character string of length 1.\n")
  }
  if (!is.character(mapFile) || length(mapFile) > 1) {
    stop("mapFile should be a character string of length 1.\n")
  }
  if (!file.exists(mapFile)) {
    stop("mapFile doesn't exist.\n")
  }
  if (!is.numeric(evalDist) || length(evalDist) > 1 || evalDist < 0) {
    stop("evalDist should be a positive numerical value.\n")
  }
  ## Get number of crosses.
  nCross <- length(markerFiles)
  ## Calculate IBD probabilities per cross.
  crossIBD <- lapply(X = seq_len(nCross), FUN = function(i) {
    if (verbose) {
      cat(paste0("calculating IBD in cross: ", crossNames[i], ".\n"))
    }
    statgenIBD::calcIBD(popType = popType,
                        markerFile = markerFiles[i],
                        mapFile = mapFile,
                        evalDist = evalDist,
                        grid = grid,
                        verbose = verbose)
  })
  ## Concatenate results.
  crossIBD <- do.call(what = `c`, args = crossIBD)
  ## Replace cross names by cross names provided in input.
  if (nCross == 1) {
    covar <- data.frame(cross = rep(x = crossNames,
                                    times = nrow(crossIBD$markers)),
                        row.names = rownames(crossIBD$markers))
  } else {
    crossDat <- data.frame(cross = paste0("cross", seq_len(nCross)),
                           crossName = crossNames)
    covar <- attr(crossIBD, "genoCross")
    covar[["cross"]] <-
      crossDat[["crossName"]][match(x = covar[["cross"]],
                                    table = crossDat[["cross"]])]
    row.names(covar) <- covar[["geno"]]
    covar <- covar[, colnames(covar) != "geno", drop = FALSE]
  }
  ## Bind phenotypic data.
  if (inherits(pheno, "list")) {
    phenoTot <- do.call(what = rbind, args = pheno)
  } else {
    phenoTot <- pheno
  }
  ## Get number of parents.
  parents <- crossIBD$parents
  nPar <- length(parents)
  ## Construct empty marker matrix.
  if (packageVersion("statgenIBD") <= "1.0.4") {
    ## Get marker names.
    markerNames <- rownames(crossIBD$markers)
    markers <- array(NA_real_, dim = c(dim(crossIBD$markers)[c(2, 1)], nPar),
                     dimnames = c(dimnames(crossIBD$markers)[c(2, 1)],
                                  list(parents)))
  } else {
    ## Get marker names.
    markerNames <- colnames(crossIBD$markers)
    markers <- array(NA_real_, dim = c(dim(crossIBD$markers)[1:2], nPar),
                     dimnames = c(dimnames(crossIBD$markers)[1:2],
                                  list(parents)))
  }
  ## Fill marker matrix.
  for (i in seq_along(markerNames)) {
    markers[, i, ] <- markers3DtoMat(markers = crossIBD$markers,
                                     parents = parents,
                                     markerSel = markerNames[i])
  }
  ## Read original map.
  mapOrig <- read.delim(mapFile, header = FALSE)
  rownames(mapOrig) <- mapOrig[[1]]
  mapOrig <- mapOrig[, -1]
  colnames(mapOrig) <- c("chr", "pos")
  ## Create gDataMPP object.
  MPPobj <- createGDataMPPInternal(geno = markers,
                                   map = crossIBD$map,
                                   pheno = phenoTot,
                                   covar = covar)
  attr(x = MPPobj, which = "popType") <- crossIBD$popType
  attr(x = MPPobj, which = "pedigree") <- crossIBD$pedigree
  attr(x = MPPobj, which = "genoCross") <- attr(x = crossIBD, which = "genoCross")
  attr(x = MPPobj, which = "mapOrig") <- mapOrig
  return(MPPobj)
}

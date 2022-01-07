#' Multi round genome scans for QTL detection
#'
#' Multi round genome scans for QTL detection.\cr\cr
#' Several rounds of QTL detection are performed. First a model is fitted
#' without cofactors. If for at least one marker the \eqn{-10log(p)} value is
#' above the threshold the marker with the lowest p-Value is added as cofactor
#' in the next round of QTL detection. This process continues until there are
#' no new markers with a \eqn{-10log(p)} value above the threshold or until
#' the maximum number of cofactors is reached.
#'
#' @param MPPobj An object of class gData, typically the output of either
#' \code{\link{calcIBDmpp}} or \code{\link{readRABBIT}}.
#' @param trait A character string indicating the trait QTL mapping is done for.
#' @param QTLwindow A numerical value indicating the window around a QTL that
#' is considered as part of that QTL.
#' @param threshold A numerical value indicating the threshold for the
#' \eqn{-10logp} value of a marker to be considered a QTL.
#' @param CIM Should Composite Interval Mapping be done? If \code{FALSE} only
#' one round of QTL mapping is done without cofactors.
#' @param maxCofactors A numerical value, the maximum number of cofactor to
#' include in the model.
#' @param verbose Should progress and intermediate plots be output?
#'
#' @export
selQTLmpp <- function(MPPobj,
                      trait = NULL,
                      QTLwindow = 10,
                      threshold = 3,
                      CIM = TRUE,
                      maxCofactors = 5,
                      verbose = FALSE) {
  if (!inherits(MPPobj, "gData")) {
    stop("MPPobj should be an object of class gData.\n")
  }
  map <- MPPobj$map
  markers <- MPPobj$markers
  pheno <- MPPobj$pheno[[1]]
  covar <- MPPobj$covar
  if (is.null(map)) {
    stop("MPP object should contain a map.\n")
  }
  if (is.null(markers)) {
    stop("MPP object should contain a marker matrix.\n")
  }
  if (length(dim(markers)) != 3) {
    stop("markers should be a 3D array with IBD probabilities.\n")
  }
  if (is.null(pheno)) {
    stop("MPP object should contain phenotypic data.\n")
  }
  if (!is.character(trait) || length(trait) > 1 ||
      !hasName(x = pheno, name = trait)) {
    stop("trait should be a character string of length one present in pheno.\n")
  }
  if (!is.numeric(QTLwindow) || length(QTLwindow) > 1 || QTLwindow < 0) {
    stop("QTLwindow should be a positive numerical value.\n")
  }
  if (!is.numeric(threshold) || length(threshold) > 1 || threshold < 0) {
    stop("threshold should be a positive numerical value.\n")
  }
  if (!is.numeric(maxCofactors) || length(maxCofactors) > 1 ||
      maxCofactors < 0) {
    stop("maxCofactors should be a positive numerical value.\n")
  }
  parents <- dimnames(markers)[[3]]
  nPar <- length(parents)
  markerNames <- dimnames(markers)[[2]]
  map <- MPPobj$map
  ## Construct model data by merging phenotypic and genotypic data.
  ## Merge phenotypic data and covar (cross).
  modDat <- merge(pheno, covar, by.x = "genotype", by.y = "row.names")
  ## Flatten markers to 2D structure.
  markers <- do.call(cbind, apply(X = markers, MARGIN = 2,
                                  FUN = I, simplify = FALSE))
  colnames(markers) <- paste0(rep(markerNames, each = nPar), "_", parents)
  ## Merge markers to modDat.
  modDat <- merge(modDat, markers, by.x = "genotype", by.y = "row.names")
  ## For Simple Interval Mapping do only 1 iteration.
  if (!CIM) {
    maxCofactors <- 0
  }
  ## Initialize parameters.
  cofactors <- NULL
  while (length(cofactors) <= maxCofactors) {
    if (verbose) {
      cat(paste0("QTL scan for trait ", trait, ", ",
                 length(cofactors), " cofactors\n"))
    }
    scanRes <- scanQTL(modDat = modDat,
                       map = map,
                       parents = parents,
                       trait = trait,
                       QTLwindow = QTLwindow,
                       cof = cofactors,
                       verbose = verbose)
    if (verbose) {
      plotIntermediateScan(scanRes,
                           threshold = threshold,
                           cofactors = cofactors,
                           trait = trait)
    }
    ## Restrict to markers outside 'known' QTLRegions.
    scanSel <- scanRes[!scanRes[["QTLRegion"]] &
                         !is.na(scanRes[["minlog10p"]]), ]
    minlog10pMax <- max(scanSel[["minlog10p"]])
    if (minlog10pMax < threshold) {
      cofactors <- sort(cofactors)
      break
    }
    ## Add new cofactor to list of cofactor for next round of scanning.
    cofactors <- c(cofactors, scanSel[which.max(scanSel[["minlog10p"]]), "snp"])
  }
  ## Construct GWAResult and signSnp
  colnames(scanRes)[colnames(scanRes) == "minlog10p"] <- "LOD"
  GWARes <- scanRes[ , colnames(scanRes) != "QTLRegion"]
  signSnp <- scanRes[scanRes[["QTLRegion"]], colnames(scanRes) != "QTLRegion"]
  signSnp[["snpStatus"]] <-
    as.factor(ifelse(signSnp[["snp"]] %in% cofactors,
                     "significant SNP",
                     "within QTL window of significant SNP"))
  ## Construct GWASInfo
  GWASInfo <- list(parents = parents)
  res <- createGWAS(GWAResult = list(pheno = GWARes),
                    signSnp = list(pheno = signSnp),
                    thr = list(pheno = setNames(threshold, trait)),
                    GWASInfo = GWASInfo)
  ## Add QTLmpp class to simplify providing generic functions.
  class(res) <- c("QTLmpp", class(res))
  return(res)
}

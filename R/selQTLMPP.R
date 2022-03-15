#' Multi round genome scans for QTL detection
#'
#' Multi round genome scans for QTL detection.\cr\cr
#' Several rounds of QTL detection are performed. First a model is fitted
#' without cofactors. If for at least one marker the \eqn{-log10(p)} value is
#' above the threshold the marker with the lowest p-Value is added as cofactor
#' in the next round of QTL detection. This process continues until there are
#' no new markers with a \eqn{-log10(p)} value above the threshold or until
#' the maximum number of cofactors is reached.
#'
#' @param MPPobj An object of class gDataMPP, typically the output of either
#' \code{\link{calcIBDMPP}} or \code{\link{readRABBIT}}.
#' @param trait A character string indicating the trait QTL mapping is done for.
#' @param QTLwindow A numerical value indicating the window around a QTL that
#' is considered as part of that QTL.
#' @param threshold A numerical value indicating the threshold for the
#' \eqn{-log10p} value of a marker to be considered a QTL.
#' @param maxCofactors A numerical value, the maximum number of cofactors to
#' include in the model. If \code{NULL} cofactors are added until no new
#' cofactors are found.
#' @param verbose Should progress and intermediate plots be output?
#'
#' @return An object of class \code{QTLMPP}
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
#' ## Simple Interval Mapping.
#' ABC_SIM <- selQTLMPP(ABC, trait = "pheno", maxCofactors = 0)
#'
#' ## Composite Interval Mapping
#' \dontrun{
#' ABC_CIM <- selQTLMPP(ABC, trait = "pheno")
#' }
#'
#' @importFrom utils head tail
#' @export
selQTLMPP <- function(MPPobj,
                      trait = NULL,
                      QTLwindow = 10,
                      threshold = 3,
                      maxCofactors = NULL,
                      verbose = FALSE) {
  if (!inherits(MPPobj, "gDataMPP")) {
    stop("MPPobj should be an object of class gDataMPP.\n")
  }
  map <- MPPobj$map
  markers <- MPPobj$markers
  pheno <- MPPobj$pheno[[1]]
  covar <- MPPobj$covar
  mapOrig <- attr(x = MPPobj, which = "mapOrig")
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
  if (!is.null(maxCofactors) &&
      (!is.numeric(maxCofactors) || length(maxCofactors) > 1 ||
       maxCofactors < 0)) {
    stop("maxCofactors should be a positive numerical value.\n")
  }
  if (is.null(maxCofactors)) {
    maxCofactors <- dim(markers)[2]
  }
  parents <- dimnames(markers)[[3]]
  nPar <- length(parents)
  markerNames <- dimnames(markers)[[2]]
  map <- MPPobj$map
  ## Construct model data by merging phenotypic and genotypic data.
  ## Merge phenotypic data and covar (cross).
  modDat <- merge(pheno, covar, by.x = "genotype", by.y = "row.names")
  ## Remove missing values for trait from modDat.
  ## Not strictly necessary, but it prevents warnings from LMMsolve later on.
  modDat <- droplevels(modDat[!is.na(modDat[[trait]]), ])
  ## Flatten markers to 2D structure.
  markers <- do.call(cbind, apply(X = markers, MARGIN = 2,
                                  FUN = I, simplify = FALSE))
  colnames(markers) <- paste0(rep(markerNames, each = nPar), "_", parents)
  ## Merge markers to modDat.
  modDat <- merge(modDat, markers, by.x = "genotype", by.y = "row.names")
  ## Initialize parameters.
  cofactors <- NULL
  mapScan <- map
  while (length(cofactors) < maxCofactors) {
    if (verbose) {
      cat(paste0("QTL scan for trait ", trait, ", ",
                 length(cofactors), " cofactors\n"))
    }
    scanRes <- scanQTL(modDat = modDat,
                       map = mapScan,
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
    ## Add new cofactor to list of cofactors for next round of scanning.
    cofactors <- c(cofactors, scanSel[which.max(scanSel[["minlog10p"]]), "snp"])
    qtlWinScan <- scanRes[scanRes[["QTLRegion"]] &
                             !scanRes[["snp"]] %in% cofactors, "snp"]
    if (length(cofactors) == 1) {
      lowPScan <- scanRes[!is.na(scanRes[["minlog10p"]]) &
                            scanRes[["minlog10p"]] < 0.31, "snp"]
    }
    mapScan <- mapScan[!rownames(mapScan) %in% c(lowPScan, qtlWinScan), ]
  }
  ## Final run with all markers and all cofactors.
  scanRes <- scanQTL(modDat = modDat,
                     map = map,
                     parents = parents,
                     trait = trait,
                     QTLwindow = QTLwindow,
                     cof = cofactors,
                     verbose = verbose)
  ## Fit model with all cofactors for computation of variance explained.
  finMod <- randomQTLmodel(modDat = modDat, map = map, parents = parents,
                           trait = trait, scanMrk = NULL, cofMrk = cofactors,
                           NULLmodel = TRUE)
  ## Get weighted residual error.
  crossN <- table(modDat[["cross"]])
  crossResErr <- tail(finMod$VarDf, length(crossN))
  resErr <- sum((crossN * crossResErr[["Variance"]])) / sum(crossN)
  ## Construct data.frame with explained variances.
  varQTL <- head(finMod$VarDf, length(cofactors))
  varAllQTLs <- sum(varQTL[["Variance"]])
  varQTL[["varExpl"]] <- varQTL[["Variance"]] / (varAllQTLs + resErr)
  varQTL <- varQTL[c("VarComp", "varExpl")]
  ## Construct GWAResult and signSnp
  colnames(scanRes)[colnames(scanRes) == "minlog10p"] <- "LOD"
  GWARes <- scanRes[ , colnames(scanRes) != "QTLRegion"]
  signSnp <- scanRes[scanRes[["QTLRegion"]], colnames(scanRes) != "QTLRegion"]
  ## Add explained variance.
  signSnp <- merge(signSnp, varQTL, by.x = "snp", by.y = "VarComp",
                   all.x = TRUE, sort = FALSE)
  ## Add nearest real marker.
  signSnp[["mrkNear"]] <- sapply(X = seq_len(nrow(signSnp)), FUN = function(i) {
    mrkChr <- signSnp[i, "chr"]
    mrkPos <- signSnp[i, "pos"]
    mapOrigChr <- mapOrig[mapOrig[["chr"]] == mrkChr, ]
    mrkNear <- mapOrigChr[which.min(abs(mapOrigChr[["pos"]] - mrkPos)), ]
    return(rownames(mrkNear))
  })
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
  ## Add QTLMPP class to simplify providing generic functions.
  class(res) <- c("QTLMPP", class(res))
  return(res)
}

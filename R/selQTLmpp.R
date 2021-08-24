#' multi-round genome scans and select cofactors
#'
#' @param MPPobj ...
#' @param trait ...
#' @param QTLwindow ...
#' @param threshold ...
#' @param CIM ...
#' @param maxIter ...
#' @param maxCofactors ...
#' @param verbose ...
#'
#' @export
selQTLmpp <- function(MPPobj,
                      trait = NULL,
                      QTLwindow = 10,
                      threshold = 3,
                      CIM = TRUE,
                      maxIter = 100,
                      maxCofactors = 5,
                      verbose = FALSE) {
  parents <- dimnames(MPPobj$markers)[[3]]
  nPar <- length(parents)
  markerNames <- dimnames(MPPobj$markers)[[2]]
  ## Construct model data by merging phenotypic and genotypic data.
  ## Merge phenotypic data and covar (cross).
  modDat <- merge(MPPobj$pheno[[1]], MPPobj$covar,
                  by.x = "genotype", by.y = "row.names")
  ## Flatten markers to 2D structure.
  markers <- do.call(cbind, apply(X = MPPobj$markers, MARGIN = 2,
                                  FUN = I, simplify = FALSE))
  colnames(markers) <- paste0(rep(markerNames, each = nPar), "_", parents)
  ## Merge markers to modDat.
  modDat <- merge(modDat, markers, by.x = "genotype", by.y = "row.names")
  ## Initialize cofactors.
  cofactors <- NULL
  ## For Simple Interval Mapping do only 1 iteration.
  if (!CIM) maxIter <- 1
  for (i in seq_len(maxIter)) {
    if (verbose) {
      cat(paste0("QTL scan for trait ", trait, ", ",
                 length(cofactors), " cofactors\n"))
    }
    scanRes <- scanQTL(modDat = modDat,
                       map = MPPobj$map,
                       parents = parents,
                       trait = trait,
                       QTLwindow = QTLwindow,
                       cof = cofactors)
    plotQTLscan(scanRes,
                threshold = threshold,
                cofactors = cofactors,
                trait = trait)
    ## Restrict to markers outside 'known' QTLRegions.
    scanSel <- scanRes[!scanRes[["QTLRegion"]] &
                             !is.na(scanRes[["minlog10p"]]), ]
    minlog10pMax <- max(scanSel[["minlog10p"]])
    if (minlog10pMax < threshold) {
      cofactors <- sort(cofactors)
      break
    }
    ## Add new cofactor to list of cofactor for next round of scanning.
    cofactors <- c(cofactors, scanSel[which.max(scanSel[["minlog10p"]]), "ndx"])
    if (length(cofactors) > maxCofactors) break
  }
  results <- list(QTLcandidates = cofactors, scanResults = scanRes)
  MPPobj[["Result"]] <- results
  return(MPPobj)
}

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
                      trait = "pheno",
                      QTLwindow = 10,
                      threshold = 3,
                      CIM = TRUE,
                      maxIter = 100,
                      maxCofactors = 5,
                      verbose = FALSE) {
  modDat <- MPPobj$calcIBDres$IBDdata
  ## Initialize cofactors.
  cofactors <- NULL
  ## For Simple Interval Mapping do only 1 iteration.
  if (!CIM) maxIter <- 1
  for (i in seq_len(maxIter)) {
    if (verbose) {
      cat(paste0("QTL scan for trait ", trait, ", ",
                 length(cofactors), " cofactors\n"))
    }
    scanRes <- scanQTL(MPPobj, QTLwindow, cof = cofactors, trait)
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

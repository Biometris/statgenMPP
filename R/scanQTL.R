#' one-round genome scan with specified cofactor(s)
#'
#' @importFrom stats coef pchisq
#' @keywords internal
scanQTL <- function(modDat,
                    map,
                    parents,
                    QTLwindow = 10,
                    cof = NULL,
                    trait = NULL,
                    maxIter = 100,
                    verbose = FALSE) {
  ## Get info from MPP object.
  nMarkers <- nrow(map)
  ## Fit NULL model without cofactors.
  fitNULLmod <- randomQTLmodel(modDat = modDat,
                               map = map,
                               parents = parents,
                               trait = trait,
                               cofMrk = NULL,
                               NULLmodel = TRUE)
  ## Initialize output.
  tst <- foreach::foreach(i = seq_len(nMarkers)) %dopar% {
    scanMrk <- rownames(map)[i]
    cofMrk <- selectCofactors(map = map,
                              marker = i,
                              cofactors = cof,
                              QTLwindow = QTLwindow)
    ## Fit model for current marker.
    fitModMrk <- randomQTLmodel(modDat = modDat,
                                map = map,
                                parents = parents,
                                trait = trait,
                                scanMrk = scanMrk,
                                cofMrk = cofMrk,
                                NULLmodel = FALSE)
    dev <- 2.0 * fitModMrk$logL - 2.0 * fitNULLmod$logL
    ## Refit NULL model only if cofactors are present for current marker.
    if (!is.null(cofMrk)) {
      fitModCof <- randomQTLmodel(modDat = modDat,
                                  map = map,
                                  parents = parents,
                                  trait = trait,
                                  cofMrk = cofMrk,
                                  NULLmodel = TRUE)
      dev <- 2.0 * fitModMrk$logL - 2.0 * fitModCof$logL
    }
    list(length(cofMrk) != length(cof), # QTLRegion
         -log10(0.5 * pchisq(dev, 1, lower.tail = FALSE)), # minlog10p
         coef(fitModMrk)[[rownames(map)[i]]]) #effects
  }
  QTLRegion <- sapply(X = tst, FUN = `[[`, 1)
  minlog10p <- sapply(X = tst, FUN = `[[`, 2)
  effects <- do.call(rbind, args = lapply(X = tst, FUN = `[[`, 3))
  dimnames(effects) <- list(rownames(map), paste0("eff_", parents))
  scanRes <- data.table::data.table(trait = trait,
                                    snp = rownames(map),
                                    map,
                                    pValue = 10 ^ -minlog10p,
                                    effects,
                                    minlog10p = minlog10p,
                                    QTLRegion = QTLRegion)
  return(scanRes)
}


#' one-round genome scan with specified cofactor(s)
#'
#' @importFrom stats coef pchisq model.matrix
#' @keywords internal
scanQTL <- function(modDat,
                    map,
                    markers,
                    parents,
                    QTLwindow = 10,
                    cof = NULL,
                    trait = NULL,
                    maxIter = 100,
                    parallel = FALSE,
                    verbose = FALSE) {
  ## Get info from input.
  nPar <- length(parents)
  nGeno <- nrow(markers)
  nCross <- length(unique(modDat[["cross"]]))
  ## Create general model matrices.
  y <- modDat[[trait]]
  X <- spam::as.spam(model.matrix(if (nCross > 1) ~cross else ~ 1,
                                  data = modDat))
  lRinv <- constructRinv(modDat, residual = ~cross, weights = 1)
  ## Fit NULL model with all cofactors.
  Z0 <- do.call(spam::cbind.spam, lapply(X = cof, FUN = function(mrk) {
    markers[, mrk, ]
  }))
  lGinv0 <- lapply(X = seq_along(cof), FUN = function(i) {
    spam::diag.spam(x = c(rep(0, nPar * (i - 1)),
                          rep(1, nPar),
                          rep(0, nPar * (length(cof) - i))))
  })
  fitModNULL <- sparseMixedModels(y = y, X = X, Z = Z0,
                                  lRinv = lRinv, lGinv = lGinv0,
                                  tolerance = 1e-3)
  ## Fit models for each marker.
  `%op%` <- getOper(parallel && foreach::getDoParRegistered())
  scanFull <- foreach::foreach(scanMrk = rownames(map)) %op% {
    cofMrk <- selectCofactors(map = map,
                              marker = scanMrk,
                              cofactors = cof,
                              QTLwindow = QTLwindow)
    selMrk <- c(scanMrk, cofMrk)
    Z <- do.call(spam::cbind.spam, lapply(X = selMrk, FUN = function(mrk) {
      markers[, mrk, ]
    }))
    lGinv <- lapply(X = seq_along(selMrk), FUN = function(i) {
      spam::diag.spam(x = c(rep(0, nPar * (i - 1)),
                            rep(1, nPar),
                            rep(0, nPar * (length(selMrk) - i))))
    })
    ## Fit model for current marker.
    fitModMrk <- sparseMixedModels(y = y, X = X, Z = Z,
                                   lRinv = lRinv, lGinv = lGinv,
                                   tolerance = 1e-3)
    ## Compute change in deviance.
    dev <- 2 * fitModMrk$logL - 2 * fitModNULL$logL
    ## Refit NULL model only if cofactors differ for current marker.
    if (length(cofMrk) != length(cof)) {
      ## Fit NULL model with restricted cofactors.
      Z1 <- do.call(spam::cbind.spam, lapply(X = cofMrk, FUN = function(mrk) {
        markers[, mrk, ]
      }))
      lGinv1 <- lapply(X = seq_along(cofMrk), FUN = function(i) {
        spam::diag.spam(x = c(rep(0, nPar * (i - 1)),
                              rep(1, nPar),
                              rep(0, nPar * (length(cofMrk) - i))))
      })
      fitModCof <- sparseMixedModels(y = y, X = X, Z = Z1,
                                     lRinv = lRinv, lGinv = lGinv1,
                                     tolerance = 1e-3)
      dev <- 2 * fitModMrk$logL - 2 * fitModCof$logL
    }
    list(length(cofMrk) != length(cof), # QTLRegion
         -log10(0.5 * pchisq(dev, 1, lower.tail = FALSE)), # minlog10p
         fitModMrk$a[(1 + nCross):(length(parents) + nCross)]) #effects
  }
  QTLRegion <- sapply(X = scanFull, FUN = `[[`, 1)
  minlog10p <- sapply(X = scanFull, FUN = `[[`, 2)
  effects <- do.call(rbind, args = lapply(X = scanFull, FUN = `[[`, 3))
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


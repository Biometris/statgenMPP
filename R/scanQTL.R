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
                    KInv = NULL,
                    trait = NULL,
                    maxIter = 100,
                    parallel = FALSE,
                    verbose = FALSE) {
  ## Get info from input.
  nPar <- length(parents)
  nGeno <- nrow(modDat)
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
  ## Add genotype.
  if (!is.null(KInv)) {
    ZGeno <- spam::diag.spam(x = 1, nrow = nGeno, ncol = nGeno)
    Z0 <- cbind(Z0, ZGeno)
  }
  lGinv0 <- lapply(X = seq_along(cof), FUN = function(i) {
    spam::diag.spam(x = c(rep(0, nPar * (i - 1)),
                          rep(1, nPar),
                          rep(0, nPar * (length(cof) - i) +
                                if (!is.null(KInv)) nGeno else 0)))
  })
  ## Fit models for each marker.
  chrs <- unique(map[["chr"]])
  `%op%` <- getOper(parallel && foreach::getDoParRegistered())
  scanFull <- foreach::foreach(i = seq_along(chrs)) %op% {
    KInvChr <- KInv[[i]]
    if (!is.null(KInvChr)) {
      lGinv1 <- spam::bdiag.spam(spam::spam(x = 0, nrow = length(cof) * nPar,
                                            ncol = length(cof) * nPar),
                                 KInvChr)
      lGinv0 <- c(lGinv0, list(lGinv1))
    }
    ## Fit NULL model with all cofactors.
    fitModNULL <- sparseMixedModels(y = y, X = X, Z = Z0,
                                    lRinv = lRinv, lGinv = lGinv0,
                                    tolerance = 1e-3)
    chrMrk <- rownames(map)[map[["chr"]] == chrs[i]]
    QTLRegion <- setNames(logical(length = length(chrMrk)), chrMrk)
    minlog10p <- setNames(numeric(length = length(chrMrk)), chrMrk)
    effects <- matrix(nrow = length(chrMrk), ncol = nPar,
                      dimnames = list(chrMrk, parents))
    for (scanMrk in chrMrk) {
      ## Get cofactors for current markers.
      cofMrk <- selectCofactors(map = map,
                                marker = scanMrk,
                                cofactors = cof,
                                QTLwindow = QTLwindow)
      selMrk <- c(scanMrk, cofMrk)
      Z <- do.call(spam::cbind.spam, lapply(X = selMrk, FUN = function(mrk) {
        markers[, mrk, ]
      }))
      if (!is.null(KInvChr)) {
        Z <- cbind(Z, ZGeno)
      }
      lGinv <- lapply(X = seq_along(selMrk), FUN = function(i) {
        spam::diag.spam(x = c(rep(0, nPar * (i - 1)),
                              rep(1, nPar),
                              rep(0, nPar * (length(selMrk) - i) +
                                    if (!is.null(KInvChr)) nGeno else 0)))
      })
      if (!is.null(KInvChr)) {
        lGinv1 <- spam::bdiag.spam(spam::spam(x = 0, nrow = length(selMrk) * nPar,
                                              ncol = length(selMrk) * nPar),
                                   KInvChr)
        lGinv <- c(lGinv, list(lGinv1))
      }
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
        if (!is.null(KInvChr)) {
          Z1 <- cbind(Z1, ZGeno)
        }
        lGinv1 <- lapply(X = seq_along(cofMrk), FUN = function(i) {
          spam::diag.spam(x = c(rep(0, nPar * (i - 1)),
                                rep(1, nPar),
                                rep(0, nPar * (length(cofMrk) - i) +
                                      if (!is.null(KInvChr)) nGeno else 0)))
        })
        if (!is.null(KInvChr)) {
          lGinv1a <- spam::bdiag.spam(spam::spam(x = 0,
                                                 nrow = length(cofMrk) * nPar,
                                                 ncol = length(cofMrk) * nPar),
                                      KInvChr)
          lGinv1 <- c(lGinv1, list(lGinv1a))
        }
        fitModCof <- sparseMixedModels(y = y, X = X, Z = Z1,
                                       lRinv = lRinv, lGinv = lGinv1,
                                       tolerance = 1e-3)
        dev <- 2 * fitModMrk$logL - 2 * fitModCof$logL
      }
      QTLRegion[scanMrk] <- length(cofMrk) != length(cof)
      minlog10p[scanMrk] <- min(-log10(0.5 * pchisq(dev, 1,
                                                    lower.tail = FALSE)), 300)
      effects[scanMrk, ] <- fitModMrk$a[(1 + nCross):(length(parents) + nCross)]
    }
    list(QTLRegion, minlog10p, effects)
  }
  QTLRegion <- do.call(c, lapply(X = scanFull, FUN = `[[`, 1))
  minlog10p <- do.call(c, lapply(X = scanFull, FUN = `[[`, 2))
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


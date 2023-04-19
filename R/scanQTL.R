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
                    Usc = NULL,
                    trait = NULL,
                    maxIter = 100,
                    se = FALSE,
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
  ## Fit models for each marker.
  chrs <- unique(map[["chr"]])
  `%op%` <- getOper(parallel && foreach::getDoParRegistered())
  scanFull <- foreach::foreach(i = seq_along(chrs)) %op% {
    UscChr <- Usc[[i]]
    nEigChr <- ncol(UscChr)
    ## Add genotype.
    lGinv0 <- lapply(X = seq_along(cof), FUN = function(i) {
      spam::diag.spam(x = c(rep(0, nPar * (i - 1)),
                            rep(1, nPar),
                            rep(0, nPar * (length(cof) - i) +
                                  if (!is.null(Usc)) nEigChr else 0)))
    })
    if (!is.null(UscChr)) {
      ZChr <- spam::cbind.spam(Z0, UscChr)
      lGinv1 <- spam::diag.spam(x = c(rep(0, times = length(cof) * nPar),
                                          rep(1, nEigChr)))
      LGinvChr <- c(lGinv0, list(lGinv1))
    } else {
      ZChr <- Z0
      LGinvChr <- lGinv0
    }
    ## Fit NULL model with all cofactors.
    fitModNULL <- sparseMixedModels(y = y, X = X, Z = ZChr,
                                    lRinv = lRinv, lGinv = LGinvChr,
                                    tolerance = 1e-3)
    chrMrk <- rownames(map)[map[["chr"]] == chrs[i]]
    QTLRegion <- setNames(logical(length = length(chrMrk)), chrMrk)
    minlog10p <- setNames(numeric(length = length(chrMrk)), chrMrk)
    effects <- se_effects <- matrix(nrow = length(chrMrk), ncol = nPar,
                                    dimnames = list(chrMrk, parents))
    for (scanMrk in chrMrk) {
      ## Get cofactors for current markers.
      cofMrk <- selectCofactors(map = map,
                                marker = scanMrk,
                                cofactors = cof,
                                QTLwindow = QTLwindow)
      selMrk <- c(scanMrk, cofMrk)
      ZMrk <- do.call(spam::cbind.spam, lapply(X = selMrk, FUN = function(mrk) {
        markers[, mrk, ]
      }))
      lGinvMrk <- lapply(X = seq_along(selMrk), FUN = function(i) {
        spam::diag.spam(x = c(rep(0, nPar * (i - 1)),
                              rep(1, nPar),
                              rep(0, nPar * (length(selMrk) - i) +
                                    if (!is.null(UscChr)) nEigChr else 0)))
      })
      if (!is.null(UscChr)) {
        ZMrk <- spam::cbind.spam(ZMrk, UscChr)
        lGinv1Mrk <- spam::diag.spam(x = c(rep(0, times = length(selMrk) * nPar),
                                        rep(1, nEigChr)))
        lGinvMrk <- c(lGinvMrk, list(lGinv1Mrk))
      }
      ## Fit model for current marker.
      fitModMrk <- sparseMixedModels(y = y, X = X, Z = ZMrk,
                                     lRinv = lRinv, lGinv = lGinvMrk,
                                     tolerance = 1e-3)
      ## Compute change in deviance.
      dev <- 2 * fitModMrk$logL - 2 * fitModNULL$logL
      ## Refit NULL model only if cofactors differ for current marker.
      if (length(cofMrk) != length(cof)) {
        ## Fit NULL model with restricted cofactors.
        ZMrk2 <- do.call(spam::cbind.spam, lapply(X = cofMrk,
                                                  FUN = function(mrk) {
                                                    markers[, mrk, ]
                                                  }))
        lGinvMrk2 <- lapply(X = seq_along(cofMrk), FUN = function(i) {
          spam::diag.spam(x = c(rep(0, nPar * (i - 1)),
                                rep(1, nPar),
                                rep(0, nPar * (length(cofMrk) - i) +
                                      if (!is.null(UscChr)) nEigChr else 0)))
        })
        if (!is.null(UscChr)) {
          ZMrk2 <- cbind(ZMrk2, UscChr)
          lGinv1Mrk2 <- spam::diag.spam(x = c(rep(0, times = length(cofMrk) * nPar),
                                             rep(1, nEigChr)))
          lGinvMrk2 <- c(lGinvMrk2, list(lGinv1Mrk2))
        }
        fitModCof <- sparseMixedModels(y = y, X = X, Z = ZMrk2,
                                       lRinv = lRinv, lGinv = lGinvMrk2,
                                       tolerance = 1e-3)
        dev <- 2 * fitModMrk$logL - 2 * fitModCof$logL
      }
      QTLRegion[scanMrk] <- length(cofMrk) != length(cof)
      minlog10p[scanMrk] <- min(-log10(0.5 * pchisq(dev, 1,
                                                    lower.tail = FALSE)), 300)
      effects[scanMrk, ] <- fitModMrk$a[(1 + nCross):(length(parents) + nCross)]
      if (se) {
        se_tot <- calcStandardErrors(C = fitModMrk$C,
                                     D = spam::diag.spam(x = 1, nrow = nrow(fitModMrk$C)))
        se_effects[scanMrk, ] <- se_tot[(1 + nCross):(length(parents) + nCross)]
      }
    }
    list(QTLRegion, minlog10p, effects, se_effects)
  }
  QTLRegion <- do.call(c, lapply(X = scanFull, FUN = `[[`, 1))
  minlog10p <- do.call(c, lapply(X = scanFull, FUN = `[[`, 2))
  effects <- do.call(rbind, args = lapply(X = scanFull, FUN = `[[`, 3))
  se_effects <- do.call(rbind, args = lapply(X = scanFull, FUN = `[[`, 4))
  dimnames(effects) <- list(rownames(map), paste0("eff_", parents))
  dimnames(se_effects) <- list(rownames(map), paste0("se_eff_", parents))
  scanRes <- data.table::data.table(trait = trait,
                                    snp = rownames(map),
                                    map,
                                    pValue = 10 ^ -minlog10p,
                                    effects,
                                    se_effects,
                                    minlog10p = minlog10p,
                                    QTLRegion = QTLRegion)
  return(scanRes)
}


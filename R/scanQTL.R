#' one-round genome scan with specified cofactor(s)
#'
#' @importFrom stats coef pchisq
#' @keywords internal
scanQTL <- function(MPPobj,
                    QTLwindow = 10,
                    cof = NULL,
                    trait = "pheno",
                    verbose = FALSE) {
  ## Get info from MPP object.
  map <- MPPobj$calcIBDres$map
  parents <- MPPobj$MPPinfo$parents
  nMarkers <- nrow(map)
  ## Fit NULL model without cofactors.
  fitNULLmod <- randomQTLmodel(MPPobj,
                               trait,
                               cofPos = NULL,
                               NULLmodel = TRUE)
  ## Initialize output.
  minlog10p <- numeric(length = nMarkers)
  QTLRegion <- rep(FALSE, nMarkers)
  effects <- matrix(nrow = nMarkers, ncol = length(parents),
                    dimnames = list(rownames(map), parents))
  for (i in seq_len(nMarkers)) {
    mrkCof <- selectCofactors(map = map,
                              marker = i,
                              cofactors = cof,
                              QTLwindow = QTLwindow)
    if (length(mrkCof) != length(cof)) {
      QTLRegion[i] <- TRUE
    }
    ## Fit model for current marker.
    fitModMrk <- randomQTLmodel(MPPobj,
                                trait,
                                scanPos = i,
                                cofPos = mrkCof,
                                NULLmodel = FALSE)
    dev <- 2.0 * fitModMrk$logL - 2.0 * fitNULLmod$logL
    ## Refit NULL model only if cofactors are present for current marker.
    if (!is.null(mrkCof)) {
      fitModCof <- randomQTLmodel(MPPobj,
                                  trait,
                                  cofPos = mrkCof,
                                  NULLmodel = TRUE)
      dev <- 2.0 * fitModMrk$logL - 2.0 * fitModCof$logL
    }
    minlog10p[i] <- -log10(0.5 * pchisq(dev, 1, lower.tail = FALSE))
    effects[i, ] <- coef(fitModMrk)[[rownames(map)[i]]]
    if (verbose && i %% 25 == 0) {
      cat(paste(i, "\n"))
    }
  }
  res <- data.frame(ndx = seq_len(nMarkers),
                    map,
                    minlog10p = minlog10p,
                    eff = effects,
                    QTLRegion = QTLRegion,
                    trait = trait)
  return(res)
}

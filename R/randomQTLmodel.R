#' one position fitted in the mixed model
#'
#' @importFrom stats as.formula
#' @keywords internal
randomQTLmodel <- function(MPPobj,
                           trait = "pheno",
                           scanPos = 1,
                           cofPos = NULL,
                           NULLmodel = FALSE) {
  modDat <- MPPobj$calcIBDres$IBDdata
  map <- MPPobj$calcIBDres$map
  parents <- MPPobj$MPPinfo$parents
  nCross <- length(unique(modDat[["cross"]]))
  nCof <- length(cofPos)
  if (nCross == 1) { # for MAGIC-type pop with one cross.
    fixed <- as.formula(paste(trait, "~1"))
  } else { # for NAM or diallel-type pops with > 1 cross
    fixed <- as.formula(paste(trait, "~cross"))
  }
  modDat[["cross"]] <- as.factor(modDat[["cross"]])
  selPos <- c(cofPos, if (!NULLmodel) scanPos)
  Lgrp <- list()
  if (length(selPos) > 0) {
    for (i in seq_along(selPos)) {
      selName <- rownames(map)[selPos[i]]
      selIBDNames <- paste0(selName, "_", parents)
      Lgrp[[selName]] <- which(colnames(modDat) %in% selIBDNames)
    }
  }
  fitMod <- LMMsolver::LMMsolve(fixed = fixed,
                                randomMatrices = if (length(Lgrp) > 0) Lgrp,
                                residualterm = "cross",
                                data = modDat,
                                eps = 1.0e-8,
                                monitor = FALSE)
  return(fitMod)
}

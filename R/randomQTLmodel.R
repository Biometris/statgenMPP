#' one position fitted in the mixed model
#'
#' @importFrom stats as.formula
#' @keywords internal
randomQTLmodel <- function(modDat,
                           map,
                           parents,
                           trait = "pheno",
                           scanPos = 1,
                           cofPos = NULL,
                           NULLmodel = FALSE) {
  nCross <- length(unique(modDat[["cross"]]))
  nCof <- length(cofPos)
  if (nCross == 1) { # for MAGIC-type pop with one cross.
    fixed <- as.formula(paste(trait, "~1"))
  } else { # for NAM or diallel-type pops with > 1 cross
    fixed <- as.formula(paste(trait, "~cross"))
  }
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

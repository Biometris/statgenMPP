#' one position fitted in the mixed model
#'
#' @importFrom stats as.formula
#' @keywords internal
randomQTLmodel <- function(modDat,
                           map,
                           parents,
                           trait = "pheno",
                           scanMrk = NULL,
                           cofMrk = NULL,
                           NULLmodel = FALSE) {
  nCross <- length(unique(modDat[["cross"]]))
  nCof <- length(cofMrk)
  if (nCross == 1) { # for MAGIC-type pop with one cross.
    fixed <- as.formula(paste(trait, "~1"))
  } else { # for NAM or diallel-type pops with > 1 cross
    fixed <- as.formula(paste(trait, "~cross"))
  }
  selMrk <- c(cofMrk, if (!NULLmodel) scanMrk)
  Lgrp <- list()
  ## MB: add grp() terms to random part of model
  ranTerm <- NULL
  if (length(selMrk) > 0) {
    for (selName in selMrk) {
      selIBDNames <- paste0(selName, "_", parents)
      Lgrp[[selName]] <- which(colnames(modDat) %in% selIBDNames)
      if (is.null(ranTerm)) {
        ranTerm <- paste0(ranTerm, "~grp(`", selName, "`)")
      } else {
        ranTerm <- paste0(ranTerm, "+grp(`", selName, "`)")
      }
    }
  }
  if (!is.null(ranTerm)) {
    ranTerm = as.formula(ranTerm)
  }
  fitMod <- LMMsolver::LMMsolve(fixed = fixed,
                                random = ranTerm,
                                group = if (length(Lgrp) > 0) Lgrp,
                                residual = ~cross,
                                data = modDat,
                                tolerance = 1.0e-3)
  return(fitMod)
}

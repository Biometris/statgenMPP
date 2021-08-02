#' one position fitted in the mixed model
#'
#' @importFrom stats as.formula
#' @keywords internal
randomQTLmodel <- function(MPPobj,
                           trait.name = "pheno",
                           scan_pos = 1,
                           cof_pos = NULL,
                           NULLmodel = FALSE) {
  data <- MPPobj$calcIBDres$IBDdata
  map <- MPPobj$calcIBDres$map
  unipar.names <- MPPobj$calcIBDres$par.names
  npos <- nrow(map)
  cross <- "cross"
  if (length(unique(data$cross)) == 1) { # for MAGIC-type pop with one cross
    fix <- as.formula(paste0(trait.name, "~1"))
  }
  if (length(unique(data$cross)) > 1) {  # for NAM or diallel-type pops with > 1 cross
    fix <- as.formula(paste0(trait.name, "~", cross))
  }
  data$cross <- as.factor(data$cross)
  if (NULLmodel) { ##NULL mode without cofactor
    if (length(cof_pos) == 0) {
      obj <- LMMsolver::LMMsolve(fixed = fix,
                                 residualterm = cross,
                                 data = data,
                                 eps = 1.0e-8,
                                 monitor = FALSE,
                                 display = FALSE)
    }  ##NULL mode with cofactor(s)
    if (length(cof_pos) != 0) {
      Lgrp <- list()
      for (i in 1:length(cof_pos)) {
        cof_name <- rownames(map)[cof_pos[i]]
        cof_IBD_name <- paste0(cof_name, "_p", unipar.names)
        Lgrp[[cof_name]] <- which(colnames(data) %in% cof_IBD_name)
      }
      obj <- LMMsolver::LMMsolve(fixed = fix,
                                 randomMatrices = Lgrp,
                                 residualterm = cross,
                                 data = data,
                                 eps = 1.0e-8,
                                 monitor = FALSE)
    }
  }
  if (!NULLmodel) { ## FULL model
    sel_pos <- c(scan_pos, cof_pos)
    Lgrp <- list()
    for (i in 1:length(sel_pos)) {
      sel_name <- rownames(map)[sel_pos[i]]
      sel_IBD_name <- paste0(sel_name, "_p", unipar.names)
      Lgrp[[sel_name]] <- which(colnames(data) %in% sel_IBD_name)
    }
    obj <- LMMsolver::LMMsolve(fixed = fix,
                               randomMatrices = Lgrp,
                               residualterm = cross,
                               data = data,
                               eps = 1.0e-8,
                               monitor = FALSE)
  }
  return(obj)
}

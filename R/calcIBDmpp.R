#' IBDs calculation and prepare the design matrix
#'
#' @param crossNames ...
#' @param locFiles ...
#' @param quaDat ...
#' @param poptypes ...
#' @param mapFile ...
#' @param evaldist ...
#'
#' @importFrom utils read.table
#' @export
calcIBDmpp <- function(crossNames,
                       locFiles,
                       quaDat,
                       poptypes,
                       mapFile,
                       evaldist,
                       verbose = FALSE) {
  ## Get number of crosses
  nCross <- length(locFiles)
  ## Calculate IBD probabilities per cross.
  crossIBD <- lapply(X = seq_len(nCross), FUN = function(i) {
    if (verbose) {
      cat(paste0("calc IBD in cross: ", crossNames[i], ".\n"))
    }
    statgenIBD::calcIBD(poptype = poptypes[i],
                        locfile = locFiles[i],
                        mapfile = mapFile,
                        evaldist = evaldist)
  })
  ## Concatenate results.
  crossIBD <- do.call(what = `c`, args = crossIBD)
  ## Construct data.frame for further analyses.
  crossProbs <- statgenIBD::getProbs(IBDprob = crossIBD,
                                     markers = rownames(crossIBD$markers),
                                     sumProbs = TRUE)
  ## Replace cross names by cross names provided in input.
  if (nCross == 1) {
    crossProbs[["cross"]] <- crossNames
  } else {
    crossDat <- data.frame(cross = paste0("cross", seq_len(nCross)),
                           crossName = crossNames)
    crossProbs[["cross"]] <-
      crossDat[["crossName"]][match(x = crossProbs[["cross"]],
                                    table = crossDat[["cross"]])]
  }
  ## Bind phenotypic data.
  phenoTot <- do.call(what = rbind, args = quaDat)
  ## Merge phenotypic data to probabilities.
  IBDdata <- merge(phenoTot, crossProbs, by.x = "ID", by.y = "geno")
  ## Construct output opbject.
  MPPobj <- structure(list(MPPinfo = list(parents = crossIBD$parents,
                                          crossNames = crossNames,
                                          locFiles = locFiles,
                                          quaFiles = quaDat,
                                          poptypes = poptypes,
                                          mapfile = mapFile,
                                          evaldist = evaldist),
                           calcIBDres = list(parents = crossIBD$parents,
                                             map = crossIBD$map,
                                             IBDdata = IBDdata)),
                      class = c("IBDProbsMPP", "list"))
  return(MPPobj)
}

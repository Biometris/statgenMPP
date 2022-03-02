#' multi-round genome scans and select cofactors
#'
#' @param MPPobj ...
#' @param QTLwindow ...
#' @param threshold ...
#' @param trait.name ...
#' @param CIM ...
#'
#' @export
selQTLMPPOrig <- function(MPPobj,
                          QTLwindow = 10,
                          threshold = 3,
                          trait.name = "pheno",
                          CIM = TRUE) {
  data <- MPPobj$calcIBDres$IBDdata
  par.names <- MPPobj$calcIBDres$par.names
  map <- MPPobj$calcIBDres$map
  nloc <- nrow(map)
  cofactors <- NULL
  for (i in 1:100) {
    print(paste0("QTL scan for trait ", trait.name, ", ",
                 length(cofactors), " cofactors"))
    df.scan <- scanQTLOrig(MPPobj, QTLwindow, cof = cofactors, trait.name)
    plotQTLscanOrig(df.scan, threshold = threshold, cofactors,
                    trait.name = trait.name)
    # remove NAs from minlog10p values
    df.scan.sel <- df.scan[!df.scan[["QTLregion"]] &
                             !is.na(df.scan[["minlog10p"]]), ]
    ord <- order(df.scan.sel$minlog10p, decreasing = TRUE)
    max_value <- df.scan.sel[ord[1], ]$minlog10p
    if (max_value < threshold) {
      if (!is.null(cofactors)) cofactors <- sort(cofactors)
    }
    if (max_value < threshold) break
    if (!CIM) break
    cofactors <- c(cofactors, df.scan.sel[ord[1], ]$ndx)
    # if (length(cofactors)==5) break # if there are too many cofactors,
    # we can set where it should stop
  }
  results <- list(QTLcandidates = cofactors, scanResults = df.scan)
  MPPobj[["Result"]] <- results
  return(MPPobj)
}

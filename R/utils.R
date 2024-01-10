#' cofactor selection based on the QTL windowsize
#'
#' @keywords internal
calcDistance <- function(map,
                         m1,
                         m2) {
  dist <- ifelse(map[m1, "chr"] == map[m2, "chr"],
                 abs(map[m1, "pos"] - map[m2, "pos"]),
                 Inf)
  return(dist)
}

#' @keywords internal
selectCofactors <- function(map,
                            marker,
                            cofactors,
                            QTLwindow) {
  if (length(cofactors) == 0) return(NULL)
  minDist <- 0.5 * QTLwindow
  dist <- sapply(cofactors, function(x) calcDistance(map, marker, x))
  if (min(dist) <= minDist) {
    cofactors <- cofactors[-which.min(dist)]
    if (length(cofactors) == 0) return(NULL)
  }
  return(cofactors)
}

#' @importFrom stats aggregate
#'
#' @keywords internal
plotIntermediateScan <- function(dat,
                                 threshold,
                                 cofactors,
                                 trait) {
  ## Construct title.
  title <- paste("QTL-profile for trait ", trait)
  if (length(cofactors) == 0) {
    title <- paste0(title, ", no cofactors")
  } else if (length(cofactors) == 1) {
    title <- paste0(title, ", one cofactor")
  } else {
    title <- paste0(title,", ", length(cofactors)," cofactors")
  }
  map <- dat[c("chr", "pos")]
  map[["chr"]] <- factor(map[["chr"]], levels = unique(map[["chr"]]))
  ## Get the boundary for each of the chromosomes.
  ## Has to be converted to numeric to avoid integer overflow in the next step.
  chrBnd <- aggregate(x = map[["pos"]], by = list(map[["chr"]]),
                      FUN = function(p) {as.numeric(max(p))})
  ## Compute cumulative positions.
  addPos <- data.frame(chr = chrBnd[, 1],
                       add = c(0, cumsum(chrBnd[, 2]))[1:nrow(chrBnd)],
                       stringsAsFactors = FALSE)
  map <- merge(map, addPos, by = "chr")
  map[["cumPos"]] <- map[["pos"]] + map[["add"]]
  manhattanPlot(xValues = map[["cumPos"]],
                yValues = dat$minlog10p,
                plotType = "lines",
                map = map,
                yThr = threshold,
                title = title)
}

#' Helper function for creating summaries that always display NA.
#'
#' @noRd
#' @keywords internal
summaryNA <- function(dat) {
  if (!any(is.na(dat))) {
    return(c(summary(dat), "NA's" = 0))
  } else{
    return(summary(dat))
  }
}

#' Helper function for accessing parallel computing functions.
#'
#' @noRd
#' @importFrom foreach `%do%` `%dopar%`
#' @keywords internal
getOper <- function(x) {
  if (x) {
    `%dopar%`
  } else {
    `%do%`
  }
}

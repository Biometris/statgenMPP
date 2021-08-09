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
  if (min(dist) > minDist) {
    return(cofactors)
  } else {
    return(cofactors[-which.min(dist)])
  }
}

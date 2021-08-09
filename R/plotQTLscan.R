#' @importFrom graphics abline axis box
#' @importFrom stats aggregate
#'
#' @keywords internal
plotQTLscan <- function(dat,
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

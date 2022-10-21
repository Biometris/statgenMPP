#' Create a matrix plot
#'
#' Create a matrix plot of effect sizes and directions per parent for
#' significant snps.
#'
#' @import grDevices graphics
#'
#' @keywords internal
effectPlot <- function(effectDat,
                       signSnp,
                       map,
                       xLab = "Chromosomes",
                       yLab = "Parents",
                       chrBoundaries,
                       title = NULL,
                       trait = NULL,
                       ...,
                       output = TRUE) {
  ## Extract central chromosome positions from map.
  ## Differentiate cases to deal with character chromosomes.
  if (is.numeric(map$chr)) {
    chrs <- as.numeric(levels(as.factor(map$chr)))
  } else {
    chrs <- levels(as.factor(map$chr))
  }
  ## Recompute cumulative positions.
  addPos <- data.frame(chr = chrBoundaries[, 1],
                       add = c(0, cumsum(chrBoundaries[, 2] + 5))[1:nrow(chrBoundaries)],
                       stringsAsFactors = FALSE)
  map <- merge(map[, !colnames(map) %in% c("add", "cumPos")],
               addPos, by = "chr")
  map[["cumPos"]] <- map[["pos"]] + map[["add"]]
  map[["width"]] <- unlist(
    aggregate(x = map[["pos"]], by = list(map[["chr"]]),
              FUN = function(x) {
                diff0 <- diff(c(0, x), lag = 2) / 2
                return(c(diff0, tail(diff0, 1)))
              }, simplify = FALSE)[, 2])
  ## Compute postions of labels for chromosomes.
  xMarks <- aggregate(x = map[["cumPos"]], by = list(map[["chr"]]),
                      FUN = function(x) {
                        min(x) + (max(x) - min(x)) / 2
                      })[, 2]
  ## Compute chromosome boundaries.
  chrBnd <- c(0, aggregate(x = map[["cumPos"]],
                           by = list(map[["chr"]]), FUN = max)[[2]])
  chrBnd <- chrBnd + c(0, aggregate(x = map[["width"]], by = list(map[["chr"]]),
                                    FUN = tail, 1)[[2]])
  ## Add cumulative position from map to effects.
  parEffData <- merge(effectDat, map[, c("snp", "cumPos", "width")],
                      by = "snp", sort = FALSE)
  ## Only plotting the effects for significant SNPs. Remove all others.
  parEffData <- parEffData[interaction(parEffData[["snp"]], parEffData[["trait"]]) %in%
                              interaction(signSnp[["snp"]], signSnp[["trait"]]), ]
  if (nrow(parEffData) > 0) {
    maxVal <- max(abs(parEffData$effect), na.rm = TRUE)
  }
  ## Create title.
  if (is.null(title)) {
    title <- paste("Parental effects at QTLs for", trait)
  }
  p <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(breaks = xMarks, labels = chrs,
                                expand = c(0, 0)) +
    ggplot2::labs(title = title, x = xLab, y = yLab) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA,
                                                        color = "black",
                                                        size = 0.5,
                                                        linetype = "solid"),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  if (nrow(parEffData) > 0) {
    p <- p +
      ggplot2::geom_tile(ggplot2::aes_string(x = "cumPos", y = "trait",
                                             fill = "effect", width = "width"),
                         height = 1, data = parEffData) +
      ggplot2::scale_fill_gradientn(colors = c("blue", "cyan", "white",
                                               "yellow","red"),
                                    values = scales::rescale(c(-1,
                                                               - sqrt(.Machine$double.eps),
                                                               0,
                                                               sqrt(.Machine$double.eps),
                                                               1)),
                                    limits = c(-maxVal, maxVal),
                                    na.value = "white") +
      ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev)
  }
  p <- p + ggplot2::geom_vline(xintercept = chrBnd, color = "grey20",
                               lty = 2, size = 0.3)
  if (output) {
    plot(p)
  }
  invisible(p)
}


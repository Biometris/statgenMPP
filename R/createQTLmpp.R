#' Plot function for the class \code{QTLmpp}
#'
#' Creates a plot of an object of S3 class \code{QTLmpp}. The following types of
#' plot can be made:
#' \itemize{
#' \item{a manhattan plot, i.e. a plot of LOD-scores per SNP}
#' \item{a plot of effect sizes and directions per parent}
#' }
#' Manhattan plots, qq plots and matrix plots are made for a single trait See
#' details for a detailed description of the plots and the plot options
#' specific to the different plots.
#'
#' @section Manhattan Plot:
#' A LOD-profile of all marker positions and corresponding LOD-scores is
#' plotted. Significant markers are highlighted with red dots. By default these
#' are taken from the result of the GWAS analysis however the LOD-threshold for
#' significant parameters may be modified using the parameter \code{yThr}. The
#' treshold is plotted as a horizontal line. If there are previously known
#' marker effect, false positives and true negatives can also be marked.\cr
#' Extra parameter options:
#' \describe{
#' \item{\code{xLab}}{A character string, the x-axis label. Default =
#' \code{"Chromosomes"}}
#' \item{\code{yLab}}{A character string, the y-axis label. Default =
#' \code{-log10(p)}}
#' \item{\code{effects}}{A character vector, indicating which SNPs correspond
#' to a real (known) effect. Used for determining true/false positives and
#' false negatives. True positives are colored green, false positives orange and
#' false negatives yellow.}
#' \item{\code{colPalette}}{A color palette used for plotting. Default
#' coloring is done by chromosome, using black and grey.}
#' \item{\code{yThr}}{A numerical value for the LOD-threshold. The value from
#' the GWAS analysis is used as default.}
#' \item{\code{signLwd}}{A numerical value giving the thickness of the
#' points that are false/true positives/negatives. Default = 0.6}
#' \item{\code{lod}}{A positive numerical value. For the SNPs with a LOD-value
#' below this value, only 5\% is plotted. The chance of a SNP being plotting is
#' proportional to its LOD-score. This option can be useful when plotting a
#' large number of SNPs.}
#' \item{\code{chr}}{A vector of chromosomes to be plotted. By default, all
#' chromosomes are plotted. Using this option allows restricting the plot to a
#' subset of chromosomes.}
#' }
#'
#' @section Parental effects Plot:
#' A plot of effect sizes for each of the parents for the significant SNPs
#' found is created.\cr
#' Extra parameter options:
#' \describe{
#' \item{\code{xLab}}{A character string, the x-axis label. Default =
#' \code{"Chromosomes"}}
#' \item{\code{yLab}}{A character string, the y-axis label. Default =
#' \code{"Parents"}}
#' }
#'
#' @param x An object of class \code{QTLmpp}.
#' @param ... further arguments to be passed on to the actual plotting
#' functions.
#' @param plotType A character string indicating the type of plot to be made.
#' One of "manhattan", and "matrix".
#' @param title A character string, the title of the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE}, only a list of ggplot objects is invisibly returned.
#'
#' @export
plot.QTLmpp <- function(x,
                        ...,
                        plotType = c("manhattan", "parEffs", "QTLRegion"),
                        title = NULL,
                        output = TRUE) {
  plotType <- match.arg(plotType)
  if (plotType == "manhattan") {
    plot.GWAS(x = x, ... = ..., plotType = plotType, trial = NULL,
              trait = NULL, title = title, type = "lines", output = output)
  } else {
    dotArgs <- list(...)
    ## Get results.
    GWAResult <- x$GWAResult[[1]]
    ## Get peaks.
    signSnp <- x$signSnp[[1]]
    signSnp <- signSnp[signSnp[["snpStatus"]] == "significant SNP", ]
    ## Compute chromosome boundaries and map.
    GWAResComp <- GWAResult[GWAResult[["trait"]] == GWAResult[["trait"]][1], ]
    chrBnd <- aggregate(x = GWAResComp$pos, by = list(GWAResComp$chr),
                        FUN = max)
    ## Compute cumulative positions.
    addPos <- data.frame(chr = chrBnd[, 1],
                         add = c(0, cumsum(chrBnd[, 2]))[1:nrow(chrBnd)],
                         stringsAsFactors = FALSE)
    map <- GWAResComp[, c("snp", "chr", "pos", "LOD")]
    map <- merge(map, addPos, by = "chr")
    map[["cumPos"]] <- map[["pos"]] + map[["add"]]
    if (plotType == "parEffs") {
      parents <- x$GWASInfo$parents
      trait <- GWAResult[["trait"]][1]
      ## Create effect data - basically a long format version of GWAResult
      effectDat <- lapply(X = parents, FUN = function(parent) {
        parDat <- GWAResult
        parDat[["trait"]] <- parent
        parDat[["effect"]] <- parDat[[paste0("eff_", parent)]]
        return(parDat)
      })
      effectDat <- do.call(rbind, effectDat)
      ## Convert signSnp to long format as well.
      signSnpLong <- lapply(X = parents, FUN = function(parent) {
        parDat <- signSnp
        parDat[["trait"]] <- parent
        parDat[["effect"]] <- parDat[[paste0("eff_", parent)]]
        return(parDat)
      })
      signSnpLong <- do.call(rbind, signSnpLong)
      ## Compute chromosome boundaries.
      GWAResult <- GWAResult[!is.na(GWAResult$pos), ]
      ## Select specific chromosome(s) for plotting.
      if (!is.null(dotArgs$chr)) {
        GWAResult <- GWAResult[GWAResult$chr %in% dotArgs$chr, ]
        if (nrow(GWAResult) == 0) {
          stop("Select at least one valid chromosome for plotting.\n")
        }
      }
      do.call(effectPlot,
              args = c(list(effectDat = effectDat,
                            signSnp = signSnpLong,
                            map = map,
                            chrBoundaries = chrBnd,
                            title = title,
                            trait = trait,
                            output = output),
                       dotArgs[!(names(dotArgs) %in% c("effectDat", "signSnp",
                                                       "map", "chrBoundaries"))]))

    } else if (plotType == "QTLRegion") {
        do.call(QTLRegionPlot,
                args = c(list(signSnp = signSnp,
                              map = map,
                              chrBoundaries = chrBnd,
                              title = title,
                              output = output),
                         dotArgs[!(names(dotArgs) %in% c("signSnp",
                                                         "chrBoundaries"))]))
    }
  }
}

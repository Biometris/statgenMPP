#' Plot function for the class \code{QTLmpp}
#'
#' Creates a plot of an object of S3 class \code{QTLmpp}. The following types of
#' plot can be made:
#' \itemize{
#' \item{QTLProfile}{ A QTL profile plot, i.e. a plot of \eqn{-10log(p)} values
#' per marker.}
#' \item{parEffs} { A plot of effect sizes and directions per parent.}
#' \item{QTLRegion}{ A plot highlighting the QTLs found on the genetic map.}
#' \item{QTLProfileExt}{ A combination of the QTL profile and parental effects
#' plotted above each other.}
#' }
#' See details for a detailed description of the plots and the plot options
#' specific to the different plots.
#'
#' @section QTL Profile plot:
#' A profile of all marker positions and corresponding \eqn{-10log(p)} values is
#' plotted. QTLs found are highlighted with red dots. The threshold is plotted
#' as a horizontal line. If there are previously known marker effects, false
#' positives and true negatives can also be marked.\cr
#' Extra parameter options:
#' \describe{
#' \item{\code{xLab}}{A character string, the x-axis label. Default =
#' \code{"Chromosomes"}}
#' \item{\code{yLab}}{A character string, the y-axis label. Default =
#' \code{-log10(p)}}
#' \item{\code{effects}}{A character vector, indicating which markers correspond
#' to a real (known) effect. Used for determining true/false positives and
#' false negatives. True positives are colored green, false positives orange and
#' false negatives yellow.}
#' \item{\code{colPalette}}{A color palette used for plotting. Default
#' coloring is done by chromosome, using black and grey.}
#' \item{\code{signLwd}}{A numerical value giving the thickness of the
#' points that are false/true positives/negatives. Default = 0.6}
#' \item{\code{chr}}{A vector of chromosomes to be plotted. By default, all
#' chromosomes are plotted. Using this option allows restricting the plot to a
#' subset of chromosomes.}
#' }
#'
#' @section Parental effects Plot:
#' A plot of effect sizes for each of the parents for the QTLs found is
#' created.\cr
#' Extra parameter options:
#' \describe{
#' \item{\code{xLab}}{A character string, the x-axis label. Default =
#' \code{"Chromosomes"}}
#' \item{\code{yLab}}{A character string, the y-axis label. Default =
#' \code{"Parents"}}
#' }
#'
#' @section QTL Region Plot:
#' A plot of effect is created on which the positions of the QTLs found in
#' the QTL detection are highlighted in red.\cr
#' No extra parameter options.
#'
#' @section Extended QTL Profile Plot:
#' An extended version of the QTL Profile Plot, in which the QLT profile plot is
#' combined with the parental effect plot to make it easier to assess the
#' effects for each specific QTL found.\cr
#' No extra parameter options.
#'
#' @param x An object of class \code{QTLmpp}.
#' @param ... further arguments to be passed on to the actual plotting
#' functions.
#' @param plotType A character string indicating the type of plot to be made.
#' One of "QTLProfile", "parEffs", and "QTLRegion".
#' @param title A character string, the title of the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE}, only a ggplot object is invisibly returned.
#'
#' @examples
#' ## Read phenotypic data.
#' pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
#'                                package = "statgenMPP"))
#' ## Rename first column to genotype.
#' colnames(pheno)[1] <- "genotype"
#'
#'
#' ## Compute IBD probabilities for simulated population - AxB, AxC.
#' ABC <- calcIBDmpp(crossNames = c("AxB", "AxC"),
#'                   markerFiles = c(system.file("extdata/multipop", "AxB.txt",
#'                                               package = "statgenMPP"),
#'                                   system.file("extdata/multipop", "AxC.txt",
#'                                               package = "statgenMPP")),
#'                   pheno = pheno,
#'                   popType = "F4DH",
#'                   mapFile = system.file("extdata/multipop", "mapfile.txt",
#'                                         package = "statgenMPP"),
#'                   evalDist = 5)
#'
#' \dontrun{
#' ## Composite Interval Mapping.
#' ABC_CIM <- selQTLmpp(ABC, trait = "pheno")
#'
#' ## QTL Profile plot.
#' plot(ABC_CIM, plotType = "QTLProfile")
#'
#' ## Plot of parental effects for QTLs found.
#' plot(ABC_CIM, plotType = "parEffs")
#'
#' ## Plot of genetic map highlighting positions of QTLs found.
#' plot(ABC_CIM, plotType = "QTLRegion")
#'
#' ## Extended QTL Profile plot.
#' plot(ABC_CIM, plotType = "QTLProfileExt")
#' }
#'
#' @export
plot.QTLmpp <- function(x,
                        ...,
                        plotType = c("QTLProfile", "parEffs", "QTLRegion",
                                     "QTLProfileExt"),
                        title = NULL,
                        output = TRUE) {
  plotType <- match.arg(plotType)
  if (plotType == "QTLProfile") {
    p <- plot.GWAS(x = x, ... = ..., plotType = "manhattan", trial = NULL,
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
    GWAResComp[["chr"]] <- factor(GWAResComp[["chr"]],
                                  levels = unique(GWAResComp[["chr"]]))
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
      p <- do.call(effectPlot,
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
      p <- do.call(geneticMapPlot,
                   args = c(list(map = map,
                                 highlight = signSnp,
                                 title = title,
                                 output = output),
                            dotArgs[!(names(dotArgs) %in% c("highlight", "map"))]))
    } else if (plotType == "QTLProfileExt") {
      ## Construct title.
      if (is.null(title)) {
        title <- paste("QTL Profile (upper panel) and parental effects",
                       "(lower panel):", GWAResult[["trait"]][1])
      }
      ## Construct data for vertical lines in QTL profile.
      ## Center line in empty space between chromosomes.
      vertDat <- addPos
      minPos <-  aggregate(x = GWAResComp$pos, by = list(GWAResComp$chr),
                           FUN = min)
      vertDat[["x"]] <- vertDat[["add"]] + c(0, minPos[-1, "x"]) / 2
      ## Add title here to assure font is the same as for other plots.
      p1 <- plot(x, plotType = "QTLProfile", title = title,
                 output = FALSE) +
        ggplot2::geom_vline(ggplot2::aes_string(xintercept = "x"),
                            linetype = "dashed", data = vertDat[-1, ]) +
        ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       plot.margin = ggplot2::unit(c(0.2, 0, 0, 0), "cm"),
                       panel.background = ggplot2::element_blank(),
                       panel.border = ggplot2::element_rect(fill = NA,
                                                            color = "black",
                                                            size = 0.5,
                                                            linetype = "solid"))
      p2 <- plot(x, plotType = "parEffs", title = "", output = FALSE) +
        ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"))
      ## Get widths.
      g1 <- ggplot2::ggplotGrob(p1)
      g2 <- ggplot2::ggplotGrob(p2)
      maxWidth = grid::unit.pmax(g1$widths[2:9], g2$widths[2:9])
      ## Set widths.
      g1$widths[2:9] <- maxWidth
      g2$widths[2:9] <- maxWidth

      p <- gridExtra::grid.arrange(g1, g2)
    }
  }
  invisible(p)
}

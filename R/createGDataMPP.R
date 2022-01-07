#' S3 Class gData
#'
#' \code{createGDataMPP} creates an object of S3 class gData with genotypic and
#' phenotypic data for usage in further analysis. All input to the function is
#' optional, however at least one input should be provided. It is possible to
#' provide an existing \code{gData} object as additional input in which case
#' data is added to this object. Existing data will be overwritten with a
#' warning.
#'
#' @param gData An optional gData object to be modified. If \code{NULL}, a new
#' gData object is created.
#' @param geno A matrix or data.frame with genotypes in the rows and markers in
#' the columns. A matrix from the \code{matrix} in the base package may be
#' provided as well as as matrix from the Matrix package.\cr
#' A three dimensional array of probabilities may be provided as well with
#' genotypes in the first, markers in the second and alleles in the third
#' dimension.\cr
#' If no row names are provided, they are taken from \code{pheno} (if supplied and
#' dimension matches). If no column names are provided, the row names
#' from \code{map} are used (if supplied and dimension matches).
#' @param map A data.frame with columns \code{chr} for chromosome and
#' \code{pos} for position. Positions can be in base pair (bp) or centimorgan (cM). They
#' should not be cumulative over the chromosomes. Other columns are ignored.
#' Marker names should be in the row names. These should match the marker names
#' in \code{geno} (if supplied).
#' @param kin A kinship matrix or list of kinship matrices with genotype in
#' rows and colums. These matrices can be from the \code{matrix} class, as
#' defined in the base package, or from the \code{dsyMatrix} class, the class
#' of symmetric matrices in the Matrix package.\cr
#' The genotypes should be identical to the genotypes in \code{geno}.\cr
#' If a list of kinship matrices is provided these are supposed to be
#' chromosome specific matrices. In that case their names should match
#' the names of the chromosomes in \code{map}. If no names are
#' provided, the number of matrices should match the number of chromosomes
#' in \code{map} in which case default names are provided.
#' @param pheno A data.frame or a list of data.frames with phenotypic data,
#' with genotypes in the first column \code{genotype} and traits in the
#' following columns. The trait columns should be numerical columns only.
#' A list of data.frames can be used for replications, i.e. different
#' environments.
#' @param covar A data.frame with extra covariates per genotype. Genotypes
#' should be in the rows.
#'
#' @return An object of class \code{gData} with the following components:
#' \item{\code{map}}{a data.frame containing map data. Map is sorted by
#' chromosome and position.}
#' \item{\code{markers}}{a sparse matrix from the Matrix package containing
#' marker information in case of numerical genotypic data, a standard matrix
#' otherwise.\cr
#' If \code{geno} is a three dimensional array, \code{markers} is a three dimensional
#' array as well.}
#' \item{\code{pheno}}{a list of data.frames containing phenotypic data.}
#' \item{\code{kinship}}{a kinship matrix of class \code{dsyMatrix} from the
#'  Matrix package.}
#' \item{\code{covar}}{a data.frame with extra covariates.}
#'
#' @seealso \code{\link{summary.gData}}
#'
#' @examples
#' set.seed(1234)
#' ## Create genotypic data.
#' geno <- matrix(sample(x = c(0, 1, 2), size = 15, replace = TRUE), nrow = 3)
#' dimnames(geno) <- list(paste0("G", 1:3), paste0("M", 1:5))
#'
#' ## Construct map.
#' map <- data.frame(chr = c(1, 1, 2, 2, 2), pos = 1:5,
#'                   row.names = paste0("M", 1:5))
#'
#' ## Create phenotypic data.
#' pheno <- data.frame(paste0("G", 1:3),
#'                     matrix(rnorm(n = 12, mean = 50, sd = 5), nrow = 3),
#'                     stringsAsFactors = FALSE)
#' dimnames(pheno) = list(paste0("G", 1:3), c("genotype", paste0("T", 1:4)))
#'
#' ## Combine all data in gData object.
#' gData <- createGDataMPP(geno = geno, map = map, pheno = pheno)
#' summary(gData)
#'
#' ## Construct covariate.
#' covar <- data.frame(C1 = c("a", "a", "b"), row.names = paste0("G", 1:3))
#'
#' ## Add covariates to previously created gData object.
#' gData2 <- createGDataMPP(gData = gData, covar = covar)
#'
#' @noRd
#' @keywords internal
createGDataMPP <- function(gData = NULL,
                           geno = NULL,
                           map = NULL,
                           kin = NULL,
                           pheno = NULL,
                           covar = NULL) {
  ## Check gData
  if (!is.null(gData) && !inherits(gData, "gData")) {
    stop("Provided gData object should be of class gData.\n")
  }
  ## Check that at least one input argument, other than gData, is provided.
  if (is.null(geno) && is.null(map) && is.null(kin) && is.null(pheno) &&
      is.null(covar)) {
    stop("At least one of geno, map, kin, pheno and covar should be",
         "provided.\n")
  }
  if (is.null(geno) || length(dim(geno)) == 2) {
    gDataNw <- statgenGWAS::createGData(gData = gData, geno = geno, map = map,
                                        kin = kin, pheno = pheno, covar = covar)
  } else {
    if (!is.null(map) || !is.null(pheno)) {
      gDataNw <- statgenGWAS::createGData(gData = gData, map = map,
                                          pheno = pheno)
      map <- gDataNw$map
      pheno <- gDataNw$pheno
    } else {
      gDataNw <- structure(list(map = NULL, markers = NULL, pheno = NULL,
                                kinship = NULL, covar = NULL), class = "gData")
    }
    ## Modify geno.
    if (!is.null(geno)) {
      ## The regular 2d case is covered by the function from statgenGWAS.
      ## Here only the 3d case has to be taken care of.
      if (!(is.array(geno) && length(dim(geno)) == 3)) {
        stop("geno should be a three-dimensional array.\n")
      }
      markers <- geno
      ## Check for row names in markers. If not available take them from pheno
      ## or use default names.
      if (all(rownames(markers) == as.character(1:nrow(markers)))) {
        if (is.null(pheno)) {
          ## Default names are constructed as g001, g002, etc. with the number
          ## of zeros dependent on the number of rows.
          rownames(markers) <-
            paste0("g", formatC(1:nrow(markers),
                                width = ceiling(log10(nrow(markers))),
                                flag = "0"))
          warning("geno contains no genotype names. Default names used.\n",
                  call. = FALSE)
        } else {
          ## Phenotypic data available. Try to copy names of genotypes from
          ## genotypic data. If dimensions don't match throw an error.
          if (nrow(pheno[[1]]) == nrow(markers)) {
            rownames(markers) <- rownames(pheno[[1]])
            warning("geno contains no genotype names. Names taken from pheno.\n",
                    call. = FALSE)
          } else {
            stop("geno contains no genotype names. Dimensions between ",
                 "geno and pheno differ.\n")
          }
        }
      } else {
        ## Sort alphabetically by genotypes.
        markers <- markers[order(rownames(markers)), , , drop = FALSE]
      }
      if (is.null(colnames(markers))) {
        ## Check for column names in markers. If not available take them from map.
        ## If map not available or dimensions don't match throw an error.
        if (is.null(map)) {
          stop("geno contains no marker names. Map not available.\n")
        }
        if (nrow(map) != ncol(markers)) {
          stop("geno contains no marker names. Dimensions between geno",
               "and map differ.\n")
        }
        colnames(markers) <- rownames(map)
        warning("geno contains no marker names. Names taken from map.\n",
                call. = FALSE)
      } else if (!is.null(map)) {
        ## Both markers and map available. Make sure markernames in markers match
        ## markernames in map. Remove non-matching markers from markers.
        ## Map may still contain markers that are not in markers.
        if (any(!colnames(markers) %in%
                rownames(map)[rownames(map) %in% colnames(markers)])) {
          warning("Not all markers in geno are in map. Extra markers ",
                  "will be removed.\n", call. = FALSE)
        }
        markers <- markers[, colnames(markers) %in% rownames(map), ,
                           drop = FALSE]
        ## Check that probabilities for each marker sum to one.
        ## Always rescale values.
        ## Throw a warning if difference from one is too large.
        genoMrk <- setNames(as.data.frame(matrix(nrow = 0, ncol = 2)),
                            c("geno", "marker"))
        for (mrk in colnames(markers)) {
          mrkProbs <- rowSums(markers[, mrk, ], na.rm = TRUE)
          if (any(abs(mrkProbs - 1) > 1e-2)) {
            genoMrk <-
              rbind(genoMrk,
                    data.frame(geno = names(mrkProbs[abs(mrkProbs - 1) > 1e-2]),
                               marker = mrk, stringsAsFactors = FALSE))
          }
          markers[, mrk, ] <- markers[, mrk, ] / mrkProbs
        }
        if (nrow(genoMrk) > 0) {
          genoMrk$genoMarker <- paste(genoMrk$geno, "\t", genoMrk$marker)
          warning("Probabilities differ from 1 for the following",
                  "combinations of genotype and markers:\n",
                  paste(genoMrk$genoMarker, collapse = "\n"), call. = FALSE)
        }
      }
      if (!is.null(gData$markers)) {
        ## gData already contained a markers object. Overwrite with a warning.
        warning("existing geno will be overwritten.\n", call. = FALSE)
      }
      gDataNw$markers <- markers
    }
    if (!is.null(kin) || !is.null(covar)) {
      gDataNw <- suppressWarnings(statgenGWAS::createGData(gData = gDataNw,
                                                           kin = kin,
                                                           covar = covar))
    }
  }
  return(gDataNw)
}


#' Plot function for the class \code{gData}
#'
#' Creates a plot of an object of S3 class \code{gData}. The following types of
#' plot can be made:
#' \itemize{
#' \item{\code{genMap}}{ A plot of the genetic map.}
#' \item{\code{allGeno}}{ A plot showing for all genotypes the IBD
#' probabilities of the parent with the highest probability per marker.}
#' \item{\code{pedigree}}{ A plot showing the structure of the pedigree of
#' the population.}
#' }
#' \cr
#' See the respective sections for more details on the plots.
#'
#' @section genMap:
#' A plot is made showing the lengths of the chromosomes and the position of
#' the markers that are present in the map. It is possible to highlight one
#' or more markers using the extra parameter \code{highlight}.
#'
#' @section allGeno:
#' A plot is made showing all genotypes and markers. Each combination of
#' genotype and marker is colored according to the parent with the highest
#' probability. A darker color indicates a higher probability.
#'
#' @section pedigree:
#' A plot is made showing the structure of the pedigree for the population in
#' the \code{gData} object.
#'
#' @param x An object of class \code{gData}.
#' @param ... Further arguments to be passed on to the actual plotting
#' functions.
#' @param plotType A character string indicating the type of plot to be made.
#' One of "genMap", "allGeno" or "pedigree".
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
#' ## Plot the genetic map.
#' plot(ABC, plotType = "genMap")
#'
#' ## Plot the genetic map and highlight marker EXT_3_30.
#' plot(ABC, plotType = "genMap", highlight = "EXT_3_30")
#'
#' ## Plot the IBD probabilities across the genome for all genotypes.
#' plot(ABC, plotType = "allGeno")
#'
#' ## Plot the pedigree.
#' plot(ABC, plotType = "pedigree")
#'
#' @return A ggplot object is invisibly returned.
#'
#' @export
plot.gData <- function(x,
                       ...,
                       plotType = c("genMap", "allGeno", "pedigree"),
                       title = NULL,
                       output = TRUE) {
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  map <- x$map
  markers <- aperm(x$markers, c(2, 1, 3))
  parents <- dimnames(markers)[[3]]
  pedigree <- attr(x = x, which = "pedigree")
  genoCross <- attr(x = x, which = "genoCross")
  popType <- attr(x = x, which = "popType")
  if (plotType == "genMap") {
    highlight <- dotArgs$highlight
    if (length(highlight) > 0) {
      highlightDat <- map[highlight, ]
    } else {
      highlightDat <- NULL
    }
    p <- geneticMapPlot(map = map, highlight = highlightDat,
                        title = title, output = FALSE)
  } else if (plotType == "allGeno") {
    p <- allGenoPlot(markers = markers, map = map, parents = parents,
                     title = title)
  } else if (plotType == "pedigree") {
    if (popType == "RABBIT") {
      stop("pedigree plot not implemented for RABBIT output.\n")
    }
    p <- pedPlot(pedigree = pedigree, offSpring = colnames(markers),
                 popType = popType,
                 multiCross = length(unique(genoCross[["cross"]])) > 1,
                 genoCross = genoCross, title = title)
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}


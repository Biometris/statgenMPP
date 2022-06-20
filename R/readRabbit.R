#' Read IBD probabilities
#'
#' Read a file with IBD probabilities computed by the RABBIT software package.
#' It is possible to additionally read the pedigree file that is also used by
#' RABBIT. Reading this file allows for plotting the pedigree. Phenotypic data
#' can be added from a data.frame.
#'
#' @param infile A character string, a link to a .csv file with IBD
#' probabilities.
#' @param pedFile A character string, a link to a .csv file with pedigree
#' information as used by RABBIT as input.
#' @param pheno A data frame with at least columns "genotype" for the
#' "genotype" and one or more numerical columns containing phenotypic
#' information. A column "cross" can be used for indicating the cross the
#' genotype comes from.
#'
#' @return A \code{gDataMPP} object with map and markers corresponding to the
#' imported information in the imported .csv file.
#'
#' @examples
#' ## Read RABBIT data for barley.
#' genoFile <- system.file("extdata/barley", "barley_magicReconstruct.zip",
#'                        package = "statgenMPP")
#' barleyMPMPP <- readRABBIT(unzip(genoFile, exdir = tempdir()))
#'
#' @references Reference to RABBIT.....
#'
#' @importFrom utils hasName unzip
#' @export
readRABBIT <- function(infile,
                       pedFile = NULL,
                       pheno = NULL) {
  if (missing(infile) || !is.character(infile) || length(infile) > 1 ||
      file.access(infile, mode = 4) == -1 || tools::file_ext(infile) != "csv") {
    stop("infile should be a character string indicating a readable .csv file")
  }
  if (!is.null(pedFile) && (!is.character(pedFile) || length(pedFile) > 1 ||
                            file.access(pedFile, mode = 4) == -1 ||
                            tools::file_ext(pedFile) != "csv")) {
    stop("pedFile should be a character string indicating a readable .csv file")
  }
  if (!is.null(pheno)) {
    if (!inherits(pheno, "data.frame")) {
      stop("pheno should be a data.frame.\n")
    }
    ## If pedFile is included cross will be read from there.
    minCols <- "genotype"
    missCols <- minCols[!sapply(X = minCols, FUN = function(minCol) {
      hasName(x = pheno, name = minCol)
    })]
    if (!is.null(pedFile) && hasName(x = pheno, name = "cross")) {
      pheno <- pheno[-which(colnames(pheno) == "cross")]
    }
    if (length(missCols) > 0) {
      stop("The following columns are missing in pheno:\n",
           paste(missCols, collapse = ", "))
    }
    trtCols <- colnames(pheno)[!colnames(pheno) %in% c(minCols, "cross")]
    nonNumCols <- trtCols[!sapply(X = pheno[trtCols], FUN = is.numeric)]
    if (length(nonNumCols) > 0) {
      stop("The following columns in pheno are not numeric:\n",
           paste(nonNumCols, collapse = ", "))
    }
  }
  ## Read first header to determine whether infile is created
  ## with julia or with mathematica
  con <- file(infile, "r")
  firstHeader <- readLines(con = con, n = 1)
  close(con)
  if (firstHeader == "MagicReconstruct,designinfo") {
    rabbitRes <- readRABBITJulia(infile)
  } else {
    rabbitRes <- readRABBITMathematica(infile)
  }
  map <- rabbitRes$map
  founderProbs <- rabbitRes$founderProbs
  genoNames <- rownames(founderProbs)
  parents <- dimnames(founderProbs)[[3]]
  if (is.null(pedFile) && is.null(rabbitRes$pedDat)) {
    covar <- NULL
    if (!is.null(pheno) && hasName(x = pheno, name = "cross")) {
      ## Split pheno in pheno and covar.
      covar <- pheno["cross"]
      rownames(covar) <- pheno[["genotype"]]
      pheno <- pheno[-which(colnames(pheno) == "cross")]
    }
    genoCross <- covar
    genoCross[["geno"]] <- rownames(covar)
    popType <- "RABBIT"
  } else {
    if (is.null(rabbitRes$pedDat)) {
      ## Read pedigree.
      pedDat <- data.table::fread(pedFile, skip = "Generation",
                                  data.table = FALSE, fill = TRUE)
    } else {
      ## Pedigree is included and read from julia input.
      pedDat <- rabbitRes$pedDat
    }
    ## Get offspring.
    offDat <- pedDat[pedDat[["Generation"]] %in% genoNames, ]
    offDat[["ID"]] <- offDat[["Generation"]]
    ## Construct genoCross.
    genoCross <- offDat[c("MemberID", "Generation")]
    colnames(genoCross) <- c("cross", "geno")
    ## Construct covar.
    covar <- genoCross["cross"]
    rownames(covar) <- genoCross[["geno"]]
    ## Remove offspring.
    pedDat <- pedDat[!is.na(pedDat[["MotherID"]]), ]
    ## Generation is read as character because of presence of 2nd table.
    ## Convert to numeric so max works as expected.
    suppressWarnings(pedDat[["Generation"]] <- as.numeric(pedDat[["Generation"]]))
    ## Split in generation0 and everything else.
    gen0 <- pedDat[pedDat[["Generation"]] == 0, ]
    gen1 <- pedDat[pedDat[["Generation"]] > 0 &
                     pedDat[["MotherID"]] != pedDat[["FatherID"]], ]
    ## Compress selfing levels at the end.
    gen2 <- pedDat[pedDat[["Generation"]] > 0 &
                     pedDat[["MotherID"]] == pedDat[["FatherID"]], ]
    ## Get number of selfing levels for popType.
    nSelf <- length(unique(gen2[["Generation"]]))
    popType <- paste0("F", nSelf)
    gen2[["MotherID"]] <- gen2[gen2[["Generation"]] == min(gen2[["Generation"]]),
                               "MotherID"]
    gen2[["FatherID"]] <- gen2[gen2[["Generation"]] == min(gen2[["Generation"]]),
                               "FatherID"]
    gen2 <- gen2[gen2[["Generation"]] == max(gen2[["Generation"]]), ]
    ## Set ID.
    gen0[["ID"]] <- parents
    gen1[["ID"]] <- paste0("H", 1:(nrow(gen1)))
    ## Set type.
    gen0[["type"]] <- "INBPAR"
    gen1[["type"]] <- paste0("HYBRID", gen1[["Generation"]])
    pedDat <- rbind(gen0, gen1)
    ## Get parents.
    pedDat[["par1"]] <- pedDat[["ID"]][match(pedDat[["MotherID"]],
                                             table = pedDat[["MemberID"]])]
    pedDat[["par2"]] <- pedDat[["ID"]][match(pedDat[["FatherID"]],
                                             table = pedDat[["MemberID"]])]
    pedDat[is.na(pedDat[["par1"]]), "par1"] <- 0
    pedDat[is.na(pedDat[["par2"]]), "par2"] <- 0
    ## Set parents for offDat through gen2.
    gen2[["par1"]] <- pedDat[["par1"]][match(gen2[["MotherID"]],
                                             table = pedDat[["MemberID"]])]
    gen2[["par2"]] <- pedDat[["par2"]][match(gen2[["FatherID"]],
                                             table = pedDat[["MemberID"]])]
    offDat[["par1"]] <- gen2[["par1"]][match(offDat[["MemberID"]],
                                             table = gen2[["MemberID"]])]
    offDat[["par2"]] <- gen2[["par2"]][match(offDat[["MemberID"]],
                                             table = gen2[["MemberID"]])]
    ## Remove last generation from pedDat.
    pedDat <- pedDat[pedDat[["Generation"]] != max(pedDat[["Generation"]]), ]
    offDat[["type"]] <- popType
    pedDat <- rbind(pedDat, offDat)
    pedDat <- pedDat[c("ID", "par1", "par2", "type")]
  }
  ## Create gDataMPP object.
  res <- createGDataMPP(geno = founderProbs, map = map, pheno = pheno,
                        covar = covar)
  attr(x = res, which = "popType") <- popType
  attr(x = res, which = "genoCross") <- genoCross
  attr(x = res, which = "mapOrig") <- map
  if (!is.null(pedDat)) {
    attr(x = res, which = "pedigree") <- pedDat
  }
  return(res)
}

#' Read infile for Mathematica
#'
#' @noRd
#' @keywords internal
readRABBITMathematica <- function(infile) {
  ## Read map and marker probabilities.
  markMap <- data.table::fread(infile, skip = "haploprob", fill = TRUE,
                               data.table = FALSE)
  ## Extract map.
  map <- data.frame(chr = as.numeric(markMap[3, -1]),
                    pos = as.numeric(markMap[4, -1]),
                    row.names = as.character(markMap[2, -1]))
  markMap <- markMap[5:(which(markMap[[2]] == "ibdprob") - 1), ]
  ## Get names of genotypes and compute number of founder alleles per genotype.
  genoNames <- unique(sapply(X = strsplit(x = markMap[, 1],
                                          split = "_haplotype"), FUN = "[[", 1))
  nAlleles = nrow(markMap) / length(genoNames)
  ## Convert markers to 3D array.
  founderProbs <- array(dim = c(length(genoNames), nrow(map), nAlleles))
  for (i in 1:nrow(map)) {
    founderProbs[, i, ] <- matrix(as.numeric(markMap[, i + 1]),
                                  ncol = nAlleles, byrow = TRUE)
  }
  ## Read parent names from file.
  parents <- data.table::fread(infile, header = FALSE, nrows = nAlleles + 2,
                               skip = "haplotypes in order",
                               data.table = FALSE, select = 3)
  parents <- as.character(parents[3:nrow(parents), ])
  ## Add dimnames to markers: genotypes x markers x parents.
  dimnames(founderProbs) <- list(genoNames, rownames(map), parents)
  res <- list(map = map, founderProbs = founderProbs)
  return(res)
}


#' Read infile for Julia
#'
#' @noRd
#' @keywords internal
readRABBITJulia <- function(infile) {
  ## Get the line numbers of the RABBIT csv table blocks.
  rabbitHeaderLineIndexes <- getRabbitHeaderLines(infile)
  ## Get the pedigree.
  pedDat <- data.table::fread(infile,
                              skip = rabbitHeaderLineIndexes$designinfo$lineNr,
                              nrows = rabbitHeaderLineIndexes$designinfo$numLines,
                              data.table = FALSE,
                              header = TRUE)
  colnames(pedDat) <- c("MemberID", "MotherID", "FatherID",
                        "Gender", "Generation")
  ## Get the founders.
  founderInfo <-
    data.table::fread(infile,
                      skip = rabbitHeaderLineIndexes$founderinfo$lineNr,
                      nrows = rabbitHeaderLineIndexes$founderinfo$numLines,
                      data.table = FALSE,
                      header = TRUE)
  ## Get the offspring.
  offspringInfo <-
    data.table::fread(infile,
                      skip = rabbitHeaderLineIndexes$offspringinfo$lineNr,
                      nrows = rabbitHeaderLineIndexes$offspringinfo$numLines,
                      data.table = FALSE,
                      header = TRUE)
  offPed <- offspringInfo #<- offspringInfo[!offspringInfo[["isoutlier"]], ]
  colnames(offPed) <- c("Generation", "MemberID", "MotherID", "FatherID",
                        "Gender")
  offPed[["MotherID"]] <- offPed[["FatherID"]] <- NA
  pedDat <- rbind(pedDat, offPed)
  ## Get the markers.
  founderGeno <-
    data.table::fread(infile,
                      skip = rabbitHeaderLineIndexes$foundergeno$lineNr,
                      nrows = rabbitHeaderLineIndexes$foundergeno$numLines,
                      data.table = FALSE,
                      header = TRUE)
  ## Extract map.
  map <- data.frame(chr = as.numeric(founderGeno[, 2]),
                    pos = as.numeric(founderGeno[, 3]),
                    row.names = as.character(founderGeno[, 1]))
  ## Get the probabilities.
  haploProb <-
    data.table::fread(infile,
                      skip = rabbitHeaderLineIndexes$haploprob$lineNr,
                      nrows = rabbitHeaderLineIndexes$haploprob$numLines,
                      data.table = FALSE,
                      header = TRUE)
  ## Prepare 3D array.
  nMarkers <- nrow(founderGeno)
  nOffspring <- nrow(offspringInfo)
  nFounders <- nrow(founderInfo)
  founderProbs <- array(dim = c(nOffspring, nMarkers, nFounders),
                        dimnames = list(offspringInfo$individual,
                                        founderGeno$marker,
                                        founderInfo$individual))
  nChrom <- nrow(haploProb) / nOffspring
  ## Fill the 3D array (offspring, markers, founders).
  chrStart <- 1
  chrLength <- 0
  for (i in 1:nrow(haploProb)) {
    markerindex <- as.numeric(stringi::stri_split_regex(haploProb[i, "markerindex"], "[|]", simplify = TRUE))
    haploindex <- as.numeric(stringi::stri_split_regex(haploProb[i, "haplotype_index"], "[|]", simplify = TRUE))
    probabilities <- as.numeric(stringi::stri_split_regex(haploProb[i, "haploprob"], "[|]", simplify = TRUE))
    mat <- Matrix::sparseMatrix(i = markerindex, j = haploindex,
                                x = probabilities, check = FALSE)
    if (i %% nOffspring == 1) {
      chrStart <- chrStart + chrLength
      chrLength <- dim(mat)[1]
      chrEnd <- chrStart + chrLength - 1
    }
    founderProbs[if (i %% nOffspring == 0) nOffspring else i %% nOffspring,
                 chrStart:chrEnd, ] <- as.numeric(mat)
  }
  res <- list(map = map, founderProbs = founderProbs, pedDat = pedDat)
  return(res)
}

#' @noRd
#' @keywords internal
getRabbitHeaderLines <- function(filepath) {
  con <- file(filepath, "r")
  magicHeaders <- list()
  i <- 0
  curHeaderLineCount <- 0
  curHeader <- NULL
  while (TRUE) {
    line <- readLines(con, n = 1)
    curHeaderLineCount <- curHeaderLineCount + 1
    if (!length(line)) {
      magicHeaders[[curHeader]]$numLines = curHeaderLineCount - 2
      break
    } else if (!is.null(line) & startsWith(line, "MagicReconstruct,")) {
      if (!is.null(curHeader)) {
        magicHeaders[[curHeader]]$numLines = curHeaderLineCount - 2
      }
      header <- gsub("MagicReconstruct,", "", line)
      curHeader <- header
      curHeaderLineCount <- 0
      magicHeaders[[header]] <- list(lineNr = i + 1, numLines = 0)
    }
    i <- i + 1
  }
  close(con)
  return(magicHeaders);
}

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
#' information. if \code{pedFile} is not specified also a column "cross"
#' indicating the cross the genotype comes from is required.
#'
#' @return A \code{gData} object with map and markers corresponding to the
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
    minCols <- c("genotype", if (is.null(pedFile)) "cross")
    missCols <- minCols[!sapply(X = minCols, FUN = function(minCol) {
      hasName(x = pheno, name = minCol)
    })]
    if (!is.null(pedFile)) {
      pheno <- pheno[-which(colnames(pheno) == "cross")]
    }
    if (length(missCols) > 0) {
      stop("The following columns are missing in pheno:\n",
           paste(missCols, collapse = ", "))
    }
    trtCols <- colnames(pheno)[!colnames(pheno) %in% minCols]
    nonNumCols <- trtCols[!sapply(X = pheno[trtCols], FUN = is.numeric)]
    if (length(nonNumCols) > 0) {
      stop("The following columns in pheno are not numeric:\n",
           paste(nonNumCols, collapse = ", "))
    }
  }
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
  markArr <- array(dim = c(length(genoNames), nrow(map), nAlleles))
  for (i in 1:nrow(map)) {
    markArr[, i, ] <- matrix(as.numeric(markMap[, i + 1]),
                             ncol = nAlleles, byrow = TRUE)
  }
  ## Read parent names from file.
  parents <- data.table::fread(infile, header = FALSE, nrows = nAlleles + 2,
                               skip = "haplotypes in order",
                               data.table = FALSE, select = 3)
  parents <- as.character(parents[3:nrow(parents), ])
  ## Add dimnames to markers: genotypes x markers x parents.
  dimnames(markArr) <- list(genoNames, rownames(map), parents)
  ## Split pheno in pheno and covar.
  if (is.null(pedFile)) {
    if (!is.null(pheno)) {
      covar <- pheno["cross"]
      rownames(covar) <- pheno[["genotype"]]
      pheno <- pheno[-which(colnames(pheno) == "cross")]
    } else {
      covar <- NULL
    }
  }
  ## Read pedigree.
  if (!is.null(pedFile)) {
    pedDat <- data.table::fread(pedFile, skip = "Generation",
                                data.table = FALSE, fill = TRUE)
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
  } else {
    genoCross <- covar
    genoCross[["geno"]] <- rownames(covar)
    popType <- "RABBIT"
  }
  ## Create gData object.
  res <- createGDataMPP(geno = markArr, map = map, pheno = pheno, covar = covar)
  attr(x = res, which = "popType") <- popType
  attr(x = res, which = "genoCross") <- genoCross
  if (!is.null(pedFile)) {
    attr(x = res, which = "pedigree") <- pedDat
  }
  return(res)
}

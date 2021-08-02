#' IBDs calculation and prepare the design matrix
#'
#' @param par.names ...
#' @param cross.names ...
#' @param loc.names ...
#' @param qua.names ...
#' @param pop.types ...
#' @param mapfile ...
#' @param evaldist ...
#'
#' @importFrom utils read.table
#' @export
calcIBDmpp <- function(par.names,
                       cross.names,
                       loc.names,
                       qua.names,
                       pop.types,
                       mapfile,
                       evaldist) {
  n.cross <- length(loc.names)
  if (n.cross == 1) { # if there are only one cross, e.g., MAGIC
    IBD.cross1 <- statgenIBD::calcIBD(poptype = pop.types[[1]],
                                      locfile = loc.names[[1]],
                                      mapfile = mapfile,
                                      evaldist = evaldist)
    map <- IBD.cross1$map
    pos.names <- rownames(IBD.cross1$map)
    unipar.names <- unique(unlist(par.names))
    for (i in 1:length(pos.names)) {
      one.pos <- pos.names[i]
      df1 <- as.data.frame(statgenIBD::getProbs(IBD.cross1, one.pos))
      one.pos.IBDs <- matrix(0, ncol = length(unipar.names), nrow = nrow(df1))
      one.pos.IBDs <- as.data.frame(one.pos.IBDs)
      colnames(one.pos.IBDs) <- paste0(one.pos, "_p", unipar.names)
      rownames(one.pos.IBDs) <- df1$geno
      IBD.Parent <- df1[, 2:(length(unipar.names) + 1)]
      IBD.Het <- df1[, (length(unipar.names) + 2):ncol(df1)]
      for (p in 1:length(unipar.names)) {
        one.parent <- unipar.names[p]
        par.ind <- grep(one.parent, names(IBD.Parent))
        # Het.names<-unlist(strsplit(names(IBD.Het), "_pHET_"))
        #
        # HET.par.names<-Het.names[-which(Het.names==one.pos)]
        HET.ind <- grep(one.parent, colnames(IBD.Het))
        one.pos.IBDs[, p] <- IBD.Parent[, par.ind] * 2 +
          apply(IBD.Het[, HET.ind], 1, sum)
      }
      if (i == 1) {
        IBDmatrix <- one.pos.IBDs
      } else {
        IBDmatrix <- cbind(IBDmatrix, one.pos.IBDs)
      }
    }
    cross.info <- data.frame(ID = rownames(IBDmatrix), cross = cross.names[[1]])
    data <- cbind(cross.info, IBDmatrix)
    pheno1 <- read.table(qua.names[[1]], header = TRUE)
    data <- merge(pheno1, data, by = "ID")
  }
  if (n.cross > 1) { # if there are multiple bi-parental crosses, e.g., NAM and diallel
    for (c in 1:n.cross) {
      print(paste0("calc IBD in cross: ", cross.names[[c]]))
      IBD.cross1 <- statgenIBD::calcIBD(poptype = pop.types[[c]],
                                        locfile = loc.names[[c]],
                                        mapfile = mapfile,
                                        evaldist = evaldist)
      map <- IBD.cross1$map
      pos.names <- rownames(IBD.cross1$map)
      unipar.names <- unique(unlist(par.names))
      for (i in 1:length(pos.names)) {
        one.pos <- pos.names[i]
        df1 <- as.data.frame(statgenIBD::getProbs(IBD.cross1, c(one.pos)))
        one.pos.IBDs <- matrix(0, ncol = length(unipar.names), nrow = nrow(df1))
        one.pos.IBDs <- as.data.frame(one.pos.IBDs)
        colnames(one.pos.IBDs) <- paste0(one.pos, "_p", unipar.names)
        rownames(one.pos.IBDs) <- df1$geno
        one.pos.IBDs[which(names(one.pos.IBDs) == names(df1)[2])] <- df1[2]*2 #+df1[4]
        one.pos.IBDs[which(names(one.pos.IBDs) == names(df1)[3])] <- df1[3]*2 #+df1[4]
        if (i == 1) {
          all.pos.IBDs <- one.pos.IBDs
        } else {
          all.pos.IBDs <- cbind(all.pos.IBDs, one.pos.IBDs)
        }
      }
      cross.info <- data.frame(ID = rownames(all.pos.IBDs),
                               cross = cross.names[[c]])
      df.data <- cbind(cross.info, all.pos.IBDs)
      pheno1 <- read.table(qua.names[[c]], header = TRUE)
      df.data <- merge(pheno1, df.data, by = "ID")
      if (c == 1){
        data <- df.data
        IBDmatrix <- all.pos.IBDs
      } else {
        data <- rbind(data, df.data)
        IBDmatrix <- rbind(IBDmatrix, all.pos.IBDs)
      }
    }
  }
  MPPobj <- list(MPPinfo = list(par.names = par.names,
                                cross.names = cross.names,
                                loc.names = loc.names,
                                qua.names = qua.names,
                                pop.types = pop.types,
                                mapfile = mapfile,
                                evaldist = evaldist),
                 calcIBDres = list(par.names = unipar.names,
                                   map = map,
                                   IBDdata = data,
                                   IBDmatrix = IBDmatrix))
}

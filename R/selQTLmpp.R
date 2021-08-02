#' multi-round genome scans and select cofactors
#'
#' @export
selQTLmpp<-function(MPPobj,
                    QTLwindow=10,
                    threshold=3,
                    trait.name='pheno',
                    CIM=TRUE){

  data<-MPPobj$calcIBDres$IBDdata
  par.names<-MPPobj$calcIBDres$par.names
  map<- MPPobj$calcIBDres$map

  nloc = nrow(map)
  cofactors = NULL
  for (i in 1:100) {
    print(paste0("QTL scan for trait ", trait.name,", ", length(cofactors), " cofactors"))
    df.scan <- scanQTL(MPPobj,QTLwindow,cof=cofactors,trait.name)
    plotQTLscan(df.scan,threshold=threshold,cofactors,trait.name=trait.name)

    df.scan.sel = filter(df.scan, QTLregion==FALSE)
    # remove NAs from minlog10p values
    df.scan.sel = filter(df.scan.sel, is.na(minlog10p) == FALSE)
    ord = order(df.scan.sel$minlog10p, decreasing = TRUE)
    max_value = df.scan.sel[ord[1],]$minlog10p

    if (max_value < threshold) {
      if (!is.null(cofactors)) cofactors <- sort(cofactors)
    }

    if (max_value < threshold) break

    if (CIM==F) break

    cofactors = c(cofactors, df.scan.sel[ord[1],]$ndx)
    # if (length(cofactors)==5) break # if there are too many cofactors, we can set where it should stop

  }

  results <- list(QTLcandidates=cofactors, scanResults=df.scan)
  MPPobj[['Result']]<-results
  return(MPPobj)
}

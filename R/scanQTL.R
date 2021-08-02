#' one-round genome scan with specified cofactor(s)
#'
#' @keywords internal
scanQTL<-function(MPPobj,
                  QTLwindow=10,
                  cof=NULL,
                  trait.name='pheno'){

  map<- MPPobj$calcIBDres$map
  unipar.names<-MPPobj$calcIBDres$par.names

  nloc = nrow(map)
  sel_markers=c(1:nloc)
  Loci=map

  minlog10p = vector(length=nloc)
  QTL_region = rep(FALSE,nloc)
  effects = matrix(nrow=nloc,ncol=length(unipar.names))
  colnames(effects) = unipar.names

  for (i in 1:nloc) {
    m = sel_markers[i]
    sel_cof = select.cofactors(Loci, m, cof, QTLwindow)
    if (length(sel_cof)!=length(cof)) {
      QTL_region[i] = TRUE
    }
    obj0<-randomQTLmodel(MPPobj, trait.name,scan_pos=m,cof_pos=sel_cof,NULLmodel=TRUE)

    obj1<-randomQTLmodel(MPPobj, trait.name,scan_pos=m,cof_pos=sel_cof,NULLmodel=FALSE)

    dev = 2.0*obj1$logL - 2.0*obj0$logL
    minlog10p[i] = -log10(0.5*pchisq(dev,1,lower.tail=FALSE))
    effects[i,] = coef(obj1)[[rownames(map)[m]]]

    if (i %% 25 == 0) {
      print(i)
    }
  }

  df = data.frame(ndx=sel_markers, Loci[sel_markers,],minlog10p=minlog10p,eff=effects, QTLregion = QTL_region, trait=trait.name)

  return(df)
}

#' cofactor selection based on the QTL windowsize
#'
#' @keywords internal
calc.distance = function(markers,m1,m2)
{
  dist = ifelse(markers[m1,]$chr==markers[m2,]$chr,
                abs(markers[m1,]$pos - markers[m2,]$pos), 1.0e5)
  dist
}

#' @keywords internal
select.cofactors = function(markers, m, cofactors, QTLwindow)
{
  if (length(cofactors)==0) return(NULL)
  min.dist = 0.5*QTLwindow
  dist = sapply(cofactors,function(x) calc.distance(markers,m,x))
  ord = order(dist)
  if (dist[ord[1]] > min.dist)
  {
    return(cofactors)
  }
  cofactors[-c(ord[1])]
}

plotQTLscan = function(data, threshold, cofactors, trait.name) {
  title = paste('QTL-profile for trait ', trait.name)
  if (length(cofactors)==0) {
    title = paste0(title,', no cofactors')
  } else if (length(cofactors) == 1) {
    title = paste0(title,', one cofactor')
  } else {
    title = paste0(title,', ',length(cofactors),' cofactors')
  }

  nchr = nlevels(factor(data$chr))
  cmin = tapply(data$pos, data$chr, min)
  cmax = tapply(data$pos, data$chr, max)
  range = (cmax-cmin)+8.0
  sumrange = cumsum(range)
  start_chr = c(0,sumrange[1:(nchr-1)])
  ver_line = sumrange[1:(nchr-1)]
  xtics = (start_chr+sumrange)*0.5
  x = start_chr[data$chr] + (data$pos-cmin[data$chr])+4

  chromF = as.factor(data$chr)
  col_chr = c("black","red","green","blue","cyan1","purple","gold","brown","chartreuse","black",
              "black","black")
  col=col_chr[chromF]
  max_y = 5.0*ceiling(max(data$minlog10p)/5.0)
  plot(x=x,y=data$minlog10p,ylim=c(0,max_y),col='blue',yaxs="i",cex = 0.55,type='l',
       axes=FALSE,xpd=FALSE,pch=16,ylab='-log10(P)',xlab='Chromosomes',main=title)
  # defines the end points of the chromosomes:
  abline(v=ver_line,col='red', lwd = 2.0,lty=3)
  abline(h=threshold,col='red',lwd = 0.7)

  if (nchr==1) {
    axis(1,at=xtics[1],labels=c(nchr))
  } else {
    axis(1,at=xtics,labels=c(1:nchr))
  }

  axis(2)
  box()
}

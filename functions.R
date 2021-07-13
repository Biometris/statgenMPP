# IBDs calculation and prepare the design matrix----

calcIBDmpp<-function(par.names,cross.names,loc.names,qua.names,pop.types,mapfile,evaldist){
  
  n.cross<- length(loc.names)


  if (n.cross==1) { # if there are only one cross, e.g., MAGIC
    IBD.cross1 <- statgenIBD::calcIBD(poptype = pop.types[[1]],
                                      locfile = loc.names[[1]],
                                      mapfile = mapfile,
                                      evaldist = evaldist)
    map<-IBD.cross1$map
    pos.names<-rownames(IBD.cross1$map)
    unipar.names<- unique(unlist(par.names))
    
    for (i in 1:length(pos.names)) {
      one.pos<-pos.names[i]
      df1 <- as.data.frame(getQTL(IBD.cross1, c(one.pos)))
      
      one.pos.IBDs<-matrix(0, ncol = length(unipar.names),nrow = nrow(df1))
      one.pos.IBDs<-as.data.frame(one.pos.IBDs)
      
      colnames(one.pos.IBDs)<-paste0(one.pos,'_p',unipar.names)
      rownames(one.pos.IBDs)<-rownames(df1)
      
      IBD.Parent<-df1[,1:length(unipar.names)]
      IBD.Het<-df1[,-c(1:length(unipar.names))]
      
      for (p in 1:length(unipar.names)) {
        one.parent<-unipar.names[p]
        par.ind<-grep(one.parent,names(IBD.Parent))
        
        Het.names<-unlist(strsplit(names(IBD.Het), "_pHET_"))
        HET.par.names<-Het.names[-which(Het.names==one.pos)]
        HET.ind<-grep(par.ind,HET.par.names)
        
        one.pos.IBDs[,p]<-IBD.Parent[,par.ind]*2+apply(IBD.Het[,HET.ind], 1, sum)
      }
      
      if (i==1){
        IBDmatrix<-one.pos.IBDs
      }else{
        IBDmatrix<-cbind(IBDmatrix,one.pos.IBDs)
      }
    }
    
    cross.info<-data.frame(ID=rownames(IBDmatrix),cross=cross.names[[1]])
    data<-cbind(cross.info,IBDmatrix)
    
    pheno1 <- read.table(qua.names[[1]],header=TRUE)
    data<-merge(pheno1,data,by='ID')
  }

  if (n.cross >1) { # if there are multiple bi-parent crosses, e.g., NAM and diallel
    for (c in 1:n.cross){
      print(paste0('calc IBD in cross: ',cross.names[[c]]))
      IBD.cross1 <- statgenIBD::calcIBD(poptype = pop.types[[c]],
                                        locfile = loc.names[[c]],
                                        mapfile = mapfile,
                                        evaldist = evaldist)
      map<-IBD.cross1$map
      pos.names<-rownames(IBD.cross1$map)
      unipar.names<- unique(unlist(par.names))
      
      for (i in 1:length(pos.names)) {
        one.pos<-pos.names[i]
        df1 <- as.data.frame(getQTL(IBD.cross1, c(one.pos)))
        
        one.pos.IBDs<-matrix(0, ncol = length(unipar.names),nrow = nrow(df1))
        one.pos.IBDs<-as.data.frame(one.pos.IBDs)
        
        colnames(one.pos.IBDs)<-paste0(one.pos,'_p',unipar.names)
        rownames(one.pos.IBDs)<-rownames(df1)
        
        one.pos.IBDs[which(names(one.pos.IBDs)==names(df1)[1])]<-df1[1]*2+df1[3]
        one.pos.IBDs[which(names(one.pos.IBDs)==names(df1)[2])]<-df1[2]*2+df1[3]
        
        if (i==1){
          all.pos.IBDs<-one.pos.IBDs
        }else{
          all.pos.IBDs<-cbind(all.pos.IBDs,one.pos.IBDs)
        }
      }
      
      cross.info<-data.frame(ID=rownames(all.pos.IBDs),cross=cross.names[[c]])
      df.data<-cbind(cross.info,all.pos.IBDs)
      
      pheno1 <- read.table(qua.names[[c]],header=TRUE)
      df.data<-merge(pheno1,df.data,by='ID')
      
      if (c==1){
        data<-df.data
        IBDmatrix<-all.pos.IBDs
      } else{
        data<-rbind(data,df.data)
        IBDmatrix<-rbind(IBDmatrix,all.pos.IBDs)
      }
    }
    
    
  }
  
  MPPobj<-list(MPPinfo=list(par.names=par.names,cross.names=cross.names,loc.names=loc.names,
                            qua.names=qua.names,pop.types=pop.types,mapfile=mapfile,evaldist=evaldist),
               calcIBDres=list(par.names=unipar.names,map=map,IBDdata=data,IBDmatrix=IBDmatrix))
  
}

#one position fitted in the mixed model----

randomQTLmodel<-function(MPPobj, trait.name='pheno',scan_pos=1,cof_pos=NULL,NULLmodel=FALSE){
  
  data<-MPPobj$calcIBDres$IBDdata
  map<- MPPobj$calcIBDres$map
  unipar.names<-MPPobj$calcIBDres$par.names
  npos<-nrow(map)
  cross<-'cross'
  
  if (length(unique(data$cross))==1) { # for MAGIC-type pop with one cross
    fix = as.formula(paste0(trait.name, "~1"))
  }
  
  if (length(unique(data$cross))>1) {  # for NAM or diallel-type pops with > 1 cross
    fix = as.formula(paste0(trait.name, "~",cross))
  }
  
  
  if (NULLmodel==TRUE) { ##NULL mode wwithout cofactor
      if (length(cof_pos)==0) {
        obj <- LMMsolve(fixed=fix, residualterm=cross, data=data, eps=1.0e-8,
                        monitor=FALSE, display= FALSE)
        
      }  ##NULL mode with cofactor(s)
      if (length(cof_pos)!=0) {
        Lgrp<-list()
        
        for (i in 1:length(cof_pos)){
          cof_name<-rownames(map)[cof_pos[i]]
          cof_IBD_name<-paste0(cof_name,'_p',unipar.names)
          Lgrp[[cof_name]]<-which(colnames(data) %in% cof_IBD_name)
        }
        
        obj <- LMMsolve(fixed=fix, randomMatrices=Lgrp, 
                        residualterm=cross, data=data, eps=1.0e-8,
                        monitor=FALSE)
        
      }
    }
    
  if (NULLmodel!=TRUE){ ## FULL model
      sel_pos<-c(scan_pos,cof_pos)
      Lgrp<-list()
      for (i in 1:length(sel_pos)){
        sel_name<-rownames(map)[sel_pos[i]]
        sel_IBD_name<-paste0(sel_name,'_p',unipar.names)
        Lgrp[[sel_name]]<-which(colnames(data) %in% sel_IBD_name)
      }
      
      obj <- LMMsolve(fixed=fix, randomMatrices=Lgrp, 
                      residualterm=cross, data=data, eps=1.0e-8,
                      monitor=FALSE)
      
    }
  
  return(obj)
}

# cofactor selection based on the QTL windowsize---- 
calc.distance = function(markers,m1,m2)
{
  dist = ifelse(markers[m1,]$chr==markers[m2,]$chr,
                abs(markers[m1,]$pos - markers[m2,]$pos), 1.0e5)
  dist
}

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

#one-round genome scan with specified cofactor(s)----
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


#multi-round genome scans and select cofactors----
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


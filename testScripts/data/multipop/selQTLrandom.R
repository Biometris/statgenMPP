library(LMMsolver)  
library(statgenIBD)  
library(dplyr)
library(asreml)

#1.pop info----
par.names<-list(par.cross1=c('A','B'),
                par.cross2=c('A','C'))

cross.names<-list(cross1.name=c('AxB'),
                cross2.name=c('AxC'))

loc.names<-list(loc.cross1="AxB.loc",
             loc.cross1="AxC.loc")

qua.names<-list(loc.cross1="AxB.qua",
                loc.cross1="AxC.qua")

pop.types<-list(poptype.cross1= "F4DH",
                poptype.cross2= "F4DH")

mapfile = "mapfile.map"

evaldist= 5


##calculate IBDs and create the data frame----

calcIBDmpp<-function(par.names,cross.names,loc.names,qua.names,pop.types,mapfile,evaldist){
  
  n.cross=length(loc.names)
  
  for (c in 1:n.cross){
    IBD.cross1 <- statgenIBD::calcIBD(poptype = pop.types[[c]],
                                      locfile = loc.names[[c]],
                                      mapfile = mapfile,
                                      evaldist = evaldist)
    map<-IBD.cross1$map
    pos.names<-rownames(IBD.cross1$map)
    unipar.names<- unique(unlist(par.names))
    # diffpar.names<-setdiff(unipar.names,par.names[[c]])
    
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
    all.pos.IBDs<-cbind(cross.info,all.pos.IBDs)
    
    pheno1 <- read.table(qua.names[[c]],header=TRUE)
    df.data<-merge(pheno1,all.pos.IBDs,by='ID')
    
    if (c==1){
      data<-df.data
    } else{
      data<-rbind(data,df.data)
    }
  }
  
  MPPobj<-list(MPPinfo=list(par.names=par.names,cross.names=cross.names,loc.names=loc.names,
                            qua.names=qua.names,pop.types=pop.types,mapfile=mapfile,evaldist=evaldist),
               calcIBDres=list(par.names=unipar.names,map=map,IBDdata=data))
  
}


MPPobj<-calcIBDmpp(par.names,cross.names,loc.names,qua.names,pop.types,mapfile,evaldist)
str(MPPobj)
head(MPPobj$calcIBDres$IBDdata)
head(MPPobj$calcIBDres$map)

# 2. Model solved by LMMsolver/Asreml----

randomQTLmodel<-function(MPPobj, trait.name='pheno',scan_pos=1,cof_pos=NULL,
                         NULLmodel=FALSE,use=c('asreml','LMMsolver')){
  
  data<-MPPobj$calcIBDres$IBDdata
  map<- MPPobj$calcIBDres$map
  unipar.names<-MPPobj$calcIBDres$par.names
  
  
  npos<-nrow(map)

  cross<-'cross'
  fix = as.formula(paste0(trait.name, "~", cross))
  
  ##use LMMsolver to solve the mixed model----
  if(use=='LMMsolver') { 
    if (NULLmodel==TRUE) {
      ##NULL mode if cof_pos=NULL
      if (length(cof_pos)==0) {
        ###using LMMsolve
        obj <- LMMsolve(fixed=fix, residualterm=cross, data=data, eps=1.0e-8,
                        monitor=FALSE, display= FALSE)

      }
      ##NULL mode if cof_pos !=NULL
      if (length(cof_pos)!=0) {
        ### using LMMsolver
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
    
    if (NULLmodel!=TRUE){
      ## FULL model
      sel_pos<-c(scan_pos,cof_pos)
      ### using LMMsolver
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
  }
  
  ##use asreml to solve the mixed model----
  if(use=='asreml') { 
    if (NULLmodel==TRUE) {
      ##NULL mode if cof_pos=NULL
      if (length(cof_pos)==0) {
        ###using LMMsolve

        obj = asreml(fixed = fix,
                     rcov =~at(cross):units,
                     # residual = ~dsum(~units|cross),
                     data=data, trace=F)
        # obj$loglik
      }
    ##NULL mode if cof_pos !=NULL
      if (length(cof_pos)!=0) {
        ### using LMMsolver
        Lgrp<-list()
        for (i in 1:length(cof_pos)){
          cof_name<-rownames(map)[cof_pos[i]]
          cof_IBD_name<-paste0(cof_name,'_p',unipar.names)
          Lgrp[[cof_name]]<-which(colnames(data) %in% cof_IBD_name)
        }

        forcofLs = paste0("grp(", rownames(map)[cof_pos], ")", collapse = "+")
        random=as.formula(paste0('~ ',forcofLs))

        obj = asreml(fixed = pheno~cross,random=random,
                     rcov =~at(cross):units,
                     # residual = ~dsum(~units|cross),
                     group = Lgrp[],
                     data=data,trace=F)
      }
    }
    
    if (NULLmodel!=TRUE){
      ## FULL model
      sel_pos<-c(scan_pos,cof_pos)
      # ### using LMMsolver
      Lgrp<-list()
      for (i in 1:length(sel_pos)){
        sel_name<-rownames(map)[sel_pos[i]]
        sel_IBD_name<-paste0(sel_name,'_p',unipar.names)
        Lgrp[[sel_name]]<-which(colnames(data) %in% sel_IBD_name)
      }
      
      ## using asreml
      forcofLs = paste0("grp(", rownames(map)[sel_pos], ")", collapse = "+")
      random=as.formula(paste0('~ ',forcofLs))

      obj= asreml(fixed = pheno~cross,
                  random=random,
                  rcov =~at(cross):units,
                  # residual = ~dsum(~units|cross),
                  group = Lgrp[],
                  data=data,trace=F)
    }
  }
return(obj)
}

obj<-randomQTLmodel(MPPobj, trait.name='pheno',scan_pos=5,cof_pos=NULL,
                         NULLmodel=FALSE,use=c('LMMsolver'))
obj$logL

obj<-randomQTLmodel(MPPobj, trait.name='pheno',scan_pos=5,cof_pos=NULL,
                    NULLmodel=FALSE,use=c('asreml'))
obj$loglik


# 3. Scan the whole genome with cof_pos(s)----

select.cofactors = function(markers, m, cofactors, QTLwindow)
{
  if (length(cofactors)==0) return(NULL)
  min.dist = 0.5*QTLwindow
  dist = sapply(cofactors,function(x) calc.distance(markers,m,x))
  ord = order(dist)
  if (dist[ord[1]] > min.dist)
  {
    # use all cofactors:
    return(cofactors)
  }
  # only remove nearest cofactor:
  cofactors[-c(ord[1])]
  #print(paste("distances",dist))
  #selcof = cofactors[dist>min.dist]
  #if (length(selcof)==0) return(NULL)
  #selcof
}
calc.distance = function(markers,m1,m2)
{
  dist = ifelse(markers[m1,]$chr==markers[m2,]$chr,
                abs(markers[m1,]$pos - markers[m2,]$pos), 1.0e5)
  dist
}

scanQTL<-function(MPPobj,
                  QTLwindow=10,
                  cof=NULL,
                  use=c('LMMsolver','asreml'),
                  trait.name='pheno'){
  
  # data<-MPPobj$calcIBDres$IBDdata
  map<- MPPobj$calcIBDres$map
  unipar.names<-MPPobj$calcIBDres$par.names
  
  nloc = nrow(map)
  sel_markers=c(1:nloc)
  Loci=map
  
  minlog10p = vector(length=nloc)
  QTL_region = rep(FALSE,nloc)
  effects = matrix(nrow=nloc,ncol=length(unipar.names))
  colnames(effects) = unipar.names
  
  for (i in 1:nloc)
  {
    m = sel_markers[i]
    sel_cof = select.cofactors(Loci, m, cof, QTLwindow)
    if (length(sel_cof)!=length(cof)) {
      QTL_region[i] = TRUE
    }
    obj0<-randomQTLmodel(MPPobj, trait.name,scan_pos=m,cof_pos=sel_cof,
                         NULLmodel=TRUE,use)
    
    obj1<-randomQTLmodel(MPPobj, trait.name,scan_pos=m,cof_pos=sel_cof,
                         NULLmodel=FALSE,use)
    if (use=='LMMsolver'){
      dev = 2.0*obj1$logL - 2.0*obj0$logL
      minlog10p[i] = -log10(0.5*pchisq(dev,1,lower.tail=FALSE))
      effects[i,] = coef(obj1)[[rownames(map)[m]]]
    }
    
    if (use=='asreml'){
      dev = 2.0*obj1$loglik - 2.0*obj0$loglik
      minlog10p[i] = -log10(0.5*pchisq(dev,1,lower.tail=FALSE))
      scan_posname<-paste0('grp(',rownames(map)[m],')')
      effects[i,] = coef(obj1,list=TRUE)[[scan_posname]]
    }
    
    
    if (i %% 25 == 0) {
      print(i)
    }
  }
  
  df = data.frame(ndx=sel_markers, Loci[sel_markers,],minlog10p=minlog10p,eff=effects, QTLregion = QTL_region, trait=trait.name)
  
  return(df)
}

scanRes.LMM<-scanQTL(MPPobj,QTLwindow=10,cof=NULL,use=c('LMMsolver','asreml')[1],trait.name='pheno')
scanRes.asr<-scanQTL(MPPobj,QTLwindow=10,cof=NULL,use=c('LMMsolver','asreml')[2],trait.name='pheno')

# 4. Select QTL candidates----

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
  axis(1,at=xtics,labels=c(1:nchr))
  axis(2)
  box()
}

selQTLMPP<-function(MPPobj,
                 QTLwindow=10,
                 threshold=3,
                 use=c('LMMsolver','asreml'),
                 trait.name='pheno',CIM=TRUE){
  
  data<-MPPobj$calcIBDres$IBDdata
  par.names<-MPPobj$calcIBDres$par.names
  map<- MPPobj$calcIBDres$map
  
  nloc = nrow(map)
  cofactors = NULL
  for (i in 1:100) {
    print(paste0("QTL scan for trait ", trait.name,", ", length(cofactors), " cofactors"))
    df.scan <- scanQTL(MPPobj,QTLwindow,cof=cofactors,use,trait.name)
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
    # if (length(cofactors)==5) break
    
  }
  
  results <- list(QTLcandidates=cofactors, scanResults=df.scan)
  MPPobj[['Result']]<-results
  return(MPPobj)
}
#5. Results----

# MPPobj<-calcIBDmpp(par.names,cross.names,loc.names,qua.names,pop.types,mapfile,evaldist)
a=Sys.time()
MPPobj.LMM<-selQTLMPP(MPPobj,QTLwindow=10,threshold=3,use=c('LMMsolver','asreml')[1],trait.name='pheno',CIM=TRUE)
b=Sys.time()
b-a

a=Sys.time()
MPPobj.asr<-selQTLMPP(MPPobj,QTLwindow=10,threshold=3,use=c('LMMsolver','asreml')[2],trait.name='pheno',CIM=TRUE)
b=Sys.time()
b-a

MPPobj.LMM$Result$QTLcandidates
MPPobj.asr$Result$QTLcandidates






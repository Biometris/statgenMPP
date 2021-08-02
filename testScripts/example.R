library(LMMsolver)  
library(statgenIBD)  
library(dplyr)
library(ggplot2)

source('functions.R')
mywd<-getwd()
# 1 specify files in lists----
## 1.1. AxBxC ----

datawd<-"data/multipop"
par.names<-list(par.cross1=c('A','B'),
                par.cross2=c('A','C')) # parents names per cross
cross.names<-list(cross1.name=c('AxB'),
                  cross2.name=c('AxC')) # cross name(s)
loc.names<-list(loc.cross1="AxB.loc",
                loc.cross2="AxC.loc") # loc files
qua.names<-list(loc.cross1="AxB.qua",
                loc.cross2="AxC.qua") # qua files
pop.types<-list(poptype.cross1= "F4DH",
                poptype.cross2= "F4DH") # cross type(s)
mapfile = "mapfile.map" #map file names
evaldist= 5 #distance for IBD calc

## 1.2. popC4S3 ----

datawd<-"data/popC4S3"
par.names<-list(par.cross1=c('BLACK','BLUE','GREEN','RED'))
cross.names<-list(cross1.name=c('popC4S3'))
loc.names<-list(loc.cross1="cross.loc")
qua.names<-list(loc.cross1="cross.qua")
pop.types<-list(poptype.cross1= "C4S3")
mapfile = "mapfile.map"
evaldist=1
  
## 1.3. popC4S3 ----

datawd<-"data/popC3S4DH"
par.names<-list(par.cross1=c('BLACK','BLUE','GREEN'))
cross.names<-list(cross1.name=c('popC4S3'))
loc.names<-list(loc.cross1="cross.loc")
qua.names<-list(loc.cross1="cross.qua")
pop.types<-list(poptype.cross1= "C4S3")
mapfile = "mapfile.map"
evaldist= 1

# 2. calculate IBDs and create the data frame----
setwd(paste0(mywd,'/',datawd))
MPPobj<-calcIBDmpp(par.names,cross.names,loc.names,qua.names,pop.types,mapfile,evaldist)
head(MPPobj$calcIBDres$IBDdata)[1:6,1:10] # pheno + design matrix for LMMsolve latter
head(MPPobj$calcIBDres$IBDmatrix)[1:6,1:6] #only IBD design matrix
setwd(mywd)

#3. genome scan for multi-QTLs----
MPPobj<-selQTLmpp(MPPobj,QTLwindow=10,threshold=3,trait.name='pheno',CIM=TRUE) #
MPPobj$Result$QTLcandidates


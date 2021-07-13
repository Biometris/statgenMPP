# Martin Boer, Biometris, 22 may 2020
rm(list=ls())

library(mblib)
library(spam)
library(calcIBD)
library(dplyr)

# using asreml4:
library(asreml)

# two crosses: AxB and AxC:

# cross AxB:
df1 <- calcIBD::calcIBD(poptype = "F4DH",
               locfile     = "AxB.loc",
               mapfile     = "mapfile.map",
               evalposfile = "eval.txt")
df1 <- cbind(cross='AxB', df1)
df1 <- select(df1,-pHET)
df1$pC <- 0.0
head(df1)

pheno1 <- read.table("AxB.qua",header=TRUE)

# cross AxC:
df2 <- calcIBD::calcIBD(poptype = "F4DH",
               locfile     = "AxC.loc",
               mapfile     = "mapfile.map",
               evalposfile = "eval.txt")
df2 <- cbind(cross='AxC', df2)
df2 <- select(df2,-pHET)
df2$pB <- 0.0
df2 <- select(df2,cross,chr, pos, ind, pA, pB, pC)

pheno2 <- read.table("AxC.qua",header=TRUE)

# combine the two crosses:
df <- rbind(df1, df2)
df$cross <- as.factor(df$cross)
pheno <- rbind(pheno1,pheno2)

# Multiple QTL model, using simulated positions:
df_QTL1 <- filter(df, chr==1 & pos == 25)
df_QTL2 <- filter(df, chr==2 & pos == 75)
df_QTL3 <- filter(df, chr==3 & pos == 75)

colnames(df_QTL1) <- c('cross','chr','pos','ID','QTL1_pA','QTL1_pB','QTL1_pC')
colnames(df_QTL2) <- c('cross','chr','pos','ID','QTL2_pA','QTL2_pB','QTL2_pC')
colnames(df_QTL3) <- c('cross','chr','pos','ID','QTL3_pA','QTL3_pB','QTL3_pC')

df <- cbind(select(df_QTL1,cross,ID,QTL1_pA,QTL1_pB, QTL1_pC), 
            select(df_QTL2,QTL2_pA,QTL2_pB, QTL2_pC),
            select(df_QTL3,QTL3_pA,QTL3_pB, QTL3_pC))
df <- left_join(df, pheno,by='ID')

# analysis using asreml:

head(df)

Lgrp <- list(QTL1=c(3:5), QTL2 = c(6:8), QTL3 = c(9:11))
obj1.asr = asreml(fixed = pheno~cross,random=~grp(QTL1)+grp(QTL2)+grp(QTL3),
                  residual = ~dsum(~units|cross),
                  group = Lgrp[],
                  data=df,trace=TRUE)
summary(obj1.asr)

coefficients(obj1.asr,list=TRUE)
obj1.asr$loglik

# analysis using mblib:

y = df$pheno
X = model.matrix(~cross, df)
Z1 <- as.matrix(df[,3:5])
Z2 <- as.matrix(df[,6:8])
Z3 <- as.matrix(df[,9:11])
Z <- cbind(Z1, Z2, Z3)
lRinv <- makeRlist(df,"cross")
lGinv <- list()
lGinv[[1]] <- bdiag.spam(diag.spam(1, 3), diag.spam(0,3), diag.spam(0,3))
lGinv[[2]] <- bdiag.spam(diag.spam(0, 3), diag.spam(1,3), diag.spam(0,3))
lGinv[[3]] <- bdiag.spam(diag.spam(0, 3), diag.spam(0,3), diag.spam(1,3))

names(lGinv) = c('QTL1', 'QTL2', 'QTL3')

obj1 <- sparseMixedModels(y, X, Z=Z, lGinv=lGinv, lRinv, eps=1.0e-8, monitor=TRUE, display= FALSE)

# compare asreml with mblib, should be equal:
obj1$logL
obj1.asr$loglik
obj1$logL-obj1.asr$loglik

# effectiv dimensions
obj1$ED

# total equal to number of observations
ncol(X) + sum(obj1$ED)


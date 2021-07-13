# Martin Boer, Biometris
library(LMMsolver)  
library(statgenIBD)  
library(dplyr)

# using asreml4:
library(asreml)

# two crosses: AxB and AxC:

# cross AxB:
res1 <- statgenIBD::calcIBD(poptype = "F4DH",
                            locfile     = "AxB.loc",
                            mapfile     = "mapfile.map",
                            evaldist    = 3.0)
head(res1$map)

# Get QTLs.
QTLs <- getQTL(res1, c("M1_1","M1_2"))
head(QTLs)

Q12 <- getQTL(res1, c("EXT8", "EXT53"))
df1 <- as.data.frame(Q12)
df1 <- cbind(cross='AxB', ind=rownames(Q12), df1)
df1 <- select(df1,-EXT8_pHet)#WHL: EXT8_pAB-->EXT8_pHet; Het is direactly deleted instead x1 and add to pA and pB
df1$EXT8_pC <- 0.0
df1$EXT53_pC <- 0.0
df1 <- select(df1, cross, ind, EXT8_pA, EXT8_pB, EXT8_pC,
       EXT53_pA, EXT53_pB, EXT53_pC)
head(df1)

pheno1 <- read.table("AxB.qua",header=TRUE)

# cross AxC:
res2 <- statgenIBD::calcIBD(poptype = "F4DH",
               locfile     = "AxC.loc",
               mapfile     = "mapfile.map",
               evaldist    = 3.0)
df2<- as.data.frame(getQTL(res2, c("EXT8","EXT53")))
df2 <- cbind(cross='AxC', ind=rownames(df2), df2)
df2 <- select(df2,-EXT8_pHet)#WHL
df2$EXT8_pB <- 0.0
df2$EXT53_pB <- 0.0
df2 <- select(df2, cross, ind, EXT8_pA, EXT8_pB, EXT8_pC,
              EXT53_pA, EXT53_pB, EXT53_pC)

pheno2 <- read.table("AxC.qua",header=TRUE)

# combine the two crosses:
df <- rbind(df1, df2)
df$cross <- as.factor(df$cross)
pheno <- rbind(pheno1,pheno2)

df <- cbind(df,pheno)
df <- select(df,-geno,-error,-ID)
head(df)

# NULL MODEL, no marker:
# obj0.asr = asreml(fixed = pheno~cross,
#                   residual = ~dsum(~units|cross),
#                   data=df, trace=TRUE)

obj0.asr = asreml(fixed = pheno~cross,
                  rcov =~ idv(units),
                  data=df, trace=TRUE) #WHL


summary(obj0.asr)$varcomp

# alternative model, using 1 marker...
Lgrp <- list(QTL=c(3:5))


# obj1.asr = asreml(fixed = pheno~cross,
#                   random=~grp(QTL),
#                   residual = ~dsum(~units|cross),
#                   group = Lgrp[],
#                   data=df,trace=TRUE)
obj1.asr = asreml(fixed = pheno~cross,
                  random=~grp(QTL),
                  rcov =~at(cross):units,
                  group = Lgrp[],
                  data=df,trace=TRUE)



summary(obj1.asr)$varcomp
dev = 2.0*obj1.asr$log - 2.0*obj0.asr$log
dev
minlog10p = -log10(0.5*pchisq(dev,1,lower.tail=FALSE))
minlog10p

# NULL mode, using LMMsolve:
obj0 <- LMMsolve(fixed=pheno~cross, residualterm='cross', data=df, eps=1.0e-8,
                          monitor=TRUE, display= FALSE)
# check with asreml:
obj0$logL
obj0.asr$loglik

# effective degrees of freedom 
df %>% group_by(cross) %>% tally()
obj0$ED

head(df)

# include QTL, using LMMsolve
# Lgrp as before for asreml:
Lgrp <- list(QTL=c(3:5))
obj1 <- LMMsolve(fixed=pheno~cross, randomMatrices=Lgrp, 
                 residualterm='cross', data=df, eps=1.0e-8,
                 monitor=TRUE, display= FALSE)

# check with asreml:
obj1$logL
obj1.asr$loglik

# effective degrees of freedom: maximum effective dimension for QTL is 2, 
# because ofsum to one constraint:
obj1$ED

# coefficients for QTL-effect:
coefficients(obj1.asr,list=TRUE)
coef(obj1)

# sum of effects equal to zero:
sum(coef(obj1)$QTL)



#WHL:With cofactor
Lgrp <- list(QTL1=c(3:5), COFACTOR1=c(6:8))
obj2 <- LMMsolve(fixed=pheno~cross, randomMatrices=Lgrp, 
                 residualterm='cross', data=df, eps=1.0e-8,
                 monitor=TRUE)
coef(obj2)$QTL1
coef(obj2)$COFACTOR1

# in terms of asreml, very similar:
# obj2.asr = asreml(fixed = pheno~cross,random=~grp(QTL1)+grp(COFACTOR1),
#                   residual = ~dsum(~units|cross),
#                   group = Lgrp[],
#                   data=df,trace=TRUE)


obj2.asr = asreml(fixed = pheno~cross,random=~grp(QTL1)+grp(COFACTOR1),
                  rcov =  ~at(cross):units,
                  group = Lgrp[],
                  data=df,trace=TRUE)

obj2$logL
obj2.asr$loglik
coef(obj2.asr)$random

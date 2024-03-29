% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selQTLMPP.R
\name{selQTLMPP}
\alias{selQTLMPP}
\title{Multi round genome scans for QTL detection}
\usage{
selQTLMPP(
  MPPobj,
  trait = NULL,
  QTLwindow = 10,
  threshold = 3,
  maxCofactors = NULL,
  K = NULL,
  computeKin = FALSE,
  parallel = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{MPPobj}{An object of class gDataMPP, typically the output of either
\code{\link{calcIBDMPP}} or \code{\link{readRABBITMPP}}.}

\item{trait}{A character string indicating the trait QTL mapping is done for.}

\item{QTLwindow}{A numerical value indicating the window around a QTL that
is considered as part of that QTL.}

\item{threshold}{A numerical value indicating the threshold for the
\eqn{-log10p} value of a marker to be considered a QTL.}

\item{maxCofactors}{A numerical value, the maximum number of cofactors to
include in the model. If \code{NULL} cofactors are added until no new
cofactors are found.}

\item{K}{A list of chromosome specific kinship matrices. If
\code{NULL} and \code{computeKin = FALSE} no kinship matrix is included in
the models.}

\item{computeKin}{Should chromosome specific kinship matrices be computed?}

\item{parallel}{Should the computation of variance components be done in
parallel? This requires a parallel back-end to be registered. See examples.}

\item{verbose}{Should progress and intermediate plots be output?}
}
\value{
An object of class \code{QTLMPP}
}
\description{
Multi round genome scans for QTL detection.\cr\cr
Several rounds of QTL detection are performed. First a model is fitted
without cofactors. If for at least one marker the \eqn{-log10(p)} value is
above the threshold the marker with the lowest p-Value is added as cofactor
in the next round of QTL detection. This process continues until there are
no new markers with a \eqn{-log10(p)} value above the threshold or until
the maximum number of cofactors is reached.
}
\details{
By default only family specific effects and residual variances and no
kinship relations are included in the model. It is possible to include
kinship relations by either specifying \code{computeKin = TRUE}. When doing
so the kinship matrix is computed by averaging \eqn{Z Z^t} over all markers,
where \eqn{Z} is the genotype x parents matrix for the marker. It is also
possible to specify a list of precomputed chromosome
specific kinship matrices in \code{K}. Note that adding a kinship matrix
to the model increases the computation time a lot, especially for large
populations.
}
\examples{
## Read phenotypic data.
pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                               package = "statgenMPP"))
## Rename first column to genotype.
colnames(pheno)[1] <- "genotype"


## Compute IBD probabilities for simulated population - AxB, AxC.
ABC <- calcIBDMPP(crossNames = c("AxB", "AxC"),
                  markerFiles = c(system.file("extdata/multipop", "AxB.txt",
                                              package = "statgenMPP"),
                                  system.file("extdata/multipop", "AxC.txt",
                                              package = "statgenMPP")),
                  pheno = pheno,
                  popType = "F4DH",
                  mapFile = system.file("extdata/multipop", "mapfile.txt",
                                        package = "statgenMPP"),
                  evalDist = 5)

## Single-QTL Mapping.
ABC_SQM <- selQTLMPP(ABC, trait = "yield", maxCofactors = 0)

## Multi-QTL Mapping.
\dontrun{
## Register parallel back-end with 2 cores.
doParallel::registerDoParallel(cores = 2)

## Run multi-QTL mapping.
ABC_MQM <- selQTLMPP(ABC, trait = "yield", parallel = TRUE)

## Run multi-QTL mapping - include kinship matrix.
ABC_MQM_kin <- selQTLMPP(ABC, trait = "yield", parallel = TRUE,
                        computeKin = TRUE)
}

}
\seealso{
\code{\link{kinshipIBD}}
}

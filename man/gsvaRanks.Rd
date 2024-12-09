% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsva.R
\name{gsvaRanks}
\alias{gsvaRanks}
\alias{gsvaRanks,gsvaParam-method}
\alias{gsvaScores}
\alias{gsvaScores,gsvaRanksParam-method}
\title{GSVA ranks and scores}
\usage{
\S4method{gsvaRanks}{gsvaParam}(param, verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose))

\S4method{gsvaScores}{gsvaRanksParam}(param, verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose))
}
\arguments{
\item{param}{A parameter object of the \code{\linkS4class{gsvaRanksParam}} class.}

\item{verbose}{Gives information about each calculation step. Default: \code{TRUE}.}

\item{BPPARAM}{An object of class \code{\link{BiocParallelParam}} specifying parameters
related to the parallel execution of some of the tasks and calculations
within this function.}
}
\value{
In the case of the \code{gsvaRanks()} method, an object of class
\code{\linkS4class{gsvaRanksParam}}.

In the case of the \code{gsvaScores()} method, a gene-set by sample matrix
of GSVA enrichment scores stored in a container object of the same type as
the input ranks data container. If
the input was a base matrix or a \code{\linkS4class{dgCMatrix}} object, then the output will
be a base matrix object with the gene sets employed in the calculations
stored in an attribute called \code{geneSets}. If the input was an
\code{\link{ExpressionSet}} object, then the output will be also an \code{\link{ExpressionSet}}
object with the gene sets employed in the calculations stored in an
attributed called \code{geneSets}. If the input was an object of one of the
classes described in \code{\link{GsvaExprData}}, such as a \code{\link{SingleCellExperiment}},
then the output will be of the same class, where enrichment scores will be
stored in an assay called \code{es} and the gene sets employed in the
calculations will be stored in the \code{rowData} slot of the object under the
column name \code{gs}.
}
\description{
Calculate GSVA scores in two steps: (1) calculate GSVA
ranks; and (2) calculate GSVA scores using the previously calculated
ranks.
}
\examples{
library(GSVA)

p <- 10 ## number of genes
n <- 30 ## number of samples
nGrp1 <- 15 ## number of samples in group 1
nGrp2 <- n - nGrp1 ## number of samples in group 2

## consider three disjoint gene sets
geneSets <- list(gset1=paste0("g", 1:3),
                 gset2=paste0("g", 4:6),
                 gset3=paste0("g", 7:10))

## sample data from a normal distribution with mean 0 and st.dev. 1
y <- matrix(rnorm(n*p), nrow=p, ncol=n,
            dimnames=list(paste("g", 1:p, sep="") , paste("s", 1:n, sep="")))

## genes in set1 are expressed at higher levels in the last 'nGrp1+1' to 'n' samples
y[geneSets$set1, (nGrp1+1):n] <- y[geneSets$set1, (nGrp1+1):n] + 2

## build GSVA parameter object
gsvapar <- gsvaParam(y, geneSets)

## calculate GSVA ranks
gsvarankspar <- gsvaRanks(gsvapar)
gsvarankspar
## calculate GSVA scores
gsva_es <- gsvaScores(gsvarankspar)
gsva_es

## calculate now GSVA scores in a single step
gsva_es1 <- gsva(gsvapar)

## both approaches give the same result with the same input gene sets
all.equal(gsva_es1, gsva_es)

## however, results will be (obviously) different with different gene sets
geneSets2 <- list(gset1=paste0("g", 3:6),
                  gset2=paste0("g", c(1, 2, 7, 8)))

## note that there is no need to calculate the GSVA ranks again
geneSets(gsvarankspar) <- geneSets2
gsvaScores(gsvarankspar)

}
\references{
Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
variation analysis for microarray and RNA-Seq data.
\emph{BMC Bioinformatics}, 14:7, 2013.
\href{https://doi.org/10.1186/1471-2105-14-7}{DOI}
}
\seealso{
\code{\linkS4class{gsvaParam}}, \code{\linkS4class{gsvaRanksParam}}, \code{\link{gsva}}
}

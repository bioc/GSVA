% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R, R/gsvaParam.R
\docType{class}
\name{gsvaParam-class}
\alias{gsvaParam-class}
\alias{gsvaRanksParam-class}
\alias{gsvaParam}
\alias{geneSets<-,gsvaRanksParam,GsvaGeneSets-method}
\alias{geneSets<-}
\title{\code{gsvaParam} class}
\usage{
gsvaParam(
  exprData,
  geneSets,
  assay = NA_character_,
  annotation = NULL,
  minSize = 1,
  maxSize = Inf,
  kcdf = c("auto", "Gaussian", "Poisson", "none"),
  kcdfNoneMinSampleSize = 200,
  tau = 1,
  maxDiff = TRUE,
  absRanking = FALSE,
  sparse = TRUE,
  checkNA = c("auto", "yes", "no"),
  use = c("everything", "all.obs", "na.rm")
)

\S4method{geneSets}{gsvaRanksParam,GsvaGeneSets}(object) <- value
}
\arguments{
\item{exprData}{The expression data set.  Must be one of the classes
supported by \code{\linkS4class{GsvaExprData}}.  For a list of these classes, see its
help page using \code{help(GsvaExprData)}.}

\item{geneSets}{The gene sets.  Must be one of the classes supported by
\code{\linkS4class{GsvaGeneSets}}.  For a list of these classes, see its help page using
\code{help(GsvaGeneSets)}.}

\item{assay}{Character vector of length 1.  The name of the assay to use in
case \code{exprData} is a multi-assay container, otherwise ignored.  By default,
the first assay is used.}

\item{annotation}{An object of class \code{\linkS4class{GeneIdentifierType}} from
package \code{GSEABase} describing the gene identifiers used as the row names of
the expression data set.  See \code{\link{GeneIdentifierType}} for help on available
gene identifier types and how to construct them.  This
information can be used to map gene identifiers occurring in the gene sets.

If the default value \code{NULL} is provided, an attempt will be made to extract
the gene identifier type from the expression data set provided as \code{exprData}
(by calling \code{\link{gsvaAnnotation}} on it).  If still not successful, the
\code{NullIdentifier()} will be used as the gene identifier type, gene identifier
mapping will be disabled and gene identifiers used in expression data set and
gene sets can only be matched directly.}

\item{minSize}{Numeric vector of length 1.  Minimum size of the resulting gene
sets after gene identifier mapping. By default, the minimum size is 1.}

\item{maxSize}{Numeric vector of length 1.  Maximum size of the resulting gene
sets after gene identifier mapping. By default, the maximum size is \code{Inf}.}

\item{kcdf}{Character vector of length 1 denoting the kernel to use during
the non-parametric estimation of the empirical cumulative distribution
function (ECDF) of expression levels across samples. The value \code{kcdf="auto"}
will allow GSVA to automatically choose one of the possible values. The
value \code{kcdf="Gaussian"} is suitable when input expression values are
continuous, such as microarray fluorescent units in logarithmic scale,
RNA-seq log-CPMs, log-RPKMs, or log-TPMs. When input expression values are
integer counts, such as those derived from RNA-seq experiments, then this
argument should be set to \code{kcdf="Poisson"}. When we do not want to use a
kernel approach for the estimation of the ECDF, then we should set
\code{kcdf="none"}.}

\item{kcdfNoneMinSampleSize}{Integer vector of length 1. When \code{kcdf="auto"},
this parameter decides at what minimum sample size \code{kcdf="none"}, i.e., the
estimation of the empirical cumulative distribution function (ECDF) of
expression levels across samples is performed directly without using a
kernel. By default, this value is set to 200; see the \code{kcdf} slot.}

\item{tau}{Numeric vector of length 1.  The exponent defining the weight of
the tail in the random walk performed by the \code{GSVA} (Hänzelmann et al.,
2013) method.  The default value is 1 as described in the paper.}

\item{maxDiff}{Logical vector of length 1 which offers two approaches to
calculate the enrichment statistic (ES) from the KS random walk statistic.
\itemize{
\item \code{FALSE}: ES is calculated as the maximum distance of the random walk
from 0. This approach produces a distribution of enrichment scores that is
bimodal, but it can give large enrichment scores to gene sets whose genes
are not concordantly activated in one direction only.
\item \code{TRUE} (the default): ES is calculated as the magnitude difference between
the largest positive and negative random walk deviations. This default value
gives larger enrichment scores to gene sets whose genes are concordantly
activated in one direction only.
}}

\item{absRanking}{Logical vector of length 1 used only when \code{maxDiff=TRUE}.
When \code{absRanking=FALSE} (default) a modified Kuiper statistic is used to
calculate enrichment scores, taking the magnitude difference between the
largest positive and negative random walk deviations. When
\code{absRanking=TRUE} the original Kuiper statistic that sums the largest
positive and negative random walk deviations is used.}

\item{sparse}{Logical vector of length 1 used only when the input expression
data in \code{exprData} is stored in a sparse matrix (e.g., a \code{dgCMatrix} or a
\code{SingleCellExperiment} object storing the expression data in a \code{dgCMatrix}).
In such a case, when \code{sparse=TRUE} (default), a sparse version of the GSVA
algorithm will be applied. Otherwise, when \code{sparse=FALSE}, the classical
version of the GSVA algorithm will be used.}

\item{checkNA}{Character vector of length 1 specifying whether the input
expression data should be checked for the presence of missing (\code{NA}) values.
This must be
one of the strings \code{"auto"} (default), \code{"yes"}, or \code{"no"}. The default value
\code{"auto"} means that the software will perform that check only when the input
expression data is provided as a base \code{\link{matrix}}, an \code{\link{ExpressionSet}} or a
\code{\link{SummarizedExperiment}} object, while every other type of input expression
data container (e.g., \code{\link{SingleCellExperiment}}, etc.) will not be checked.
If \code{checkNA="yes"}, then the input expression data will be checked for
missing values irrespective of the object class of the data container, and
if \code{checkNA="no"}, then that check will not be performed.}

\item{use}{Character vector of length 1 specifying a policy for dealing with
missing values (\code{NA}s) in the input expression data argument \code{exprData}. It
only applies when either \code{checkNA="yes"}, or \code{checkNA="auto"} (see the
\code{checkNA} parameter. The argument value must be one of the strings
\code{"everything"} (default), \code{"all.obs"}, or \code{"na.rm"}. The policy of the
default value \code{"everything"} consists of propagating \code{NA}s so that the
resulting enrichment score will be \code{NA}, whenever one or more of its
contributing values is \code{NA}, giving a warning when that happens. When
\code{use="all.obs"}, the presence of \code{NA}s in the input expression data will
produce an error. Finally, when \code{use="na.rm"}, \code{NA} values in the input
expression data will be removed from calculations, giving a warning when that
happens, and giving an error if no values are left after removing the \code{NA}
values.}

\item{object}{For the replacement method, an object of class
\code{\linkS4class{gsvaRanksParam}}.}

\item{value}{For the replacement method, an object of the classes supported by
\code{\linkS4class{GsvaGeneSets}}.}
}
\value{
A new \code{\linkS4class{gsvaParam}} object.
}
\description{
S4 class for GSVA method parameter objects.

Objects of class \code{gsvaParam} contain the parameters for running
the \code{GSVA} method.
}
\details{
In addition to the common parameter slots inherited from \verb{[GsvaMethodParam]},
this class has slots for the six method-specific parameters of the GSVA
method described below.

In addition to a number of parameters shared with all methods
implemented by package GSVA, \code{GSVA} takes six method-specific parameters.
All of these parameters are described in detail below.
}
\section{Slots}{

\describe{
\item{\code{kcdf}}{Character vector of length 1 denoting the kernel to use during
the non-parametric estimation of the empirical cumulative distribution
function (ECDF) of expression levels across samples. The value \code{kcdf="auto"}
will allow GSVA to automatically choose one of the possible values. The
value \code{kcdf="Gaussian"} is suitable when input expression values are
continuous, such as microarray fluorescent units in logarithmic scale,
RNA-seq log-CPMs, log-RPKMs, or log-TPMs. When input expression values are
integer counts, such as those derived from RNA-seq experiments, then this
argument should be set to \code{kcdf="Poisson"}. When we do not want to use a
kernel approach for the estimation of the ECDF, then we should set
\code{kcdf="none"}.}

\item{\code{kcdfNoneMinSampleSize}}{Integer vector of length 1. When \code{kcdf="auto"},
this parameter decides at what minimum sample size \code{kcdf="none"}, i.e., the
estimation of the empirical cumulative distribution function (ECDF) of
expression levels across samples is performed directly without using a
kernel; see the \code{kcdf} slot.}

\item{\code{tau}}{Numeric vector of length 1.  The exponent defining the weight of
the tail in the random walk performed by the GSVA (Hänzelmann et al., 2013)
method.}

\item{\code{maxDiff}}{Logical vector of length 1 which offers two approaches to
calculate the enrichment statistic (ES) from the KS random walk statistic.
\itemize{
\item \code{FALSE}: ES is calculated as the maximum distance of the random walk from 0.
\item \code{TRUE}: ES is calculated as the magnitude difference between
the largest positive and negative random walk deviations.
}}

\item{\code{absRanking}}{Logical vector of length 1 used only when \code{maxDiff=TRUE}.
When \code{absRanking=FALSE} a modified Kuiper statistic is used to calculate
enrichment scores, taking the magnitude difference between the largest
positive and negative random walk deviations. When \code{absRanking=TRUE} the
original Kuiper statistic that sums the largest positive and negative
random walk deviations, is used. In this latter case, gene sets with genes
enriched on either extreme (high or low) will be regarded as ’highly’
activated.}

\item{\code{sparse}}{Logical vector of length 1 used only when the input expression
data in \code{exprData} is stored in a sparse matrix (e.g., a \code{dgCMatrix} or a
container object, such as a \code{SingleCellExperiment}, storing the expression
data in a \code{dgCMatrix}).
In such a case, when \code{sparse=TRUE}, a sparse version of the GSVA algorithm
will be applied. Otherwise, when \code{sparse=FALSE}, the classical version of
the GSVA algorithm will be used.}

\item{\code{checkNA}}{Character vector of length 1. One of the strings \code{"auto"}
(default), \code{"yes"}, or \code{"no"}, which refer to whether the input expression
data should be checked for the presence of missing (\code{NA}) values.}

\item{\code{didCheckNA}}{Logical vector of length 1, indicating whether the input
expression data was checked for the presence of missing (\code{NA}) values.}

\item{\code{anyNA}}{Logical vector of length 1, indicating whether the input
expression data contains missing (\code{NA}) values.}

\item{\code{use}}{Character vector of length 1. One of the strings \code{"everything"}
(default), \code{"all.obs"}, or \code{"na.rm"}, which refer to three different policies
to apply in the presence of missing values in the input expression data; see
\code{\link{ssgseaParam}}.}
}}

\examples{
library(GSVA)
library(GSVAdata)

data(leukemia)
data(c2BroadSets)

## for simplicity, use only a subset of the sample data
ses <- leukemia_eset[1:1000, ]
gsc <- c2BroadSets[1:100]
gp1 <- gsvaParam(ses, gsc)
gp1


}
\references{
Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
variation analysis for microarray and RNA-Seq data.
\emph{BMC Bioinformatics}, 14:7, 2013.
\href{https://doi.org/10.1186/1471-2105-14-7}{DOI}
}
\seealso{
\code{\linkS4class{GsvaExprData}},
\code{\linkS4class{GsvaGeneSets}},
\code{\linkS4class{GsvaMethodParam}},
\code{\linkS4class{plageParam}},
\code{\linkS4class{zscoreParam}},
\code{\linkS4class{ssgseaParam}}
}

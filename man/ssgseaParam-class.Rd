% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R, R/ssgseaParam.R
\docType{class}
\name{ssgseaParam-class}
\alias{ssgseaParam-class}
\alias{ssgseaParam}
\alias{anyNA,ssgseaParam-method}
\title{\code{ssgseaParam} class}
\usage{
ssgseaParam(
  exprData,
  geneSets,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE,
  checkNA = c("auto", "yes", "no"),
  use = c("everything", "all.obs", "na.rm")
)

\S4method{anyNA}{ssgseaParam}(x, recursive = FALSE)
}
\arguments{
\item{exprData}{The expression data.  Must be one of the classes
supported by \code{\linkS4class{GsvaExprData}}. Type \code{help(GsvaExprData)} to consult
the available classes.}

\item{geneSets}{The gene sets.  Must be one of the classes supported by
\code{\linkS4class{GsvaGeneSets}}.}

\item{assay}{The name of the assay to use in case \code{exprData} is a multi-assay
container, otherwise ignored.  By default, the first assay is used.}

\item{annotation}{The name of a Bioconductor annotation package for the gene
identifiers occurring in the row names of the expression data matrix.  This
can be used to map gene identifiers occurring in the gene sets if those are
provided in a \code{\link{GeneSetCollection}}.  By default gene identifiers used in
expression data matrix and gene sets are matched directly.}

\item{minSize}{Minimum size of the resulting gene sets after gene identifier
mapping. By default, the minimum size is 1.}

\item{maxSize}{Maximum size of the resulting gene sets after gene identifier
mapping. By default, the maximum size is \code{Inf}.}

\item{alpha}{Numeric vector of length 1.  The exponent defining the
weight of the tail in the random walk performed by the \code{ssGSEA} (Barbie et
al., 2009) method.  The default value is 0.25 as described in the paper.}

\item{normalize}{Logical vector of length 1; if \code{TRUE} runs the \code{ssGSEA}
method from Barbie et al. (2009) normalizing the scores by the absolute
difference between the minimum and the maximum, as described in their paper.
Otherwise this last normalization step is skipped.}

\item{checkNA}{Character string specifying whether the input expression data
should be checked for the presence of missing (\code{NA}) values. This must be
one of the strings \code{"auto"} (default), \code{"yes"}, or \code{"no"}. The default value
\code{"auto"} means that the software will perform that check only when the input
expression data is provided as a base \code{\link{matrix}}, an \code{\link{ExpressionSet}} or a
\code{\link{SummarizedExperiment}} object, while every other type of input expression
data container (e.g., \code{\link{SingleCellExperiment}}, etc.) will not be checked.
If \code{checkNA="yes"}, then the input expression data will be checked for
missing values irrespective of the object class of the data container, and
if \code{checkNA="no"}, then that check will not be performed.}

\item{use}{Character string specifying a policy for dealing with missing
values (\code{NA}s) in the input expression data argument \code{exprData}. It only
applies when either \code{checkNA="yes"}, or \code{checkNA="auto"} (see the \code{checkNA}
parameter. The argument value must be one of the strings \code{"everything"}
(default), \code{"all.obs"}, or \code{"na.rm"}. The policy of the default value
\code{"everything"} consists of propagating \code{NA}s so that the resulting enrichment
scores will be \code{NA}, whenever one or more of its contributing values is \code{NA},
giving a warning when that happens. When \code{use="all.obs"}, the presence of
\code{NA}s in the input expression data will produce an error. Finally, when
\code{use="na.rm"}, \code{NA} values in the input expression data will be removed from
calculations, giving a warning when that happens, and giving an error if no
values are left after removing the \code{NA} values.}

\item{x}{An object of class \code{\link{ssgseaParam}}.}

\item{recursive}{Not used with \code{x} being an object of class \code{\link{ssgseaParam}}.}
}
\value{
A new \code{\linkS4class{ssgseaParam}} object.
}
\description{
Method-specific parameters for the ssGSEA method.

Objects of class \code{ssgseaParam} contain the parameters for running
the \code{ssGSEA} method.
}
\details{
In addition to the two common parameter slots inherited from
\verb{[GsvaMethodParam]}, this class has slots for the two method-specific
parameters of the \code{ssGSEA} method described below.

In addition to an expression data set and a collection of
gene sets, \code{ssGSEA} takes two method-specific parameters as described below.
}
\section{Slots}{

\describe{
\item{\code{alpha}}{Numeric vector of length 1.  The exponent defining the
weight of the tail in the random walk performed by the ssGSEA (Barbie et
al., 2009) method.}

\item{\code{normalize}}{Logical vector of length 1.  If \code{TRUE} runs the ssGSEA
method from Barbie et al. (2009) normalizing the scores by the absolute
difference between the minimum and the maximum, as described in their paper.
Otherwise this last normalization step is skipped.}

\item{\code{checkNA}}{Character string. One of the strings \code{"auto"} (default),
\code{"yes"}, or \code{"no"}, which refer to whether the input expression data should
be checked for the presence of missing (\code{NA}) values.}

\item{\code{didCheckNA}}{Logical vector of length 1, indicating whether the input
expression data was checked for the presence of missing (\code{NA}) values.}

\item{\code{anyNA}}{Logical vector of length 1, indicating whether the input
expression data contains missing (\code{NA}) values.}

\item{\code{use}}{Character string. One of the strings \code{"everything"} (default),
\code{"all.obs"}, or \code{"na.rm"}, which refer to three different policies to apply
in the presence of missing values in the input expression data; see
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
sp1 <- ssgseaParam(ses, gsc)
sp1

}
\references{
Barbie, D.A. et al. Systematic RNA interference reveals that
oncogenic KRAS-driven cancers require TBK1.
\emph{Nature}, 462(5):108-112, 2009.
\href{https://doi.org/10.1038/nature08460}{DOI}
}
\seealso{
\code{\linkS4class{GsvaExprData}},
\code{\linkS4class{GsvaGeneSets}},
\code{\linkS4class{GsvaMethodParam}},
\code{\linkS4class{plageParam}},
\code{\linkS4class{zscoreParam}},
\code{\linkS4class{gsvaParam}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsvaNewAPI.R
\name{geneSets}
\alias{geneSets}
\alias{geneSetSizes}
\alias{geneSets,GsvaMethodParam-method}
\alias{geneSets,SummarizedExperiment-method}
\alias{geneSets,SingleCellExperiment-method}
\alias{geneSets,SpatialExperiment-method}
\alias{geneSets,GsvaExprData-method}
\alias{geneSetSizes,GsvaMethodParam-method}
\alias{geneSetSizes,GsvaExprData-method}
\title{Retrieve or Determine Gene Sets}
\usage{
\S4method{geneSets}{GsvaMethodParam}(obj)

\S4method{geneSets}{SummarizedExperiment}(obj)

\S4method{geneSets}{SingleCellExperiment}(obj)

\S4method{geneSets}{SpatialExperiment}(obj)

\S4method{geneSets}{GsvaExprData}(obj)

\S4method{geneSetSizes}{GsvaMethodParam}(obj)

\S4method{geneSetSizes}{GsvaExprData}(obj)
}
\arguments{
\item{obj}{An object of one of the following classes:
\itemize{
\item An expression data object of one of the classes described in
\code{\linkS4class{GsvaExprData}} that is the return value of a call to \code{gsva()}.
\item A parameter object of one of the classes described in
\code{\linkS4class{GsvaMethodParam}} that could be used in a call to \code{gsva()}.
}}
}
\value{
The \code{geneSets()} methods return a named list of character vectors
where each character vector contains the gene IDs of a gene set.
The \code{geneSetSizes()} methods return a named integer vector of gene set sizes.
}
\description{
Retrieves or determines the gene sets that have been used
or would be used in a \code{gsva()} gene set analysis.  These are not necessarily
the same as the input gene sets.  See Details.
}
\details{
The gene sets used in a \code{gsva()} gene set analysis, or just their
sizes, may be a valuable input to subsequent analyses.  However, they are not
necessarily the same as the original input gene sets, or their sizes: based
on user choices, the gene annotation used, or presence/absence of genes in
gene sets and expression data set, \code{gsva()} may have to modify them during
the preparation of an analysis run.
In order to make use of these gene sets or their sizes, you can either
\itemize{
\item retrieve them from the object returned by \code{gsva()} by passing this object
to \code{geneSets()} or \code{geneSetSizes()}, or
\item predict them by calling \code{geneSets()} or \code{geneSetSizes()} on the parameter
object that would also be passed to \code{gsva()}.  This is much slower and should
only be done if you do not intend to run an actual gene set analysis.
}

\code{geneSetSizes()} is a convenience wrapper running \code{lengths()} on the list of
gene sets returned by \code{geneSets()}.
}

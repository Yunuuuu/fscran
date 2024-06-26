% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runTSNE.R
\name{runTSNE}
\alias{runTSNE}
\alias{runTSNE.SingleCellExperiment}
\alias{runTSNE.Seurat}
\alias{runTSNE.default}
\title{Compute the t-stochastic neighbor embedding with \code{scran.chan}}
\usage{
runTSNE(object, ...)

\method{runTSNE}{SingleCellExperiment}(
  object,
  ...,
  dimred = "PCA",
  n_dimred = NULL,
  assay = NULL,
  name = "TSNE"
)

\method{runTSNE}{Seurat}(
  object,
  ...,
  dimred = "PCA",
  n_dimred = NULL,
  assay = NULL,
  layer = NULL,
  name = "TSNE"
)

\method{runTSNE}{default}(
  object,
  perplexity = 30L,
  interpolate = NULL,
  max_depth = 7L,
  max_iter = 500L,
  ...,
  approximate = TRUE,
  seed = 1234L,
  threads = NULL
)
}
\arguments{
\item{object}{A matrix-like object. Rows are features and columns are cells.}

\item{...}{\itemize{
\item \code{default} method: cannot be used.
\item \code{SingleCellExperiment} and \code{Seurat} method: additional arguments passed to
\code{default} methods.
}}

\item{dimred}{String or integer scalar specifying the existing dimensionality
reduction results to use.}

\item{n_dimred}{Integer scalar or vector specifying the dimensions to use if
dimred is specified.}

\item{assay}{Integer scalar or string indicating which assay of x
contains the expression values.}

\item{name}{String specifying the name to be used to store the result in the
\link[SingleCellExperiment:reducedDims]{reducedDims} or
\link[SeuratObject:Seurat-class]{reductions} of the output.}

\item{layer}{Name of the layer to get from the assay data.}

\item{perplexity}{Numeric scalar specifying the perplexity to use in the
t-SNE algorithm.}

\item{interpolate}{Integer scalar specifying the grid resolution for
interpolating repulsive forces (larger is slower but more accurate). A value
of zero disables interpolation, while a value of NULL only uses interpolation
for large datasets.}

\item{max_depth}{Integer scalar specifying the maximum depth of the
Barnes-Hut quad trees (larger is slower but more accurate). Set to a large
integer (e.g., 1000) to eliminate depth restrictions.}

\item{max_iter}{Integer scalar specifying the maximum number of iterations to
perform.}

\item{approximate}{Logical scalar specifying whether to perform an approximate neighbor search.}

\item{seed}{Integer scalar specifying the seed to use.}

\item{threads}{Integer scalar specifying the number of threads to use. If
\code{NULL}, all detected threads will be used. See
\link[parallel:detectCores]{detectCores}.}
}
\value{
\itemize{
\item \code{default} method: A numeric matrix where rows are cells and columns are
the two dimensions of the embedding.
\item \code{SingleCellExperiment} method: embedding was added into
\link[SingleCellExperiment:reducedDims]{reducedDims} named as \code{name}.
\item \code{Seurat} method: embedding was added into
\link[SeuratObject:Seurat-class]{reductions} named as \code{name}.
}
}
\description{
Compute the t-stochastic neighbor embedding with \code{scran.chan}
}
\seealso{
\link[scran.chan:runTSNE.chan]{runTSNE.chan}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRscores.r
\name{smoothScoreLR}
\alias{smoothScoreLR}
\title{Smooth existing LR-scores using spatial nearest neighbors}
\usage{
smoothScoreLR(object, score.use = NULL, eps = NULL, lrscore.name = NULL)
}
\arguments{
\item{object}{A CytoSignal object with existing LR-scores available in
"lrscore" slot.}

\item{score.use}{Name of LR-scores to be smoothed. Default use the most
recently calculated LR-scores.}

\item{eps}{The radius of the epsilon ball. Default use parameter inferred
with \code{\link{inferEpsParams}}. See Description as well.}

\item{lrscore.name}{Name of the smoothed LR-scores to be stored in the object.
Default appending \code{"_smoothed"} after \code{score.use}.}
}
\value{
Returns the original object with result updated in "lrscore" slot.
An "lrScores" object is created with its "score" slot filled a cell by
interaction sparse matrix indicating the smoothed LR-score. When a new GauEps
NN graph is inferred, it will be stored in the "imputation" slot. See
\code{\link{findNNGauEB}} for detail.
}
\description{
This function takes an existing LR-score out and smooth it using the spatial
nearest neighbor graph inferred with Gaussian-Epsilon method. If the radius
queried (\code{eps}) is the same as any pre-calculated NN graph during
imputation stage, this function will directly use it. However, a different
\code{eps} can also be used. In this case, a new NN graph will be inferred
and stored in the object.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRscores.r
\name{inferNullScoreLR}
\alias{inferNullScoreLR}
\title{Permute LR score for specific ligand-receptor imputation obj pairs}
\usage{
inferNullScoreLR(object, slot.use = NULL)
}
\arguments{
\item{object}{A Cytosignal object}

\item{slot.use}{Which LR score to use. Use the name specified with \code{tag}
when running \code{\link{inferLRScore}}.}
}
\value{
A Cytosignal object
}
\description{
This function is a follow-up function of inferScoreLR. It computes the NULL LR-scores
using the NULL imputation results and stores the results in the LR score object.
The null distribution of the LR scores can be used to test the significance of the LR scores.
}

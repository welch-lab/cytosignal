% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRscores.r
\name{inferScoreLR}
\alias{inferScoreLR}
\title{Compute the LR score for specific ligand-receptor imputation pairs}
\usage{
inferScoreLR(
  object,
  lig.slot,
  recep.slot,
  intr.db.name = c("diff_dep", "cont_dep"),
  tag = paste0(lig.slot, "-", recep.slot)
)
}
\arguments{
\item{object}{A Cytosignal object}

\item{lig.slot}{The ligand slot to use}

\item{recep.slot}{The receptor slot to use}

\item{intr.db.name}{The intr database name to use}

\item{tag}{Name of the result to be stored in object.}
}
\value{
A Cytosignal object
}
\description{
Compute the LR score for specific ligand-receptor imputation pairs
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_gene.R
\name{heatmap_GO}
\alias{heatmap_GO}
\title{Show significant genes across top GO term hits with coefficients from
regression analysis of an interaction}
\usage{
heatmap_GO(
  intrDEG,
  GO,
  intr,
  description.col = "description",
  pval.col = "pval",
  gene.col = "genes",
  gene.split.fun = function(x) unlist(strsplit(x, ",")),
  term.topN = 20,
  gene.topN = 20,
  binary_sign = FALSE,
  text.size = 10
)
}
\arguments{
\item{intrDEG}{A \code{CytosignalIntrDEG} object, output from
\code{\link{inferIntrDEG}}.}

\item{GO}{A data.frame object for GO enrichment analysis result.}

\item{intr}{A single interaction ID (starts with "CPI") or its numeric index
within the range of \code{intrDEG}.}

\item{description.col, pval.col, gene.col}{The column names of the data.frame
\code{GO} that contains the term description, p-value, and gene hit string.
Default \code{"description"}, \code{"pval"}, \code{"genes"}.}

\item{term.topN}{Use this number of top GO terms, ranked by p-values. Default
\code{20}.}

\item{gene.topN}{Use this number of top genes, ranked by absolute value of
coefficients. Default \code{20}.}

\item{binary_sign}{Whether to convert coefficient value to binary sign value.
Default \code{FALSE}.}

\item{color_num}{Number of colors in the heatmap. Can only use \code{2} or
\code{3}, Default \code{2} use white-red color palette. \code{3} use scaled
blue-white-red color palette.}
}
\description{
Create a heatmap for an interaction on GO terms by significant genes, colored
by the coefficients of the genes in the regression model returned by
\code{\link{inferIntrDEG}}. The GO term enrichment can be done with any tools
available. The input data.frame \code{GO} must contain fields for 1. term
description, character, pointed to by \code{description.col}, 2. p-value,
numeric, pointed to by \code{pval.col}, 3. gene hit string, character,
pointed to by \code{gene.col}. The gene hit string for each term must be form
in a way that function \code{gene.split.fun} can split it into a character
vector of gene names. For example, if a gene hit string is "gene1, gene2, gene3",
then \code{gene.split.fun} should be \code{function(x) unlist(strsplit(x, ", "))},
so that a split result \code{c("gene1", "gene2", "gene3")} can be correctly
obtained.
}

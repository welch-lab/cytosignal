% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.r
\name{rankIntrSpatialVar}
\alias{rankIntrSpatialVar}
\title{Rank the inferred high-quality interactions by their spatial variability}
\usage{
rankIntrSpatialVar(object, slot.use = NULL, numCores = 1, verbose = FALSE)
}
\arguments{
\item{object}{A \linkS4class{CytoSignal} object, with
\code{\link{inferSignif}} already run.}

\item{slot.use}{Which LR score to use. Use the name specified with \code{tag}
when running \code{\link{inferLRScore}}.}

\item{numCores}{SPARK::sparkx parameter, an integer specifying the number of
threads.}

\item{verbose}{SPARK::sparkx parameter, a logical value indicating whether to
print details for debug purpose}
}
\value{
The input \linkS4class{CytoSignal} object with the spatially variable
high-quality interaction list updated at
\code{object@lrscore[[slot.use]]@res.list$result.spx}
}
\description{
This function utilizes SPARK package to calculate the spatial variability of
the high-quality interactions, using their LR scores. Please refer to
\href{https://doi.org/10.1186/s13059-021-02404-0}{Jiaqiang Zhu, et al., 2021, Genome Biology}
for more details of the method.
}
\examples{
\dontrun{
object <- findNN(object)
object <- imputeLR(object)
object <- inferScoreLR(object, lig.slot = "GauEps", recep.slot = "Raw",
                       intr.db.name = "diff_dep")
object <- permuteLR(object)
object <- inferNullScoreLR(object)
object <- inferSignif(object)
object <- rankIntrSpatialVar(object)
}
}

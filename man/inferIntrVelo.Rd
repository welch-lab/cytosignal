% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{inferIntrVelo}
\alias{inferIntrVelo}
\title{Calculate the interaction velocity from the imputed time derivative of
\eqn{L} or \eqn{R} values}
\usage{
inferIntrVelo(object, recep.smooth = FALSE, norm.method = "scanpy")
}
\arguments{
\item{object}{A \linkS4class{CytoSignal} object, with
\code{\link{imputeVeloLR}} already run.}

\item{recep.smooth}{A logical scalar. Whether to use the smoothed \eqn{R}
values which is imputed with DT method. Default \code{FALSE}.}

\item{norm.method}{The normalization method to apply to the velocity data,
need to be consistent with the normalization method used when generating the
input RNA velocity. Please consult the tool used for it, e.g. scVelo,
VeloVAE. Default is \code{"scanpy"}. Can choose from \code{"scanpy"},
\code{"cpm"}, \code{"default"} or \code{"none"}.}
}
\value{
A \linkS4class{CytoSignal} object updated. Entries in
\code{object@lrvelo} slot will be updated with the velocity scores. When
\code{recep.smooth} is by default \code{FALSE}, the \code{object@lrvelo}
slot will be updated with \code{object@lrvelo$`GauEps-Raw`} and
\code{object@lrvelo$`DT-Raw`}. When \code{recep.smooth} is \code{TRUE}, the
\code{object@lrvelo} slot will be updated with
\code{object@lrvelo$`GauEps-DT`} and \code{object@lrvelo$`DT-DT`}.
}
\description{
After running \code{\link{imputeVeloLR}}, we can calculate the interaction
velocity for each location. Please refer to the manualscript for detail of
the calculation.
}
\examples{
\dontrun{
object <- addVelo(object, velo.s, velo.u)
object <- findNN(object)
object <- imputeVeloLR(object)
object <- inferIntrVelo(object)
}
}

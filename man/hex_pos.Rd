% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_utility.r
\name{hex_pos}
\alias{hex_pos}
\title{Find centre of closest hexagon.}
\usage{
hex_pos(x, y, width, height)
}
\arguments{
\item{x}{numeric x position}

\item{y}{numeric y position}

\item{width}{of hexagon}

\item{height}{of hexagon}
}
\value{
matrix giving position of closest hexagon center
}
\description{
Find centre of closest hexagon.
}
\examples{
x <- runif(1e4)
y <- runif(1e4)
res <- hex_pos(x, y, 0.5, 0.5)
plot(x, y, type = "n")
segments(x, y, res[, 1], res[, 2], col = "grey80")
points(unique(res), pch = 20, cex = 2)
}
\keyword{internal}

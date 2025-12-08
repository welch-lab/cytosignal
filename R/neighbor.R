#' Set parameters for correctly scaling the spatial coordinates
#' @description
#' Set parameters required for modeling the neighborhood of LR interactions. A
#' plot is by default shown to help users to visually check if the parameters
#' look appropriate.
#' @param object A \linkS4class{cytosignal2} object.
#' @param micronPerUnit Numeric, 1 spatial coordinate unit equals to how many
#' microns.
#' @param ballRadiusMicron Numeric, radius of the Epsilon ball in microns.
#' Default \code{200}.
#' @param thresh Numeric, threshold for the Gaussian kernel. Default
#' \code{0.001}.
#' @param plot Logical, whether to plot the scale bar after setting the
#' parameters. Default \code{TRUE}.
#' @return The input \linkS4class{cytosignal2} object with updated parameters
#' in \code{object@parameters}.
#' \itemize{
#'    \item{\code{micronPerUnit}: 1 spatial coordinate unit equals to how many
#'    microns.}
#'    \item{\code{ballRadiusMicron}: radius of the Epsilon ball in microns.}
#'    \item{\code{ballRadiusCoord}: radius of the Epsilon ball in spatial
#'    coordinate units.}
#'    \item{\code{sigma}: sigma used for Gaussian kernel in modeling the
#'    diffusion-dependent neighborhood.}
#' }
#' @export
setParams <- function(
        object,
        micronPerUnit,
        ballRadiusMicron = 200,
        thresh = 0.001,
        plot = TRUE
) {
    object@parameters[['micronPerUnit']] <- micronPerUnit
    object@parameters[['ballRadiusMicron']] <- ballRadiusMicron
    ballRadiusCoord <- ballRadiusMicron / micronPerUnit
    object@parameters[['ballRadiusCoord']] <- ballRadiusCoord
    # sigma used for Gaussian kernal in tech resolution
    object@parameters[['sigma']] <- ballRadiusCoord / sqrt(-2 * log(thresh))

    if (isTRUE(plot)) print(plotScale(object))
    return(object)
}


#' Build weighted neighborhood graph that models ligand received by spots
#' @rdname findNeighbor
#' @export
#' @seealso [setParams()]
#' @seealso [plotNeighbor()]
#' @description
#' \code{findNeighbor()} is a wrapper function that calls
#' \code{findNeighborGauEB2()} and \code{findNeighborDT2()}.
#'
#' Parameters \bold{must} be set with \code{\link{setParams}} \bold{before}
#' running these functions, in order to appropriately scale the spatial
#' coordinates and assign the distance limits.
#'
#' \code{findNeighborGauEB2()} uses Gaussian Epsilon ball method to detect
#' neighbors within the limited distance and assign weights based on a Gaussian
#' kernel. This models the diffusion-dependent ligand-receptor interactions.
#'
#' \code{findNeighborDT2()} uses Delaunay triangulation to detect neighbors that
#' are directly contacting each other within the limited distance and assigns
#' equal weights. This models the contact-dependent ligand-receptor
#' interactions.
#' @param object A `cytosignal2` object. With \code{\link{setParams}} called.
#' @return
#' A `cytosignal2` object with \code{neighborDiff} and/or \code{neighborCont}
#' slots filled with a sparse square matrix representing the neighborhood graph.
findNeighbor <- function(object) {
    object <- findNeighborGauEB2(object)
    object <- findNeighborDT2(object)
    return(object)
}

#' @rdname findNeighbor
#' @export
findNeighborGauEB2 <- function(
        object
) {
    spatial <- object@spatial
    eps <- object@parameters$ballRadiusCoord
    sigma <- object@parameters$sigma
    if (is.null(eps) || is.null(sigma)) {
        cli::cli_abort('Set parameters first with {.fn setParams}.')
    }
    distance <- select_EB_rcpp2(spatial, eps = eps)
    nNeighbor <- diff(distance@p)
    hasNeighborIdx <- nNeighbor > 0
    gauss_vec_inplace_cpp(distance@x, sigma)
    Matrix::diag(distance)[hasNeighborIdx] <- gauss_vec_cpp(1e-9, sigma) * 5
    distance <- normalizeSparse_cpp(distance)
    object@neighborDiff <- distance
    return(object)
}

#' @rdname findNeighbor
#' @export
findNeighborDT2 <- function(
        object
) {
    maxDist <- object@parameters$ballRadiusCoord
    spatial <- object@spatial
    if (is.null(maxDist)) {
        cli::cli_abort('Set parameters first with {.fn setParams}.')
    }
    DT <- RTriangle::triangulate(RTriangle::pslg(P = spatial))
    edges <- DT$E
    node1Loc <- spatial[edges[, 1],]
    node2Loc <- spatial[edges[, 2],]
    distance <- euclidean_elementwise_cpp(node1Loc, node2Loc)
    edgeValidIdx <- distance <= maxDist
    edges <- edges[edgeValidIdx, , drop = FALSE]

    neighbors <- Matrix::sparseMatrix(
        i = c(edges[, 1], edges[, 2]),
        j = c(edges[, 2], edges[, 1]),
        x = rep(1, 2*nrow(edges)),
        dims = c(nrow(spatial), nrow(spatial)),
        dimnames = list(NULL, rownames(spatial))
    )
    hasNeighborIdx <- diff(neighbors@p) > 0
    neighbors <- normalizeSparse_cpp(neighbors)
    Matrix::diag(neighbors)[hasNeighborIdx] <- 1
    object@neighborCont <- neighbors
    return(object)
}

#' Visualize neighbor weights of a spot
#' @rdname plotNeighbor
#' @export
#' @seealso [findNeighbor()]
#' @description
#' This function makes a spatial coordinate plot colored by the weights of a
#' specified spot in the selected type of neighborhood graph. By default,
#' the spot with the most neighbors is selected.
#' @param object A `cytosignal2` object. With \code{\link{findNeighbor}} called.
#' @param type Type of neighborhood graph to visualize. Either
#' \code{'diffusion'} or \code{'contact'}.
#' @param index Index of the spot to visualize. Default \code{NULL} selects
#' the spot with the most neighbors. Can be an integer index, character ID, or
#' logical selector vector.
#' @return
#' A ggplot2 object showing the spatial coordinates colored by the neighbor
#' weights of the selected spot.
plotNeighbor <- function(
        object,
        type = c('diffusion', 'contact'),
        index = NULL
) {
    type <- match.arg(arg = type, several.ok = FALSE)
    graph <- switch(
        type,
        'diffusion' = object@neighborDiff,
        'contact' = object@neighborCont
    )
    if (is.null(graph)) {
        cli::cli_abort('Run {.fn findNeighbor} first.')
    }
    if (is.null(index)) {
        nNeighbors <- diff(graph@p)
        index <- which.max(nNeighbors)
        spotText <- 'the spot with the most neighbors'
    } else {
        index <- .checkValid.Index(index, colnames(object@rawData), 1)
        spotText <- colnames(object@rawData)[index]
    }
    title <- sprintf('%s-dependent neighbor weights of\n%s', type, spotText)
    weights <- graph[, index]
    object$weights <- weights
    plotSpatial(object, colorBy = 'weights', title = title, viridisOption = 'D')
}

#' @rdname plotNeighbor
#' @export
plotNeighborDiff <- function(
        object,
        index = NULL
) {
    plotNeighbor(object, type = 'diffusion', index = index)
}

#' @rdname plotNeighbor
#' @export
plotNeighborCont <- function(
        object,
        index = NULL
) {
    plotNeighbor(object, type = 'contact', index = index)
}

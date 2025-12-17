countEdges <- function(
        object,
        lrscore.use = NULL,
        intr.use = NULL
) {
    lrslotAvail <- names(object@lrscore)
    lrslotAvail <- lrslotAvail[lrslotAvail != 'default']
    if (!is.null(lrscore.use)) {
        lrscore.use <- rlang::arg_match(lrscore.use, lrslotAvail, multiple = TRUE)
    } else {
        if (any(grepl('_smooth$', lrslotAvail))) {
            lrscore.use <- lrslotAvail[grepl('_smooth$', lrslotAvail)]
        } else {
            warning("Neither is `lrscore.use` specified, nor are there any smoothed score. Using all available inference.")
            lrscore.use <- lrslotAvail
        }
    }
    clusterVar <- object@clusters
    if (is.null(clusterVar)) {
        stop("No cluster set at `object@clusters`")
    }
    mtxList <- list()
    for (i in seq_along(lrscore.use)) {
        scorename <- lrscore.use[i]
        res.list <- object@lrscore[[scorename]]@res.list[["result.hq"]]
        lig.slot <- object@lrscore[[scorename]]@lig.slot
        recep.slot <- object@lrscore[[scorename]]@recep.slot
        intr.list <- names(res.list)
        nn.graph <- object@imputation[[lig.slot]]@nn.graph

        barcodes <- colnames(object@raw.counts)
        # if no intr.use provided, use all intrs in the db
        intrs <- intr.use %||% intr.list

        # all intrs in intr.use must be in intr.list
        intrs <- intrs[intrs %in% intr.list]
        if (length(intrs) == 0) {
            warning("No interactions in `intr.use` found for ", scorename)
            next
        }

        # Now building one-hot matrix noting which receiver cell has which intr
        all_counts <- lapply(seq_along(res.list), function(i) {
            intr.name <- names(res.list)[i]
            res.df <- data.frame(receiver = res.list[[i]], intr = intr.name, count = 1)
        })
        all_counts <- Reduce(rbind, all_counts)
        all_counts$receiver <- match(all_counts$receiver, barcodes)
        all_counts$intr <- match(all_counts$intr, intr.list)
        receiver.intr.mtx <- sparseMatrix(
            i = all_counts$receiver,
            j = all_counts$intr,
            x = all_counts$count,
            dims = c(length(barcodes), length(intr.list)),
            dimnames = list(barcodes, intr.list)
        )

        nIntrPerReceiver <- rowSums(receiver.intr.mtx)
        # One-hot encoding of cluster labels
        cluster.cell.mtx <- fac2sparse(clusterVar)
        # One-hot encoding of sender-receiver graph. A one at graph[i,j] means cell
        # i sends to cell j
        sender.receiver.mtx <- nn.graph
        sender.receiver.mtx@x <- rep(1, length(sender.receiver.mtx@x))
        # Amount of edge to count from each sending cluster to each receiving cell
        cluster.receiver.mtx <- cluster.cell.mtx %*% sender.receiver.mtx
        # Multiply each row by nIntr per receiver cell, then summarize by
        # receiving clusters
        cluster.cluster.NIntr <- sweep(
            x = cluster.receiver.mtx,
            MARGIN = 2,
            STATS = nIntrPerReceiver,
            FUN = "*"
        ) %*%
            t(cluster.cell.mtx)
        # Now compute total number of edges between clusters if all interactions
        # are valid.
        cluster.cluster.totalIntr <- colSums(sender.receiver.mtx %*% t(cluster.cell.mtx))*length(intr.list)
        # Finally, get the fraction of actual edges drawn over total possible edges
        mtxList[[scorename]] <- sweep(cluster.cluster.NIntr, 2, cluster.cluster.totalIntr, '/')
    }
    return(mtxList)
}

#' Make Circos plot summarising number of interactions between clusters
#' @description
#' This function counts all the edges from the neighborhood graph connecting
#' between cluster pairs, weighted by the number of significant
#' interactions detected in the receiver cells. This is then normalized by the
#' total number of possible edges between the cluster pairs as is all
#' interactions are significant.
#' Users can make multiple circos plots for different types of LR scores (e.g.
#' diffusion or contact dependent) at the same time, combined by grid.
#'
#' This function behaves a bit differently when showing a single selected type
#' of interactions versus multiple types. When multiple plots are intended, by
#' either explicitly selecting multiple \code{lrscore.use} or by default when
#' multiple types of inference are available, the plots are automatically
#' arranged in a grid layout. In single-plot mode, the plot is filled into the
#' current graphics device. That means users can globally make layout settings
#' before calling and combine with other plots. See examples.
#' @param object A \code{\linkS4class{cytosignal}} object with significant
#' interactions inferred with \code{\link{inferSignif}}.
#' @param lrscore.use Which LR score(s) to use. Use the name specified with
#' \code{lrscore.name} when running \code{\link{inferScoreLR}}. Default
#' \code{NULL} uses all available smoothed scores.
#' @param intr.use Which interactions to use focus on. Default \code{NULL} uses
#' all significant interactions inferred for the selected \code{lrscore.use}.
#' @param colors.list A named list of colors for the clusters. Default
#' \code{NULL} uses built-in color palette.
#' @param title.size,label.size Font size of the title and cluster labels,
#' respectivey. Default \code{12} and \code{8}.
#' @param circle.margin A positive numeric vector of 4 values for the left,
#' right, bottom, top margin sizes of the main circos panel. Usually the top
#' (4th) value should be larger to accommodate the title. Default
#' \code{c(0.01, 0.01, 0.01, 0.2)}.
#' @return NULL value returned. Circos plots are drawn.
#' @export
#' @examples
#' \dontrun{
#' # Assuming `cs` is a cytosignal object after calling `inferSignif()`
#' # Simply show the count of diffusion-dependent interactions between cluster
#' # pairs.
#' plotCircosNIntr(cs, 'diffusion-Raw_smooth')
#'
#' # Show all the count of all types of interactions. By default, an additional
#' # inference called "Raw-Raw" is also included for a contact-dependent
#' # interactions
#' plotCircosNIntr(cs)
#'
#' # Show selected types of interactions only
#' plotCircosNIntr(cs, lrscore.use = c('diffusion-Raw_smooth',
#'                                    'diffusion-Env_smooth'))
#'
#' # Customized layout for combining with other R graphics or re-arranging the
#' # plots for multiple types of interactions
#' par(mfrow = c(1,2))
#' plotCircosNIntr(cs, 'contact-Raw_smooth')
#' plotCircosNIntr(cs, 'diffusion-Raw_smooth')
#' }
#'
plotCircosNIntr <- function(
        object,
        lrscore.use = NULL,
        intr.use = NULL,
        colors.list = NULL,
        title.size = 12,
        label.size = 8,
        circle.margin = c(0.01, 0.01, 0.01, 0.2)
) {
    if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required for this function. Please install it with \ninstall.packages('circlize')",
            call. = FALSE
        )
    }
    mtxList <- countEdges(
        object = object,
        lrscore.use = lrscore.use,
        intr.use = intr.use
    )
    colors.list <- colors.list %||%
        setNames(
            uniqueColors(nlevels(object@clusters)),
            levels(object@clusters)
        )
    nPlots <- length(mtxList)
    if (nPlots == 0) {
        stop("Failed to generate any plot. Please check warning messages")
    }

    nr <- ceiling(sqrt(nPlots))
    nc <- ceiling(nPlots/nr)

    userpar <- graphics::par(no.readonly = TRUE)
    if (nPlots > 1) {
        graphics::par(mfrow = c(nr, nc))
        on.exit(graphics::par(mfrow = userpar$mfrow), add = TRUE)
    }
    graphics::par(
        cex = 1,
        cex.main = title.size/12
    )
    # Only reset what we changed. Otherwise it'll trigger a reset of layout
    # setting which is not desired.
    on.exit(graphics::par(
        cex = userpar$cex,
        cex.main = userpar$cex.main
    ), add = TRUE)
    for (i in seq_along(mtxList)) {
        scorename <- names(mtxList)[i]
        edge.mtx <- as.matrix(mtxList[[i]])
        colors <- colors.list[c(rownames(edge.mtx), colnames(edge.mtx))]
        rownames(edge.mtx) <- paste0("s-", rownames(edge.mtx))
        colnames(edge.mtx) <- paste0("r-", colnames(edge.mtx))
        names(colors) <- c(rownames(edge.mtx), colnames(edge.mtx))

        # set the order based on the rowSums and colSums
        col.order <- colnames(edge.mtx)[order(colSums(edge.mtx), decreasing = T)]
        row.order <- rownames(edge.mtx)[order(rowSums(edge.mtx), decreasing = F)]
        use.order <- c(row.order, col.order)
        circlize::circos.par(
            "track.height" = 0.8,
            circle.margin = circle.margin # L, R, B, T
        )
        circlize::chordDiagram(
            x = edge.mtx,
            big.gap = 15,
            grid.col = colors,
            annotationTrack = "grid",
            order = use.order,
            annotationTrackHeight = c(0.03, 0.01),
            preAllocateTracks = list(
                track.height = max(graphics::strwidth(unlist(dimnames(edge.mtx))))/2
            )
        )
        circlize::circos.track(
            track.index = 1,
            panel.fun = function(x, y) {
                circlize::circos.text(
                    circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1],
                    circlize::CELL_META$sector.index, facing = "clockwise",
                    niceFacing = TRUE, adj = c(0, 0.5), cex = label.size/12
                )
            },
            bg.border = NA
        )
        graphics::title(main = scorename, line = -1)
        circlize::circos.clear()
    }
    return(invisible(NULL))
}

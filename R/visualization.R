countEdges2 <- function(
        object,
        type = NULL,
        splitType = FALSE,
        csdb = NULL,
        fdrThresh = 0.05
) {
    csdb <- subsetCSDB(csdb, type = type)
    dblist <- split(csdb, csdb$type)
    dblist <- dblist[sapply(dblist, nrow) > 0]
    # Sanity check for cluster variable
    clusterVarname <- object@parameters$cluster
    if (is.null(clusterVarname)) {
        cli::cli_abort(c(
            x = 'No cluster variable set.',
            i = 'See {.code ?setCluster}.'
        ))
    }
    clusterVar <- object@metadata[[clusterVarname]]
    if (is.null(clusterVar)) {
        cli::cli_abort("Selected cluster variable {.val {clusterVarname}} is not available in metadata")
    }

    # Pull data to be reused
    fdr <- object@significance$spatialFDR
    if (is.null(fdr)) {
        cli::cli_abort('Spatial FDR not found. Please run {.fn inferLRScore} first.')
    }
    # One-hot encoding of cluster labels
    clusterCellMtx <- fac2sparse(clusterVar)
    # Count number of edges by interaction type, since different graph will be used
    mtxList <- list()
    for (i in seq_along(dblist)) {
        type <- names(dblist)[i]
        subdb <- dblist[[i]]
        if (nrow(subdb) == 0) next
        if (type == 'contact') senderReceiverMtx <- object@neighborCont
        else if (type == 'diffusion') senderReceiverMtx <- object@neighborDiff
        else {
            cli::cli_warn("Unknown interaction type {.val {type}}. Ask maintainer to fix it.")
            next
        }

        # Now building one-hot matrix noting which receiver cell has which intr
        receiverIntrMtx <- as(fdr[, subdb$interactors, drop = FALSE] < fdrThresh, "dgCMatrix")
        nIntrPerReceiver <- rowSums(receiverIntrMtx, na.rm = TRUE)
        # One-hot encoding of sender-receiver graph. A one at graph[i,j] means cell
        # i sends to cell j
        senderReceiverMtx@x <- rep(1, length(senderReceiverMtx@x))
        # Amount of edge to count from each sending cluster to each receiving cell
        clusterReceiverMtx <- clusterCellMtx %*% senderReceiverMtx
        # Multiply each row by nIntr per receiver cell, then summarize by
        # receiving clusters
        clusterClusterNIntr <- sweep(
            x = clusterReceiverMtx,
            MARGIN = 2,
            STATS = nIntrPerReceiver,
            FUN = "*"
        ) %*%
            t(clusterCellMtx)
        # Now compute total number of edges between clusters if all interactions
        # are valid.
        if (isTRUE(splitType)) {
            nIntr <- nrow(subdb)
        } else {
            nIntr <- nrow(csdb)
        }
        clusterClusterTotalIntr <- colSums(senderReceiverMtx %*% t(clusterCellMtx))*nIntr
        # Finally, get the fraction of actual edges drawn over total possible edges
        mtxList[[type]] <- as.matrix(sweep(clusterClusterNIntr, 2, clusterClusterTotalIntr, '/'))
    }

    if (isTRUE(splitType)) {
        return(mtxList)
    } else {
        if (length(dblist) == 1) return(mtxList)
        else return(list(all = Reduce('+', mtxList)))
    }
}

#' Make Circos plot summarising number of interactions between clusters
#' @description
#' This function counts all the edges from the neighborhood graph connecting
#' between cluster pairs, weighted by the number of significant
#' interactions detected in the receiver cells. This is then normalized by the
#' total number of possible edges between the cluster pairs as if all
#' interactions are significant.
#' @details
#' \code{splitType} controls if all interactions are counted together in one
#' plot or counted separately and shown in multiple plots. When splitting by
#' type, the count is normalized by total number of interactions of that type
#' only.
#'
#' When making a single plot, users can make use of R's graphic layout system
#' to combine with other plots or arrange the plot in a customized way. See
#' examples. When making multiple plot in one call, the plots are arranged by
#' default and the function ignores any global layout settings.
#' @param object A \code{\linkS4class{cytosignal}} object with significant
#' interactions inferred with \code{\link{inferSignif}}.
#' @param type Selecting the type of interaction to show. Default \code{NULL}
#' considers both options of \code{'contact'} and \code{'diffusion'}.
#' @param splitType Whether to split the interactions by type and show multiple
#' plots. Default \code{FALSE}. See details.
#' @param interaction Which interactions to use focus on. Default \code{NULL}
#' uses all significant interactions selected by \code{\link{getTopIntr}}.
#' Otherwise, users can provide a vector of interaction selection (name or
#' indices applicable to \code{intrDB(object)}), or a pre-filtered \code{csdb}
#' object (e.g. non-default \code{getTopIntr} result).
#' @param fdrThresh,minExp,minSpot Parameters passed to \code{\link{getTopIntr}}
#' for filtering significant interactions when \code{interaction = NULL}.
#' \code{fdrThresh} is additionally applied to claim significant interactions in
#' each spot.
#' @param colors A vector of colors for the clusters. If \code{NULL}, a built-in
#' color palette is used. Default \code{NULL}.
#' @param titles Customized title texts for the subplots. Default \code{NULL}
#' decides the titles automatically. When multiple plots are generated, the
#' lenght of \code{titles} must match the number of plots.
#' @param titleTextSize,label.size Font size of the title and cluster labels,
#' respectivey. Default \code{12} and \code{8}.
#' @param circleMargin A positive numeric vector of 4 values for the left,
#' right, bottom, top margin sizes of the main circos panel. Usually the top
#' (4th) value should be larger to accommodate the title. Default
#' \code{c(0.01, 0.01, 0.01, 0.2)}.
#' @return NULL value returned. Circos plots are drawn.
#' @export
#' @examples
#' \dontrun{
#' # Assuming `cs` is a cytosignal2 object after calling `inferLRScore()`
#' # Summarize the count of all types of interactions in one plot.
#' plotCircosNIntr2(cs)
#'
#' # Show the count of only diffusion-dependent interactions between cluster
#' # pairs.
#' plotCircosNIntr2(cs, 'diffusion')
#'
#' # Show the type of both types of interactions in separate panels
#' plotCircosNIntr2(cs, splitType = TRUE)
#'
#' # Customized layout for combining with other R graphics or re-arranging the
#' # plots for multiple types of interactions
#' par(mfrow = c(1,2))
#' plotCircosNIntr2(cs, 'diffusion')
#' plot(1:10, 1:10)
#' }
#'
plotCircosNIntr2 <- function(
        object,
        type = NULL,
        splitType = FALSE,
        interaction = NULL,
        fdrThresh = 0.05,
        minExp = 100,
        minSpot = 100,
        colors = NULL,
        titles = NULL,
        titleTextSize = 12,
        labelTextSize = 8,
        circleMargin = c(0.01, 0.01, 0.01, 0.2)
) {
    if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required for this function. Please install it with \ninstall.packages('circlize')",
             call. = FALSE
        )
    }
    csdb <- .checkIntrSelection(
        object = object,
        interaction = interaction,
        fdrThresh = fdrThresh,
        minExp = minExp,
        minSpot = minSpot,
        error = TRUE
    )
    mtxList <- countEdges2(
        object = object,
        type = type,
        splitType = splitType,
        csdb = csdb,
        fdrThresh = fdrThresh
    )
    clusterVarname <- object@parameters$cluster
    clusterVar <- object@metadata[[clusterVarname]]
    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
    nCluster <- nlevels(clusterVar)
    if (is.null(colors)) {
        colors <- csColors[seq_len(nCluster)]
    } else {
        if (length(colors) < nCluster) {
            cli::cli_abort(
                "Need {nCluster} colors while {length(colors)} {?was/were} provided."
            )
        }
        colors <- colors[seq_len(nCluster)]
    }
    names(colors) <- levels(clusterVar)

    nPlots <- length(mtxList)
    if (nPlots == 0) {
        stop("Failed to generate any plot. Please check warning messages")
    }
    if (is.null(titles)) {
        if (nPlots > 1) titles <- names(mtxList)
        if (nPlots == 1) {
            if (is.null(interaction) && isFALSE(splitType)) {
                titles <- "All interactions"
            } else if (nrow(csdb) == 1) {
                titles <- csdb$interactors[1]
            } else {
                titles <- "Selected interactions"
            }
        }
    } else {
        if (length(titles) < nPlots) {
            cli::cli_abort(
                "{nPlots} subplot{?s} are generated while only {length(titles)} title{?s} {?is/are} provided."
            )
        }
        titles <- titles[seq_len(nPlots)]
    }

    nc <- ceiling(sqrt(nPlots))
    nr <- ceiling(nPlots/nc)

    userpar <- graphics::par(no.readonly = TRUE)
    if (nPlots > 1) {
        graphics::par(mfrow = c(nr, nc))
        on.exit(graphics::par(mfrow = userpar$mfrow), add = TRUE)
    }
    graphics::par(
        cex = 1,
        cex.main = titleTextSize/12
    )
    # Only reset what we changed. Otherwise it'll trigger a reset of layout
    # setting which is not desired.
    on.exit(graphics::par(
        cex = userpar$cex,
        cex.main = userpar$cex.main
    ), add = TRUE)
    for (i in seq_along(mtxList)) {
        scorename <- names(mtxList)[i]
        edgeMtx <- as.matrix(mtxList[[i]])
        colorsExtend <- colors[c(rownames(edgeMtx), colnames(edgeMtx))]
        rownames(edgeMtx) <- paste0("s-", rownames(edgeMtx))
        colnames(edgeMtx) <- paste0("r-", colnames(edgeMtx))
        names(colorsExtend) <- c(rownames(edgeMtx), colnames(edgeMtx))

        # set the order based on the rowSums and colSums
        col.order <- colnames(edgeMtx)[order(colSums(edgeMtx), decreasing = T)]
        row.order <- rownames(edgeMtx)[order(rowSums(edgeMtx), decreasing = F)]
        use.order <- c(row.order, col.order)
        circlize::circos.par(
            "track.height" = 0.8,
            circle.margin = circleMargin # L, R, B, T
        )
        circlize::chordDiagram(
            x = edgeMtx,
            big.gap = 15,
            grid.col = colorsExtend,
            annotationTrack = "grid",
            order = use.order,
            annotationTrackHeight = c(0.03, 0.01),
            preAllocateTracks = list(
                track.height = max(graphics::strwidth(unlist(dimnames(edgeMtx))))/2
            )
        )
        circlize::circos.track(
            track.index = 1,
            panel.fun = function(x, y) {
                circlize::circos.text(
                    circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1],
                    circlize::CELL_META$sector.index, facing = "clockwise",
                    niceFacing = TRUE, adj = c(0, 0.5), cex = labelTextSize/12
                )
            },
            bg.border = NA
        )
        graphics::title(main = titles[i], line = -1)
        circlize::circos.clear()
    }
    return(invisible(NULL))
}



#' Combined plots of all information for selected interactions
#' @description
#' This function generates combined plots showing spatial expression of
#' ligand and receptor genes, spatial LR score, and cluster annotation for
#' selected interactions.
#' @inheritParams generalParam
#' @param edge,velo Not implemented yet
#' @param legendNCol,legendNRow Controls the layout of cluster legend. When only
#' two genes are to be shown (one ligand and one receptor), only
#' \code{legendNCol} is used to set the number of columns in the legend. When
#' more than two genes are to be shown (complex ligand/receptor), only
#' \code{legendNRow} is used to set the number of rows in the legend.
#' @inheritParams plotSpatial
#' @return A \code{\link[patchwork]{patchwork}} object when one interaction is
#' provided; A list of \code{\link[patchwork]{patchwork}} objects when multiple
#' interactions are provided.
#' @export
plotIntrSummary <- function(
        object,
        interaction,
        fdrThresh = 0.05,
        edge = FALSE,
        velo = FALSE,
        legendNCol = NULL,
        legendNRow = NULL,
        colors = NULL,
        legendDotSize = 4,
        dotSize = 0.5,
        dotAlpha = 1,
        paletteOption = 'C',
        paletteDirection = -1,
        zeroAsNA = TRUE,
        naColor = 'grey80',
        titleTextSize = 12,
        subtitleTextSize = 10,
        legendTextSize = 8
) {
    if (is.null(interaction)) {
        cli::cli_abort('NULL selection is not allowed.')
    }
    csdb <- .checkIntrSelection(object, interaction, error = FALSE)

    ligands <- strsplit(csdb$ligands, ';')
    receptors <- strsplit(csdb$receptors, ';')
    plist <- list()
    nSignif <- nSignifSpots(object@significance$spatialFDR, fdrThresh)[csdb$interactors]
    for (i in seq_len(nrow(csdb))) {
        mode <- 0
        intr <- csdb$interactors[i]
        ligGenes <- ligands[[i]]
        pLigs <- lapply(ligGenes, function(gene) {
            plotSpatialGene(
                object,
                gene = gene,
                dotSize = dotSize,
                dotAlpha = dotAlpha,
                paletteOption = paletteOption,
                paletteDirection = paletteDirection,
                zeroAsNA = zeroAsNA,
                naColor = naColor,
                titleTextSize = titleTextSize,
                subtitleTextSize = subtitleTextSize,
                legendTextSize = legendTextSize
            ) + labs(title = paste0('Ligand gene: ', gene))
        })
        recGenes <- receptors[[i]]
        pRecep <- lapply(recGenes, function(gene) {
            plotSpatialGene(
                object,
                gene = gene,
                dotSize = dotSize,
                dotAlpha = dotAlpha,
                paletteOption = paletteOption,
                paletteDirection = paletteDirection,
                zeroAsNA = zeroAsNA,
                naColor = naColor,
                titleTextSize = titleTextSize,
                subtitleTextSize = subtitleTextSize,
                legendTextSize = legendTextSize
            ) + labs(title = paste0('Receptor gene: ', gene))
        })
        pGenes <- c(pLigs, pRecep)
        nGene <- length(pGenes)
        pIntr <- plotSpatialLRScore(
            object = object,
            interaction = intr,
            dotSize = dotSize,
            dotAlpha = dotAlpha,
            paletteOption = paletteOption,
            paletteDirection = paletteDirection,
            zeroAsNA = zeroAsNA,
            naColor = naColor,
            titleTextSize = titleTextSize,
            subtitleTextSize = subtitleTextSize,
            legendTextSize = legendTextSize
        )
        pIntr$labels$subtitle <- paste0(
            pIntr$labels$subtitle,
            sprintf('\nSignificant in %d spots (FDR < %.3f)', nSignif[i], fdrThresh)
        )
        if (nGene == 2) {
            pCluster <- plotSpatialMetadata(
                object,
                colors = colors,
                legendDotSize = legendDotSize,
                dotSize = dotSize,
                dotAlpha = dotAlpha,
                paletteOption = paletteOption,
                paletteDirection = paletteDirection,
                zeroAsNA = zeroAsNA,
                naColor = naColor,
                titleTextSize = titleTextSize,
                subtitleTextSize = subtitleTextSize,
                legendTextSize = legendTextSize
            ) +
                labs(title = NULL) +
                guides(colour = guide_legend(
                    title = object@parameters$cluster,
                    override.aes = list(size = legendDotSize),
                    ncol = legendNCol %||% 1,
                    theme = theme(legend.title = element_text(
                        size = titleTextSize,
                        face = 'bold'
                    ))
                ))
            lgd <- cowplot::get_legend(pCluster, legend = 'right')
            pCluster <- pCluster + theme(legend.position = 'none')
            pComb <- ((pIntr | pCluster) / (pGenes[[1]] | pGenes[[2]]) | lgd) +
                patchwork::plot_layout(widths = c(4, 1))
        } else {
            pCluster <- plotSpatialMetadata(
                object = object,
                colors = colors,
                legendDotSize = legendDotSize,
                dotSize = dotSize,
                dotAlpha = dotAlpha,
                paletteOption = paletteOption,
                paletteDirection = paletteDirection,
                zeroAsNA = zeroAsNA,
                naColor = naColor,
                titleTextSize = titleTextSize,
                subtitleTextSize = subtitleTextSize,
                legendTextSize = legendTextSize
            ) +
                labs(title = NULL) +
                guides(colour = guide_legend(
                    title = object@parameters$cluster,
                    override.aes = list(size = legendDotSize),
                    nrow = legendNRow %||% 8,
                    theme = theme(legend.title = element_text(
                        size = titleTextSize,
                        face = 'bold'
                    ))
                ))
            lgd <- cowplot::get_legend(pCluster, legend = 'right')
            pCluster <- pCluster + theme(legend.position = 'none')
            topList <- list(pIntr, pCluster, lgd)
            topList <- c(topList, rep(list(patchwork::plot_spacer()), nGene - 3))
            pComb <- patchwork::wrap_plots(c(topList, pGenes), nrow = 2)
        }
        plist[[intr]] <- pComb

    }
    if (length(plist) == 1) return(plist[[1]])
    else return(plist)
}


#' 3D plot showing edges connecting sender and receiver spots for selected interactions
#' @description
#' This function selects the receiver spots with significant LRscores for
#' the selected interactions, and connects them to corresponding neighbor sender
#' spots with positive ligand expression.
#'
#' A spatial coordinate panel colored by cluster is shown at bottom and top
#' plains of the 3D plot, with solid dots representing selected sender/receiver
#' spots, and translucent dots representing other spots. The edges connecting
#' the sender and receiver spots are shown as segments between the two plains.
#' The edges are downsampled for clearer visualization.
#' @inheritParams generalParam
#' @param showTitle,showSubTitle,showLegend Logical, whether to show title,
#' subtitle, and legend in the plot. Default \code{TRUE}.
#' @param showBox Logical, whether to show the box frame lines. Axis arrows and
#' labels are also hidden when set to \code{FALSE}. Default \code{TRUE}.
#' @param edgeLocalMinSize Minimum number of edges to sample from each local
#' spatial cluster of receiver spots. This balances the downsampling when
#' satellite interaction activities are present. Default \code{10}.
#' @param edgeSampleRate Proportion of edges to sample from major interaction
#' activity areas. Default \code{0.05}.
#' @param theta 3D plot horizontal rotation angle, positive value rotate
#' clockwise when looking from top. Default \code{-17}.
#' @param phi 3D plot vertical rotation angle, positive value rotate clockwise
#' when looking from left side, or say the side closer to users goes downward.
#' Default \code{30}.
#' @param boxHeight Height of the 3D box, a value between 0 and 1. Default
#' \code{0.7}.
#' @param perspTrans A value which can be used to vary the strength of the
#' perspective transformation. Values greater than 1 will lessen the perspective
#' effect and values less and 1 will exaggerate it. Default \code{5}.
#' @param dotSizeSignif,dotAlphaSignif Dot size and transparency for spots
#' to be highlighted. Default \code{0.3} and \code{1} (non-transparent).
#' @param dotSizeOther,dotAlphaOther Dot size and transparency for other spots.
#' Default \code{0.16} and \code{0.3} (mostly transparent).
#' @param legendTitleSize Font size for the legend title. Default \code{10}.
#' @param labelTextSize Font size for the z-axis label. Default \code{8}.
#' @param edgeColor Color of the edges. Default \code{'grey'}.
#' @param edgeLinewidth Line width of the edges. Default \code{0.2}.
#' @param seed Random seed for edge sampling for reproducibility. Default
#' \code{NULL} does not set seed. Setting a seed does not interfere with R's
#' global random state.
#' @inheritParams plotSpatial
#' @return A \code{plot3D} plist object when one interaction is provided; A list
#' of \code{plot3D} plist objects when multiple interactions are provided. A
#' plist object can be printed in console and shown in viewer, or saved as RDS
#' file for later use.
#' @export
plotEdge2 <- function(
        object,
        interaction,
        colors = NULL,
        fdrThresh = 0.05,
        showTitle = TRUE,
        showSubTitle = TRUE,
        showLegend = TRUE,
        showBox = TRUE,
        edgeLocalMinSize = 10,
        edgeSampleRate = 0.05,
        theta = -17,
        phi = 30,
        boxHeight = 0.7,
        perspTrans = 5,
        dotSizeSignif = 0.3,
        dotAlphaSignif = 1,
        dotSizeOther = 0.16,
        dotAlphaOther = 0.3,
        edgeColor = 'grey',
        edgeLinewidth = 0.2,
        titleTextSize = 12,
        subtitleTextSize = 10,
        legendTitleSize = 10,
        legendTextSize = 8,
        legendNCol = NULL,
        labelTextSize = 8,
        seed = NULL
) {
    csdb <- .checkIntrSelection(
        object,
        interaction = interaction,
        fdrThresh = fdrThresh,
        error = FALSE
    )
    if (nrow(csdb) == 0) {
        cli::cli_abort('No valid interaction selected')
    }

    fdr <- object@significance$spatialFDR
    if (is.null(fdr)) {
        cli::cli_abort('Spatial FDR not available. Please run {.fn inferLRScore} first.')
    }

    clusterVarname <- object@parameters$cluster
    if (is.null(clusterVarname)) {
        cli::cli_alert_warning(c(
            x = 'No cluster variable set. Setting all spots to the same color.',
            i = 'See {.code ?setCluster}.'
        ))
        object@metadata$.__all__ <- factor(rep('All', nrow(object@metadata)))
        clusterVarname <- '.__all__'
        showLegend <- FALSE
    }
    if (!clusterVarname %in% colnames(object@metadata)) {
        cli::cli_abort("Selected cluster variable {.val {clusterVarname}} is not available in metadata")
    }
    clusterVar <- object@metadata[[clusterVarname]]
    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
    colors <- colors %||% csColors
    uniqColors <- colors[seq_len(nlevels(clusterVar))]
    colors <- colors[as.integer(clusterVar)]

    spatial <- object@spatial
    xlim <- range(spatial[,1])
    xlim <- c(xlim[1] - 0.05 * diff(xlim), xlim[2] + 0.05 * diff(xlim))
    ylim <- range(spatial[,2])
    ylim <- c(ylim[1] - 0.05 * diff(ylim), ylim[2] + 0.05 * diff(ylim))

    plist <- list()

    legendParams <- list(
        cex = legendTextSize/12,
        title.cex = legendTitleSize/12,
        text = levels(clusterVar),
        col = uniqColors,
        title = clusterVarname,
        ncol = legendNCol %||% ceiling(length(uniqColors)/25)
    )
    for (i in seq_len(nrow(csdb))) {
        intr <- csdb$interactors[i]
        intrType <- as.character(csdb$type[i])
        receiverIdx <- fdr[,intr] < fdrThresh
        receiverIdx[is.na(receiverIdx)] <- FALSE
        if (csdb$type[i] == 'contact') graph <- object@neighborCont
        else if (csdb$type[i] == 'diffusion') graph <- object@neighborDiff
        else {
            cli::cli_warn("Unknown interaction type {.val {csdb$type[i]}} for {.val {intr}}. Ask maintainer to fix it.")
            next
        }
        dimnames(graph) <- list(colnames(object@rawData), colnames(object@rawData))
        senderIdxUniq <- unique(graph[,receiverIdx, drop = FALSE]@i + 1)
        ligGenes <- strsplit(csdb$ligands[i], ';')[[1]]
        ligGenesExp <- object@rawData[ligGenes, senderIdxUniq, drop = FALSE]
        if (intrType == 'diffusion') {
            # Require at least one ligand component exists in neighborhood
            senderIdx <- colnames(ligGenesExp)[colSums(ligGenesExp) > 0]
        } else (
            # Require ligand complex members co-express in a sender spot
            senderIdx <- colnames(ligGenesExp)[colSums(ligGenesExp > 0) == nrow(ligGenesExp)]
        )

        withr::with_seed(
            seed = seed,
            {
                edgeDF <- .takeBalanceEdges(
                    graph = graph,
                    spatial = spatial,
                    senderIdx = senderIdx,
                    receiverIdx = receiverIdx,
                    eps = object@parameters$ballRadiusCoord/2,
                    minSize = edgeLocalMinSize,
                    sampleRate = edgeSampleRate
                )
            }
        )

        if (isTRUE(showTitle)) {
            title <- intr
        } else {
            title <- NULL
        }
        if (isTRUE(showSubTitle)) {
            intrTypeStr <- csdb[csdb$interactors == intr, 'type'] %>%
                pull(.data[['type']]) %>% as.character() %>% toupper() %>% lowwords() %>%
                paste0('-dependent')
            nSignifStr <- sprintf(', significant in %d spots', sum(receiverIdx))
            altNameStr <- csdb[csdb$interactors == intr, 'alt_name']
            altNameStr <- ifelse(nchar(altNameStr) > 0, sprintf('\nAlt. name: %s', altNameStr), '')
            subtitle <- paste0(intrTypeStr, nSignifStr, altNameStr)
        } else {
            subtitle <- NULL
        }

        # Clear up device
        grDevices::pdf(nullfile())
        # Bottom panel, solid color for highlighted sender spots,
        # translucent for others
        senderIdxBool <- colnames(object@rawData) %in% senderIdx
        plot3D::points3D(
            x = spatial[senderIdxBool, 1], y = spatial[senderIdxBool, 2], z = rep(0, sum(senderIdxBool)),
            xlim = xlim, ylim = ylim, zlim = c(-0.01,1.01),
            expand = boxHeight, d = perspTrans,
            col = colors[senderIdxBool], colvar = NULL,
            pch = 16, cex = dotSizeSignif, alpha = dotAlphaSignif,
            theta = theta, phi = phi,
            main = title, cex.main = titleTextSize/12,
            sub = subtitle, cex.sub = subtitleTextSize/12,
            xlab = "", ylab = "", zlab = 'Receiver <- Sender',
            cex.lab = labelTextSize/12,
            box = isTRUE(showBox),
            add = FALSE, plot = FALSE
        )
        plot3D::points3D(
            x = spatial[!senderIdxBool, 1], y = spatial[!senderIdxBool, 2], z = rep(0, sum(!senderIdxBool)),
            col = colors[!senderIdxBool], colvar = NULL,
            pch = 16, cex = dotSizeOther, alpha = dotAlphaOther,
            add = TRUE, plot = FALSE
        )
        # Segments showing sampled ones from true edges
        plot3D::segments3D(
            x0 = edgeDF$sender_x, y0 = edgeDF$sender_y, z0 = rep(0, nrow(edgeDF)),
            x1 = edgeDF$receiver_x, y1 = edgeDF$receiver_y, z1 = rep(1, nrow(edgeDF)),
            col = edgeColor, lwd = edgeLinewidth, lty = 1,
            add = TRUE, plot = FALSE
        )
        # Top panel, solid color for highlighted receiver spots,
        # translucent for others
        plot3D::points3D(
            x = spatial[receiverIdx,1], y = spatial[receiverIdx, 2], z = rep(1, sum(receiverIdx)),
            col = colors[receiverIdx], colvar = NULL,
            pch = 16, cex = dotSizeSignif, alpha = dotAlphaSignif,
            add = TRUE, plot = FALSE
        )
        plot3D::points3D(
            x = spatial[!receiverIdx,1], y = spatial[!receiverIdx, 2], z = rep(1, sum(!receiverIdx)),
            col = colors[!receiverIdx], colvar = NULL,
            pch = 16, cex = dotSizeOther, alpha = dotAlphaOther,
            add = TRUE, plot = FALSE
        )
        plist[[intr]] <- plot3D::getplist()
        plist[[intr]]$legend <- if (isTRUE(showLegend)) legendParams else NULL
        # line is eventually passed to graphics::title, controls how far main
        # and subtitles are from the plot
        plist[[intr]]$dot$line <- 0.5
        # plot3D works in a way that it creates the plotting region from a fixed
        # location (fraction) of the figure device region, using a 4-item vector
        # plist$plt$main <- c(x0, x1, y0, y1).
        # x0 won't work for anything after test.
        plist[[intr]]$plt$main[1] <- 0
        # x1 being 1 sets the right-most end of the plotting region to the
        # right edge of the screen. This is needed to fix the legend box to the
        # right side without hiding any long-text label out of the screen.
        plist[[intr]]$plt$main[2] <- 1
        # y0 and y1 sets the bottom and top edges so they are adapted to keep or
        # remove the area for subtitle and title, respectively.
        if (!isTRUE(showSubTitle)) {
            plist[[intr]]$plt$main[3] <- 0
        } else {
            plist[[intr]]$plt$main[3] <- 0.1
        }
        if (!isTRUE(showTitle)) {
            plist[[intr]]$plt$main[4] <- 1
        } else {
            plist[[intr]]$plt$main[4] <- 0.9
        }
        grDevices::dev.off()
    }
    if (length(plist) == 1) return(plist[[1]])
    else return(plist)
}

.takeBalanceEdges <- function(
        graph,
        spatial,
        senderIdx,
        receiverIdx,
        eps,
        minSize = 5,
        sampleRate = 0.1
) {
    spatialSub <- spatial[receiverIdx, , drop = FALSE]
    dbscanRes <- dbscan::dbscan(x = spatialSub, eps = eps, minPts = 1)
    receiverSpatCluster <- stats::setNames(factor(dbscanRes$cluster), rownames(spatialSub))
    edgeDFList <- list()
    for (i in seq_along(levels(receiverSpatCluster))) {
        receiverIdxSub <- names(receiverSpatCluster)[receiverSpatCluster == levels(receiverSpatCluster)[i]]
        graphSub <- graph[senderIdx, receiverIdxSub, drop = FALSE]
        edgeDFSub <- data.frame(
            sender = rownames(graphSub)[graphSub@i + 1],
            receiver = colnames(graphSub)[rep(seq_len(ncol(graphSub)), diff(graphSub@p))]
        ) %>%
            mutate(
                sender_x = spatial[.data[['sender']], 1],
                sender_y = spatial[.data[['sender']], 2],
                receiver_x = spatial[.data[['receiver']], 1],
                receiver_y = spatial[.data[['receiver']], 2]
            ) %>%
            select(.data[['sender_x']], .data[['sender_y']],
                   .data[['receiver_x']], .data[['receiver_y']])
        nEdgesSub <- nrow(edgeDFSub)
        if (nEdgesSub <= minSize) {
            edgeDFList[[i]] <- edgeDFSub
        } else if (nEdgesSub*sampleRate <= minSize) {
            edgeDFList[[i]] <- edgeDFSub[sample(nEdgesSub, minSize, replace = FALSE),]
        } else {
            edgeDFList[[i]] <- edgeDFSub[sample(nEdgesSub, floor(nEdgesSub * sampleRate), replace = FALSE),]
        }
    }
    do.call(rbind, edgeDFList)
}

#' Dot plot showing cluster-wise counts of edges between spots for each interaction
#' @description
#' This function collects all valid sender and receiver spot pairs for each
#' interaction, and summarizes the counts of such pairs (edges) by sender
#' cluster and receiver cluster. This function further calculates the total
#' expected possible edges between clusters based on the neighbor graph
#' corresponding to the type of each interaction. A proportion is then derived
#' by dividing the observed edge counts by the expected possible counts.
#' @inheritParams generalParam
#' @param senders Character vector of sender cluster names to only shown in the
#' plot. Default \code{NULL} shows all clusters that send valid signal.
#' @param receivers Character vector of receiver cluster names to only shown in
#' the plot. Default \code{NULL} shows all clusters that receive valid signal.
#' @param sumIntrs Logical, whether to sum the edges from all selected
#' interactions into one plot. Default \code{FALSE} generates separate plots for
#' each interaction.
#' @return A \code{ggplot} object when \code{sumIntrs = FALSE} and only one
#' interaction is selected; A list of \code{ggplot} objects when
#' \code{sumIntrs = FALSE} and multiple interactions are selected; A single
#' \code{ggplot} object when \code{sumIntrs = TRUE}.
#' @export
plotEdgeDot <- function(
        object,
        interaction,
        senders = NULL,
        receivers = NULL,
        sumIntrs = FALSE,
        fdrThresh = 0.05
) {
    csdb <- .checkIntrSelection(
        object,
        interaction = interaction,
        error = FALSE
    )
    if (nrow(csdb) == 0) {
        cli::cli_abort('No valid interaction selected')
    }

    clusterVarname <- object@parameters$cluster
    if (is.null(clusterVarname)) {
        cli::cli_abort(c(
            x = 'No cluster variable set. Setting all spots to the same color.',
            i = 'See {.code ?setCluster}.'
        ))
    }
    if (!clusterVarname %in% colnames(object@metadata)) {
        cli::cli_abort("Selected cluster variable {.val {clusterVarname}} is not available in metadata")
    }
    clusterVar <- object@metadata[[clusterVarname]]
    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)

    edges <- getAllEdges(object, csdb, fdrThresh) %>%
        mutate(
            sender_cluster = clusterVar[.data[['senders']]],
            receiver_cluster = clusterVar[.data[['receiver']]]
        )
    csdb <- csdb[csdb$interactors %in% edges$intr, , drop = FALSE]
    if (!is.null(senders)) {
        edges <- edges %>%
            filter(.data[['sender_cluster']] %in% senders)
    }
    if (!is.null(receivers)) {
        edges <- edges %>%
            filter(.data[['receiver_cluster']] %in% receivers)
    }

    clusterCellMtx <- fac2sparse(cs$cluster)
    nnCont <- object@neighborCont
    nnCont@x <- rep(1, length(nnCont@x))
    contClusterClusterExpect <- clusterCellMtx %*% nnCont %*% t(clusterCellMtx)
    nnDiff <- object@neighborDiff
    nnDiff@x <- rep(1, length(nnDiff@x))
    diffClusterClusterExpect <- clusterCellMtx %*% nnDiff %*% t(clusterCellMtx)
    plotlist <- list()
    if (!isTRUE(sumIntrs)) {
        for (i in seq_len(nrow(csdb))) {
            intrName <- csdb$interactors[i]
            intrType <- as.character(csdb$type[i])
            intrAltName <- csdb$alt_name[i]
            subtitle <- sprintf('%s-dependent interactions', lowwords(intrType))
            if (nchar(intrAltName) > 0) {
                subtitle <- paste0(subtitle, sprintf('\nAlt. name: %s', intrAltName))
            }
            if (intrType == 'contact') expect <- contClusterClusterExpect
            else if (intrType == 'diffusion') expect <- diffClusterClusterExpect

            edges_sub <- edges %>%
                filter(.data[['intr']] == intrName)
            tmp <- edges_sub %>%
                group_by(.data[['sender_cluster']], .data[['receiver_cluster']]) %>%
                summarise(Count = n())
            edges_sub <- edges_sub %>%
                group_by(.data[['sender_cluster']], .data[['receiver_cluster']]) %>%
                summarise(Count = n()) %>%
                mutate(
                    Prop = .data[['Count']] / diag(expect[.data[['sender_cluster']], .data[['receiver_cluster']], drop = FALSE])
                )
            plotlist[[intrName]] <- ggplot(
                data = edges_sub,
                mapping = aes(
                    x = .data[['receiver_cluster']],
                    y = .data[['sender_cluster']],
                    size = .data[['Count']],
                    fill = .data[['Prop']]
                )
            ) +
                geom_point(shape = 21, stroke = 0.3) +
                scale_size(name = '# Edges', transform = 'log10') +
                labs(
                    title = intrName,
                    subtitle = subtitle,
                    x = 'Receiver',
                    y = 'Sender',
                ) +
                scale_fill_viridis_c(direction = -1, labels = scales::percent) +
                theme_bw() +
                theme(
                    panel.grid = element_line(color = 'grey90', size = 0.1),
                    axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5),
                    plot.title = element_text(face = 'bold')
                ) +
                coord_fixed()
        }
        if (length(plotlist) == 1) return(plotlist[[1]])
        else return(plotlist)
    } else {

    }
}

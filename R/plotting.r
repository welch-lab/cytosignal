#' Plot significant interactions ranked by the user-specified metric
#'
#' @param object A cytosignal object
#' @param num.plot Number of interactions to plot
#' @param res_dir Directory to save the plots
#' @param plot.details Whether to plot NULL imputed values and scores
#' @param slot.use The LRscore slot to use for plotting
#' @param signif.use The metric used to rank the interactions, by default "result.hq.pear"
#' @param plot.clusters Whether to plot the clusters
#' @param plot.velo Whether to plot the velocity
#' @param colors.list A list of colors to use for plotting
#' @param pt.size Size of the points
#' @param pt.stroke Stroke of the points
#' @param u_width Width of the plot
#' @param u_hgt Height of the plot
#' @param set.res Resolution of the plot
#' @param return.plot Whether to return the plot
#'
#' @return A plot if return.plot is TRUE. Otherwise, plots are saved to the specified directory.
#'
#' @export
plotSignif <- function(object, num.plot = NULL, res_dir, plot.details = T, slot.use = NULL, signif.use = NULL, plot.clusters = T,
                    plot.velo = F, colors.list = NULL, pt.size=0.1, pt.stroke = 0.2, u_width = 6, u_hgt = 5, set.res = 200,
                    return.plot = F
){
    if (length(object@lrscore) == 0){
        stop("No score matrix found!\n")
    }

    if (is.null(slot.use)){
        slot.use = object@lrscore[["default"]]
    }

    score.obj = object@lrscore[[slot.use]]
    intr.db <- object@intr.valid[[score.obj@intr.slot]]

    if (is.null(signif.use)){
        signif.use = "result.hq.pear"
    }

    if (!signif.use %in% names(score.obj@res.list)){
        stop("No such significance level: ", signif.use, " found!\n")
    }

    if (is.null(num.plot)) {
        index.len <- length(object@lrscore[[slot.use]]@res.list[[signif.use]])
        num.plot <- intersect(union(1:20, (index.len-10):index.len), 1:index.len)
    }

    dge.raw <- object@counts
    clusters <- object@clusters
    lig.slot <- score.obj@lig.slot
    recep.slot <- score.obj@recep.slot
    dge.lig <- object@imputation[[lig.slot]]@imp.data
    dge.recep <- object@imputation[[recep.slot]]@imp.data

    # sample null values for plotting
    sample.index <- sample(ncol(score.obj@lig.null), ncol(dge.lig))
    null.dge.gau <- score.obj@lig.null[, sample.index]
    null.dge.dt <- score.obj@recep.null[, sample.index]
    colnames(null.dge.gau) <- colnames(null.dge.dt) <- colnames(dge.lig)
    rownames(null.dge.gau) <- rownames(null.dge.dt) <- rownames(dge.lig)

    # null.dge.dt <- dge.recep

    # sample null scores for plotting
    score.mtx <- score.obj@score
    null.lrscore.mtx <- score.obj@score.null
    null.lrscore.mtx <- null.lrscore.mtx[sample(nrow(null.lrscore.mtx), nrow(score.mtx)), ]
    rownames(null.lrscore.mtx) <- rownames(score.mtx)

    res.list <- score.obj@res.list[[signif.use]]
    res.intr.list = names(res.list)
    cells.loc = as.data.frame(object@cells.loc)
    cells.loc = cells.loc[colnames(dge.lig), ]

    if (is.null(colors.list)){
        levels.use <- levels(clusters)
        colors.list <- uniqueColors(length(levels.use))
        # as.character(paletteer::paletteer_d("ggsci::default_igv",
                        #n = length(levels.use)))
        names(colors.list) = levels.use
    }

    if (nrow(cells.loc) != ncol(dge.lig)){
        stop("Incorrect beads positions!\n")
    }

    cat("Now plotting INTRs...\nRank No.")

    plots.list <- lapply(num.plot, function(i){ # for each intr in the res.list
        cat(i, ", ", sep = "")
        sub.df = cells.loc
        null.sub.df = cells.loc

        pair.index = which(object@intr.valid$intr.index$id_cp_interaction == res.intr.list[i])

        ligand.name = object@intr.valid$intr.index[pair.index, 4]
        if (ligand.name == ""){ # if complex
            ligand.name = object@intr.valid$intr.index[pair.index, 2]
        }
        ligand.name = gsub("_HUMAN", "", ligand.name)

        receptor.name = object@intr.valid$intr.index[pair.index, 5]
        if (receptor.name == ""){
            receptor.name = object@intr.valid$intr.index[pair.index, 3]
        }
        receptor.name = gsub("_HUMAN", "", receptor.name)

        # generate interaction names
        # intr.name = paste0(gsub("_HUMAN", "", ligand.name), "-", gsub("_HUMAN", "", receptor.name))
        intr.name = paste0(ligand.name, "-", receptor.name)

        colnames(sub.df) = c("x", "y")
        colnames(null.sub.df) = c("x", "y")

        ligands = dge.lig[names( intr.db[[2]][intr.db[[2]] == res.intr.list[i]] ), ]
        ligands.ori = dge.raw[names(intr.db[[2]][intr.db[[2]] == res.intr.list[i]]), ]
        ligands.null = null.dge.gau[names(intr.db[[2]][intr.db[[2]] == res.intr.list[i]]), ]

        receps = dge.recep[names(intr.db[[3]][intr.db[[3]] == res.intr.list[i]]), ]
        receps.ori = dge.raw[names(intr.db[[3]][intr.db[[3]] == res.intr.list[i]]), ]
        receps.null = null.dge.dt[names(intr.db[[3]][intr.db[[3]] == res.intr.list[i]]), ]

        if (!is.null(nrow(ligands))){ligands = Matrix::colSums(ligands)}
        if (!is.null(nrow(receps))){receps = Matrix::colSums(receps)}
        if (!is.null(nrow(ligands.ori))){ligands.ori = Matrix::colSums(ligands.ori)}
        if (!is.null(nrow(ligands.null))){ligands.null = Matrix::colSums(ligands.null)}
        if (!is.null(nrow(receps.ori))){receps.ori = Matrix::colSums(receps.ori)}
        if (!is.null(nrow(receps.null))){receps.null = Matrix::colSums(receps.null)}

        sub.df$ligands = ligands
        sub.df$ligands_ori = ligands.ori
        sub.df$receptors = receps
        sub.df$receptors_ori = receps.ori
        null.sub.df$ligands_null = ligands.null
        null.sub.df$receptors_null = receps.null

        p.lig = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = ligands))+
            ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
            ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            cowplot::theme_cowplot(12)+
            labs(color = paste0("Ligand \n", ligand.name), x = NULL, y = NULL)

        p.lig.ori = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = ligands_ori))+
            ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
            ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            cowplot::theme_cowplot(12)+
            labs(color = paste0("Ligand-ori \n", ligand.name), x = NULL, y = NULL)

        p.lig.null = ggplot2::ggplot(null.sub.df, ggplot2::aes(x, y, color = ligands_null))+
            ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
            ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            cowplot::theme_cowplot(12)+
            labs(color = paste0("Ligand-null \n", ligand.name), x = NULL, y = NULL)

        p.recep = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = receptors))+
            ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
            ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            cowplot::theme_cowplot(12)+
            labs(color = paste0("Receptor \n", receptor.name), x = NULL, y = NULL)

        p.recep.ori = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = receptors_ori))+
            ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
            ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            cowplot::theme_cowplot(12)+
            labs(color = paste0("Receptor-ori \n", receptor.name), x = NULL, y = NULL)

        p.recep.null = ggplot2::ggplot(null.sub.df, ggplot2::aes(x, y, color = receptors_null))+
            ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
            ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            cowplot::theme_cowplot(12)+
            labs(color = paste0("Receptor-null \n", receptor.name), x = NULL, y = NULL)

        p.list = list(p.lig, p.lig.ori, p.lig.null, p.recep, p.recep.ori, p.recep.null)
        names(p.list) = c("lig", "lig.ori", "lig.null", "recep", "recep.ori", "recep.null")

        ##### plot scores

        if (ncol(dge.lig) != nrow(score.mtx)){
            cat("Setting LR-scores of cells with no neighbors to empty.\n")
            sub.df$scores = 0
            sub.df[Matrix::rownames(score.mtx), "scores"] = score.mtx[, res.intr.list[i]]

            null.sub.df$null_scores = 0
            null.sub.df[Matrix::rownames(null.lrscore.mtx), "null_scores"] = null.lrscore.mtx[, res.intr.list[i]]
        } else{
            sub.df$scores = score.mtx[, res.intr.list[i]]
            null.sub.df$null_scores = null.lrscore.mtx[, res.intr.list[i]]

        }

        p.score = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = scores))+
            ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
            ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            cowplot::theme_cowplot(12)+
            labs(x = NULL, y = NULL)

        p.list[["score"]] = p.score

        p.score.null = ggplot2::ggplot(null.sub.df, ggplot2::aes(x, y, color = null_scores))+
            ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
            ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            cowplot::theme_cowplot(12)+
            labs(x = NULL, y = NULL)

        p.list[["score.null"]] = p.score.null

        ##### plot significant beads

        sub.df$group = "beads"
        # sub.df$group[res.list[[i]]] = "significant"
        sub.df[res.list[[i]], "group"] = "significant"

        p.sig = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = group))+
            ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
            # geom_segment(ggplot2::aes(xend = xend, yend = yend), data = edges.df, color = "red", size = 0.2)+
            scale_colour_manual(values = c("beads" = "grey", "significant" = "red"))+
            cowplot::theme_cowplot(12)+
            labs(x = NULL, y = NULL)

        p.list[["sig"]] = p.sig

        if (plot.clusters){
            clusters = clusters[rownames(sub.df)]
            sub.df$clusters = as.character(clusters)

            sig.sub.df = sub.df[sub.df$group == "significant", ]

            p.clust = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = clusters))+
                ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
                # scale_colour_viridis_d(direction = -1)+
                # ggsci::scale_color_nejm()+
                # ggsci::scale_color_igv(palette = "default")+
                scale_colour_discrete(type = colors.list, breaks = levels.use) +
                ggplot2::geom_point(data = sig.sub.df, ggplot2::aes(x, y), color = "red", size = pt.size*3, alpha = 0.15)+
                cowplot::theme_cowplot(12)+
                guides(color = guide_legend(override.aes = list(size = 4)))+
                labs(x = NULL, y = NULL)

            p.list[["clust"]] = p.clust
        }

        if (plot.velo){
            sub.df$velo = velo.mtx[, res.intr.list[i]]

            p.velo = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = velo))+
                ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
                ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
                cowplot::theme_cowplot(12)+
                labs(x = NULL, y = NULL)

            p.list[["velo"]] = p.velo
        }

        if (return.plot){ return(p.list) }

        plot.len = length(p.list)
        # cat(paste0("Length of plot list: ", plot.len, "\n"))
        intr.name = paste0(intr.name, "-n_", length(res.list[[i]]))

        if (plot.details){
            p.list <- p.list[c(1, 4, 7, 9, 3, 2, 5, 8, 10, 6)]
            png(paste0(res_dir, "/Rank_", i, "_", intr.name, ".png"),
                        res = set.res, width = u_width*4, height = u_hgt*2, unit="in")
            print(cowplot::plot_grid(plotlist = p.list, nrow = 2))
            dev.off()

        } else{
            p.list <- p.list[c(1, 4, 7, 9, 2, 5, 8, 10)]
            png(paste0(res_dir, "/Rank_", i, "_", intr.name, ".png"),
                        res = set.res, width = u_width*4, height = u_hgt*2, unit="in")
            print(cowplot::plot_grid(plotlist = p.list, nrow = 2))
            dev.off()
        }

    })

    cat("End.\nFinished!\n")

    if (return.plot){return(plots.list)}
}

#' @export
plotCluster <- function(
        object,
        colors.list = NULL,
        pt.size = 0.5,
        pt.stroke = 0.2,
        legendNCol = 2,
        legendNRow = NULL
) {
    colors.list <- .checkColorList(object, colors.list)
    plotDF <- as.data.frame(object@cells.loc)
    colnames(plotDF) <- c("x", "y")
    clusters <- object@clusters
    clusters <- clusters[rownames(plotDF)]
    plotDF$clusters = clusters
    xRange <- c(min(plotDF$x), max(plotDF$x))
    yRange <- c(min(plotDF$y), max(plotDF$y))
    ggplot2::ggplot(plotDF,
                    ggplot2::aes(.data[["x"]], .data[["y"]],
                                 color = .data[["clusters"]])) +
        ggplot2::geom_point(size = pt.size, stroke = pt.stroke) +
        ggplot2::scale_colour_manual(values = levels(colors.list)) +
        ggplot2::guides(
            color = ggplot2::guide_legend(override.aes = list(size = 4),
                                          ncol = legendNCol, nrow = legendNRow)
        ) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.line.x = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(fill = NA)
        ) +
        ggplot2::coord_fixed(xlim = xRange, ylim = yRange)
}

#' @importFrom rlang .data
#' @export
plotIntrValue <- function(
        object,
        intr = NULL,
        type = c("ligand", "ligand_ori", "ligand_null",
                 "receptor", "receptor_ori", "receptor_null",
                 "score", "score_null"),
        slot.use = NULL,
        signif.use = NULL,
        color.flrs = FALSE,
        pt.size = 0.5,
        pt.stroke = 0.2,
        raster = NULL
) {
    values <- getIntrValue(object, intr = intr, type = type,
                           slot.use = slot.use, signif.use = signif.use)
    if (is.null(raster)) {
        # Automatically determine whether to rasterize
        if (utils::object.size(values) >= 1e7) {
            warning("Too much information queried, rastering the output plot ",
                    "automatically. To force vectorized plotting please set ",
                    "`raster = FALSE`.", immediate. = TRUE)
            raster <- TRUE
        }
        else raster <- FALSE
    }
    if (isTRUE(raster)) {
        if (!requireNamespace("scattermore", quietly = TRUE)) {
            stop("Package \"scattermore\" required for rastering the plots.\n",
                 "Please install with command:\n",
                 "install.packages(\"scattermore\")")
        }
    }

    outputList <- list()
    for (i in seq_along(values)) {
        plotDF <- values[[i]]
        intrName <- names(values)[i]
        ligRec <- strsplit(intrName, "-")[[1]]
        ligName <- ligRec[1]
        recName <- ligRec[2]
        xRange <- c(min(plotDF$x), max(plotDF$x))
        yRange <- c(min(plotDF$y), max(plotDF$y))
        plotList <- lapply(type, function(info) {
            legendTitle <- info
            if (startsWith(info, "ligand")) {
                title <- ligName
                #title <- paste0(info, "\n", ligName)
            } else if (startsWith(info, "receptor")) {
                title <- recName
                #title <- paste0(info, "\n", recName)
            } else {
                title <- intrName
                #title <- info
            }

            if (isTRUE(color.flrs)) {
                plotDF[[info]] <- sqrt(plotDF[[info]])
                p <- ggplot2::ggplot(plotDF, ggplot2::aes(x = .data[["x"]],
                                                          y = .data[["y"]]))
                if (isFALSE(raster)) {
                    p <- p +
                        ggplot2::geom_point(size = pt.size, stroke = pt.stroke,
                                            color = "#323C6E", alpha = .9) +
                        ggplot2::geom_point(mapping = ggplot2::aes(alpha = .data[[info]]),
                                            size = pt.size, stroke = pt.stroke,
                                            color = "#FA1E1E")
                } else (
                    p <- p + scattermore::geom_scattermore(pointsize = pt.size + 1,
                                                           stroke = pt.stroke,
                                                           color = "#323C6E") +
                        scattermore::geom_scattermore(mapping = ggplot2::aes(alpha = .data[[info]]),
                                                      pointsize = pt.size + 1.5, stroke = pt.stroke,
                                                      color = "#FA1E1E")
                )
                p <- p + ggplot2::scale_alpha_continuous(
                    range = c(0, 1), breaks = round(seq(from = 0, to = floor(max(plotDF[[info]])), length = 3), 2)
                )
            } else {
                p <- ggplot2::ggplot(plotDF, ggplot2::aes(x = .data[["x"]],
                                                          y = .data[["y"]],
                                                          colour = .data[[info]]))

                if (isFALSE(raster)) {
                    p <- p + ggplot2::geom_point(size = pt.size,
                                                 stroke = pt.stroke)
                } else (
                    p <- p + scattermore::geom_scattermore(pointsize = pt.size + 1.4,
                                                           stroke = pt.stroke)
                )

                p <- p +
                    ggplot2::scale_color_viridis_c(option = "plasma",
                                                   direction = -1,
                                                   na.value = '#F5F5F5')
            }

            p <- p +
                ggplot2::labs(color = legendTitle, title = title,
                              x = NULL, y = NULL) +
                ggplot2::theme_classic() +
                ggplot2::theme(
                    axis.line.x = ggplot2::element_blank(),
                    axis.line.y = ggplot2::element_blank(),
                    axis.ticks.x = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank(),
                    axis.text.x = ggplot2::element_blank(),
                    axis.text.y = ggplot2::element_blank(),
                    panel.border = ggplot2::element_rect(fill = NA),
                    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                    legend.title.align = 1,
                    legend.title = ggplot2::element_text(vjust = 0.8),
                    legend.position = "bottom",
                    legend.justification = 0.9,
                    legend.margin = ggplot2::margin(t = 0)
                ) +
                ggplot2::coord_fixed(xlim = xRange, ylim = yRange)
            if (isTRUE(color.flrs)) {
                p <- p +
                    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black"))
            }
            return(p)
        })
        names(plotList) <- type
        outputList[[intrName]] <- plotList
    }
    return(outputList)
}

#' Plot significant interactions ranked by the user-specified metric
#' @param object A cytosignal object
#' @param intr Specify interactions to be plotted. A vector of either the unique
#' ID of interactions or the numeric indices. Available IDs can be shown with
#' \code{showIntr(object)}. Availability of an interaction depends on the
#' LRscore slot to be used as well as the significance metric to be used.
#' @param edge,velo Logical, whether to plot edge or velocity, respectively.
#' @param slot.use The LRscore slot to use for plotting
#' @param signif.use The metric used to rank the interactions, by default "result.hq.pear"
#' @param colors.list A list of colors to use for plotting
#' @param plot_dir Directory to save the plots
#' @param pt.size Size of the points
#' @param pt.stroke Stroke of the points
#' @param return.plot Whether to return the plot
#' @param plot.fmt Format of output file. "png", "pdf", or "svg".
#' @param resolution Resolution of the output figure.
#' @param verbose Logical, whether to show progress.
#'
#' @return A plot if return.plot is TRUE. Otherwise, plots are saved to the specified directory.
#'
#' @export
plotSignif2 <- function(
        object,
        intr,
        edge = FALSE,
        velo = FALSE,
        slot.use = NULL,
        signif.use = NULL,
        colors.list = NULL,
        pt.size = 0.1,
        pt.stroke = 0.2,
        return.plot = FALSE,
        plot_dir = "csSignifPlot/",
        plot.fmt = c("png", "pdf", "svg"),
        raster = NULL,
        resolution = 300,
        verbose = FALSE
) {
    if (is.numeric(intr)) {
        intrRanks <- intr
        intr <- getCPIs(object, intrRanks, slot.use = slot.use,
                        signif.use = signif.use)
    } else if (is.character(intr)) {
        intrRanks <- getIntrRanks(object, intr, slot.use = slot.use,
                                  signif.use = signif.use)
    }
    names(intrRanks) <- intr

    # check before plotting to avoid wasting time
    if (isTRUE(velo)) {
        slot.use <- .checkSlotUse(object, slot.use, velo = TRUE)
    } else {
        slot.use <- .checkSlotUse(object, slot.use)
    }
    signif.use <- .checkSignifUse(object, signif.use, slot.use)
    intr <- .checkIntrAvail(object, intr, slot.use, signif.use)
    col.fac <- .checkColorList(object, colors.list)
    intrNames <- getIntrNames(object, intr)
    plot.fmt <- match.arg(plot.fmt)

    ggs <- plotIntrValue(object, intr = intr, slot.use = slot.use,
                         signif.use = signif.use, pt.size = pt.size,
                         pt.stroke = pt.stroke, raster = raster)
    edgePaths <- NULL
    veloPaths <- NULL
    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    if (isTRUE(edge)) {
        # check possible existence of edge plots
        edgePaths <- file.path(plot_dir, "Edge_plots")

        if (dir.exists(edgePaths)) {
            testFiles <- list.files(edgePaths, pattern = "Edge-.*.png")
            trueFiles <- paste0("Edge-", intrNames, ".png")
            if (all(trueFiles %in% testFiles)) {
                message("Edge plots already exist. Skip plotting edges.")
                do.plotEdge <- FALSE
            } else {
                message("Edge plots exist but not for all interactions. ",
                        "Replotting edges.")
                do.plotEdge <- TRUE
            }
        } else {
            do.plotEdge <- TRUE
        }
        if (isTRUE(do.plotEdge)) {
            plotEdge(object, intr, slot.use = slot.use, plot_dir = edgePaths,
            signif.use = signif.use, return.plot = FALSE,
            colors.list = levels(col.fac), width = 6, height = 6,
            pt.size = pt.size, pt.stroke = pt.stroke, verbose = FALSE)
        }

        edgePaths <- file.path(edgePaths, paste0("Edge-", intrNames, ".png"))
        names(edgePaths) <- intr
    }
    if (isTRUE(velo)) {
        veloPaths <- file.path(plot_dir, "Velo_plots")

        if (dir.exists(veloPaths)) {
            testFiles <- list.files(veloPaths, pattern = "Velo-.*.png")
            trueFiles <- paste0("Velo-", intrNames, ".png")
            if (all(trueFiles %in% testFiles)) {
                message("Velocity plots already exist. Skip plotting velocities.")
                do.plotVelo <- FALSE
            } else {
                message("Velocity plots exist but not for all interactions. ",
                        "Replotting velocities.")
                do.plotVelo <- TRUE
            }
        } else {
            do.plotVelo <- TRUE
        }

        if (isTRUE(do.plotVelo)) {
            plotVelo(object, intr, slot.use = slot.use, plot_dir = veloPaths,
                     signif.use = signif.use, return.plot = FALSE,
                     colors.list = levels(col.fac), width = 6, height = 6,
                     pt.size = pt.size, pt.stroke = pt.stroke, verbose = FALSE)
        }

        veloPaths <- file.path(veloPaths, paste0("Velo-", intrNames, ".png"))
        names(veloPaths) <- intr
    }
    plotList <- list()
    # if (isFALSE(edge) && isFALSE(velo)) legendNCol <- 1 else legendNCol <- 2
    pclust <- plotCluster(object, colors.list = levels(col.fac),
                          legendNCol = 1, pt.size = pt.size,
                          pt.stroke = pt.stroke)
    clustLegend <- cowplot::get_legend(pclust)
    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    message("Now plotting for interactions: ", appendLF = FALSE)
    for (intrx in intr) {
        intrName <- getIntrNames(object, intrx)
        message(intrName, ", ", appendLF = FALSE)

        gglist <- ggs[[intrName]]
        main <- cowplot::plot_grid(
            grid::textGrob(getLigandNames(object, intrx),
                           gp = grid::gpar(fontface = "bold")),
            gglist$ligand + ggplot2::ggtitle(NULL),
            gglist$ligand_ori + ggplot2::ggtitle(NULL),
            grid::textGrob(getReceptorNames(object, intrx),
                           gp = grid::gpar(fontface = "bold")),
            gglist$receptor + ggplot2::ggtitle(NULL),
            gglist$receptor_ori + ggplot2::ggtitle(NULL),
            grid::textGrob(intrName,
                           gp = grid::gpar(fontface = "bold")),
            gglist$score + ggplot2::ggtitle(NULL),
            pclust + ggplot2::theme(legend.position = "none"),
            nrow = 3, byrow = FALSE,
            align = "hv", axis = "lrtb",
            rel_heights = c(0.12, 1, 1)
        )
        main <- cowplot::plot_grid(main,
                                      clustLegend,
                                      rel_widths = c(4, 1))
        if (isFALSE(edge) && isFALSE(velo)) {
            combine <- main
        } #else {
            # add another sanity check
            # main <- cowplot::plot_grid(
            #     grid::textGrob(getLigandNames(object, intrx),
            #                    gp = grid::gpar(fontface = "bold")),
            #     gglist$ligand + ggplot2::ggtitle(NULL),
            #     gglist$ligand_ori + ggplot2::ggtitle(NULL),
            #     grid::textGrob(getReceptorNames(object, intrx),
            #                    gp = grid::gpar(fontface = "bold")),
            #     gglist$receptor + ggplot2::ggtitle(NULL),
            #     gglist$receptor_ori + ggplot2::ggtitle(NULL),
            #     grid::textGrob(intrName,
            #                    gp = grid::gpar(fontface = "bold")),
            #     gglist$score + ggplot2::ggtitle(NULL),
            #     clustLegend,
            #     nrow = 3, byrow = FALSE,
            #     align = "hv", axis = "lrtb",
            #     rel_heights = c(0.12, 1, 1)
            # )
        if (isTRUE(edge) && isTRUE(velo)) {
            edgeGrob <- png_as_grob(edgePaths[intrx])
            veloGrob <- png_as_grob(veloPaths[intrx])

            lower <- cowplot::plot_grid(
                edgeGrob, veloGrob, nrow = 1
            )
            combine <- cowplot::plot_grid(main, lower, nrow = 2,
                                          rel_heights = c(2, 1.2))
        } else if (isTRUE(edge) && isFALSE(velo)) {
            edgeGrob <- png_as_grob(edgePaths[intrx])
            combine <- cowplot::plot_grid(main, edgeGrob,
                                          nrow = 1,
                                          rel_widths = c(3, 2))
        } else if (isFALSE(edge) && isTRUE(velo)) {
            veloGrob <- png_as_grob(veloPaths[intrx])
            combine <- cowplot::plot_grid(main, veloGrob,
                                          nrow = 1,
                                          rel_widths = c(3, 2))
        }
        # }
        if (isTRUE(return.plot)) {
            plotList[[intrName]] <- combine
        } else {
            filenameIntr <- paste0(plot_dir, "/Rank_", intrRanks[intrx], "_", intrName)
            if (isFALSE(edge) && isFALSE(velo)) {
                width <- 12
                height <- 8
                filenameIntr <- paste0(filenameIntr, ".", plot.fmt)
            } else if (isTRUE(edge) && isTRUE(velo)) {
                width <- 12
                height <- 13
                filenameIntr <- paste0(filenameIntr, "_edge_velo.", plot.fmt)
            } else {
                width <- 22
                height <- 8
                filenameIntr <- paste0(filenameIntr, "_",
                                       ifelse(edge, "edge", "velo"),
                                       ".", plot.fmt)
            }
            if (plot.fmt == "png") {
                png(filenameIntr, width = width, height = height,
                    res = resolution, units = "in")
            } else if (plot.fmt == "pdf") {
                pdf(filenameIntr, width = width, height = height)
            } else if (plot.fmt == "svg") {
                svg(filenameIntr, width = width, height = height)
            }
            print(combine)
            dev.off()
            if (isTRUE(verbose)) {
                message("Combined figure saved at: ", normalizePath(filenameIntr))
            }
        }
    }
    message("Finished!")
    if (isFALSE(return.plot)) return(invisible())
    else return(plotList)
}

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
                    plot.velo = F, colors.list = NULL, pt.size=0.5, pt.stroke = 0.2, u_width = 6, u_hgt = 5, set.res = 200,
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
        colors.list = as.character(paletteer::paletteer_d("ggsci::default_igv",
                        n = length(levels.use)))
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

        p.lig = ggplot(sub.df, aes(x, y, color = ligands))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Ligand \n", ligand.name), x = NULL, y = NULL)

        p.lig.ori = ggplot(sub.df, aes(x, y, color = ligands_ori))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Ligand-ori \n", ligand.name), x = NULL, y = NULL)

        p.lig.null = ggplot(null.sub.df, aes(x, y, color = ligands_null))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Ligand-null \n", ligand.name), x = NULL, y = NULL)

        p.recep = ggplot(sub.df, aes(x, y, color = receptors))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Receptor \n", receptor.name), x = NULL, y = NULL)

        p.recep.ori = ggplot(sub.df, aes(x, y, color = receptors_ori))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Receptor-ori \n", receptor.name), x = NULL, y = NULL)

        p.recep.null = ggplot(null.sub.df, aes(x, y, color = receptors_null))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
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

        p.score = ggplot(sub.df, aes(x, y, color = scores))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(x = NULL, y = NULL)

        p.list[["score"]] = p.score

        p.score.null = ggplot(null.sub.df, aes(x, y, color = null_scores))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(x = NULL, y = NULL)

        p.list[["score.null"]] = p.score.null

        ##### plot significant beads

        sub.df$group = "beads"
        # sub.df$group[res.list[[i]]] = "significant"
        sub.df[res.list[[i]], "group"] = "significant"

        p.sig = ggplot(sub.df, aes(x, y, color = group))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            # geom_segment(aes(xend = xend, yend = yend), data = edges.df, color = "red", size = 0.2)+
            scale_colour_manual(values = c("beads" = "grey", "significant" = "red"))+
            theme_cowplot(12)+
            labs(x = NULL, y = NULL)

        p.list[["sig"]] = p.sig

        if (plot.clusters){
            clusters = clusters[rownames(sub.df)]
            sub.df$clusters = as.character(clusters)

            sig.sub.df = sub.df[sub.df$group == "significant", ]

            p.clust = ggplot(sub.df, aes(x, y, color = clusters))+
                geom_point(size = pt.size, stroke = pt.stroke)+
                # scale_colour_viridis_d(direction = -1)+
                # ggsci::scale_color_nejm()+
                # ggsci::scale_color_igv(palette = "default")+
                scale_colour_discrete(type = colors.list, breaks = levels.use) +
                geom_point(data = sig.sub.df, aes(x, y), color = "red", size = pt.size*3, alpha = 0.15)+
                theme_cowplot(12)+
                guides(color = guide_legend(override.aes = list(size = 4)))+
                labs(x = NULL, y = NULL)

            p.list[["clust"]] = p.clust
        }

        if (plot.velo){
            sub.df$velo = velo.mtx[, res.intr.list[i]]

            p.velo = ggplot(sub.df, aes(x, y, color = velo))+
                geom_point(size = pt.size, stroke = pt.stroke)+
                scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
                theme_cowplot(12)+
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



#' Plot 3D LR-velo ranked by the user-specified metric
#'
#' @param object A cytosignal object
#' @param num.plot Number of interactions to plot
#' @param res_dir Directory to save the plots
#' @param slot.use The LRscore slot to use for plotting
#' @param signif.use The metric used to rank the interactions, by default "result.hq.pear"
#' @param plot.clusters Whether to plot the clusters
#' @param plot.velo Whether to plot the velocity
#' @param colors.list A list of colors to use for plotting
#' @param pt.size Size of the points
#' @param pt.stroke Stroke of the points
#' @param u_width Width of the plot
#' @param u_hgt Height of the plot
#' @param z.scaler Scaling factor for the z-axis
#' @param use.cex Point size
#' @param use.shape Point shape
#' @param use_xbins Number of bins for the x-axis
#' @param use_ybins Number of bins for the y-axis
#' @param arrow.line.width Width of the arrow line
#' @param arrow.width Width of the arrow
#' @param use.phi Set view angel: phi the colatitude
#' @param use.theta Set view angel: theta gives the azimuthal direction
#' @param set.res Resolution of the plot
#' @param return.plot Whether to return the plot
#'
#' @return A plot if return.plot is TRUE. Otherwise, plots are saved to the specified directory.
#'
#' @export

plotVelo <- function(
    object,
    num.plot = NULL,
    plot_dir,
    plot.fmt,
    slot.use = NULL,
    signif.use = NULL,
    plot.clusters = TRUE,
    colors.list = NULL,
    z.scaler = 0.03,
    use.cex = 0.1,
    use.shape = 16,
    use_xbins = 15,
    use_ybins = 15,
    arrow.line.width = 0.6,
    arrow.width = 0.06,
    width = 3,
    height = 3,
    use.phi = 45,
    use.theta = -17
) {
    if (is.null(slot.use)){
        slot.use <- object@lrvelo[["default"]]
    }

    if (!slot.use %in% names(object@lrvelo)){
        stop("The slot.use is not in the lrvelo slot.\n")
    }

    if (!slot.use %in% names(object@lrscore)){
        stop("The slot.use is not in the lrscore slot.\n")
    }

    velo.obj <- object@lrvelo[[slot.use]]

    if (is.null(colors.list)){
        levels.use <- levels(object@clusters)
        colors.list = as.character(paletteer::paletteer_d("ggsci::default_igv",
                        n = length(levels.use)))
        names(colors.list) = levels.use
    }

    if (length(colors.list) != length(levels(object@clusters))){
        stop("The length of colors.list is not equal to the length of levels(clusters).\n")
    }

    col.fac = object@clusters
    levels(col.fac) = colors.list

    if (is.null(signif.use)){
        signif.use <- "result.hq.pear"
        # use.res.list <- object@lrscore[[signif.use]]@res.list[["result.hq.pear"]]
    }

    if (!signif.use %in% names(object@lrscore[[slot.use]]@res.list)) {
        stop("The rank is not in the corresponded object lrscore slot.\n")
    }

    use.res.list <- object@lrscore[[slot.use]]@res.list[[signif.use]]

    # intrs that are in the lrscore@result.list
    # velo may lack some genes so the intrs may be different

    velo.intr.index <- which(names(use.res.list) %in% colnames(velo.obj@intr.velo))

    if (!plot.fmt %in% c("png", "pdf")){
        stop("The plot.fmt is not supported.\n")
    }

    if (is.null(num.plot)) {
        index.len <- length(object@lrscore[[slot.use]]@res.list[[signif.use]])
        num.plot <- intersect(union(1:20, (index.len-10):index.len), 1:index.len)
    }

    cat("Plotting the results...\n")

    plot.list <- lapply(num.plot, function(i){
        intr.rank <- velo.intr.index[i]
        cat("No.", intr.rank, ", ", sep = "")
        use.intr = names(use.res.list)[intr.rank]

        # get value for the ranked #1 intr
        plot.df = as.data.frame(object@cells.loc)
        plot.df$velo = velo.obj@intr.velo[rownames(plot.df), use.intr]

        pt.df = plot.df[, c(1,2)]
        pt.df$z = 0

        pt.df$col = as.character(col.fac[rownames(pt.df)])
        x.scale = max(pt.df$x) - min(pt.df$x)
        y.scale = max(pt.df$y) - min(pt.df$y)

        #------------------------------------- plotting df for arrows -------------------------------------#
        # weight -> cell counts in each bin; var4 -> average velo in each bin
        bin = hex_bin(plot.df$x, plot.df$y, var4=plot.df$velo, var4.to.color = T, xbins = use_xbins, ybins = use_ybins)

        arrows.df = as.data.frame(bin[, c(1,2,3)])
        colnames(arrows.df) = c("x", "y", "bin_velo")

        ars.pos = arrows.df[arrows.df$bin_velo > 0, ]
        ars.neg = arrows.df[arrows.df$bin_velo < 0, ]
        ars.zero = arrows.df[arrows.df$bin_velo == 0, ]
        # use.scale = max(abs(arrows.df$bin_velo))
        z.scale = max(arrows.df$bin_velo) - min(arrows.df$bin_velo)

        # scale the length of each arrow by the maximum abs(bin_velo)
        # overall, scale by each intr

        z.hold = z.scale*z.scaler # set the interval between arrows and points

        if (nrow(ars.pos) > 0){
            ars.pos$length = abs(ars.pos$bin_velo)/z.scale
            ars.pos.plot.df = data.frame(
                x0 = ars.pos$x, y0 = ars.pos$y, z0 = z.hold,
                x1 = ars.pos$x, y1 = ars.pos$y, z1 = ars.pos$length+z.hold
            )
        }

        if (nrow(ars.neg) > 0){
            ars.neg$length = abs(ars.neg$bin_velo)/z.scale
            ars.neg.plot.df = data.frame(
                x0 = ars.neg$x, y0 = ars.neg$y, z0 = ars.neg$length+z.hold,
                x1 = ars.neg$x, y1 = ars.neg$y, z1 = z.hold
            )
        }

        if (nrow(ars.zero) > 0) {
            ars.zero.plot.df = data.frame(
                x0 = ars.zero$x, y0 = ars.zero$y, z0 = z.hold,
                col = "#f7f7f7"
            )
        }

        intr.name = getIntrNames(use.intr, object@intr.valid$intr.index)

        sig.len <- lengths(use.res.list)[intr.rank]

        if (plot.fmt == "png") {
            png(paste0(plot_dir, "/", "Rank_", intr.rank, "_", intr.name, "-n_", sig.len, "-3D.png"), width = width, height = height, units = "in", res = 400)
        } else if (plot.fmt == "pdf") {
            pdf(paste0(plot_dir, "/", "Rank_", intr.rank, "_", intr.name, "-n_", sig.len, "-3D.pdf"), width = width, height = height)
        } else {
            stop("Plotting format not supported.\n")
        }

        # cex: control size of points
        plot3D::points3D(pt.df$x, pt.df$y, pt.df$z,
            xlim = c(min(pt.df$x) - x.scale*0.025, max(pt.df$x) + x.scale*0.025),
            ylim = c(min(pt.df$y) - y.scale*0.025, max(pt.df$y) + y.scale*0.025),
            zlim = c(-0.01, 1.1), expand = 0.3,
            theta = use.theta, phi = use.phi, d = 2,
            colvar = NULL, col = pt.df$col,
            colkey = FALSE, pch = use.shape, cex = use.cex,
            main = "CytoSignalVelo", zlab = "velocity",
            xlab = "", ylab = ""
        )

        # plot points with velo = 0 to be grey points
        if (nrow(ars.zero) > 0) {
            ars.zero.plot.df = data.frame(
                x0 = ars.zero$x, y0 = ars.zero$y, z0 = z.hold
            )
            plot3D::points3D(ars.zero.plot.df$x0, ars.zero.plot.df$y0, ars.zero.plot.df$z0,
                colvar = NULL, col = "#f7f7f7",
                colkey = FALSE, pch = 19, cex = use.cex, add = T)
        }

        ### parameters for plotting arrows
        ### lwd: line width; length: arrow edge length; angle: arrow angle

        # positive arrows
        if (nrow(ars.pos) > 0 ){
            plot3D::arrows3D(
                x0 = ars.pos.plot.df$x0, y0 = ars.pos.plot.df$y0, z0 = ars.pos.plot.df$z0,
                x1 = ars.pos.plot.df$x1, y1 = ars.pos.plot.df$y1, z1 = ars.pos.plot.df$z1,
                colvar = NULL, col = "#d73027", lwd = arrow.line.width, length = arrow.width,
                clab = intr.name, d = 3, add = T
            )
        }

        # negative arrows
        if (nrow(ars.neg) > 0){
            plot3D::arrows3D(
                x0 = ars.neg.plot.df$x0, y0 = ars.neg.plot.df$y0, z0 = ars.neg.plot.df$z0,
                x1 = ars.neg.plot.df$x1, y1 = ars.neg.plot.df$y1, z1 = ars.neg.plot.df$z1,
                colvar = NULL, col = "#3c24f1", lwd = arrow.line.width, length = arrow.width,
                clab = intr.name, d = 3, add = T
            )
        }
        return(list(p = pt.df, a1 = ars.pos.plot.df, a2 = ars.neg.plot.df))
        dev.off()

    })
    return(plot.list)
    cat("End.\nFinished!\n")
}

#' @importFrom rlang .data
#' @export
plotIntrValue <- function(
        object,
        intr.use = NULL,
        type = c("ligand", "ligand_ori", "ligand_null",
                 "receptor", "receptor_ori", "receptor_null",
                 "score", "score_null"),
        slot.use = NULL,
        signif.use = NULL,
        pt.size = 0.5,
        pt.stroke = 0.2,
        raster = NULL
) {
    values <- getIntrValue(object, intr.use = intr.use, type = type,
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
        plotList <- lapply(type, function(info) {
            if (startsWith(info, "ligand")) {
                title <- paste0(info, "\n", ligName)
            } else if (startsWith(info, "receptor")) {
                title <- paste0(info, "\n", recName)
            } else {
                title <- info
            }
            p <- ggplot2::ggplot(plotDF, ggplot2::aes(x = .data[["x"]],
                                                      y = .data[["y"]],
                                                      colour = .data[[info]]))

            if (isFALSE(raster)) {
                p <- p + ggplot2::geom_point(size = pt.size, stroke = pt.stroke)
            } else (
                p <- p + scattermore::geom_scattermore(pointsize = pt.size,
                                                       stroke = pt.stroke)
            )

            p <- p +
                ggplot2::scale_color_viridis_c(option = "plasma", direction = -1,
                                               na.value = '#F5F5F5') +
                ggplot2::labs(color = title, x = NULL, y = NULL) +
                ggplot2::theme_classic() +
                ggplot2::theme(
                    axis.line.x = ggplot2::element_blank(),
                    axis.line.y = ggplot2::element_blank(),
                    axis.ticks.x = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank(),
                    axis.text.x = ggplot2::element_blank(),
                    axis.text.y = ggplot2::element_blank(),
                    panel.border = ggplot2::element_rect(fill = NA)
                )
        })
        names(plotList) <- type
        outputList[[intrName]] <- plotList
    }
    return(outputList)
}


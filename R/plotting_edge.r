#' Plotting edge for a given interaction from a CytoSignal object
#'
#' @param object \code{\linkS4class{CytoSignal}} object.
#' @param intr Interaction to plot. See available options with
#' \code{\link{showIntr}}.
#' @param type Type of plot, either \code{"sender"} or \code{"receiver"}.
#' Default \code{"sender"}.
#' @param slot.use Slot to use for plotting
#' @param signif.use Significance level to use for plotting
#' @param colors.list List of colors to use for plotting
#' @param return.plot Whether to return "plist" object for figure organization.
#' \code{TRUE} returns "plist" which can be shown on current display device with
#' \code{plot()}. \code{FALSE} saves the figure to disk. Default \code{TRUE}.
#' @param plot_dir Path where the figure will be saved when
#' \code{return.plot = FALSE}. Default \code{"csEdgePlot/"} (under current
#' working directory).
#' @param filename Filename of the figure inside \code{plot_dir}. Please match
#' extension name with \code{plot.fmt} and do not include path. Default
#' \code{NULL} and the exact filename will be determined by interaction
#' information and \code{plot.fmt}.
#' @param plot.fmt Format of plot, either \code{"png"} or \code{"pdf"}.
#' @param title Title of plot
#' @param use.cex Size of points
#' @param use.shape Shape of points
#' @param line.width Width of lines
#' @param use.phi Angle of 3D plot
#' @param use.theta Angle of 3D plot
#' @param z.scaler Scaling factor for z-axis
#' @param z.pt.interval Interval of points on z-axis
#' @param pt.stroke Width of points
#' @param u_width Width of plot
#' @param u_hgt Height of plot
#' @param set.res Resolution of plot
#' @return \code{plist} object when \code{return.plot = TRUE}, no in memory
#' object is returned when \code{return.plot = FALSE} but the figure will be
#' saved on disk.
#' @export
plotEdge <- function(
        object, intr, type = c("receiver", "sender"), slot.use = NULL,
        signif.use = NULL, colors.list = NULL,
        return.plot = TRUE, plot_dir = "csEdgePlot/", filename = NULL,
        plot.fmt = c("png", "pdf", "svg"),
        title = NULL, edge.size = 500, use.shape = 16,
        line.width = 0.01, use.phi = 30, use.theta = -17, z.scaler = 0.03,
        z.pt.interval = 1, pt.size = 0.1, pt.stroke = 0.2, width = 5, height = 5,
        set.res = 300, verbose = TRUE
){
    type <- match.arg(type)
    plot.fmt <- match.arg(plot.fmt)
    slot.use <- .checkSlotUse(object, slot.use)
    signif.use <- .checkSignifUse(object, signif.use, slot.use)
    intr <- .checkIntrAvail(object, intr, slot.use, signif.use)
    filename <- .checkArgLen(filename, length(intr))
    title <- .checkArgLen(title, length(intr))
    score.obj <- object@lrscore[[slot.use]]
    res.list <- score.obj@res.list[[signif.use]]

    lig.slot <- score.obj@lig.slot
    #recep.slot <- score.obj@recep.slot

    cells.loc <- as.data.frame(object@cells.loc)
    # nn.id <- object@imputation[[lig.slot]]@nn.id
    nn.graph <- object@imputation[[lig.slot]]@nn.graph

    col.fac <- .checkColorList(object, colors.list)
    intr.names <- getIntrNames(object, intr)
    plotList <- list()
    for (i in seq_along(intr)) {
        intrx <- intr[i]
        receiver.cells <- res.list[[intrx]]
        receiver.idx <- sort(match(receiver.cells, rownames(cells.loc)))
        # nn.id.sig <- nn.id[nn.id %in% as.character(receiver.idx), drop = T]
        nn.graph.sig <- nn.graph[, receiver.idx]
        if (is.null(title)) {
            titleIntr <- paste0("Edge", "-", intr.names[intrx])
        } else {
            titleIntr <- title[i]
        }
        if (isTRUE(verbose)) {
            message("Now plotting edges for interaction: ",
                    intr.names[intrx], " ...")
        }

        p1 <- .plotEdgeMatrix(
            cells.loc = cells.loc, type = type, edge.size = edge.size,
            nn.graph.sig = nn.graph.sig, receiver.idx = receiver.idx,
            col.fac = col.fac, res.list = res.list, intr = intr,
            title = titleIntr,
            pt.size = pt.size, use.shape = use.shape, line.width = line.width,
            use.phi = use.phi, use.theta = use.theta, z.scaler = z.scaler,
            z.pt.interval = z.pt.interval, pt.stroke = pt.stroke
        )

        if (isTRUE(return.plot)) {
            plotList[[intrx]] <- p1
        } else {
            if (!dir.exists(plot_dir)) dir.create(plot_dir)
            if (is.null(filename))
                filenameIntr <- paste0(plot_dir, "/", titleIntr, ".", plot.fmt)
            else
                filenameIntr <- file.path(plot_dir, filename[i])

            if (plot.fmt == "png") {
                png(filenameIntr, width = width, height = height, units = "in",
                    res = set.res)
            } else if (plot.fmt == "pdf") {
                pdf(filenameIntr, width = width, height = height)
            } else if (plot.fmt == "svg") {
                svg(filenameIntr, width = width, height = height)
            }

            graphics::par(mar = c(2, 2, 4, 2), oma = c(0, 0, 0, 0),
                          mgp = c(0, 0, 0), xpd = TRUE)
            plot(p1)
            dev.off()
            if (isTRUE(verbose)) {
                message("Interaction edge plot saved at: ",
                        normalizePath(filenameIntr))
            }

        }
    }
    if (isTRUE(return.plot)) {
        if (length(plotList) == 1) return(plotList[[1]])
        else return(plotList)
    }
    else invisible()
}

.plotEdgeMatrix <- function(
        cells.loc, nn.graph.sig, receiver.idx, col.fac, res.list, intr,
        type = c("sender", "receiver"), edge.size = 500,
        title = NULL, pt.size = 0.1, use.shape = 16, line.width = 0.01,
        use.phi = 30, use.theta = -17, z.scaler = 0.03, z.pt.interval = 1,
        pt.stroke = 0.2
){
    type <- match.arg(type)
    pt.df <- as.data.frame(cells.loc)
    pt.df$z <- 0
    colnames(pt.df) <- c("x", "y", "z")
    pt.df$col <- as.character(col.fac[rownames(pt.df)])
    # pt.df$group <- "beads"
    # pt.df[res.list[[intr]], "group"] <- "sig"

    x.scale <- max(pt.df$x) - min(pt.df$x)
    y.scale <- max(pt.df$y) - min(pt.df$y)

    xlim <- c(min(pt.df$x) - x.scale*0.025, max(pt.df$x) + x.scale*0.025)
    ylim <- c(min(pt.df$y) - y.scale*0.025, max(pt.df$y) + y.scale*0.025)

    # find sender
    senders <- unique(nn.graph.sig@i) + 1
    sender.idx <- nn.graph.sig@i + 1

    # pt.df.send <- pt.df[senders, ]
    # pt.df.n_send <- pt.df[-senders, ]

    # sender.x <- cells.loc[sender.idx, "x"]
    # sender.y <- cells.loc[sender.idx, "y"]

    # pt.df.rec <- pt.df[receiver.idx, ]
    # pt.df.n_rec <- pt.df[-receiver.idx, ]

    # recep.x <- rep(cells.loc[receiver.idx, "x"], times = diff(nn.graph.sig@p))
    # recep.y <- rep(cells.loc[receiver.idx, "y"], times = diff(nn.graph.sig@p))

    if (type == "sender") { # sender cells on the top panel

        up.df <- pt.df[senders, ]
        up.df$z <- up.df$z + z.pt.interval
        up.df.transp <- pt.df[-senders, ]
        up.df.transp$z <- up.df.transp$z + z.pt.interval

        down.df <- pt.df[receiver.idx, ]
        down.df.transp <- pt.df[-receiver.idx, ]

        seg.up.x <- cells.loc[sender.idx, "x"]
        seg.up.y <- cells.loc[sender.idx, "y"]
        seg.up.z <- rep(z.pt.interval, times = length(sender.idx)) - 0.01

        seg.down.x <- rep(cells.loc[receiver.idx, "x"], times = diff(nn.graph.sig@p))
        seg.down.y <- rep(cells.loc[receiver.idx, "y"], times = diff(nn.graph.sig@p))
        seg.down.z <- rep(0, times = length(sender.idx))

        # recep.z <- rep(0, times = length(sender.idx))
        # sender.z <- rep(z.pt.interval, times = length(sender.idx)) - 0.01
        zlab.use <- "Sender  ->  Receiver"

    } else { # receiver cells on the top panel

        up.df <- pt.df[receiver.idx, ]
        up.df$z <- up.df$z + z.pt.interval
        up.df.transp <- pt.df[-receiver.idx, ]
        up.df.transp$z <- up.df.transp$z + z.pt.interval

        down.df <- pt.df[senders, ]
        down.df.transp <- pt.df[-senders, ]

        seg.up.x <- rep(cells.loc[receiver.idx, "x"], times = diff(nn.graph.sig@p))
        seg.up.y <- rep(cells.loc[receiver.idx, "y"], times = diff(nn.graph.sig@p))
        seg.up.z <- rep(z.pt.interval, times = length(sender.idx)) - 0.01

        seg.down.x <- cells.loc[sender.idx, "x"]
        seg.down.y <- cells.loc[sender.idx, "y"]
        seg.down.z <- rep(0, times = length(sender.idx))

        # recep.z <- rep(z.pt.interval, times = length(sender.idx)) - 0.01
        # sender.z <- rep(0, times = length(sender.idx))
        zlab.use <- "Receiver  <-  Sender"
    }

    sample.idx <- sample(1:length(sender.idx),
                        size = min(edge.size, length(sender.idx)))

    grDevices::pdf(nullfile())

    # plot down panel, involved cells
    plot3D::points3D(
        down.df$x, down.df$y, down.df$z,
        xlim = xlim, ylim = ylim, zlim = c(-0.01, z.pt.interval + 0.1),
        expand = 0.7, theta = use.theta, phi = use.phi, d = 5,
        colvar = NULL, col = down.df$col, alpha = 1,
        colkey = FALSE, pch = use.shape, cex = pt.size,
        main = title, zlab = zlab.use,
        xlab = "", ylab = "", plot = FALSE
    )

    # plot down panel, non-involved cells
    plot3D::points3D(down.df.transp$x, down.df.transp$y, down.df.transp$z,
        colvar = NULL, col = down.df.transp$col, alpha = 0.2,
        colkey = FALSE, pch = use.shape, cex = pt.size,
        plot = FALSE, add = TRUE
    )

    plot3D::segments3D(
        x0 = seg.down.x[sample.idx], y0 = seg.down.y[sample.idx],
        z0 = seg.down.z[sample.idx],
        x1 = seg.up.x[sample.idx], y1 = seg.up.y[sample.idx],
        z1 = seg.up.z[sample.idx],
        col = "grey", lwd = 0.1, lty = 1,
        plot = FALSE, add = TRUE
    )

    # plot receiver cells
    plot3D::points3D(up.df$x, up.df$y, up.df$z,
                     colvar = NULL, col = up.df$col, alpha = 1,
                     colkey = FALSE, pch = use.shape, cex = pt.size,
                     plot = FALSE, add = TRUE)

    # plot non-receiver cells
    plot3D::points3D(x = up.df.transp$x, y = up.df.transp$y,
                     z = up.df.transp$z,
                     colvar = NULL, col = up.df.transp$col, alpha = 0.2,
                     colkey = FALSE, pch = use.shape, cex = pt.size,
                     plot = FALSE, add = TRUE)
    #pr <- grDevices::recordPlot()
    p1 <- plot3D::getplist()
    dev.off()
    return(p1)
}



#' Plot edges by each cluster as sender or receiver cells
#'
#' By default ranked by the number of significant cells in each cluster
#'
#' @param object A \code{CytoSignal} object
#' @param plot_dir Directory to save plots
#' @param cluster.list A list of clusters to plot
#' @param intr.num Number of interactions to plot
#' @param type Either "sender" or "receiver"
#' @param slot.use Slot to use for plotting
#' @param signif.use Significance threshold to use for plotting
#' @param colors.list A list of colors to use for each cluster
#' @param plot.fmt Plot format
#' @param edge.size Number of edges to plot
#' @param all.in.one Plot all clusters in one plot
#' @param plot.all.sig Plot all significant edges
#' @param pt.size Size of points
#' @param use.shape Shape of points
#' @param line.width Width of lines
#' @param use.phi Angle of view
#' @param use.theta Angle of view
#' @param z.scaler Scale of z-axis
#' @param z.pt.interval Interval of z-axis
#' @param pt.stroke Stroke of points
#' @param u_width Width of plot
#' @param u_hgt Height of plot
#' @param set.res Resolution of plot
#' @param return.plot Return plot object
#' @param ... Other arguments
#'
#' @return A list of plots
#'
#' @export
#'
plotSigCluster <- function(object, plot_dir, cluster.list = NULL, intr.num = 10, type = c("sender", "receiver"), slot.use = NULL, signif.use = NULL,
                    colors.list = NULL, plot.fmt = "png", edge.size = 2000, all.in.one = T, plot.all.sig = F,
                    use.cex = 0.1, use.shape = 16, line.width = 0.02,
                    use.phi = 30, use.theta = -17,
                    z.scaler = 0.03, z.pt.interval = 1,
                    pt.stroke = 0.2, u_width = 6, u_hgt = 6, set.res = 400,
                    return.plot = F
){
    if (length(type)!=1 || !(type %in% c("sender", "receiver"))){
        stop("type must be either 'sender' or 'receiver'!\n")
    }

    if (length(object@lrscore) == 0){
        stop("No score matrix found!\n")
    }

    if (is.null(slot.use)){
        slot.use = object@lrscore[["default"]]
    }

    score.obj = object@lrscore[[slot.use]]
    # intr.db <- object@intr.valid[[score.obj@intr.slot]]

    if (is.null(signif.use)){
        signif.use = "result.hq.pear"
    }

    if (!signif.use %in% names(score.obj@res.list)){
        stop("No such significance level: ", signif.use, " found!\n")
    }

    clusters <- object@clusters
    if (is.null(colors.list)){
        levels.use <- levels(clusters)
        colors.list <- uniqueColors(length(levels.use))
        #as.character(paletteer::paletteer_d("ggsci::default_igv",
        #                n = length(levels.use)))
        names(colors.list) = levels.use
    }
    col.fac = clusters
    levels(col.fac) = colors.list

    if (is.null(cluster.list)){
        cluster.list = names(sort(table(clusters), decreasing = T))[1:5]
    }

    if (length(intersect(cluster.list, levels(clusters))) == 0){
        stop("No such clusters in: ", cluster, " found!\n")
    }

    # initialize all needed variables
    lig.slot <- score.obj@lig.slot
    recep.slot <- score.obj@recep.slot
    res.list <- score.obj@res.list[[signif.use]]
    intr.list <- names(res.list)
    cells.loc = as.data.frame(object@cells.loc)
    # nn.id = object@imputation[[lig.slot]]@nn.id
    nn.graph = object@imputation[[lig.slot]]@nn.graph

    # make a df ranking by cell type, n(intr) X n(cell type)
    if (type == "sender") {
        intr.rank.clust <- lapply(intr.list, function(intr) {
            receiver.cells = res.list[[intr]]
            receiver.idx = sort(match(receiver.cells, rownames(cells.loc)))
            nn.graph.sig = nn.graph[, receiver.idx]
            senders = unique(nn.graph.sig@i)+1
            sender.cells = rownames(cells.loc)[senders]
            sender.clust = clusters[rownames(cells.loc)[senders]]
            sender.clust.vec = as.numeric(table(sender.clust)[levels(clusters)])
        })
    } else if (type == "receiver") {
        intr.rank.clust <- lapply(intr.list, function(intr) {
            receiver.cells = res.list[[intr]]
            receiver.idx = sort(match(receiver.cells, rownames(cells.loc)))
            rec.cells = rownames(cells.loc)[receiver.idx]
            rec.clust = clusters[rownames(cells.loc)[receiver.idx]]
            rec.clust.vec = as.numeric(table(rec.clust)[levels(clusters)])
        })
    } else {
        stop("type must be either 'sender' or 'receiver'!\n")
    }

    intr.rank.clust = Reduce(cbind, intr.rank.clust)
    dimnames(intr.rank.clust) = list(levels(clusters), intr.list)

    intr.symbols = lapply(intr.list, function(intr) {
        pair.index = which(object@intr.valid$intr.index$id_cp_interaction == intr)

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

        intr.name = paste0(ligand.name, "-", receptor.name)

        return(intr.name)
    })
    intr.symbols = unlist(intr.symbols)
    names(intr.symbols) = intr.list

    p.list = lapply(cluster.list, function(clust) {

        clust.idx = match(names(clusters[clusters == clust]), rownames(cells.loc))
        intr.order = sort(intr.rank.clust[clust, ], decreasing = T)

        p.intr = lapply(1:intr.num, function(n) {
            # get intr symbols
            sig.num = intr.order[n]
            intr = names(intr.order)[n]
            intr.name = intr.symbols[intr]
            title = paste0(clust, " #", n, " - ", intr.name, " (", sig.num, " cells)")

            if (type == "sender") {
                # get all intr sig cells and their index
                receiver.cells = res.list[[intr]]
                receiver.idx = sort(match(receiver.cells, rownames(cells.loc)))
                nn.graph.sig = nn.graph[, receiver.idx]

                # filter edges by sender cells only in this cluster
                sender.idx = nn.graph.sig@i + 1
                sender.idx.clust = which(sender.idx %in% clust.idx)
                nn.graph.sig = as(nn.graph.sig, "TsparseMatrix")
                nn.graph.sig@i = nn.graph.sig@i[sender.idx.clust]
                nn.graph.sig@j = nn.graph.sig@j[sender.idx.clust]
                nn.graph.sig@x = nn.graph.sig@x[sender.idx.clust]
                nn.graph.sig = as(nn.graph.sig, "CsparseMatrix")

                # remove empty columns
                empty.cols = which(diff(nn.graph.sig@p) == 0)
                nn.graph.sig = nn.graph.sig[, -empty.cols]

                if (!plot.all.sig){
                    receiver.idx = receiver.idx[-empty.cols]
                }

                p1 = plotEdge.matrix_like(
                    cells.loc, "sender", edge.size, nn.graph.sig, receiver.idx, col.fac, res.list, intr,
                    title, use.cex, use.shape, line.width, use.phi, use.theta,
                    z.scaler, z.pt.interval, pt.stroke, u_width, u_hgt, set.res
                )

            } else if (type == "receiver") {
                # get all intr sig cells and their index
                receiver.cells = res.list[[intr]]
                receiver.idx = match(receiver.cells, rownames(cells.loc))

                # filter receiver cells only in this cluster
                receiver.idx = sort(intersect(receiver.idx, clust.idx))
                nn.graph.sig = nn.graph[, receiver.idx]

                p1 = plotEdge.matrix_like(
                    cells.loc, "receiver", edge.size, nn.graph.sig, receiver.idx, col.fac, res.list, intr,
                    title, use.cex, use.shape, line.width, use.phi, use.theta,
                    z.scaler, z.pt.interval, pt.stroke, u_width, u_hgt, set.res
                )

            } else {
                stop("Type not supported.\n")
            }

            return(p1)
        })

        # rename plots with intr symbols
        names(p.intr) = intr.symbols[names(intr.order)[1:intr.num]]
        return(p.intr)
    })

    # rename plots with cluster names
    names(p.list) = cluster.list

    if (return.plot) {
        return(p.list)
    }

    # make plot dir for each cluster
    if (!dir.exists(plot_dir)){
        dir.create(plot_dir)
    }

    if (all.in.one) {
        clust.plot.dir = paste0(plot_dir, "/Cluster_In_One-as_", type, "/")
        if (!dir.exists(clust.plot.dir)){
            dir.create(clust.plot.dir, showWarnings = F)
        }

        do.plot <- lapply(names(p.list), function(clust) {
            message("Now plotting cluster: ", clust, " ...")

            plot.list = p.list[[clust]]

            # plot by row, maximum 4 plots per row
            if (length(plot.list) <= 4) {
                nrow = 1
                ncol = length(plot.list)
            } else {
                nrow = ceiling(length(plot.list)/4)
                ncol = 4
            }

            use.title = paste0(clust, " - ", type)

            if (plot.fmt == "png") {
                png(paste0(clust.plot.dir, "/", use.title, ".png"), width = u_width*ncol, height = u_hgt*nrow, units = "in", res = set.res)
            } else if (plot.fmt == "pdf") {
                pdf(paste0(clust.plot.dir, "/", use.title, ".pdf"), width = u_width*ncol, height = u_hgt*nrow)
            } else {
                stop("Plotting format not supported.\n")
            }

            # par(mfrow = c(nrow, ncol), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 0), mgp = c(0, 0, 0), xpd = T)
            par(mfrow = c(nrow, ncol), mar = c(1, 1, 1 ,1), oma = c(0, 0, 0, 0), mgp = c(0, 0, 0), xpd = T)
            for (i in 1:length(plot.list)) {
                plot(plot.list[[i]]); gc()
            }

            dev.off()

        })
    } else {
        plot_dir = paste0(plot_dir, "/Cluster_Separate-as_", type, "/")
        if (!dir.exists(plot_dir)){
            dir.create(plot_dir, showWarnings = F)
        }

        do.plot <- lapply(names(p.list), function(clust) {
            message("Now plotting cluster: ", clust, " ...")

            clust.plot.dir = paste0(plot_dir, "/", clust, "/")
            if (!dir.exists(clust.plot.dir)){
                dir.create(clust.plot.dir, showWarnings = F)
            }

            plot.list = p.list[[clust]]

            for (i in 1:length(plot.list)) {
                intr = names(plot.list)[i]
                p1 = plot.list[[i]]
                rank = i

                if (plot.fmt == "png") {
                    png(paste0(clust.plot.dir, "/#", rank, " - ", intr, ".png"), width = u_width, height = u_hgt, units = "in", res = set.res)
                } else if (plot.fmt == "pdf") {
                    pdf(paste0(clust.plot.dir, "/#", rank, " - ", intr, ".pdf"), width = u_width, height = u_hgt)
                } else {
                    stop("Plotting format not supported.\n")
                }
                par(mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0), mgp = c(0, 0, 0), xpd = T)
                plot(p1)
                dev.off()
                gc()
            }
        })
    }



}


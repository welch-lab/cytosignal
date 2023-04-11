#' Plotting edge for a given interaction from a CytoSignal object
#' 
#' @param object CytoSignal object
#' @param plot_dir Directory to save plots
#' @param intr Interaction to plot
#' @param type Type of plot, either "sender" or "receiver"
#' @param slot.use Slot to use for plotting
#' @param signif.use Significance level to use for plotting
#' @param colors.list List of colors to use for plotting
#' @param plot.fmt Format of plot, either "png" or "pdf"
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
#' @param return.plot Whether to return plot
#' 
#' @return Plot of edge
#' 
#' @export
#' 
plotEdge <- function(
    object,
    ...
) {
    UseMethod(generic = "plotEdge", object = object)
}

#' Plotting edge for a given interaction from a CytoSignal object
#' 
#' @return Plot of edge
#' 
#' @export
#' 
plotEdge.CytoSignal <- function(object, plot_dir, intr, type = c("sender", "receiver"), edge.size = 1000,
                    slot.use = NULL, signif.use = NULL,
                    colors.list = NULL, plot.fmt = "png", title = NULL,
                    use.cex = 0.1, use.shape = 16, line.width = 0.01,
                    use.phi = 30, use.theta = -17,
                    z.scaler = 0.03, z.pt.interval = 1,
                    pt.stroke = 0.2, u_width = 6, u_hgt = 5, set.res = 400,
                    return.plot = T
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
    intr.db <- object@intr.valid[[score.obj@intr.slot]]

    if (is.null(signif.use)){
        signif.use = "result.hq.pear"
    }

    if (!signif.use %in% names(score.obj@res.list)){
        stop("No such significance level: ", signif.use, " found!\n")
    }

    if (!intr %in% names(score.obj@res.list[[signif.use]])){
        stop("No such interaction: ", intr, " found!\n")
    }

    clusters <- object@clusters
    lig.slot <- score.obj@lig.slot
    recep.slot <- score.obj@recep.slot
    res.list <- score.obj@res.list[[signif.use]]
    cells.loc = as.data.frame(object@cells.loc)
    # nn.id = object@imputation[[lig.slot]]@nn.id
    nn.graph = object@imputation[[lig.slot]]@nn.graph

    # res.intr.list = names(res.list)
    receiver.cells = res.list[[intr]]
    receiver.idx = sort(match(receiver.cells, rownames(cells.loc)))
    # nn.id.sig = nn.id[nn.id %in% as.character(receiver.idx), drop = T]
    nn.graph.sig = nn.graph[, receiver.idx]

    if (is.null(colors.list)){
        levels.use <- levels(clusters)
        colors.list = as.character(paletteer::paletteer_d("ggsci::default_igv",
                        n = length(levels.use)))
        names(colors.list) = levels.use
    }
    col.fac = clusters
    levels(col.fac) = colors.list

    # get intr symbols
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

    message("Now plotting INTR: ", intr.name, " ...")

    if (is.null(title)){
        title = paste0("Edge", "-", intr.name)
    } else {
        title = paste0(title, "Edge", "-", intr.name)
    }

    if (type == "sender") {
        p1 = plotEdge.matrix_like(
            cells.loc, type = type, edge.size = edge.size, nn.graph.sig, receiver.idx, col.fac, res.list, intr,
            title, use.cex, use.shape, line.width, use.phi, use.theta,
            z.scaler, z.pt.interval, pt.stroke, u_width, u_hgt, set.res
        )
    } else if (type == "receiver") {
        p1 = plotEdge.matrix_like(
            cells.loc, type = type, edge.size = edge.size, nn.graph.sig, receiver.idx, col.fac, res.list, intr,
            title, use.cex, use.shape, line.width, use.phi, use.theta,
            z.scaler, z.pt.interval, pt.stroke, u_width, u_hgt, set.res
        )
    } else {
        stop("Type not supported.\n")
    }

    if (return.plot){
        return(p1)
    }

    if (plot.fmt == "png") {
        png(paste0(plot_dir, "/", title, ".png"), width = u_width, height = u_hgt, units = "in", res = set.res)
    } else if (plot.fmt == "pdf") {
        pdf(paste0(plot_dir, "/", title, ".pdf"), width = u_width, height = u_hgt)
    } else {
        stop("Plotting format not supported.\n")
    }
    
    par(mfrow = c(1, 1), mar = c(1, 1, 1 ,1), oma = c(0, 0, 0, 0), mgp = c(0, 0, 0), xpd = T)
    plot(p1)

    dev.off()
    

}

#' Sub function for plotEdge
#' 
#' @return a plist object
#' 
#' @export
#' 
plotEdge.matrix_like <- function(
        cells.loc, type = NULL, edge.size = 2000, nn.graph.sig, receiver.idx, col.fac, res.list, intr,
        title, use.cex, use.shape, line.width, use.phi, use.theta,
        z.scaler, z.pt.interval, pt.stroke, u_width, u_hgt, set.res
){
    pt.df = as.data.frame(cells.loc)
    pt.df$z = 0
    colnames(pt.df) = c("x", "y", "z")
    pt.df$col = as.character(col.fac[rownames(pt.df)])
    # pt.df$group = "beads"
    # pt.df[res.list[[intr]], "group"] = "sig"

    x.scale = max(pt.df$x) - min(pt.df$x)
    y.scale = max(pt.df$y) - min(pt.df$y)

    # find sender
    # senders = as.numeric(unique(names(nn.id.sig)))
    senders = unique(nn.graph.sig@i)+1
    sender.idx = nn.graph.sig@i+1

    # cex: control size of points
    pt.df.send = pt.df[senders, ]
    pt.df.n_send = pt.df[-senders, ]

    pt.df.rec = pt.df[receiver.idx, ]
    pt.df.n_rec = pt.df[-receiver.idx, ]

    xlim = c(min(pt.df$x) - x.scale*0.025, max(pt.df$x) + x.scale*0.025)
    ylim = c(min(pt.df$y) - y.scale*0.025, max(pt.df$y) + y.scale*0.025)

    recep.x = rep(cells.loc[receiver.idx, "x"], times = diff(nn.graph.sig@p))
    recep.y = rep(cells.loc[receiver.idx, "y"], times = diff(nn.graph.sig@p))
    # recep.z = rep(z.pt.interval, times = length(sender.idx)) - 0.01
    # recep.z <- rep(0, times = length(sender.idx))

    sender.x = cells.loc[sender.idx, "x"]
    sender.y = cells.loc[sender.idx, "y"]
    # sender.z = rep(0, times = length(sender.idx))
    # sender.z = rep(z.pt.interval, times = length(sender.idx)) - 0.01

    if (type == "sender") {
        recep.z <- rep(0, times = length(sender.idx))
        sender.z <- rep(z.pt.interval, times = length(sender.idx)) - 0.01
    } else if (type == "receiver") {
        recep.z <- rep(z.pt.interval, times = length(sender.idx)) - 0.01
        sender.z <- rep(0, times = length(sender.idx))
    } else {
        stop("Type not supported.\n")
    }

    sample.idx = sample(1:length(sender.x), size = min(edge.size, length(sender.x)))

    # plot sender cells
    plot3D::points3D(pt.df.rec$x, pt.df.rec$y, pt.df.rec$z, 
        xlim = xlim, ylim = ylim, zlim = c(-0.01, z.pt.interval+0.1), expand = 0.7,
        theta = use.theta, phi = use.phi, d = 5, 
        colvar = NULL, col = pt.df.rec$col, alpha = 1,
        colkey = FALSE, pch = use.shape, cex = use.cex,
        main = title, zlab = "Inferred CCC",
        xlab = "", ylab = "", plot = F
    )

    # plot non-sender cells
    plot3D::points3D(pt.df.n_rec$x, pt.df.n_rec$y, pt.df.n_rec$z, 
        colvar = NULL, col = pt.df.n_rec$col, alpha = 0.2,
        colkey = FALSE, pch = use.shape, cex = use.cex, 
        plot = F, add = T
    )

    plot3D::segments3D(
        x0 = sender.x[sample.idx], y0 = sender.y[sample.idx], z0 = sender.z[sample.idx],
        x1 = recep.x[sample.idx], y1 = recep.y[sample.idx], z1 = recep.z[sample.idx],
        col = "grey", lwd = 0.1, lty = 1,
        plot = F, add = T
    )

    # plot receiver cells
    plot3D::points3D(pt.df.send$x, pt.df.send$y, pt.df.send$z + z.pt.interval, 
        colvar = NULL, col = pt.df.send$col, alpha = 1,
        colkey = FALSE, pch = use.shape, cex = use.cex, 
        plot = F, add = T
    )

    # plot non-receiver cells
    plot3D::points3D(pt.df.n_send$x, pt.df.n_send$y, pt.df.n_send$z + z.pt.interval, 
        colvar = NULL, col = pt.df.n_send$col, alpha = 0.2,
        colkey = FALSE, pch = use.shape, cex = use.cex, 
        plot = F, add = T
    )

    p1 = plot3D::getplist()

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
#' @param use.cex Size of points
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
    intr.db <- object@intr.valid[[score.obj@intr.slot]]

    if (is.null(signif.use)){
        signif.use = "result.hq.pear"
    }

    if (!signif.use %in% names(score.obj@res.list)){
        stop("No such significance level: ", signif.use, " found!\n")
    }

    clusters <- object@clusters
    if (is.null(colors.list)){
        levels.use <- levels(clusters)
        colors.list = as.character(paletteer::paletteer_d("ggsci::default_igv",
                        n = length(levels.use)))
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


#' Plot 3D LR-velo ranked by the user-specified metric
#'
#' @param object A cytosignal object
#' @param num.plot Number of interactions to plot
#' @param res_dir Directory to save the plots
#' @param slot.use The LRscore slot to use for plotting
#' @param signif.use The metric used to rank the interactions, by default "result.hq.pear"
#' @param use.clusters Plot only selected clusters. Default \code{NULL} for all clusters.
#' @param plot.velo Whether to plot the velocity
#' @param colors.list A list of colors to use for plotting
#' @param pt.size Size of the points
#' @param pt.stroke Stroke of the points
#' @param width Width of the plot. Default 6
#' @param height Height of the plot. Default 6
#' @param z.scaler Scaling factor for the z-axis
#' @param pt.size Point size. Default 0.1
#' @param use.shape Point shape. Default 16
#' @param use_xbins Number of bins for the x-axis. Default 15
#' @param use_ybins Number of bins for the y-axis. Default 15
#' @param arrow.line.width Width of the arrow line. Default 0.6
#' @param arrow.width Width of the arrow. Default 0.06
#' @param use.phi Set view angel: phi the colatitude. Default 30
#' @param use.theta Set view angel: theta gives the azimuthal direction. Default -17
#' @param box Whether to show a box panel for the 3D region. Default TRUE
#' @param axis.arrow.len When box = FALSE, set the length of the axis arrows.
#' @param set.res Resolution of the plot. Default 300.
#' @param return.plot Whether to return the plot
#'
#' @return A plot if return.plot is TRUE. Otherwise, plots are saved to the specified directory.
#'
#' @export
plotVelo <- function(
    object,
    intr,
    return.plot = TRUE,
    plot_dir = "csVeloPlot/",
    filename = NULL,
    plot.fmt = c("png", "pdf"),
    slot.use = NULL,
    signif.use = NULL,
    use.clusters = NULL,
    colors.list = NULL,
    z.scaler = 0.03,
    title = NULL,
    pt.size = 0.1,
    use.shape = 16,
    use_xbins = 15,
    use_ybins = 15,
    arrow.line.width = 0.6,
    arrow.width = 0.06,
    pt.stroke = 0.2,
    use.phi = 30,
    use.theta = -17,
    box = TRUE,
    axis.arrow.len = 1,
    width = 6,
    height = 6,
    set.res = 300,
    verbose = TRUE
) {
  plot.fmt <- match.arg(plot.fmt)
  slot.use <- .checkSlotUse(object, slot.use, velo = TRUE)
  signif.use <- .checkSignifUse(object, signif.use, slot.use)
  intr <- .checkIntrAvail(object, intr, slot.use, signif.use)
  filename <- .checkArgLen(filename, length(intr))
  title <- .checkArgLen(title, length(intr))
  cells.loc <- as.data.frame(object@cells.loc)

  velo.obj <- object@lrvelo[[slot.use]]

  col.fac <- .checkColorList(object, colors.list)

  intr.names <- getIntrNames(object, intr)

  plotList <- list()
  for (i in seq_along(intr)) {
    intrx <- intr[i]
    velo <- velo.obj@intr.velo[rownames(cells.loc), intrx]
    if (is.null(title)) {
      titleIntr <- paste0("Velo", "-", intr.names[intrx])
    } else {
      titleIntr <- title[i]
    }
    if (isTRUE(verbose)) {
      message("Now plotting velocity for interaction: ",
              intr.names[intrx], " ...")
    }
    if (!is.null(use.clusters)) {
      if (any(!use.clusters %in% object@clusters)) {
        stop("Specified clusters not found: ",
             paste(use.clusters[!use.clusters %in% object@clusters],
                   collapse = ", "))
      }
      idx <- object@clusters %in% use.clusters
      cells.loc <- cells.loc[idx,]
      col.fac <- droplevels(col.fac[idx])
      velo <- velo[idx]
    }
    p1 <- .plotVeloMatrix(
      cells.loc = cells.loc, intr = intrx, velo = velo, col.fac = col.fac,
      use_xbins = use_xbins, use_ybins = use_ybins, title = titleIntr,
      pt.size = pt.size, use.shape = use.shape,
      arrow.line.width = arrow.line.width, arrow.width = arrow.width,
      use.phi = use.phi, use.theta = use.theta, z.scaler = z.scaler,
      pt.stroke = pt.stroke, box = box, axis.arrow.len = axis.arrow.len
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
        message("Interaction velocity plot saved at: ",
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

.plotVeloMatrix <- function(
    cells.loc, intr, velo, col.fac, use_xbins = 15, use_ybins = 15,
    title = "CytoSignalVelo", pt.size = 0.1, use.shape = 16,
    arrow.line.width = 0.6, arrow.width = 0.06,
    use.phi = 30, use.theta = -17, z.scaler = 0.03,
    pt.stroke = 0.2, box = TRUE, axis.arrow.len = 1
) {
  # get value for the ranked #1 intr
  plot.df = as.data.frame(cells.loc)
  plot.df$velo = velo

  pt.df = plot.df[, c(1,2)]
  pt.df$z = 0
  pt.df$col = as.character(col.fac[rownames(pt.df)])
  x.scale = max(pt.df$x) - min(pt.df$x)
  y.scale = max(pt.df$y) - min(pt.df$y)

  # plotting df for arrows %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # weight -> cell counts in each bin; var4 -> average velo in each bin
  bin = hex_bin(plot.df$x, plot.df$y, var4 = plot.df$velo,
                var4.to.color = TRUE, xbins = use_xbins, ybins = use_ybins)

  arrows.df = as.data.frame(bin[, c(1,2,3)])
  colnames(arrows.df) = c("x", "y", "bin_velo")

  ars.pos = arrows.df[arrows.df$bin_velo > 0, ]
  ars.neg = arrows.df[arrows.df$bin_velo < 0, ]
  ars.zero = arrows.df[arrows.df$bin_velo == 0, ]
  # use.scale = max(abs(arrows.df$bin_velo))
  z.scale = max(arrows.df$bin_velo) - min(arrows.df$bin_velo)

  # scale the length of each arrow by the maximum abs(bin_velo)
  # overall, scale by each intr

  # z.hold = z.scale*z.scaler # set the interval between arrows and points
  z.hold = 0.04

  if (nrow(ars.pos) > 0) {
    ars.pos$length = abs(ars.pos$bin_velo)/z.scale
    ars.pos.plot.df = data.frame(
      x0 = ars.pos$x, y0 = ars.pos$y, z0 = z.hold,
      x1 = ars.pos$x, y1 = ars.pos$y, z1 = ars.pos$length + z.hold
    )
  }

  if (nrow(ars.neg) > 0) {
    ars.neg$length = abs(ars.neg$bin_velo)/z.scale
    ars.neg.plot.df = data.frame(
      x0 = ars.neg$x, y0 = ars.neg$y, z0 = ars.neg$length + z.hold,
      x1 = ars.neg$x, y1 = ars.neg$y, z1 = z.hold
    )
  }

  if (nrow(ars.zero) > 0) {
    ars.zero.plot.df = data.frame(
      x0 = ars.zero$x, y0 = ars.zero$y, z0 = z.hold,
      col = "#f7f7f7"
    )
  }
  pdf(nullfile())
  # cex: control size of points
  xlim <- c(min(pt.df$x) - x.scale*0.025,
            max(pt.df$x) + x.scale*0.025)
  ylim <- c(min(pt.df$y) - y.scale*0.025,
            max(pt.df$y) + y.scale*0.025)
  plot3D::points3D(pt.df$x, pt.df$y, pt.df$z,
                   xlim = xlim, ylim = ylim,
                   zlim = c(-0.01, 1.1), expand = 0.3, pt.stroke = pt.stroke,
                   theta = use.theta, phi = use.phi, d = 2,
                   colvar = NULL, col = pt.df$col,
                   colkey = FALSE, pch = use.shape, cex = pt.size,
                   main = title, zlab = "velocity",
                   xlab = "", ylab = "", plot = FALSE, box = box
  )

  # plot points with velo = 0 to be grey points
  if (nrow(ars.zero) > 0) {
    ars.zero.plot.df = data.frame(
      x0 = ars.zero$x, y0 = ars.zero$y, z0 = z.hold
    )
    plot3D::points3D(ars.zero.plot.df$x0,
                     ars.zero.plot.df$y0,
                     ars.zero.plot.df$z0,
                     colvar = NULL, col = "#f7f7f7",
                     colkey = FALSE, pch = 19, cex = pt.size,
                     add = TRUE, plot = FALSE)
  }

  ### parameters for plotting arrows
  ### lwd: line width; length: arrow edge length; angle: arrow angle

  # positive arrows
  if (nrow(ars.pos) > 0 ) {
    plot3D::arrows3D(
      x0 = ars.pos.plot.df$x0, y0 = ars.pos.plot.df$y0, z0 = ars.pos.plot.df$z0,
      x1 = ars.pos.plot.df$x1, y1 = ars.pos.plot.df$y1, z1 = ars.pos.plot.df$z1,
      colvar = NULL, col = "#d73027", lwd = arrow.line.width,
      length = arrow.width, clab = intr, d = 3, add = TRUE, plot = FALSE
    )
  }

  # negative arrows
  if (nrow(ars.neg) > 0) {
    plot3D::arrows3D(
      x0 = ars.neg.plot.df$x0, y0 = ars.neg.plot.df$y0, z0 = ars.neg.plot.df$z0,
      x1 = ars.neg.plot.df$x1, y1 = ars.neg.plot.df$y1, z1 = ars.neg.plot.df$z1,
      colvar = NULL, col = "#3c24f1", lwd = arrow.line.width,
      length = arrow.width, clab = intr, d = 3, add = TRUE, plot = FALSE
    )
  }

  # if (isFALSE(box)) {
  #   # Manually add axis arrows when no default box (w/ those arrows)
  #   axisCenter <- c(xlim[1], ylim[2], -0.1)
  #   z.axis.len <- max(ars.pos.plot.df$z1) * axis.arrow.len
  #   x.axis.len <- .2*diff(xlim) * axis.arrow.len
  #   y.axis.len <- .2*diff(ylim) * axis.arrow.len
  #
  #   # The z-axis arrow, pointing upwards
  #   plot3D::arrows3D(x0 = axisCenter[1], y0 = axisCenter[2], z0 = axisCenter[3],
  #                    x1 = axisCenter[1], y1 = axisCenter[2], z1 = axisCenter[3] + z.axis.len,
  #                    colvar = NULL, col = "black",
  #                    xlim = xlim, ylim = ylim, plot = FALSE, add = TRUE)
  #   # The x-axis arrow, pointing to the right
  #   plot3D::arrows3D(x0 = axisCenter[1], y0 = axisCenter[2], z0 = axisCenter[3],
  #                    x1 = axisCenter[1] + x.axis.len, y1 = axisCenter[2], z1 = axisCenter[3],
  #                    colvar = NULL, col = "black",
  #                    xlim = xlim, ylim = ylim, plot = FALSE, add = TRUE)
  #   # The y-axis arrow, pointing to the right
  #   plot3D::arrows3D(x0 = axisCenter[1], y0 = axisCenter[2], z0 = axisCenter[3],
  #                    x1 = axisCenter[1], y1 = axisCenter[2] - y.axis.len, z1 = axisCenter[3],
  #                    colvar = NULL, col = "black",
  #                    xlim = xlim, ylim = ylim, plot = FALSE, add = TRUE)
  #   # plot3D::text3D(axisCenter[1], axisCenter[2], axisCenter[3], labels = "velocity",
  #   #                add = TRUE, plot = FALSE)
  # }

  p <- plot3D::getplist()
  dev.off()
  return(p)
}

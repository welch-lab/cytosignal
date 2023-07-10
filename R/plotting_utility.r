#' GGPLOT2 FUNCTIONALITY FOR MAPPING TO HEXAGON SIZE AND COLOUR AESTHETICS
#' by Robin Edwards, 2013 (geotheory.co.uk, @geotheory)
#' This has been adapted from the ggplot bin_hex.R script that underpins geom_hex, etc
#' (see https://github.com/hadley/densityvis/blob/master/R/bin-hex.r).
#'
#' These functions implement aesthetic mapping to hexagon size (area), in addition to the existing
#' colour-mapping functionality.  The key change is the addition of a new fourth variable (var4)
#' to hex_bin(), which complements the inbuilt hexagon binning functionality.  The 'frequency.to.area'
#' argument enables the default mappings of binned data to colour and var4 to size to be interchanged.
#' The hmin/hmax arguments [0,1] set area mapping constraints (hmax can exceed 1).
#' xlim/xlat enable hexagon tesselation to be constrained independently of data range.
#' There may be some bugs in the implementation. A legend for hexagon size has not been implemented.

# require(ggplot2)
# require(plyr)

#' Bin data into hexagons (2d).
#'
#' @param x a numeric vector of x positions
#' @param y a numeric vector of y positions
#' @param weight \code{NULL} or a numeric vector providing weights for each
#'   observation, replace `counts`, mapped to color
#' @param var4 \code{NULL} or a numeric vector providing weights for each
#'   observation, averaged within each bin, mapped to hex bin size
#' @param height height of each hexagon, if \code{NULL} computed from ybins
#' @param width width of each hexagon, if \code{NULL} computed from ybins
#' @param xbins number of horizontal bins, if \code{width} unspecified
#' @param ybins number of vertical bins, if \code{height} unspecified
#' @param na.rm If \code{TRUE} missing values will be silently removed,
#'   otherwise they will be removed with a warning.
#' @export
#' @seealso \code{\link{hex_pos}} for algorithm that finds hexagon center
#'   closest to each point and \code{\link{hex_coord}} that generates
#'   coordinates of each hexagon.
#' @return A data frame with columns \code{x}, \code{y} and \code{freq},
#'   and attributes \code{width} and \code{height}.

#' @importFrom plyr count
#'
#' @examples
#' plot(hex_bin(runif(1e4), runif(1e4)))
#' plot(hex_bin(rnorm(1e4), rnorm(1e4)))
#'
#' data(baseball, package = "plyr")
#' bin <- hex_bin(baseball$g, baseball$ab)
hex_bin <- function(x, y, weight = NULL, var4 = NULL, width = NULL, height = NULL,
                    xbins = 20, ybins = 20, var4.to.color = FALSE, na.rm = FALSE,
                    hmin = 0, hmax = 1, xlim = NULL, ylim = NULL, ...) {
  if(hmax > 1) warning("hmax > 1 is likely to result in some hexagon overplotting")

  cleaned <- clean_xy(x, y, weight, var4, xlim=xlim, ylim=ylim)

  if (is.null(xlim)) xlim <- c(min(cleaned$x), max(cleaned$x))
  if (is.null(ylim)) ylim <- c(min(cleaned$y), max(cleaned$y))
  if (is.null(width))  width  <- diff(xlim) / xbins
  if (is.null(height)) height <- diff(ylim) / ybins
  height <- height * sqrt(3)

  pos <- hex_pos(cleaned$x, cleaned$y, width, height)
  cleaned$x <- pos[,1]; cleaned$y <- pos[,2]

  # bin values by hexagon
  binned <- plyr::count(cleaned, c("x", "y"), "weight")

  #### var4: another variable to map to size, and will be averaged within each hex bin ####
  var4_sum <- aggregate(cleaned$var4, by=list(cleaned$x, cleaned$y), FUN=mean)
  names(var4_sum) = c('x','y','var4')
  # cols must match order of binned
  binned$var4 <- var4_sum$var4[match(paste(binned$x, binned$y), paste(var4_sum$x, var4_sum$y))]

  # swap the mappings of weight/var4 from color/size to size/color
  names(binned) <- c('x','y','col','size')
  if(var4.to.color) binned <- transform(binned, size=col, col=size)

  # 'size' field now definitely maps to hex area
  if(!is.null(var4) & min(binned$size)<0) warning("size vector cannot include negative values")

  # scale size variable and min-max parameters from hexagon area to side
  binned$size = hex_side(binned$size)
  hmax_a = hex_side(hmax)/hex_side(1)
  hmin_a = hex_side(hmin)/hex_side(1)
  hrange = hmax_a - hmin_a

  # normalise and rescale to custom min-max parameters
  binned$size <- binned$size * hrange / max(binned$size) + hmin_a

  structure(
    binned,
    width = width,
    height = height,
    class = c("bin_hex", "data.frame")
  )
}



clean_xy <- function(x, y, weight=NULL, var4=NULL, na.rm=TRUE, xlim=NULL, ylim=NULL) {
  # If !na.rm, remove missing values with a warning.
  # Otherwise just remove them
  missing <- !is.finite(x) | !is.finite(y)
  nmissing <- sum(missing)

  if (na.rm && nmissing > 0) {
    warning("Removing ", nmissing, " missing values")
  }

  # Check weights, and throw out missing values and zero-weight observations
  if (is.null(weight)) {
    weight <- rep.int(1, length(x))
  } else {
    weight[is.na(weight)] <- 0
  }

  # Check sizes, and throw out missing values and zero-weight observations
  if (is.null(var4)) {
    var4 <- rep.int(1, length(x))
  } else {
    var4[is.na(var4)] <- 0
  }

#   ok <- !missing & weight > 0 & var4 > 0
  ok <- !missing

  # if (all(ok)) { data.frame(x = x, y = y, weight = weight, var4 = var4)
  # } else {
  #   x <- x[ok]
  #   y <- y[ok]
  #   var4 <- var4[ok]
  #   weight <- weight[ok]
  # }

  if (!all(ok)){
    x <- x[ok]
    y <- y[ok]
    var4 <- var4[ok]
    weight <- weight[ok]
  }

  out <- data.frame(x = x, y = y, weight = weight, var4 = var4)

  if (!is.null(xlim)){
    out <- out[out$x >= min(xlim) & out$x <= max(xlim), ]
  }

  if (!is.null(ylim)){
    out <- out[out$y >= min(ylim) & out$y <= max(ylim), ]
  }

  return(out)
}


# plot.bin_hex <- function(x, ...) {
#   if (!require("scales")) {
#     message("Scales package required for plotting 2d densities")
#     return()
#   }

#   if(all(x$col == x$col[1])) { # no variation in colour variable
#     col <- rep("black", nrow(x))
#   } else col <- cscale(x$col, seq_gradient_pal(low = "grey70", high = "red"))

#   hexes <- hex_coord(x=x$x, y=x$y, width=attr(x, "width"), height=attr(x, "height"),
#                      size=x$size, ...)

#   plot(hexes[,1], hexes[,2], type = "n", ...)
#   polygon(hexes, col = col, border = NA)
# }


# Binning algorithms are available for various lattices in dimensions 2-8
# (Conway and Sloane 1982). The following subroutine is a fast FORTRAN
# implementation of hexagonal binning. The key observation is that hexagon
# centers fall on two staggered lattices whose cells are rectangles. Presuming
# the long side of the rectangle is in the y direction, dividing the y
# coordinates by square root (3) [SQRT(3)] makes the cells square. Thus the
# algorithm uses two lattices with square cells. The first lattice has points
# at the integers with [0, 0] as the lower left value. The second lattice is
# shifted so that the lower left value is at [.5 , .5]. The x and y vectors
# are scaled into [0, SIZE] and [0, SIZE / SQRT(3)], respectively. SIZE
# determines the portions of the lattices that are used. For each data point,
# binning consists of finding one candidate lattice point from each lattice
# and then selecting the nearest of the two candidates.

#' Find centre of closest hexagon.
#'
#' @param x numeric x position
#' @param y numeric y position
#' @param width of hexagon
#' @param height of hexagon
#' @return matrix giving position of closest hexagon center
#' @keywords internal
#' @export
#' @examples
#' x <- runif(1e4)
#' y <- runif(1e4)
#' res <- hex_pos(x, y, 0.5, 0.5)
#' plot(x, y, type = "n")
#' segments(x, y, res[, 1], res[, 2], col = "grey80")
#' points(unique(res), pch = 20, cex = 2)


hex_pos <- function(x, y, width, height) {
  height <- height / sqrt(3)

  minx <- min(x, na.rm = TRUE)
  miny <- min(y, na.rm = TRUE)

  # Scale to [0, nrows/ncols]
  sx <- (x - minx) / width
  sy <- (y - miny) / height

  # Find closest center: [0, 0] or [0.5, 0.5]?
  fx <- round(sx)
  fy <- round(sy)

  dist_0 <- 3 * (sx - fx)^2 + (sy - fy)^2
  dist_1 <- 3 * (sx - fx + 0.5)^2 + (sy - fy + 0.5)^2
  dist_2 <- 3 * (sx - fx + 0.5)^2 + (sy - fy - 0.5)^2
  dist_3 <- 3 * (sx - fx - 0.5)^2 + (sy - fy + 0.5)^2
  dist_4 <- 3 * (sx - fx - 0.5)^2 + (sy - fy - 0.5)^2
  dist_smallest <- pmin(dist_0, dist_1, dist_2, dist_3, dist_4)

  x_offset <- rep(0, length(x))
  x_offset[dist_smallest == dist_1] <- +0.5
  x_offset[dist_smallest == dist_2] <- +0.5
  x_offset[dist_smallest == dist_3] <- -0.5
  x_offset[dist_smallest == dist_4] <- -0.5

  y_offset <- rep(0, length(y))
  y_offset[dist_smallest == dist_1] <- +0.5
  y_offset[dist_smallest == dist_2] <- -0.5
  y_offset[dist_smallest == dist_3] <- +0.5
  y_offset[dist_smallest == dist_4] <- -0.5

  # Transform back to original coordinates
  cbind(x = (fx - x_offset) * width + minx, y = (fy - y_offset) * height + miny)
}

#' Generate hexagon coordinates.
#'
#' Long axis is horizontal. Edges clock-wise from far-left, separated by
#' row of missing values.
#'
#' @param x horizontal position of center
#' @param y vertical position of center
#' @param width hex width
#' @param height hex height
#' @export
#' @keywords internal
#' @return A two column matrix with 7 times as many rows as input.
#' @examples
#' x <- runif(1000)
#' y <- runif(1000)
#' res <- unique(hex_pos(x, y, 0.5, 0.5))
#' hexes <- hex_coord(res[, 1], res[, 2], 0.6, 0.5)
#'
#' hexes <- hex_coord(res[, 1], res[, 2], rnorm(1000,.5,.3), rnorm(1000,.5,.3))
#'
#' plot(hexes, type = "n")
#' polygon(hexes)
#' points(res)


hex_coord <- function(x, y, width, height, size = 1) {
  dx <- size * width / 6
  dy <- size * height / 2 / sqrt(3)

  hex_x <- rbind(x - 2 * dx, x - dx, x + dx, x + 2 * dx, x + dx, x - dx, NA)
  hex_y <- rbind(y, y + dy, y + dy, y, y - dy, y - dy, NA)

  cbind(as.vector(hex_x), as.vector(hex_y))
}


hex_coord_df <- function(x, y, width, height, size = 1) {
  # like hex_coord but returns a dataframe of vertices grouped by an id variable
  dx <- size * width / 6
  dy <- size * height / 2 / sqrt(3)

  hex_x <- rbind(x - 2 * dx, x - dx, x + dx, x + 2 * dx, x + dx, x - dx)
  hex_y <- rbind(y, y + dy, y + dy, y, y - dy, y - dy)
  id    <- rep(1:length(x), each=6)

  data.frame(cbind(x=as.vector(hex_x), y=as.vector(hex_y), id))
}


## Functions for calculating hexagon geometries

# hexagon side from area
hex_side <- function(area) sqrt(2 * area / (sqrt(3)*3))

# hexagon area from side (not used)
hex_area <- function(side) side^2 * sqrt(3) * 3/2

# quick function for gg-plotting with 1 line (limited functionality)
qplothex <- function(x, y, var4 = NULL, f.to.a = FALSE, ...){
  bin <- hex_bin(x=x, y=y, var4=var4, frequency.to.area=f.to.a, ...)
  hexes <- hex_coord_df(x=bin$x, y=bin$y, width=attr(bin, "width"), height=attr(bin, "height"), size=bin$size)
  hexes$col <- rep(bin$col, each=6)
  ggplot(hexes, aes(x=x, y=y)) + geom_polygon(aes(fill=col, group=id))
}

png_as_grob <- function(filename) {
    image <- png::readPNG(filename)
    grid::rasterGrob(image, interpolate = TRUE)
}


#' Adaptive palette (discrete)
#'
#' Create a discrete palette that will use the first `n` colors from
#' the supplied color values when the palette has enough colors.
#' Otherwise, use an interpolated color palette.
#'
#' @param values Color values.
pal_ramp <- function(values, force.interp = FALSE) {
  force(values)
  color.func <- function(n) {
    if (force.interp) {
      grDevices:::colorRampPalette(values, alpha = TRUE)(n)
    } else {
        if (n <= length(values)) {
            values[seq_len(n)]
        } else {
            grDevices:::colorRampPalette(values, alpha = TRUE)(n)
        }
    }

  }
  return(color.func)
}


#' Adaptive color palette generator
#'
#' Adaptive color palette generator for ggsci color palettes using `pal_ramp()`.
#'
#' @param name Color palette name in ggsci
#' @param palette Color palette type in ggsci
#' @param alpha Transparency level, a real number in (0, 1].
#'
#' @details See `names(ggsci:::ggsci_db)` for all color palette names in ggsci.
#' See `names(ggsci:::ggsci_db$"pal")` for available palette types under
#' the palette `pal`.
pal_adaptive <- function(name, palette, alpha = 1, force.interp = FALSE) {
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")

  raw_cols <- ggsci:::ggsci_db[[name]][[palette]]
  raw_cols_rgb <- grDevices:::col2rgb(raw_cols)
  alpha_cols <- grDevices:::rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )

  pal_ramp(unname(alpha_cols), force.interp)
}


.checkColorList <- function(object, colors.list = NULL) {
  if (is.null(colors.list)) {
    levels.use <- levels(object@clusters)
    color.func <- pal_adaptive("igv", "default")
    # color.func <- pal_adaptive("nejm", "default")

    color.list <- color.func(length(levels.use))
    names(color.list) <- levels.use
    # colors.list <- as.character(paletteer::paletteer_d("ggsci::default_igv",
    #                                                    n = length(levels.use)))
    # names(colors.list) <- levels.use
    return(colors.list)
  }

  col.len <- length(unique(colors.list))
  
  # if (length(colors.list) != length(levels(object@clusters))) {
  if (!identical(col.len, length(levels(object@clusters)))) {
      stop("The length of `colors.list` is not equal to the number of clusters.")
  }

  col.fac <- object@clusters
  levels(col.fac) <- colors.list

  return(col.fac)
}


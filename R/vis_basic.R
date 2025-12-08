theme_csqc <- function(
        base_size = 11,
        base_family = "",
        base_line_size = base_size/22,
        base_rect_size = base_size/22
) {
    ggplot2::theme(
        panel.border = ggplot2::element_rect(fill = NA, color = 'black', linewidth = 0.5),
        panel.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank()
    )
}

theme_cs_spatial <- function() {
    ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, color = 'black', linewidth = 0.5),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = 'bold')
    )
}

plotTotalCountQC <- function(object, base_size = 11, statsTextSize = base_size/3) {
    if (!'total_counts' %in% colnames(object@metadata)) {
        cli::cli_abort('Variable {.field total_counts} is missing from metadata.')
    }
    log1p_total_counts <- log1p(object$total_counts)
    bpstats <- grDevices::boxplot.stats(log1p_total_counts)
    label_df <- data.frame(
        stats = bpstats$stats,
        label = round(exp(bpstats$stats) - 1, digits = 1)
    )
    data.frame(
        log1p_total_counts = log1p_total_counts
    ) %>%
        ggplot2::ggplot(
            mapping = ggplot2::aes(x = 'Sample', y = .data[['log1p_total_counts']])
        ) +
        ggplot2::geom_violin() +
        ggplot2::geom_boxplot(fill = NA, staplewidth = 0.1) +
        ggplot2::geom_text(
            data = label_df,
            mapping = ggplot2::aes(
                x = 'Sample', y = .data[['stats']], label = .data[['label']]
            ),
            size = statsTextSize,
            hjust = -0.2, vjust = -0.2
        ) +
        ggplot2::scale_y_continuous(
            transform = 'log1p'
        ) +
        ggplot2::labs(
            y = 'log(Total counts + 1)'
        ) +
        theme_csqc(base_size = base_size)
}

plotGeneDetectedQC <- function(object, base_size = 11, statsTextSize = base_size/3) {
    if (!'gene_detected' %in% colnames(object@metadata)) {
        cli::cli_abort('Variable {.field gene_detected} is missing from metadata.')
    }
    log1p_gene_detected <- log1p(object$gene_detected)
    bpstats <- grDevices::boxplot.stats(log1p_gene_detected)
    label_df <- data.frame(
        stats = bpstats$stats,
        label = round(exp(bpstats$stats) - 1, digits = 1)
    )
    data.frame(
        log1p_gene_detected = log1p_gene_detected
    ) %>%
        ggplot2::ggplot(
            mapping = ggplot2::aes(x = 'Sample', y = .data[['log1p_gene_detected']])
        ) +
        ggplot2::geom_violin() +
        ggplot2::geom_boxplot(fill = NA, staplewidth = 0.1) +
        ggplot2::geom_text(
            data = label_df,
            mapping = ggplot2::aes(
                x = 'Sample', y = .data[['stats']], label = .data[['label']]
            ),
            size = statsTextSize,
            hjust = -1, vjust = -0.2
        ) +
        ggplot2::labs(
            y = 'log(Number of genes detected + 1)'
        ) +
        theme_csqc(base_size = base_size)
}



plotSpatial <- function(
        object,
        colorBy = NULL,
        dotSize = 0.5,
        title = NULL,
        zeroAsNA = TRUE,
        legendDotSize = 4,
        viridisOption = 'C',
        colors = NULL,
        naColor = 'grey70'
) {
    colorBy <- colorBy %||% object@parameters$defaultColorBy
    df <- data.frame(
        x = object@spatial[, 1],
        y = object@spatial[, 2]
    )
    titleText <- title
    if (!is.null(colorBy)) {
        colorBy <- rlang::arg_match(
            arg = colorBy,
            values = colnames(object@metadata),
            multiple = FALSE
        )
        titleText <- title %||% colorBy
        showLegendTitle <- ifelse(is.null(title), FALSE, TRUE)
        colorByVal <- object@metadata[[colorBy]]
        if (isTRUE(zeroAsNA) && is.numeric(colorByVal)) {
            colorByVal[colorByVal == 0] <- NA
        }
        df[[colorBy]] <- colorByVal
        p <- ggplot2::ggplot(
            data = df,
            mapping = ggplot2::aes(
                x = .data[['x']],
                y = .data[['y']],
                color = .data[[colorBy]]
            )
        )
    } else {
        showLegendTitle <- FALSE
        p <- ggplot2::ggplot(
            data = df,
            mapping = ggplot2::aes(
                x = .data[['x']],
                y = .data[['y']]
            )
        )
    }
    p <- p +
        ggplot2::geom_point(size = dotSize) +
        ggplot2::ggtitle(titleText) +
        ggplot2::coord_fixed() +
        theme_cs_spatial() +
        ggplot2::theme(
            legend.title = switch(
                as.character(showLegendTitle),
                'TRUE' = ggplot2::element_text(),
                'FALSE' = ggplot2::element_blank()
            )
        )

    if (is.null(colorBy)) {
        return(p)
    } else if (is.logical(colorByVal)) {
        colors <- colors %||% c('TRUE' = 'red', 'FALSE' = 'grey')
        p <- p + ggplot2::scale_color_manual(
            values = colors, na.value = naColor
        ) +
            ggplot2::guides(
                color = ggplot2::guide_legend(override.aes = list(size = legendDotSize))
            )
    } else if (is.factor(colorByVal) || is.character(colorByVal)) {
        colors <- colors %||% csColors
        p <- p + ggplot2::scale_color_manual(
            values = colors,
            na.value = naColor
        ) +
            ggplot2::guides(
                color = ggplot2::guide_legend(override.aes = list(size = legendDotSize))
            )
    } else if (is.numeric(colorByVal)) {
        p <- p + ggplot2::scale_color_viridis_c(direction = -1, na.value = naColor, option = viridisOption)
    }
    return(p)
}

#' Check if parameters look right with scale bar on spatial coordinates
#' @description
#' Make a spatial coordinate plot with a scale bar added. The length of the
#' scale bar represents the radius of the Epsilon ball, which is used later for
#' modeling the neighborhood for diffusion-dependent LR interactions. Users can
#' determine if the scale bar looks appropriate based on the experimental setup
#' of the sliced tissue.
#' @param object A `cytosignal2` object.
#' @param position Position of the scale bar. One of \code{'bl'} (bottom left),
#' \code{'br'} (bottom right), \code{'tl'} (top left), or \code{'tr'} (top
#' right). Default \code{'bl'}.
#' @param direction Direction of the scale bar. Either \code{'horizontal'} or
#' \code{'vertical'}. Default \code{'horizontal'}.
#' @param labelSize Size of the labeled text. Default \code{5}.
#' @return A ggplot2 object showing the spatial coordinates with a scale bar.
#' Error is thrown if the parameters have not been set.
#' @export
plotScale <- function(
        object,
        position = c('bl', 'br', 'tl', 'tr'),
        direction = c('horizontal', 'vertical'),
        labelSize = 5
) {
    position <- match.arg(position)
    direction <- match.arg(direction)

    xrange <- range(object@spatial[, 1])
    yrange <- range(object@spatial[, 2])

    ballRaiusCoord <- object@parameters$ballRadiusCoord
    ballRadiusMicron <- object@parameters$ballRadiusMicron
    if (is.null(ballRaiusCoord) || is.null(ballRadiusMicron)) {
        cli::cli_abort('Set parameters first with {.fn setParams}.')
    }
    xstart <- switch(
        position,
        'bl' = xrange[1] + 0.1*diff(xrange),
        'br' = xrange[2] - 0.1*diff(xrange),
        'tl' = xrange[1] + 0.1*diff(xrange),
        'tr' = xrange[2] - 0.1*diff(xrange)
    )
    ystart <- switch(
        position,
        'bl' = yrange[1] + 0.1*diff(yrange),
        'br' = yrange[1] + 0.1*diff(yrange),
        'tl' = yrange[2] - 0.1*diff(yrange),
        'tr' = yrange[2] - 0.1*diff(yrange)
    )

    if (direction == 'horizontal') {
        xend <- switch(
            position,
            'bl' = xstart + ballRaiusCoord,
            'br' = xstart - ballRaiusCoord,
            'tl' = xstart + ballRaiusCoord,
            'tr' = xstart - ballRaiusCoord
        )
        yend <- ystart
        labelx <- switch(
            position,
            'bl' = xstart + ballRaiusCoord / 2,
            'br' = xstart - ballRaiusCoord / 2,
            'tl' = xstart + ballRaiusCoord / 2,
            'tr' = xstart - ballRaiusCoord / 2
        )
        labely <- ystart
        labelhjust <- 0.5
        labelvjust <- 1.4
        labelangle <- 0
    } else {
        xend <- xstart
        yend <- switch(
            position,
            'bl' = ystart + ballRaiusCoord,
            'br' = ystart + ballRaiusCoord,
            'tl' = ystart - ballRaiusCoord,
            'tr' = ystart - ballRaiusCoord
        )
        labelx <- xstart
        labely <- switch(
            position,
            'bl' = ystart + ballRaiusCoord / 2,
            'br' = ystart + ballRaiusCoord / 2,
            'tl' = ystart - ballRaiusCoord / 2,
            'tr' = ystart - ballRaiusCoord / 2
        )
        labelhjust <- 0.5
        labelvjust <- 1.4
        labelangle <- 90
    }

    plotSpatial(object) +
        ggplot2::annotate(
            geom = 'segment',
            x = xstart, y = ystart,
            xend = xend, yend = yend,
            color = 'black', linewidth = 1
        ) +
        ggplot2::annotate(
            geom = 'text',
            x = labelx, y = labely,
            label = paste0(ballRadiusMicron, ' µm'),
            hjust = labelhjust,
            vjust = labelvjust,
            angle = labelangle,
            size = labelSize
        )
}

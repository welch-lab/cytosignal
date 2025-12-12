theme_csqc <- function(
        base_size = 11,
        base_family = "",
        base_line_size = base_size/22,
        base_rect_size = base_size/22
) {
    theme(
        panel.border = element_rect(fill = NA, color = 'black', linewidth = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
    )
}

theme_cs_spatial <- function(
        titleTextSize = 12,
        subtitleTextSize = 10,
        legendTextSize = 8
) {
    theme(
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black', linewidth = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(face = 'bold', size = titleTextSize),
        plot.subtitle = element_text(size = subtitleTextSize),
        legend.text = element_text(size = legendTextSize)
    )
}

#' Quality control plotting functions
#' @description
#' These are the functions for plotting quality control (QC) metrics via violin
#' plots, with boxplots overlaid and key statistics labeled.
#' @param object A \code{\linkS4class{cytosignal2}} object.
#' @param base_size Base font size. Default \code{11}.
#' @param statsTextSize Font size of the labeled statistics. Default is
#' \code{base_size/3}.
#' @return A ggplot object showing the QC metric.
#' @rdname plotQC
#' @export
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
        ggplot(
            mapping = aes(x = 'Sample', y = .data[['log1p_total_counts']])
        ) +
        geom_violin() +
        geom_boxplot(fill = NA, staplewidth = 0.1) +
        geom_text(
            data = label_df,
            mapping = aes(
                x = 'Sample', y = .data[['stats']], label = .data[['label']]
            ),
            size = statsTextSize,
            hjust = -0.2, vjust = -0.2
        ) +
        scale_y_continuous(
            transform = 'log1p'
        ) +
        labs(
            y = 'log(Total counts + 1)'
        ) +
        theme_csqc(base_size = base_size)
}

#' @rdname plotQC
#' @export
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
        ggplot(
            mapping = aes(x = 'Sample', y = .data[['log1p_gene_detected']])
        ) +
        geom_violin() +
        geom_boxplot(fill = NA, staplewidth = 0.1) +
        geom_text(
            data = label_df,
            mapping = aes(
                x = 'Sample', y = .data[['stats']], label = .data[['label']]
            ),
            size = statsTextSize,
            hjust = -1, vjust = -0.2
        ) +
        labs(
            y = 'log(Number of genes detected + 1)'
        ) +
        theme_csqc(base_size = base_size)
}

#' Spatial coordinate plotting functions
#' @description
#' These are the functions for plotting any potential values stored in a
#' \code{\linkS4class{cytosignal2}} object that can be applied on to the spatial
#' coordinates.
#' @param object A \code{\linkS4class{cytosignal2}} object.
#' @param interaction Interaction name(s) as in the column names of
#' \code{object@LRScore}.
#' @param gene Gene name(s) as in the row names of \code{object@rawData}.
#' @param metadataCol Metadata column name(s) as in the column names of
#' \code{object@metadata}.
#' @param dotSize Size of the dots. Default \code{0.5}.
#' @param dotAlpha Transparency of the dots in the scatter plot. Default
#' \code{1} for no transparency. 0 for fully transparent.
#' @param paletteOption Color palette option for continuous variables. One of
#' letter from A to H, or exact name from \code{magma}, \code{inferno},
#' \code{plasma}, \code{viridis}, \code{cividis}, \code{rocket}, \code{mako},
#' or \code{turbo}.
#' @param paletteDirection Direction of the color palette for continuous
#' variables. Default \code{-1} for dark-solid colors representing high values.
#' Use \code{1} for light colors representing high values.
#' @param zeroAsNA Logical, whether to treat zero values as NA for continuous
#' variables. Works together with argument \code{naColor} to distinguish
#' zero value from low positive values. Default \code{TRUE}.
#' @param naColor Color to use for NA values, and NA values converted from
#' zeros. Default \code{'grey70'}.
#' @param colors Color vector to use for categorical or logical variables.
#' Default \code{NULL} adopts: when categorical, built-in colors \code{csColors}
#' that allows maximum 51 categories; when logical, red TRUE and grey FALSE.
#' @param legendDotSize Dot size in the legend for categorical or logical
#' variables. Default \code{4}.
#' @param legendNCol,legendNRow Number of columns and rows in the legend for
#' categorical or logical variables. Default \code{NULL} tries to keep at most
#' 15 rows.
#' @param titleTextSize,subtitleTextSize,legendTextSize Font sizes for the plot
#' title, subtitle, and legend text. Default \code{12}, \code{10}, and \code{8},
#' respectively.
#' @param ... Graphic parameters passed down to specific functions in this same
#' family.
#' @return If only one value is given to plot, a ggplot object is returned. If
#' multiple values are given, a patchwork object combining all plots is
#' returned. A patchwork wrap can be directly printed for combined view of all
#' plots. Users can also extract individual plots from the patchwork wrap by
#' e.g. \code{plotlist[[1]]}.
#' @rdname plotSpatial
#' @export
plotSpatial <- function(
        object,
        metadataCol = NULL,
        gene = NULL,
        interaction = NULL,
        ...
) {
    if (is.null(metadataCol) &&
        is.null(gene) &&
        is.null(interaction)) {
        # Nothing given, do the default colorBy
        return(plotSpatialMetadata(object = object, ...))
    }
    plotlist <- list()
    if (!is.null(metadataCol)) {
        pmeta <- plotSpatialMetadata(
            object = object,
            metadataCol = metadataCol,
            ...
        )
        if (length(metadataCol) == 1) plotlist[[metadataCol]] <- pmeta
        else {
            for (i in seq_along(metadataCol)) {
                plotlist[[metadataCol[i]]] <- pmeta[[i]]
            }
        }
    }
    if (!is.null(gene)) {
        pgene <- plotSpatialGene(
            object = object,
            gene = gene,
            ...
        )
        if (length(gene) == 1) plotlist[[gene]] <- pgene
        else {
            for (i in seq_along(gene)) {
                plotlist[[gene[i]]] <- pgene[[i]]
            }
        }
    }
    if (!is.null(interaction)) {
        plrscore <- plotSpatialLRScore(
            object = object,
            interaction = interaction,
            ...
        )
        if (length(interaction) == 1) plotlist[[interaction]] <- plrscore
        else {
            for (i in seq_along(interaction)) {
                plotlist[[interaction[i]]] <- plrscore[[i]]
            }
        }
    }
    if (length(plotlist) == 1) {
        return(plotlist[[1]])
    } else {
        return(patchwork::wrap_plots(plotlist))
    }
}


#' @rdname plotSpatial
#' @export
plotSpatialLRScore <- function(
        object,
        interaction,
        dotSize = 0.5,
        dotAlpha = 1,
        paletteOption = 'C',
        paletteDirection = -1,
        zeroAsNA = TRUE,
        naColor = 'grey70',
        titleTextSize = 12,
        subtitleTextSize = 10,
        legendTextSize = 8
) {
    paletteOption <- rlang::arg_match(
        arg = paletteOption,
        values = c(LETTERS[1:8], 'magma', 'inferno', 'plasma', 'viridis',
                   'cividis', 'rocket', 'mako', 'turbo'),
        multiple = FALSE
    )

    # Make data frame for plotting first
    spatial <- object@spatial
    plotDF <- data.frame(
        x = spatial[, 1],
        y = spatial[, 2]
    )
    allVars <- character()
    varType <- character()
    subtitles <- character()

    lrscore <- object@LRScore
    if (is.null(lrscore)) {
        cli::cli_abort('LRScore not yet inferred. Run {.fn inferLRScore} first.')
    }
    interaction <- rlang::arg_match(
        arg = interaction,
        values = colnames(lrscore),
        multiple = TRUE
    )
    values <- as.data.frame(as.matrix(lrscore[, interaction, drop = FALSE]))
    if (isTRUE(zeroAsNA)) values[values == 0] <- NA
    plotDF <- cbind(plotDF, values)

    # Additional annotation for each interaction
    db <- object@intrDB
    intrTypeStr <- db[match(interaction, db$interactors), ] %>%
        pull(.data[['type']]) %>%
        as.character() %>%
        toupper() %>%
        lowwords() %>%
        paste0('-dependent LRScore')
    altNameStr <- db[match(interaction, db$interactors), ] %>%
        mutate(alt_name = ifelse(
            nchar(.data[['alt_name']]) > 0,
            paste0('\nAlt. name: ', .data[['alt_name']]),
            ''
        )) %>%
        pull(.data[['alt_name']])
    subtitles <- c(subtitles, paste0(intrTypeStr, altNameStr))

    .plotSpatialFromDF(
        plotDF = plotDF,
        varNames = interaction,
        subtitles = subtitles,
        varTypes = rep('continuous', length(interaction)),
        dotSize = dotSize,
        dotAlpha = dotAlpha,
        paletteOption = paletteOption,
        paletteDirection = paletteDirection,
        naColor = naColor,
        titleTextSize = titleTextSize,
        subtitleTextSize = subtitleTextSize,
        legendTextSize = legendTextSize
    )
}

#' @rdname plotSpatial
#' @export
plotSpatialGene <- function(
        object,
        gene,
        dotSize = 0.5,
        dotAlpha = 1,
        paletteOption = 'C',
        paletteDirection = -1,
        zeroAsNA = TRUE,
        naColor = 'grey70',
        titleTextSize = 12,
        subtitleTextSize = 10,
        legendTextSize = 8
) {
    paletteOption <- rlang::arg_match(
        arg = paletteOption,
        values = c(LETTERS[1:8], 'magma', 'inferno', 'plasma', 'viridis',
                   'cividis', 'rocket', 'mako', 'turbo'),
        multiple = FALSE
    )

    # Make data frame for plotting first
    spatial <- object@spatial
    plotDF <- data.frame(
        x = spatial[, 1],
        y = spatial[, 2]
    )
    rawdata <- object@rawData
    if (any(!gene %in% rownames(rawdata))) {
        missingGenes <- gene[!gene %in% rownames(rawdata)]
        cli::cli_warn('{.val {missingGenes}} not found.')
        gene <- gene[gene %in% rownames(rawdata)]
    }
    values <- rawdata[gene, , drop = FALSE]
    # Normalize gene expression to log1p(CPM)
    libSize <- object$total_counts
    values@x <- values@x / rep.int(libSize, diff(values@p))
    values@x <- log1p(values@x * 1e6)
    values <- as.data.frame(t(as.matrix(values)))
    if (isTRUE(zeroAsNA)) values[values == 0] <- NA
    plotDF <- cbind(plotDF, values)
    subtitles <- rep('Gene expression in log(CPM+1)', length(gene))

    .plotSpatialFromDF(
        plotDF = plotDF,
        varNames = gene,
        subtitles = subtitles,
        varTypes = rep('continuous', length(gene)),
        dotSize = dotSize,
        dotAlpha = dotAlpha,
        paletteOption = paletteOption,
        paletteDirection = paletteDirection,
        naColor = naColor,
        titleTextSize = titleTextSize,
        subtitleTextSize = subtitleTextSize,
        legendTextSize = legendTextSize
    )
}

#' @rdname plotSpatial
#' @export
plotSpatialMetadata <- function(
        object,
        metadataCol = NULL,
        dotSize = 0.5,
        dotAlpha = 1,
        # Categorical/logical settings
        colors = NULL,
        legendDotSize = 4,
        legendNCol = NULL,
        legendNRow = NULL,
        # Continuous settings
        paletteOption = 'D',
        paletteDirection = -1,
        zeroAsNA = TRUE,
        naColor = 'grey70',
        # Text settings
        titleTextSize = 12,
        subtitleTextSize = 10,
        legendTextSize = 8
) {
    paletteOption <- rlang::arg_match(
        arg = paletteOption,
        values = c(LETTERS[1:8], 'magma', 'inferno', 'plasma', 'viridis',
                   'cividis', 'rocket', 'mako', 'turbo'),
        multiple = FALSE
    )
    metadataCol <- metadataCol %||% object@parameters$cluster
    if (length(metadataCol) > 0) {
        metadataCol <- rlang::arg_match(
            arg = metadataCol,
            values = colnames(object@metadata),
            multiple = TRUE
        )
    }

    # Make data frame for plotting first
    spatial <- object@spatial
    plotDF <- data.frame(
        x = spatial[, 1],
        y = spatial[, 2]
    )
    varTypes <- character()
    values <- object@metadata[, metadataCol, drop = FALSE]
    for (varName in metadataCol) {
        value <- values[[varName]]
        if (is.numeric(value)) {
            varTypes <- c(varTypes, 'continuous')
            if (isTRUE(zeroAsNA)) {
                value[value == 0] <- NA
                values[[varName]] <- value
            }
        } else if (is.logical(value)) {
            varTypes <- c(varTypes, 'logical')
        } else if (is.factor(value) || is.character(value)) {
            varTypes <- c(varTypes, 'categorical')
        } else {
            cli::cli_abort('Metadata column {.field {varName}} is {.cls {class(value)[1]}}. Expecting numeric, logical, factor, or character.')
        }
    }
    plotDF <- cbind(plotDF, values)
    .plotSpatialFromDF(
        plotDF = plotDF,
        varNames = metadataCol,
        subtitles = NULL,
        varTypes = varTypes,
        dotSize = dotSize,
        dotAlpha = dotAlpha,
        colors = colors,
        legendDotSize = legendDotSize,
        legendNCol = legendNCol,
        legendNRow = legendNRow,
        paletteOption = paletteOption,
        paletteDirection = paletteDirection,
        naColor = naColor,
        titleTextSize = titleTextSize,
        subtitleTextSize = subtitleTextSize,
        legendTextSize = legendTextSize
    )
}

.plotSpatialFromDF <- function(
        plotDF,
        varNames,
        varTypes,
        subtitles,
        dotSize,
        dotAlpha,
        # Categorical/logical settings
        colors = NULL,
        legendDotSize = 4,
        legendNCol = NULL,
        legendNRow = NULL,
        # Continuous settings
        paletteOption = 'C',
        paletteDirection = -1,
        naColor = 'grey70',
        # Text settings
        titleTextSize = 12,
        subtitleTextSize = 10,
        legendTextSize = 8
) {
    # Make each plot and fill into a list
    plotlist <- list()
    if (length(varNames) == 0) {
        plotlist[[1]] <- ggplot(
            data = plotDF,
            mapping = aes(x = .data[['x']], y = .data[['y']])
        ) +
            geom_point(size = dotSize, alpha = dotAlpha) +
            coord_fixed() +
            theme_cs_spatial(titleTextSize = titleTextSize,
                             subtitleTextSize = subtitleTextSize,
                             legendTextSize = legendTextSize)
    }
    for (i in seq_along(varNames)) {
        varName <- varNames[i]
        type <- varTypes[i]
        value <- plotDF[[varName]]
        subtitle <- subtitles[i]
        if (type %in% c('categorical')) {
            df <- plotDF[sample(nrow(plotDF)), c('x', 'y', varName)]
        } else {
            # For continuous variables, plot high values on top
            dfna <- plotDF[is.na(plotDF[[varName]]), c('x', 'y', varName)]
            dfpos <- plotDF[!is.na(plotDF[[varName]]), c('x', 'y', varName)]
            dfpos <- dfpos[order(dfpos[[varName]], decreasing = FALSE), ]
            df <- rbind(dfna, dfpos)
        }
        p <- ggplot(
            data = df,
            mapping = aes(
                x = .data[['x']],
                y = .data[['y']],
                colour = .data[[varName]]
            )
        ) +
            geom_point(size = dotSize, alpha = dotAlpha) +
            labs(title = varName) +
            coord_fixed() +
            theme_cs_spatial(titleTextSize = titleTextSize,
                             subtitleTextSize = subtitleTextSize,
                             legendTextSize = legendTextSize)
        if (!is.null(subtitle) && !is.na(subtitle))
            p <- p + labs(subtitle = subtitle)
        if (type == 'categorical') {
            if (is.factor(value)) {
                nkeys <- nlevels(value)
            } else if (is.character(value)) {
                nkeys <- length(unique(value))
            }
            if (is.null(legendNCol) && is.null(legendNRow)) {
                legendNColUse <- ceiling(nkeys / 15)
                legendNColUse <- NULL
                legendNRowUse <- NULL
            } else {
                legendNColUse <- legendNCol
                legendNRowUse <- legendNRow
            }
            if (is.null(colors)) {
                if (nkeys <= length(csColors)) {
                    colorScale <- scale_colour_manual(values = csColors, na.value = naColor)
                } else {
                    colorScale <- scale_colour_hue(na.value = naColor)
                }
            } else {
                colorScale <- scale_colour_manual(values = colors, na.value = naColor)
            }
            p <- p + colorScale +
                guides(
                    colour = guide_legend(
                        title = NULL,
                        override.aes = list(size = legendDotSize),
                        ncol = legendNColUse, nrow = legendNRowUse
                    )
                )
        } else if (type == 'logical') {
            useColors <- colors %||% c('grey', 'red')
            p <- p +
                scale_colour_manual(
                    values = useColors,
                    na.value = naColor
                ) +
                guides(
                    colour = guide_legend(
                        title = NULL,
                        override.aes = list(size = legendDotSize),
                        ncol = legendNCol, nrow = legendNRow
                    )
                )
        } else if (type == 'continuous') {
            p <- p +
                scale_colour_viridis_c(
                    option = paletteOption,
                    direction = paletteDirection,
                    na.value = naColor
                ) +
                guides(
                    colour = guide_colorbar(
                        title = NULL,
                        theme = theme(
                            legend.key.height = unit(1, "null"),
                            legend.key.width = unit(0.4, "cm"),
                            legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
                        )
                    )
                )
        }
        plotlist[[varName]] <- p
    }
    if (length(plotlist) == 1) {
        return(plotlist[[1]])
    } else {
        return(patchwork::wrap_plots(plotlist))
    }
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
        annotate(
            geom = 'segment',
            x = xstart, y = ystart,
            xend = xend, yend = yend,
            color = 'black', linewidth = 1
        ) +
        annotate(
            geom = 'text',
            x = labelx, y = labely,
            label = paste0(ballRadiusMicron, ' µm'),
            hjust = labelhjust,
            vjust = labelvjust,
            angle = labelangle,
            size = labelSize
        )
}

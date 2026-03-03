#' Test gene expression association with interaction significance
#' @description
#' This function performs a one-sided Wilcoxon's rank-sum test on each gene,
#' between the cells that receive significant interaction signal and the cells
#' that do not. The test is done per interaction and per cluster (if available).
#' For each interaction, the function collects the genes that show significant
#' up-regulation when spots receive significant interaction signal. Genes are
#' eventually ranked by an effect size metric that combines the average log
#' fold-change and the number of clusters where the gene is significantly
#' up-regulated.
#' @inheritParams generalParam
#' @param intrFdrThresh FDR threshold to claim that an interaction is
#' significantly observed in a spot. Default \code{0.05}.
#' @param minExp,minSpot Quality control parameters for top interaction
#' selection if \code{interaction} is not provided. See \code{\link{getTopIntr}}.
#' @param degFdrThresh FDR threshold to claim that a gene is differentially
#' expressed in the Wilcoxon's rank-sum test. Default \code{0.05}.
#' @param topN Maximum number of top DEGs to append per interaction per cluster.
#' Default \code{30}.
#' @return A list of character vectors, each containing the DEGs associated
#' with the significance of an interaction. If only one interaction is tested,
#' a single character vector is returned.
#' @export
#' @seealso \code{\link{refineIntrDEG}}, \code{\link{plotIntrDEGHeatmap}}
findIntrDEG <- function(
        object,
        interaction = NULL,
        intrFdrThresh = 0.05,
        minExp = 100,
        minSpot = 100,
        degFdrThresh = 0.05
) {
    intrFDR <- object@significance$spatialFDR
    if (is.null(intrFDR)) {
        cli::cli_abort('Spatial FDR not found. Please run {.fn inferLRScore} first.')
    }

    csdb <- .checkIntrSelection(
        object = object,
        interaction = interaction,
        signifOnly = TRUE,
        fdrThresh = intrFdrThresh,
        minExp = minExp,
        minSpot = minSpot
    )
    nIntr <- nrow(csdb)

    lrscore <- object@LRScore

    clusterVarname <- object@parameters$cluster
    if (is.null(clusterVarname)) {
        cli::cli_alert_warning(c(
            x = 'No cluster variable set. Testing all cells at once',
            i = 'See {.code ?setCluster}.'
        ))
        object@metadata$.__all__ <- factor(rep('All', nrow(object@metadata)))
        clusterVarname <- '.__all__'
    }
    if (!clusterVarname %in% colnames(object@metadata)) {
        cli::cli_abort("Selected cluster variable {.val {clusterVarname}} is not available in metadata")
    }
    clusterVar <- object@metadata[[clusterVarname]]
    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
    nCluster <- nlevels(clusterVar)

    dge <- object@rawData
    libSize <- object$total_counts
    if (is.null(libSize)) {
        cli::cli_abort("Field {.field total_counts} is missing from metadata.")
    }
    dge@x <- dge@x / rep.int(libSize, diff(dge@p))
    dge <- log1p(dge*1e4)

    cli::cli_progress_bar(name = "Detecting DEGs per cluster per interaction", total = nIntr*nCluster)
    resList <- list()
    for (i in seq_len(nIntr)) {
        intrName <- csdb$interactors[i]
        selection <- matrix(
            data = FALSE, nrow = nrow(dge), ncol = nCluster,
            dimnames = list(rownames(dge), levels(clusterVar))
        )
        lfcMat <- matrix(
            data = NA, nrow = nrow(dge), ncol = nCluster,
            dimnames = list(rownames(dge), levels(clusterVar))
        )
        isSignif <- intrFDR[, intrName] < intrFdrThresh
        isSignif[is.na(isSignif)] <- FALSE
        for (j in seq_len(nCluster)) {
            clusterName <- levels(clusterVar)[j]
            inCluster <- clusterVar == clusterName
            if (sum(inCluster & isSignif) < 10 || sum(inCluster & !isSignif) < 10) {
                cli::cli_progress_update()
                next
            }
            label <- factor(isSignif[inCluster])
            dgeSub <- dge[, inCluster, drop = FALSE]
            res <- wilcoxauc(dgeSub, label, alternative = "greater") %>%
                na.exclude() %>%
                filter(
                    .data[['group']] == 'TRUE',
                    .data[['padj']] < degFdrThresh
                    # .data[['logFC']] > logFCThresh
                )
            selection[res$feature, j] <- TRUE
            lfcMat[res$feature, j] <- res$logFC
            cli::cli_progress_update()
        }
        res <- data.frame(
            gene = rownames(dge),
            nCluster = rowSums(selection),
            meanLFC = rowMeans(lfcMat, na.rm = TRUE),
            effectSize = rowMeans(lfcMat, na.rm = TRUE) * sqrt(rowSums(selection))
        ) %>%
            filter(
                .data[['nCluster']] > 0,
            ) %>%
            arrange(-.data[['effectSize']])

        resList[[intrName]] <- res
    }
    cli::cli_progress_done()
    res <- list(
        geneSignif = resList,
        fdrThresh = intrFdrThresh,
        clusterVar = clusterVar
    )
    class(res) <- 'csIntrDEG'
    return(res)
}

#' Refine interaction-associated DEG with elastic net regression
#' @description
#' After detecting DEGs associated with interaction significance using
#' \code{\link{findIntrDEG}}, this function further refines the DEG list
#' by fitting an elastic net regression model per interaction. The further
#' model refines the LRScore by predicting with the model using the expression
#' of the detected DEGs and the cluster labels. The genes that are selected by
#' the model are kept as the final DEG list.
#' @inheritParams generalParam
#' @param intrDEG A \code{csIntrDEG} object returned by
#' \code{\link{findIntrDEG}}.
#' @param minSpot Minimum number of spots a gene must be expressed in to be
#' considered for modeling. Default \code{50}.
#' @param effectSizeThresh Minimum effect size a DEG must have to be considered
#' for modeling. Default \code{1}.
#' @param alphaTest A numeric vector of alpha values to test in the elastic
#' net regression, the model with the lowest loss will be considered. Default
#' \code{seq(0.5, 1, 0.1)} (i.e. 0.5, 0.6, 0.7, ..., 1).
#' @param nFold Number of folds to use in cross-validation. Default \code{5}.
#' @param seed Random seed for reproducibility. Default \code{NULL} runs under
#' user global random states.
#' @return A \code{csIntrDEG} object with refined DEG lists.
#' @export
#' @seealso \code{\link{findIntrDEG}}, \code{\link{plotIntrDEGHeatmap}}
refineIntrDEG <- function(
        object,
        intrDEG,
        minSpot = 50,
        effectSizeThresh = 1,
        alphaTest = seq(0.5, 1, 0.1),
        nFold = 5,
        seed = NULL
) {
    csdb <- .checkIntrSelection(
        object,
        interaction = names(intrDEG$geneSignif)
    )
    clusterVar <- intrDEG$clusterVar
    dge <- object@rawData
    dge <- dge[rowSums(dge > 0) > minSpot, , drop = FALSE]
    intrDEG$models <- list()
    intrDEG$clusterSignif <- list()
    for (i in seq_along(csdb$interactors)) {
        intrName <- csdb$interactors[i]
        score <- object@LRScore[, intrName]
        # Not looping over selection list in case there's a mismatch
        # Together with the check above, this assures that only interactions
        # present in both object database and the result list are used.
        selection <- intrDEG$geneSignif[[intrName]]
        # Only use genes that are shown as DEG + not the ligands/receptors +
        # of high quality (at least in `minSpot` cells)
        sigGenes <- selection$gene[selection$effectSize > effectSizeThresh]
        ligGenes <- unlist(strsplit(csdb$ligands[i], split = ';'))
        recGenes <- unlist(strsplit(csdb$receptors[i], split = ';'))
        geneUse <- setdiff(sigGenes, c(ligGenes, recGenes))
        geneUse <- geneUse[geneUse %in% rownames(dge)]
        matSub <- dge[geneUse, , drop = FALSE]
        cli::cli_process_start("Building model for {.val {intrName}} using {.val {length(geneUse)}} genes")
        # Now build the model as: y ~ gene1 + gene2 + ... + cluster1 + cluster2 + ...
        modelData <- as.data.frame(t(as.matrix(matSub)))
        modelData <- cbind(modelData, cluster = clusterVar)
        # Backticks are important for cluster levels with special characters such as spaces and dashes
        form <- paste0("`", colnames(modelData), "`", collapse = " + ")
        form <- paste0("y ~ ", form)
        form <- stats::formula(form)
        modelData <- cbind(modelData, y = score)
        model <- stats::model.matrix(form, data = modelData)
        # Remove the first '(intercept)' column
        x <- model[, seq(from = 2, to = ncol(model))]
        cli::cli_process_done()

        # Now start fitting
        bestLoss <- Inf
        bestModel <- NULL
        bestAlpha <- NULL
        withr::with_seed(
            seed = seed,
            {
                for (alpha in alphaTest) {
                    cli::cli_process_start("Fitting model for {.val {intrName}} with alpha = {.val {alpha}}; nfolds = {.val {nFold}}")
                    glmnetTMP <- glmnet::cv.glmnet(
                        x = x, y = score,
                        nfolds = nFold, alpha = alpha
                    )

                    if (glmnetTMP$lambda.1se < bestLoss) {
                        bestLoss <- glmnetTMP$lambda.1se
                        bestModel <- glmnetTMP
                        bestAlpha <- alpha
                    }
                    cli::cli_process_done()
                }
                cli::cli_alert_success('Best alpha: {.val {bestAlpha}}')
            }
        )
        lambda <- bestModel[['lambda.1se']]
        # scoreRefine <- glmnet::predict.glmnet(bestModel$glmnet.fit, newx = x, s = lambda)
        modelCoef <- stats::coef(bestModel)
        # extract model selected genes
        featureSignif <- rownames(modelCoef)[modelCoef[,1] != 0]
        featureSignif <- featureSignif[featureSignif != '(Intercept)']
        featureSignif <- trimbackstick(featureSignif)
        geneSignif <- featureSignif[featureSignif %in% geneUse]
        clusterSignif <- setdiff(featureSignif, geneSignif)
        clusterSignif <- sub("^cluster", "", clusterSignif)
        intrDEG$models[[intrName]] <- bestModel
        intrDEG$clusterSignif[[intrName]] <- clusterSignif
        intrDEG$geneSignif[[intrName]] <- selection[match(geneSignif, selection$gene), , drop = FALSE]
    }
    return(intrDEG)
}

#' csIntrDEG - Result holder for interaction-associated DEG analysis
#' @name csIntrDEG
#' @description
#' An S3 class to store the results of interaction-associated DEG detection.
#' It is basically a list object, with class attribute \code{csIntrDEG}.
#'
#' It contains the following fields:
#' \describe{
#'   \item{geneSignif}{A list of data frames, each containing the DEGs
#'   associated with an interaction. Each data frame has the following columns:
#'   \itemize{
#'     \item \code{gene}: Gene name.
#'     \item \code{nCluster}: Number of clusters where the gene is significantly
#'     up-regulated when the interaction is significant.
#'     \item \code{meanLFC}: Mean log fold-change of the gene across the
#'     significant clusters.
#'     \item \code{effectSize}: Effect size metric combining mean log
#'     fold-change and number of clusters.
#'   }}
#'   \item{fdrThresh}{FDR threshold used to define if an interaction is
#'   significantly observed at each spot.}
#'   \item{clusterVar}{A factor vector of cluster labels for each spot.}
#'   \item{models}{Only available after \code{\link{refineIntrDEG}}.
#'
#'   A list of fitted elastic net regression models per interaction.}
#'   \item{clusterSignif}{Only available after \code{\link{refineIntrDEG}}.
#'
#'   A list of character vectors, each containing the clusters that are selected
#'   by the elastic net regression model per interaction.}
#' }
#' @seealso \code{\link{findIntrDEG}}, \code{\link{refineIntrDEG}}, \code{\link{plotIntrDEGHeatmap}}
NULL

#' Print method for csIntrDEG objects
#' @param x A \code{csIntrDEG} object.
#' @param ... Additional arguments (not used).
#' @export
#' @method print csIntrDEG
#' @return Invisible NULL value.
print.csIntrDEG <- function(x, ...) {
    clusterVar <- x$clusterVar
    cat("CytoSignal2 intraction-associated DEGs\n")
    cat("---------------------------------------\n")
    uniqClusters <- levels(clusterVar)
    uniqClusters <- cli::cli_vec(uniqClusters, style = list(`vec-trunc` = 5))
    cat(cli::format_inline(
        '- Clusters ({nlevels(clusterVar)}): {.val {uniqClusters}}\n'
    ))
    cat(cli::format_inline(
        '- FDR threshold: {x$fdrThresh}\n'
    ))
    dfShow <- data.frame(
        Interaction = names(x$geneSignif),
        nDEGs = sapply(x$geneSignif, function(df) sum(df$effectSize > 1)),
        Top = sapply(x$geneSignif, function(df) {
            topDEG <- df$gene[df$effectSize > 1]
            topDEG <- head(topDEG, 5)
            paste(topDEG, collapse = ', ')
        }),
        row.names = 'Interaction'
    )
    print(dfShow)
    return(invisible(NULL))
}


#' Heatmap visualization of interaction-associated DEG expression
#' @description
#' This function makes heatmaps to visualize the expression of
#' interaction-associated DEGs detected by \code{\link{findIntrDEG}} and/or
#' \code{\link{refineIntrDEG}}. For each interaction, a heatmap is made to show
#' the expression of the ligand and receptor genes, and the DEGs. The spots
#' are grouped by their cluster labels and whether they receive significant
#' interaction signal.
#' @details
#' The top \code{topN} DEGs are directly sliced from
#' \code{intrDEG$$geneSignif[["interaction_name"]]}, as the table is pre-sorted.
#' When using a \code{csIntrDEG} object returned by \code{\link{refineIntrDEG}},
#' the clusters that show association with the interaction will be highlighted
#' in black text at the bottom of the heatmap, while other clusters are colored
#' in grey.
#' @inheritParams generalParam
#' @param intrDEG A \code{csIntrDEG} object returned by
#' \code{\link{findIntrDEG}} or \code{\link{refineIntrDEG}}.
#' @param topN Maximum number of top DEGs to show per interaction. Default
#' \code{20}.
#' @param sampleSize Number of spots to randomly sample for visualization.
#' Default \code{5000}.
#' @param clusterLabelSize Font size for cluster labels at the bottom. Default
#' \code{8}.
#' @param geneNameSize Font size for gene names on the right. Default \code{8}.
#' @param rowTitleSize Font size for row titles on the left, such as "ligand",
#' "receptor" and "DEG". Default \code{10}.
#' @param annNameSize Font size for annotation names on the right of the top
#' annotation bar. Default \code{8}.
#' @param titleSize Font size for the heatmap title. Default \code{12}.
#' @param legendTitleSize Font size for the legend title. Default \code{10}.
#' @param legendTextSize Font size for the legend text. Default \code{8}.
#' @return A \code{HeatmapList} object powered by package "ComplexHeatmap" if
#' only one interaction is plotted. A list of \code{HeatmapList} objects if
#' multiple interactions are plotted.
#' @export
#' @seealso \code{\link{findIntrDEG}}, \code{\link{refineIntrDEG}}
plotIntrDEGHeatmap <- function(
        object,
        intrDEG,
        interaction = NULL,
        topN = 20,
        sampleSize = 5000,
        clusterLabelSize = 8,
        geneNameSize = 8,
        rowTitleSize = 10,
        annNameSize = 8,
        titleSize = 12,
        legendTitleSize = 10,
        legendTextSize = 8
) {
    if (!requireNamespace('ComplexHeatmap', quietly = TRUE)) {
        cli::cli_abort(c(
            x = 'Package {.pkg ComplexHeatmap} is required to use this function.',
            i = "Please install it first via {.code BiocManager::install('ComplexHeatmap')}."
        ))
    }
    if (!is.numeric(trim) ||
        length(trim) != 2 ||
        trim[1] >= trim[2]) {
        cli::cli_abort("Argument {.field trim} must be a numeric vector of length 2 with increasing values.")
    }
    selectionList <- intrDEG$geneSignif
    interaction <- interaction %||% names(selectionList)
    csdb <- .checkIntrSelection(object = object, interaction = interaction)
    cli::cli_process_start("Making {nrow(csdb)} heatmap{?s}")
    ligGenes <- sapply(csdb$ligands, strsplit, split = ';')
    recGenes <- sapply(csdb$receptors, strsplit, split = ';')
    fdr <- cs@significance$spatialFDR
    fdrThresh <- intrDEG$fdrThresh
    clusterVar <- intrDEG$clusterVar
    mat <- object@rawData
    mat@x <- mat@x / rep.int(object$total_counts, diff(mat@p))
    mat <- log1p(mat*1e6)
    spotSample <- sample(ncol(mat), min(sampleSize, ncol(mat)))
    fdr <- fdr[spotSample, , drop = FALSE]
    clusterVar <- clusterVar[spotSample]
    colors <- colors %||% stats::setNames(csColors[seq_along(levels(clusterVar))], levels(clusterVar))
    mat <- mat[, spotSample, drop = FALSE]
    allHmList <- list()
    for (i in seq_len(nrow(csdb))) {
        intrName <- csdb$interactors[i]
        intrType <- csdb$type[i]
        titleText <- sprintf('DEGs associated with %s (%s-dependent)', intrName, intrType)
        ligGene <- ligGenes[[i]]
        recGene <- recGenes[[i]]
        highlightCluster <- intrDEG$clusterSignif[[intrName]]
        highlightCluster <- highlightCluster %||% levels(clusterVar)
        # Do not direct access with `i` since intrSelection might remove some
        genes <- selectionList[[intrName]]$gene
        if (length(genes) > topN) genes <- head(genes, n = topN)
        genes <- setdiff(genes, c(ligGene, recGene))
        allGenes <- c(ligGene, recGene, genes)
        geneLabel <- c(rep('Ligand', length(ligGene)),
                       rep('Receptor', length(recGene)),
                       rep('DEG', length(genes)))
        isSignif <- fdr[, intrName] < fdrThresh
        matSub <- mat[allGenes, , drop = FALSE]
        annCol <- data.frame(
            Significance = factor(ifelse(isSignif, 'Significant', 'N/S'))
        )
        rownames(annCol) <- colnames(matSub)
        annRow <- data.frame(
            Gene_Type = factor(geneLabel, levels = c('Ligand', 'Receptor', 'DEG'))
        )
        rownames(annRow) <- rownames(matSub)
        hmList <- lapply(seq_along(levels(clusterVar)), function(i) {
            index <- clusterVar == levels(clusterVar)[i]
            matSubSub <- matSub[, index, drop = FALSE]
            annColSub <- annCol[index, , drop = FALSE]
            matSubSub <- as.matrix(matSubSub)
            hm <- ComplexHeatmap::Heatmap(
                matrix = matSubSub, name = 'Scaled Expression',
                col = circlize::colorRamp2(
                    breaks = seq(from = 0, to = max(matSubSub), length.out = 10),
                    colors = viridis::magma(n = 10, direction = -1)
                ),
                cluster_row_slices = FALSE, cluster_rows = FALSE,
                cluster_column_slices = FALSE, cluster_columns = TRUE,
                column_gap = grid::unit(0, 'mm'),
                top_annotation = ComplexHeatmap::HeatmapAnnotation(
                    df = annColSub,
                    col = list(
                        Significance = c(
                            'Significant' = '#4DAF4A',
                            'N/S' = 'grey70'
                        )
                    ),
                    show_annotation_name = i == nlevels(clusterVar),
                    annotation_name_gp = grid::gpar(fontsize = annNameSize)
                ),

                show_row_names = TRUE, show_column_names = FALSE,
                row_title_rot = 0,
                row_title_gp = grid::gpar(fontsize = rowTitleSize),
                row_names_gp = grid::gpar(fontsize = geneNameSize),

                column_title = levels(clusterVar)[i],
                column_title_side = 'bottom',
                column_title_rot = 45,
                column_title_gp = grid::gpar(
                    fontsize = clusterLabelSize,
                    col = ifelse(
                        levels(clusterVar)[i] %in% highlightCluster,
                        'black',
                        'grey50'
                    )
                ),

                row_split = annRow$Gene_Type,
                column_split = annColSub$Significance
            )
            return(hm)
        })
        hmList <- Reduce(`+`, hmList)
        grDevices::pdf(nullfile())
        hmList <- ComplexHeatmap::draw(
            hmList,
            column_title = titleText,
            column_title_gp = grid::gpar(fontsize = titleSize, fontface = 'bold'),
            legend_title_gp = grid::gpar(fontsize = legendTitleSize, fontface = 'bold'),
            legend_labels_gp = grid::gpar(fontsize = legendTextSize),
            merge_legend = TRUE,
            ht_gap = grid::unit(1, 'mm')
        )
        dev.off()
        allHmList[[intrName]] <- hmList
    }
    cli::cli_process_done()
    if (length(allHmList) == 1) {
        return(allHmList[[1]])
    } else {
        return(allHmList)
    }
}

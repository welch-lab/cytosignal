#' Infer significant genes for each interaction
#' @description
#' This function performs wilcoxon one-sided test for each cell type, between
#' cells enriched with an interaction against other cells. The top significant
#' genes are selected for each interaction. The expression profile of the
#' selected genes, together with the cell type variable, are then feeded into a
#' regression model for further refinement. This analysis infers the significant
#' genes for each interaction, and refines the LRScore using the selected genes.
#' @param object A \code{Cytosignal} object with \code{\link{inferIntrScore}}
#' already run.
#' @param intr Specify interactions to be used. A vector of either the unique
#' ID of interactions or the numeric rank indices. Available IDs can be shown
#' with \code{showIntr(object)}. Availability of an interaction depends on the
#' LRscore slot to be used as well as the significance metric to be used.
#' @param slot.use The LRscore type to be used. See vignette for explanation.
#' @param signif.use The significance metric to be used. See vignette for
#' explanation.
#' @param num.per.choose The maximum number of genes to be chosen for each
#' wilcoxon test. A test happens for each interaction and each cell type.
#' Default \code{50}.
#' @param alpha.test A sequence of alpha values to be tested for the elastic net
#' regression model. The best alpha value will be chosen based on the smallest
#' loss. Larger alpha values result in less final selection. Default
#' \code{seq(0.5, 1, 0.1)}.
#' @param seed Random seed for controlling the cross validation. Default
#' \code{1}.
#' @param minCell Minimum number of cells where a gene is expressed to be
#' considered. Default \code{50}.
#' @param verbose Whether to print progress messages. Default \code{TRUE}.
#' @return list object, each element is the result for each interaction. Each
#' interaction result is a list object containing the following entries:
#' \itemize{
#' \item{\code{glmnet_model} - A "cv.glmnet" object, the regression model.}
#' \item{\code{score_refine} - A 1-column matrix of refined LR scores of this
#' interaction in each cell.}
#' \item{\code{coef} - A matrix of coefficients of the regression model.}
#' \item{\code{sign_clusters} - A character vector of significant clusters in
#' selected by the model.}
#' \item{\code{sign_genes} - A character vector of significant genes selected by
#' the model.}
#' \item{\code{alpha} - The alpha value used for the final chosen model.}
#' \item{\code{slot.use} - The LRScore type used for this analysis.}
#' \item{\code{signif.use} - The significance metric used for this analysis.}
#' }
#' @export
inferIntrDEG <- function(
    object,
    intr = NULL,
    slot.use = NULL,
    signif.use = NULL,
    num.per.choose = 50,
    alpha.test = seq(0.5, 1, 0.1),
    seed = 1,
    minCell = 50,
    verbose = TRUE
) {
  slot.use <- .checkSlotUse(object, slot.use)
  signif.use <- .checkSignifUse(object, signif.use = signif.use,
                                slot.use = slot.use)
  intr_selected_gene <- gene_select_group_by_permute(
    object, intr = intr, slot.use = slot.use, signif.use = signif.use,
    num_per_choose = num.per.choose, verbose = verbose
  )
  res <- refine_score(
    object,
    intr_selected_gene = intr_selected_gene,
    alpha.test = alpha.test, slot.use = slot.use, minCell = minCell,
    seed = seed, verbose = verbose
  )
  res <- lapply(res, function(x) {
    x$slot.use <- slot.use
    x$signif.use <- signif.use
    return(x)
  })
  class(res) <- "CytosignalIntrDEG"
  return(res)
}


# For each significant interaction, iterate through each cluster label and
# perform wilcoxon one-sided (greater) test on cells in the cluster and is
# significant against cells in the cluster and is insignificant. Gather the
# top differentially expressed gene sets for each interaction.
gene_select_group_by_permute <- function(
    object,
    intr = NULL,
    slot.use = NULL,
    signif.use = NULL,
    num_per_choose = 30,
    verbose = TRUE
) {
  expression <- object@raw.counts
  expression <- normCounts.dgCMatrix(expression, scale.fac = Matrix::colSums(expression), method = "default")
  permute_res <- object@lrscore[[slot.use]]@res.list[[signif.use]]
  if (is.null(intr)) {
    intr <- names(permute_res)
  }
  if (is.numeric(intr)) {
    intrRanks <- intr
    interaction <- getCPIs(object, intrRanks, slot.use = slot.use,
                    signif.use = signif.use)
  } else if (is.character(intr)) {
    interaction <- intr
    intrRanks <- getIntrRanks(object, interaction, slot.use = slot.use,
                              signif.use = signif.use)
  }
  clusters <- object@clusters
  expression = expression[rowSums(expression) != 0, ]
  if (!is.factor(clusters)) clusters <- factor(clusters)
  # cluster_type = levels(droplevels(clusters))
  cluster_type <- unique(clusters)
  intr_gene_dict = list()
  if (isTRUE(verbose))
    message("Performing wilcoxon tests to identify significant genes for ",
            "each interaction...")
  totalIter<- length(interaction) * length(cluster_type)
  currIter <- 0
  if (isTRUE(verbose))
    pb <- utils::txtProgressBar(min = 0, max = totalIter, style = 3,
                                file = stderr())
  for (intr in interaction) {
    # for each interaction, reset the chosen genes
    # This is for recording the gene has been chosen and next iteration for
    # cluster level just use the rest of genes
    chosen_genes <- rep(FALSE, nrow(expression))
    chosen_genes_names <- rownames(expression)
    names(chosen_genes) <- chosen_genes_names
    permute_sign = permute_res[[intr]]
    # get the spatial significant cells and spatial insignificant cells
    permute_insig <- colnames(expression)[!colnames(expression) %in% permute_sign]
    for (cell_clust in cluster_type) {
      currIter <- currIter + 1
      cell_idx = names(clusters[clusters == cell_clust])
      if (sum(permute_sign %in% cell_idx) < 10 ||
          sum(permute_insig %in% cell_idx) < 10) {
        if (isTRUE(verbose)) utils::setTxtProgressBar(pb, currIter)
        next
      }
      sig_use <- permute_sign[permute_sign %in% cell_idx]
      insig_use <- permute_insig[permute_insig %in% cell_idx]
      expr_test <- expression[!chosen_genes, c(sig_use, insig_use)]
      grouping <- factor(c(rep(1, length(sig_use)), rep(2, length(insig_use))))
      test.res <- wilcoxauc(expr_test, grouping, alternative = "greater")
      test.res <- na.exclude(test.res)
      test.res <- test.res[test.res$group == 1,]
      test.res <- test.res[order(test.res$padj),]
      test.res <- test.res[seq(num_per_choose),]
      test.res <- test.res[test.res$padj < 0.05,]
      chosen_genes[test.res$feature] <- TRUE
      if (isTRUE(verbose)) utils::setTxtProgressBar(pb, currIter)
    }
    intr_gene_dict[[intr]] = chosen_genes_names[chosen_genes]
  }
  if (isTRUE(verbose)) message()
  return(intr_gene_dict)
}


# x: dgCMatrix - feature x observation
# clusterVar: vector/factor object, grouping assignment
# alternative: character, one of "two.sided", "greater", "less", can be
# partially matched.
# Rcpp source code located in src/wilcoxon.cpp
# Credit to GitHub repo immunogenomics/presto, with massive modification so this
# version is minimized to only work for dgCMatrix and only produces p-value and
# adjusted p-value. Alternative assumption is added so it can be aware of the
# test that our analysis needs.
wilcoxauc <- function(
    x,
    clusterVar,
    alternative = c("two.sided", "greater", "less")
) {
  alternative <- match.arg(alternative)
  if (methods::is(x, 'dgTMatrix')) x <- methods::as(x, 'CsparseMatrix')
  if (methods::is(x, 'TsparseMatrix')) x <- methods::as(x, 'CsparseMatrix')
  if (is.null(row.names(x))) {
    rownames(x) <- paste0('Feature', seq(nrow(x)))
  }
  if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
  clusterVar <- droplevels(clusterVar)
  groupSize <- as.numeric(table(clusterVar))

  ## Compute primary statistics
  n1n2 <- groupSize * (ncol(x) - groupSize)
  # rankRes - list(X_ranked, ties), where X_ranked is obs x feature
  xRanked <- Matrix::t(x)
  # This computes the ranking of non-zero values and the ties
  ties <- cpp_rank_matrix_dgc(xRanked@x, xRanked@p,
                              nrow(xRanked), ncol(xRanked))
  ustat <- computeUstat(xRanked, clusterVar, n1n2, groupSize)
  pvals <- computePval(ustat, ties, ncol(x), n1n2, alternative = alternative)
  fdr <- apply(pvals, 2, function(p) stats::p.adjust(p, 'BH'))

  data.frame(
    feature = rep(row.names(x), times = length(levels(clusterVar))),
    group = factor(rep(levels(clusterVar), each = nrow(x)),
                   levels = levels(clusterVar)),
    # statistics = as.numeric(t(ustat)),
    # pval = as.numeric(pvals),
    padj = as.numeric(fdr)
  )
}

computeUstat <- function(Xr, cols, n1n2, groupSize) {
  grs <- rowAggregateSum_sparse(Xr, as.integer(cols) - 1, length(unique(cols)))
  # With the ranking of only non-zero features, here the tie-ranking of
  # zeros need to be added.
  nnz <- rowNNZAggr_sparse(Xr, as.integer(cols) - 1, length(unique(cols)))
  gnz <- groupSize - nnz
  zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
  ustat <- t((t(gnz) * zero.ranks)) + grs - groupSize*(groupSize + 1)/2
  return(ustat)
}

computePval <- function(ustat, ties, N, n1n2, alternative) {
  z <- ustat - .5 * n1n2
  CORRECTION <- switch(alternative,
                       "two.sided" = sign(z) * 0.5,
                       "greater" = 0.5,
                       "less" = -0.5)
  .x1 <- N ^ 3 - N
  .x2 <- 1 / (12 * (N ^ 2 - N))
  rhs <- unlist(lapply(ties, function(tvals) {
    (.x1 - sum(tvals ^ 3 - tvals)) * .x2
  }))
  usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
  z <- t((z - CORRECTION) / usigma)
  pvals <- switch(alternative,
    "greater" = matrix(stats::pnorm(as.numeric(z), lower.tail = FALSE), ncol = ncol(z)),
    "two.sided" = matrix(2 * stats::pnorm(-abs(as.numeric(z))), ncol = ncol(z)),
    "less" = matrix(stats::pnorm(as.numeric(z), lower.tail = TRUE), ncol = ncol(z))
  )

  return(pvals)
}



refine_score <- function(
    object,
    intr_selected_gene,
    alpha.test = seq(0.5, 1, 0.1),
    slot.use = NULL,
    minCell = 50,
    seed = 1,
    verbose = TRUE
) {
  slot.use <- .checkSlotUse(object, slot.use)
  # signif.use <- .checkSignifUse(object, signif.use = signif.use, slot.use = slot.use)
  expression <- object@raw.counts
  clusters <- object@clusters
  expression <- normCounts.dgCMatrix(expression,
                                     scale.fac = Matrix::colSums(expression),
                                     method = "default")
  # clean the gene names
  gene_to_uniprot <- object@intr.valid$gene_to_uniprot
  rownames(gene_to_uniprot) <- gene_to_uniprot$uniprot
  # gaupes_uni <- object@intr.valid$diff_dep$combined
  diff_uni <- object@intr.valid$diff_dep$combined
  cont_uni <- object@intr.valid$cont_dep$combined
  # remove genes expressed less than 50 cells
  genes_sum <- rowSums(expression > 0)
  genes_lowqual <- names(which(genes_sum <= minCell))
  # expression <- expression[!rownames(expression) %in% genes_rm,]
  scores <- object@lrscore[[slot.use]]@score
  res.list <- list()
  if (isTRUE(verbose)) {
    message("Traning regression model to further test gene significance...")
    pb <- utils::txtProgressBar(
      min = 0,
      max = length(intr_selected_gene) * length(alpha.test),
      style = 3, file = stderr()
    )
  }

  curr <- 0
  for (test_intr in names(intr_selected_gene)) {
    selected_genes <- intr_selected_gene[[test_intr]]
    exclude_uni <- names(cont_uni[cont_uni == test_intr])
    exclude_gene_upper <- gene_to_uniprot[exclude_uni,"gene_name"]
    # exlcude genes that make the interaction
    exclude_names<- lowwords(exclude_gene_upper)
    selected_genes <- selected_genes[!selected_genes %in% exclude_names]
    selected_genes <- selected_genes[!selected_genes %in% genes_lowqual]

    model_data <- as.data.frame(t(as.matrix(expression[selected_genes, , drop = FALSE])))
    model_data <- cbind(model_data, clusters)
    y <- scores[, test_intr]
    # fomurla of the linear model, clusters are also included as covariates, y ~ g1 + g2 + .... + c1 + c2....
    form <- paste0("`", colnames(model_data), "`", collapse = " + ")
    form <- paste0("y ~ ", form)
    form <- stats::formula(form)
    model_data <- cbind(model_data,y)
    model_X <- stats::model.matrix(form, data = model_data)
    x <- model_X[, seq(2, ncol(model_X))]
    loss_best <- Inf
    model_best <- NULL
    alpha_best <- NULL
    set.seed(seed)
    # tunning the elastic net
    for (alpha in seq(0.5, 1, 0.1)) {
      curr <- curr + 1
      glmnet_tmp <- glmnet::cv.glmnet(x = x, y = y, nfolds = 5,
                                     seed = seed, alpha = alpha)
      if (glmnet_tmp$lambda.1se < loss_best) {
        loss_best <- glmnet_tmp$lambda.1se
        model_best <- glmnet_tmp
        alpha_best <- alpha
      }
      if (isTRUE(verbose)) utils::setTxtProgressBar(pb, curr)
    }
    score_refine <- glmnet:::predict.cv.glmnet(model_best, newx = x)
    tmp_coef <- stats::coef(model_best)
    # extract model selected genes
    sign_features = rownames(tmp_coef)[tmp_coef[,1] != 0]
    sign_features = sign_features[2:length(sign_features)]
    sign_features <- trimbackstick(sign_features)
    sign_genes = sign_features[sign_features %in% selected_genes]
    sign_clusters <- setdiff(sign_features, sign_genes)
    sign_clusters <- sub("^clusters", "", sign_clusters)
    res.list[[test_intr]] <- list(
      glmnet_model = model_best,
      score_refine = score_refine,
      coef = tmp_coef,
      sign_clusters = sign_clusters,
      sign_genes = sign_genes,
      alpha = alpha_best
    )
  }
  if (isTRUE(verbose)) message()
  return(res.list)
}

#' Plot the refined score of each interaction after regression model refinement
#' @param object A \code{Cytosignal} object with \code{\link{inferIntrScore}}
#' already run.
#' @param intrDEGRes The direct output object of \code{\link{inferIntrDEG}}.
#' @param intr A vector of unique interaction IDs that are available in
#' \code{intrDEGRes}, or numerical index within its range. Default \code{NULL}
#' use all the results.
#' @param pt.size Size of the points in the plot. Default \code{0.5}.
#' @param pt.stroke Stroke size of the points in the plot. Default \code{0.1}.
#' @return List of ggplot objects, each shows the refined LR score of each
#' selected interaction.
plotRefinedScore <- function(
    object,
    intrDEGRes,
    intr = NULL,
    pt.size = 0.5,
    pt.stroke = 0.1
) {
  plotDF <- as.data.frame(object@cells.loc)
  if (is.null(intr)) intr <- names(intrDEGRes)
  if (is.numeric(intr) || is.logical(intr)) intr <- names(intrDEGRes)[intr]
  plotList <- list()
  xRange <- c(min(plotDF$x), max(plotDF$x))
  yRange <- c(min(plotDF$y), max(plotDF$y))
  for (interaction in intr) {
    intrName <- getIntrNames(object, interaction)
    score_refine <- intrDEGRes[[interaction]]$score_refine[rownames(plotDF),]
    plotDF$score <- score_refine
    plotList[[interaction]] <-
      ggplot2::ggplot(plotDF, ggplot2::aes(x = .data[["x"]],
                                           y = .data[["y"]],
                                           color = .data[["score"]])) +
      ggplot2::geom_point(size = pt.size, stroke = 0.2) +
      ggplot2::scale_color_viridis_c(option = "plasma", direction = -1,
                                     na.value = '#F5F5F5') +
      ggplot2::labs(title = intrName, color = "Refined Score",
                    x = NULL, y = NULL) +
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
  return(plotList)
}

#' Show significant genes across top GO term hits with coefficients from
#' regression analysis of an interaction
#' @description
#' Create a heatmap for an interaction on GO terms by significant genes, colored
#' by the coefficients of the genes in the regression model returned by
#' \code{\link{inferIntrDEG}}. The GO term enrichment can be done with any tools
#' available. The input data.frame \code{GO} must contain fields for 1. term
#' description, character, pointed to by \code{description.col}, 2. p-value,
#' numeric, pointed to by \code{pval.col}, 3. gene hit string, character,
#' pointed to by \code{gene.col}. The gene hit string for each term must be form
#' in a way that function \code{gene.split.fun} can split it into a character
#' vector of gene names. For example, if a gene hit string is "gene1, gene2, gene3",
#' then \code{gene.split.fun} should be \code{function(x) unlist(strsplit(x, ", "))},
#' so that a split result \code{c("gene1", "gene2", "gene3")} can be correctly
#' obtained.
#'
#' @param intrDEG A \code{CytosignalIntrDEG} object, output from
#' \code{\link{inferIntrDEG}}.
#' @param GO A data.frame object for GO enrichment analysis result.
#' @param intr A single interaction ID (starts with "CPI") or its numeric index
#' within the range of \code{intrDEG}.
#' @param description.col,pval.col,gene.col The column names of the data.frame
#' \code{GO} that contains the term description, p-value, and gene hit string.
#' Default \code{"description"}, \code{"pval"}, \code{"genes"}.
#' @param term.topN Use this number of top GO terms, ranked by p-values. Default
#' \code{20}.
#' @param gene.topN Use this number of top genes, ranked by absolute value of
#' coefficients. Default \code{20}.
#' @param color_num Number of colors in the heatmap. Can only use \code{2} or
#' \code{3}, Default \code{2} use white-red color palette. \code{3} use scaled
#' blue-white-red color palette.
#' @param binary_sign Whether to convert coefficient value to binary sign value.
#' Default \code{FALSE}.
#' @export
heatmap_GO <- function(
    intrDEG,
    GO,
    intr,
    description.col = "description",
    pval.col = "pval",
    gene.col = "genes",
    gene.split.fun = function(x) unlist(strsplit(x, ",")),
    term.topN = 20,
    gene.topN = 20,
    binary_sign = FALSE,
    text.size = 10
) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package pheatmap is required for this function to work. ",
         "Please install it with command: `install.packages('pheatmap')`.")
  }
  if ((is.numeric(intr) && (intr < 1 || intr > length(intrDEG))) ||
      (is.character(intr) && !intr %in% names(intrDEG))) {
    stop("Invalid interaction ID or index.")
  }
  data <- intrDEG[[intr]]
  tmp_coef <- data$coef
  rownames(tmp_coef) <- trimbackstick(rownames(tmp_coef))
  # remove genes that are filtered out by the elastic net model
  sign_features <- rownames(tmp_coef)[tmp_coef[,1] != 0]
  sign_features_weights = tmp_coef[tmp_coef[,1] != 0,]
  enet_df = data.frame(sign_features = sign_features, weights = sign_features_weights)
  # Remove the "(Intercept)"
  enet_df = enet_df[2:nrow(enet_df),]
  # sort the data by the absolute value of the weights
  enet_df = enet_df[order(abs(enet_df$weights), decreasing = TRUE), ]
  enet_df = enet_df[!startsWith(enet_df$sign_features, "cluster"), , drop = FALSE]
  enet_df <- enet_df[seq_len(min(nrow(enet_df), gene.topN)),]

  GO <- GO[order(GO[[pval.col]]),]
  GO <- GO[seq_len(min(nrow(GO), term.topN)),]

  # signif genes x go term
  occur_matrix = matrix(0, nrow = nrow(enet_df), ncol = nrow(GO),
                        dimnames = list(enet_df$sign_features,
                                        GO[[description.col]]))

  for (j in 1:ncol(occur_matrix)) {
    genes_string = GO[j, gene.col]
    genes_vec = gene.split.fun(genes_string)
    new_values = rep(0, nrow(occur_matrix))
    # fill the values with weights if gene hit the GO term
    new_values[enet_df$sign_features %in% genes_vec] <- enet_df$weights[enet_df$sign_features %in% genes_vec]
    occur_matrix[,j] <- if (isTRUE(binary_sign)) sign(new_values) else new_values
  }

  # filter GO terms without any selected genes hit
  occur_matrix <- occur_matrix[rowSums(occur_matrix != 0 ) > 0, colSums(occur_matrix != 0) > 0]

  color_palette <- grDevices::colorRampPalette(c("#2166AC","white", "firebrick"))(99)
  posnegmax <- max(abs(occur_matrix))
  breaks <- seq(from = -posnegmax, to = posnegmax, length.out = 100)

  pheatmap::pheatmap(t(occur_matrix),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     color = color_palette,
                     breaks = breaks,
                     fontsize = text.size,
                     scale = "none")
}

# e.g. "`hello`" becomes "hello"
trimbackstick <- function(x) {
  ifelse(startsWith(x, "`") & endsWith(x, "`"),
         sub("`$", "", sub("^`", "", x)),
         x)
}


#' Format a CytosignalIntrDEG object to string
#' @param x A \code{CytosignalIntrDEG} object.
#' @param ... Passed to other methods
#' @return A string representation of the object
#' @export
#' @method format CytosignalIntrDEG
format.CytosignalIntrDEG <- function(x, ...) {
  title <- "Cytosignal::inferIntrDEG() result"
  nIntr <- length(x)
  signif.use <- x[[1]]$signif.use
  slot.use <- x[[1]]$slot.use
  msg <- paste0(
    title, "\n",
    "Number of interactions: ", nIntr, "\n",
    "LRScore type (`slot.use`): ", slot.use, "\n",
    "Significance metric (`signif.use`): ", signif.use
  )
  return(msg)
}

#' Print the CytosignalIntrDEG object representation to screen
#' @param x A \code{CytosignalIntrDEG} object.
#' @param ... Passed to other methods
#' @export
#' @method print CytosignalIntrDEG
print.CytosignalIntrDEG <- function(x, ...) {
  cat(format(x), "\n")
}

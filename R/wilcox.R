# X: matrix of data to be tested
# y: grouping label of columns of X
# Rcpp source code located in src/wilcoxon.cpp
wilcoxauc <- function(x, clusterVar) {
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
  # ranksRes <- list(X_ranked = xT, ties = ties)

  # rankRes <- colRanking(x)
  ustat <- computeUstat(xRanked, clusterVar, n1n2, groupSize)
  auc <- t(ustat / n1n2)
  pvals <- computePval(ustat, ties, ncol(x), n1n2)
  fdr <- apply(pvals, 2, function(p) stats::p.adjust(p, 'BH'))

  ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
  groupSums <- colAggregateSum_sparse(x, as.integer(clusterVar) - 1, length(unique(clusterVar)))
  # groupSums <- colAggregateSum(x, clusterVar)
  group_nnz <- colNNZAggr_sparse(x, as.integer(clusterVar) - 1, length(unique(clusterVar)))
  # group_nnz <- colNNZAggr(x, clusterVar)
  group_pct <- t(sweep(group_nnz, 1, as.numeric(table(clusterVar)), "/"))

  group_pct_out <- sweep(-group_nnz, 2, colSums(group_nnz), "+")
  group_pct_out <- sweep(group_pct_out, 1,
                         as.numeric(length(clusterVar) - table(clusterVar)),
                         "/")
  group_pct_out <- t(group_pct_out)

  groupMeans <- t(sweep(groupSums, 1, as.numeric(table(clusterVar)), "/"))

  cs <- colSums(groupSums)
  gs <- as.numeric(table(clusterVar))
  lfc <- Reduce(cbind, lapply(seq_along(levels(clusterVar)), function(g) {
    groupMeans[, g] - (cs - groupSums[g, ])/(length(clusterVar) - gs[g])
  }))

  data.frame(
    feature = rep(row.names(x), times = length(levels(clusterVar))),
    group = factor(rep(levels(clusterVar), each = nrow(x)),
                   levels = levels(clusterVar)),
    avgExpr = as.numeric(groupMeans),
    logFC = as.numeric(lfc),
    statistic = as.numeric(t(ustat)),
    auc = as.numeric(auc),
    pval = as.numeric(pvals),
    padj = as.numeric(fdr),
    pct_in = as.numeric(100 * group_pct),
    pct_out = as.numeric(100 * group_pct_out)
  )
}

computeUstat <- function(Xr, cols, n1n2, groupSize) {
  grs <- rowAggregateSum_sparse(Xr, as.integer(cols) - 1, length(unique(cols)))
  # grs <- rowAggregateSum(Xr, cols)

  # if (inherits(Xr, 'dgCMatrix')) {
  # With the ranking of only non-zero features, here the tie-ranking of
  # zeros need to be added.
  nnz <- rowNNZAggr_sparse(Xr, as.integer(cols) - 1, length(unique(cols)))
  gnz <- groupSize - nnz
  zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
  ustat <- t((t(gnz) * zero.ranks)) + grs - groupSize*(groupSize + 1)/2
  # } else {
  #     ustat <- grs - groupSize * (groupSize + 1) / 2
  # }
  return(ustat)
}

computePval <- function(ustat, ties, N, n1n2) {
  z <- ustat - .5 * n1n2
  z <- z - sign(z) * .5
  .x1 <- N ^ 3 - N
  .x2 <- 1 / (12 * (N ^ 2 - N))
  rhs <- unlist(lapply(ties, function(tvals) {
    (.x1 - sum(tvals ^ 3 - tvals)) * .x2
  }))
  usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
  z <- t(z / usigma)
  pvals <- matrix(2 * stats::pnorm(-abs(as.numeric(z))), ncol = ncol(z))
  return(pvals)
}

validGenes <- function(object) {
    rawData <- object@rawData
    dbGenes <- .uniqGeneInDB(object@intrDB)
    geneUse <- intersect(rownames(rawData), dbGenes)
    sort(geneUse)
}

imputeNiche2 <- function(
        object,
        type = c('diffusion', 'contact')
) {
    type <- match.arg(type)

    if (type == 'diffusion') {
        graph <- object@neighborDiff
    } else if (type == 'contact') {
        graph <- object@neighborCont
    }

    if (is.null(graph)) {
        cli::cli_abort('{type}-dependent neighbor graph not built. Please run {.fn findNeighbor} first.')
    }

    # if (is.null(object$total_counts)) {
    #     object$total_counts <- Matrix::colSums(object@rawData)
    # }
    geneUse <- validGenes(object)
    rawData <- object@rawData[geneUse, , drop = FALSE]

    imputed <- rawData %*% graph;

    # weighted sum of scale factors as well
    # libSize <- matrix(object$total_counts, nrow = 1)
    # libSizeImputed <- libSize %*% graph
    # libSizeImputed <- as.numeric(libSizeImputed)

    return(imputed)
}

weightedLibSize <- function(
        object,
        type = c('diffusion', 'contact')
) {
    type <- match.arg(type)

    if (type == 'diffusion') {
        graph <- object@neighborDiff
    } else if (type == 'contact') {
        graph <- object@neighborCont
    }

    if (is.null(graph)) {
        cli::cli_abort('{type}-dependent neighbor graph not built. Please run {.fn findNeighbor} first.')
    }

    if (is.null(object$total_counts)) {
        object$total_counts <- Matrix::colSums(object@rawData)
    }
    libSize <- matrix(object$total_counts, nrow = 1)
    libSizeImputed <- libSize %*% graph
    libSizeImputed <- as.numeric(libSizeImputed)

    return(libSizeImputed)
}

#' Infer the LR score with significance tested
#' @description
#' This function imputes the niche expression of potential ligand and receptor
#' genes using the neighbor graphs, and computes the final LR score. Optionally,
#' with \code{smoothR = TRUE}, this function allows smoothing the receptor
#' expression of a spot using the neighbor graph of contact-dependent
#' neighborhood. This is useful when the data is of high spatial resolution
#' and/or high sparsity. This function further performs a permutation test to
#' infer the significance of each LR interaction in each spot.
#' @param object A \code{\linkS4class{cytosignal2}} object, with neighbor graphs
#' constructed.
#' @param smoothR Logical, whether to smooth the receptor expression using the
#' contact-dependent neighbor graph. Default is \code{FALSE}.
#' @param permSize Integer, sample size to reach in the permutation test.
#' Default \code{1e5}.
#' @return A \code{\linkS4class{cytosignal2}} object with LR scores inferred and
#' stored in the \code{lrScore} slot.
#' @export
inferLRScore <- function(
        object,
        smoothR = FALSE,
        permSize = 1e5
) {
    geneUse <- validGenes(object)
    cli::cli_process_start(
        "Calculating LRscore using {length(geneUse)} genes found in database"
    )
    # Start with diffusion-dependent interactions
    diffImpL <- imputeNiche2(object, 'diffusion')
    diffLibSize <- weightedLibSize(object, 'diffusion')
    diffImpL@x <- diffImpL@x / rep.int(diffLibSize, diff(diffImpL@p))
    rm(diffLibSize)
    diffImpL@x <- log1p(diffImpL@x * 1e4)

    contImpL <- imputeNiche2(object, 'contact')
    contLibSize <- weightedLibSize(object, 'contact')
    contImpL@x <- contImpL@x / rep.int(contLibSize, diff(contImpL@p))
    rm(contLibSize)
    contImpL@x <- log1p(contImpL@x * 1e4)

    raw <- object@rawData[geneUse, , drop = FALSE]
    object@parameters$smoothR <- smoothR
    if (isTRUE(smoothR)) {
        recepImp <- contImpL
    } else {
        recepImp <- raw
        recepImp@x <- recepImp@x / rep.int(object$total_counts, diff(recepImp@p))
        recepImp@x <- log1p(recepImp@x * 1e4)
    }

    # Locally smooth the L and R imputation
    DTgraph <- object@neighborCont
    dtAvgGraph <- to_mean(DTgraph)
    diffImpL <- diffImpL %*% dtAvgGraph
    contImpL <- contImpL %*% dtAvgGraph
    recepImp <- recepImp %*% dtAvgGraph
    rm(DTgraph, dtAvgGraph)
    colnames(diffImpL) <- colnames(object@rawData)

    diffLigMap <- intrGeneMap(object, type = 'diffusion', component = 'ligands')
    diffRecMap <- intrGeneMap(object, type = 'diffusion', component = 'receptors')
    LRdiff <- multiply_lr_cpp(diffImpL, diffLigMap, recepImp, diffRecMap)
    dimnames(LRdiff) <- list(
        colnames(object@rawData),
        colnames(diffLigMap)
    )
    rm(diffImpL, diffLigMap, diffRecMap)
    contLigMap <- intrGeneMap(object, type = 'contact', component = 'ligands')
    contRecMap <- intrGeneMap(object, type = 'contact', component = 'receptors')
    LRcont <- multiply_lr_cpp(contImpL, contLigMap, recepImp, contRecMap)
    dimnames(LRcont) <- list(
        colnames(object@rawData),
        colnames(contLigMap)
    )
    rm(contImpL, contLigMap, contRecMap, recepImp)

    lrscore <- cbind(LRdiff, LRcont)
    lrscore <- lrscore[, object@intrDB$interactors]
    # Last round of smoothing
    lrscore <- t(object@neighborDiff) %*% lrscore
    rownames(lrscore) <- colnames(object@rawData)
    object@LRScore <- lrscore
    cli::cli_process_done()

    cli::cli_alert_info(
        'Permuting the whole dataset to test the significance of LRscores'
    )
    object <- permuteTest(object, permSize = permSize)

    cli::cli_alert_info(
        "Adjusting p-values spatially with neighbor weighting"
    )
    hasContNeighbor <- which(diff(object@neighborCont@p) > 0)
    pval <- object@significance$pval[hasContNeighbor, , drop = FALSE]
    contGraph <- object@neighborCont[hasContNeighbor, hasContNeighbor, drop = FALSE]

    fdr <- matrix(NA, nrow = nrow(lrscore), ncol = ncol(lrscore))
    dimnames(fdr) <- dimnames(lrscore)
    fdrSub <- spatialGraphFDR_cpp(
        pval = pval,
        intrType = as.integer(intrDB(object)$type) - 1,
        contGraph = contGraph
    )
    fdr[hasContNeighbor, ] <- fdrSub
    object@significance$spatialFDR <- fdr
    return(object)
}

#' Calculates ECDF of null LR scores via permutation test
#' @description
#' This function performs a permutation test to derive the empirical proportion
#' of null scores not exceeding the observed score, also called an empirical
#' cumulative distribution function (ECDF). This allows the inference of
#' one-sided permutation p-values for each LR interaction in each spot by
#' subtracting the ECDF value from 1.
#' @details
#' The permutation test is performed by shuffling the spatial neighbor
#' relationships to impute null \eqn{L} and \eqn{R} values. The same procedure
#' as in \code{\link{inferLRScore}} is applied to compute null LR scores.
#' @param object A \code{\linkS4class{cytosignal2}} object with LR scores
#' inferred with \code{\link{inferLRScore}}.
#' @param permSize Integer, sample size to reach in the permutation test.
#' Default \code{1e5}.
#' @return The input object with slot \code{nullECDF}
#' @noRd
permuteTest <- function(
        object,
        permSize = 1e5
) {
    set.seed(1)
    lrscore <- object@LRScore
    if (is.null(lrscore)) {
        cli::cli_abort("LR score not found. Please run {.fn inferLRScore} first.")
    }
    recepUseDT <- object@parameters$smoothR
    if (!is.logical(recepUseDT)) {
        cli::cli_abort("Parameter {.field smoothR} corrupted. Please re-run {.fn inferLRScore}.")
    }

    if (!is.numeric(permSize) ||
        permSize < 1 ||
        permSize != round(permSize)) {
        cli::cli_abort("{.field permSize} must be a positive integer.")
    }

    n <- nrow(lrscore)
    if (permSize < n){
        cli::cli_warn("Permutation size smaller than number of spots ({n}), using {n} instead.")
        permSize <- n
    }

    raw <- object@rawData[validGenes(object), , drop = FALSE]
    libSize <- object$total_counts
    times <- ceiling(permSize / n)

    eachSize <- ceiling(permSize / times)

    # For reproducing old implementation, we keep the random number generation
    # in the same order as normally it would happen in inferIntrScore().
    diffGraph <- object@neighborDiff
    contGraph <- object@neighborCont
    dtAvgGraph <- to_mean(contGraph)

    diffPermIdxMat <- sapply(
        X = seq_len(times),
        FUN = function(i) sample(n)
    )
    diff_lig_NGL <- lapply(
        X = seq_len(times),
        FUN = function(i) t(diffGraph[diffPermIdxMat[,i], , drop = FALSE])
    )
    diff_recep_NGL <- lapply(
        X = seq_len(times),
        FUN = function(i) {
            if (recepUseDT) {
                t(contGraph[diffPermIdxMat[,i], , drop = FALSE])
            } else {
                methods::as(Matrix::Diagonal(n), 'CsparseMatrix')[, diffPermIdxMat[,i], drop = FALSE]
            }
        }
    )
    diffSampleIdx <- sample(n, eachSize)
    diff_dtAvg_NGL <- lapply(seq_len(times), function(i) {
        t(dtAvgGraph[diffPermIdxMat[, i], diffSampleIdx, drop = FALSE])
    })
    # Smoothing always use the GauEps graph
    diff_smooth_NGL <- lapply(seq_len(times), function(i) {
        t(diffGraph[sample(eachSize), sample(eachSize), drop = FALSE])
    })

    contPermIdxMat <- sapply(
        X = seq_len(times),
        FUN = function(i) sample(n)
    )
    cont_lig_NGL <- lapply(
        X = seq_len(times),
        FUN = function(i) t(contGraph[contPermIdxMat[,i], , drop = FALSE])
    )
    cont_recep_NGL <- lapply(
        X = seq_len(times),
        FUN = function(i) {
            if (recepUseDT) {
                t(contGraph[contPermIdxMat[,i], , drop = FALSE])
            } else {
                methods::as(Matrix::Diagonal(n), 'CsparseMatrix')[, contPermIdxMat[,i], drop = FALSE]
            }
        }
    )
    contSampleIdx <- sample(n, eachSize)
    cont_dtAvg_NGL <- lapply(seq_len(times), function(i) {
        t(dtAvgGraph[contPermIdxMat[, i], contSampleIdx, drop = FALSE])
    })
    # Smoothing always use the GauEps graph
    cont_smooth_NGL <- lapply(seq_len(times), function(i) {
        t(diffGraph[sample(eachSize), sample(eachSize), drop = FALSE])
    })

    intrType <- as.integer(object@intrDB$type) - 1
    Lmap <- intrGeneMap(object, type = c('diff', 'cont'), component = 'lig')
    Rmap <- intrGeneMap(object, type = c('diff', 'cont'), component = 'rec')


    ecdf <- perm_test_Rcpp(
        raw = t(raw),
        libSize = libSize,
        intrType = intrType,
        diff_lig_NGL = diff_lig_NGL,
        diff_recep_NGL = diff_recep_NGL,
        diff_dtAvg_NGL = diff_dtAvg_NGL,
        diff_smooth_NGL = diff_smooth_NGL,
        cont_lig_NGL = cont_lig_NGL,
        cont_recep_NGL = cont_recep_NGL,
        cont_dtAvg_NGL = cont_dtAvg_NGL,
        cont_smooth_NGL = cont_smooth_NGL,
        Lmap = Lmap,
        Rmap = Rmap,
        lrscore = lrscore
    )
    pval <- as.matrix(1 - ecdf)
    dimnames(pval) <- dimnames(lrscore)
    object@significance$pval <- pval
    return(object)
}



#' Computes LR score with sender-receiver directionality at cluster level
#' @description
#' While \code{\link{inferLRScore}} computes label-free LRscores at single-spot
#' level, this function calculates LRscores for each pair of sender and receiver
#' clusters. This is done by imputing ligand amount (\eqn{L}) only from the
#' neighbor spots in a specific sender cluster. A cluster assignment variable
#' must be preset with \code{\link{setCluster}}.
#' @param object A \code{\linkS4class{cytosignal2}} object with neighbor graphs
#' constructed.
#' @param smoothR Logical, whether to smooth the receptor expression using the
#' contact-dependent neighbor graph. Default is \code{FALSE}.
#' @return A 3D array of dimension (nCluster, nCluster, nInteractors), \eqn{CLR}
#' for cluster-wise LRscore. A value \code{CLR[s, r, i]} in the array means the
#' LRscore for interaction \code{i} when the sender cluster is \code{s} and the
#' receiver cluster is \code{r}.
#' dimension is the sender cluster, the second dimension is the
inferClusterLRScore <- function(
        object,
        smoothR = FALSE
) {
    useCluster <- object@parameters$cluster
    if (is.null(useCluster)) {
        cli::cli_abort("No cluster variable preset. Please set with {.fn setCluster}.")
    }
    clusterVar <-object@metadata[[useCluster]]
    if (is.character(clusterVar) || is.logical(clusterVar)) {
        clusterVar <- factor(clusterVar)
    } else if (is.numeric(clusterVar)) {
        nUniq <- length(unique(clusterVar))
        if (nUniq > 100) {
            cli::cli_abort("Numeric cluster variable has {nUniq} unique values, which does not look right.")
        }
        clusterVar <- factor(clusterVar)
    } else if (is.factor(clusterVar)) {
        # do nothing
    } else {
        cli::cli_abort("Cluster variable must be of atomic type.")
    }
    clusterLevels <- levels(clusterVar)
    clusterInt <- as.integer(clusterVar)
    Lmap <- intrGeneMap(object, component = 'lig')
    Rmap <- intrGeneMap(object, component = 'recep')
    raw <- object@rawData[validGenes(object), , drop = FALSE]
    diffLigGraph <- object@neighborDiff
    contLigGraph <- object@neighborCont
    if (isTRUE(smoothR)) {
        recepGraph <- contGraph
    } else {
        recepGraph <- methods::as(Matrix::Diagonal(ncol(raw)), 'CsparseMatrix')
    }
    dtAvgGraph <- to_mean(object@neighborCont)
    clusterWiseLRscore <- clusterWiseLRscore_cpp(
        raw = raw,
        libSize = object$total_counts,
        clusterInt = clusterInt - 1,
        nCluster = length(clusterLevels),
        intrType = as.integer(object@intrDB$type) - 1,
        diff_lig_graph = diffLigGraph,
        cont_lig_graph = contLigGraph,
        recep_graph = recepGraph,
        dtAvg_graph = dtAvgGraph,
        Lmap = Lmap,
        Rmap = Rmap
    )
    dimnames(clusterWiseLRscore) <- list(
        clusterLevels,
        clusterLevels,
        object@intrDB$interactors
    )
    object@clusterLRScore <- clusterWiseLRscore
    return(object)
}


plotClusterLRHeatmap <- function(
        clr,
        interaction,
        ...
) {
    # TODO
}

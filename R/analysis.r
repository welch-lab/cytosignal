#' Infer the parameters of the Gaussian kernel
#'
#' @param object A Cytosignal object
#' @param scale.facor 1 spatial coord unit equals to how many µm
#' @param r.eps.real The radius of the epsilon ball in tech resolution in um, default 200 um
#' @param thresh The total signal out of the epsilon ball
#'
#' @return A Cytosignal object
#' @export

inferEpsParams <- function(object, scale.factor = NULL, r.eps.real = 200, thresh = 0.001){
  r.eps.scale = r.eps.real / scale.factor # the radius of the epsilon ball in scaled resolution
  sigma.scale = r.eps.scale / sqrt(-2 * log(thresh)) # sigma in tech resolution

  object@parameters[["r.diffuse.scale"]] <- r.eps.scale
  object@parameters[["sigma.scale"]] <- sigma.scale

  return(object)
}

#' Find the neighbors of each cell in the Epsilon Ball
#' @description
#' Find neighbors within a range from each cell, and determine the neighbors'
#' weights using a Gaussian kernel.
#' @param object A CytoSignal object with cell location matrix available.
#' @param eps The radius of the epsilon ball. For CytoSignal object, by default
#' use parameter inferred with \code{\link{inferEpsParams}}.
#' @param sigma The sigma of the Gaussian kernel. By default use parameter
#' inferred with \code{\link{inferEpsParams}}.
#' @param self.weight Weight of the index cell. Use a number between 0-1, or
#' choose from \code{"auto"} or \code{"sum_1"} for automatic self-weighting.
#' @param imp.name Name prefix of the analysis.
#' @return Returns the original object with result updated in "imputation" slot.
#' An "ImpData" object is created with its "nn.graph" slot filled a NxN sparse
#' matrix indicating the NN relationship and weights. A column is considered as
#' a center of the epsilon ball and the rows containing non-zero values are
#' its neighbors.
#' @export
findNNGauEB <- function(
    object,
    eps = NULL,
    sigma = NULL,
    self.weight = "auto",
    imp.name = "GauEps"
) {
  cells.loc <- object@cells.loc

  if (is.null(imp.name)) {imp.name <- "GauEps"}

  if (is.null(eps)) {
    eps <- object@parameters$r.diffuse.scale
  } else{
    imp.name <- paste0(imp.name, "_eps-", eps)
  }

  if (is.null(sigma)) {
    sigma <- object@parameters$sigma.scale
  } else{
    imp.name <- paste0(imp.name, "_sigma-", sigma)
  }

  message("Finding neighbors in epsilon circle with imp.name: ", imp.name, "...")

  if (is.null(sigma)) stop("`sigma` has to be specified for matrix method.")

  dis_mat <- select_EB_rcpp2(cells.loc, eps = eps)
  gauss_vec_inplace_cpp(dis_mat@x, sigma)

  if (is.numeric(self.weight)) {
    if (self.weight < 0 || self.weight > 1) {
      stop("Self weight should be within (0,1)!")
    }
    message("Using manual self weight: ", self.weight, "...")
    # Gaussian kernel for neighbors, normalize neighbors for each center, and hard-set self-weight
    dis_mat <- normCounts.dgCMatrix(dis_mat, scale.fac = Matrix::colSums(dis_mat), method = "none")
    Matrix::diag(dis_mat) <- self.weight
  } else if (self.weight == "auto") {
    message("Determining self weight automatically...")
    # Gaussian kernel with 0-self-dist, multiply self-weight by 5, normalize for each center
    Matrix::diag(dis_mat) <- gauss_vec_cpp(1e-9, sigma) * 5
    dis_mat <- normCounts.dgCMatrix(dis_mat, scale.fac = Matrix::colSums(dis_mat), method = "none")
  } else if (self.weight == "sum_1") {
    message("Using self weight: all NB weights sum to 1...")
    # Gaussian kernel with 0-self-dist, normalize for each center
    Matrix::diag(dis_mat) <- gauss_vec_cpp(1e-9, sigma)
    dis_mat <- normCounts.dgCMatrix(dis_mat, scale.fac = Matrix::colSums(dis_mat), method = "none")
  } else {
    stop("Please set `self.weight` to a number within (0, 1) ",
         "or use \"auto\" or \"sum_1\".")
  }

  colnames(dis_mat) <- rownames(cells.loc)
  nn.obj <- methods::new(
    "ImpData",
    method = "GauEps",
    imp.data = new("dgCMatrix"),
    imp.velo = new("matrix"),
    nn.graph = dis_mat,
    scale.fac = new("numeric"),
    scale.fac.velo = new("numeric"),
    log = list(
      "Parameters" = paste0("eps: ", eps, ", sigma: ", sigma)
    )
  )

  object@imputation[[imp.name]] <- nn.obj
  object@imputation[["default"]] <- imp.name

  return(object)
}

#' Find the direct connected neighbor of each cell, using Delaunay triangulation
#' @description
#' Find the direct connected neighbors of each cell, and assigning an averaged
#' weight to each neighbor.
#' @param object A CytoSignal object or a matrix of cell location.
#' @param max.r The maximum radius of the edges, by default is the r.diffuse.scale,
#' or the value inferred by \code{\link{inferEpsParams}}.
#' @param weight.sum The sum of the weights
#' @param imp.name The name of the analysis.
#' @return Returns the original object with result updated in "imputation" slot.
#' An "ImpData" object is created with its "nn.graph" slot filled a NxN sparse
#' matrix indicating the NN relationship and weights. A column is considered as
#' a center of the epsilon ball and the rows containing non-zero values are
#' its neighbors.
#' @export
findNNDT <- function(
    object,
    max.r = NULL,
    weight.sum = 2,
    imp.name = "DT"
) {
  if (is.null(imp.name)) imp.name <- "DT"

  message("Finding neighbors using DT with imp.name: ", imp.name, "...")

  cells.loc <- object@cells.loc

  # filter out the edges that are over given radius
  if (is.null(max.r)) {
    max.r <- object@parameters$r.diffuse.scale
    if (is.null(max.r)) {
      stop("No default `max.r` specified or inferred.")
    }
  } else if (!is.numeric(max.r)) {
    stop("max.r should be a numeric value.")
  }

  cells.dt <- RTriangle::triangulate(RTriangle::pslg(P = cells.loc))
  edges <- cells.dt$E
  node1.loc <- cells.loc[edges[, 1],]
  node2.loc <- cells.loc[edges[, 2],]
  # Check Euclidean distance for dropping edges over max.r
  dist <- euclidean_elementwise_cpp(node1.loc, node2.loc)
  nn.valid <- dist <= max.r
  edges <- edges[nn.valid, , drop = FALSE]

  # Build neighbor graph temporarily with weights as 1, then normalize by column
  # so that the colSums is the niche size.
  weight.mtx <- Matrix::sparseMatrix(
    i = c(edges[, 1], edges[, 2]),
    j = c(edges[, 2], edges[, 1]),
    x = rep(1, 2*nrow(edges)),
    dims = c(nrow(cells.loc), nrow(cells.loc)),
    dimnames = list(NULL, rownames(cells.loc))
  )

  norm.param <- Matrix::colSums(weight.mtx)
  if (weight.sum == 1) norm.param = norm.param + 1
  # weight.mtx@x <- weight.mtx@x / rep.int(norm.param, diff(weight.mtx@p))
  weight.mtx <- normCounts.dgCMatrix(weight.mtx, scale.fac = norm.param, method = "none")
  Matrix::diag(weight.mtx) <- 1

  nn.obj <- methods::new(
    "ImpData",
    method = "DT",
    imp.data = new("dgCMatrix"),
    # imp.data.null = new("dgCMatrix"),
    # imp.data.null = list(),
    # intr.data = new("dgCMatrix"),
    imp.velo = new("matrix"),
    nn.graph = weight.mtx,
    scale.fac = new("numeric"),
    scale.fac.velo = new("numeric"),
    log = list(
      "Parameters" = "Delauany Triangulation",
      "Num of neighbors" = paste0("Mean: ", mean(norm.param), ", Median: ", median(norm.param))
    )
  )

  object@imputation[[imp.name]] <- nn.obj
  object@imputation[["default"]] <- imp.name

  return(object)
}


#' Create a ImpData object using raw data without imputation
#' @description
#' This function do not attempt to find neighbors but simply creates internal
#' data structure which is helpful for downstream analysis.
#' @param object A CytoSignal object
#' @return Returns the original object with result updated in "imputation" slot.
#' An "ImpData" object is created with its "nn.graph" taken from the result of
#' \code{\link{findNNDT}}.
#' @export
findNNRaw <- function(
    object
) {
  imp.name <- "Raw"

  message("Setting original expression as raw imputation...")

  if ("DT" %in% names(object@imputation)){
    message("DT has been done before, taking the same neighbors.")
    # nn <- new("list")
    # nn$id <- object@imputation[["DT"]]@nn.id
    # nn$dist <- object@imputation[["DT"]]@nn.dist
  } else {
    object <- findNNDT(object)
    # cells.loc <- object@cells.loc
    # nn <- findNNDT.matrix(cells.loc)
  }
  nn.graph <- object@imputation[["DT"]]@nn.graph

  nn.obj <- methods::new(
    "ImpData",
    method = "Raw",
    imp.data = object@counts,
    # imp.data.null = new("dgCMatrix"),
    # imp.norm.data = new("list"),
    imp.velo = new("matrix"),
    # intr.data = new("dgCMatrix"),
    nn.graph = nn.graph,
    # nn.id = nn$id,
    # nn.dist = nn$dist,
    scale.fac = object@parameters[["lib.size"]],
    scale.fac.velo = new("numeric"),
    log = list(
      "Parameters" = "Raw data without imputation"#,
      # "Num of neighbors" = paste0("Mean: ", mean(table(nn$id)), ", Median: ", median(table(nn$id)))
    )
  )

  object@imputation[[imp.name]] <- nn.obj
  object@imputation[["default"]] <- imp.name

  return(object)
}

#' Impute the data
#' @description
#' Impute the quantitative value of ligand received or receptor expressed at
#' each bead basing on the NN weighting previously determined with
#' \code{findNN*} functions.
#' @param object A Cytosignal object with \code{findNN*} performed.
#' @param imp.use The type of neighbors. The \code{imp.name} used when finding
#' neighbors. Default \code{NULL} uses the last NN type computed.
#' @param weights Weighting strategy. Default \code{"none"} use weights computed
#' in \code{findNN*}. Alternatively, choose from \code{"mean"}, \code{"counts"},
#' \code{"dist"}.
#' @return Returns the original object with imputation result updated. For the
#' selected "ImpData" specified by \code{imp.use}, the following slots are
#' updated:
#' \itemize{
#' \item{"imp.data"}{The imputed data.}
#' \item{"nn.graph"}{The neighbor graph, updated if \code{"weights"} is not
#' "none".}
#' \item{"scale.fac"}{The scale factor of the imputed data.}
#' }
#' @export
imputeNiche <- function(
    object,
    imp.use = NULL,
    weights = c("none", "mean", "counts", "dist")
) {
  weights <- match.arg(weights)
  imp.use <- .checkImpUse(object, imp.use)
  if (imp.use == "Raw") {
    warning("Raw imputation data is the original counts and is already available. Returning input object without doing anything.")
    return(object)
  }
  message("Imputing using ", imp.use, "...")

  dge.raw <- object@counts

  weight.mtx <- object@imputation[[imp.use]]@nn.graph

  if (ncol(dge.raw) != nrow(weight.mtx) ||
      ncol(dge.raw) != ncol(weight.mtx)) {
    stop("Number of beads in DGE does not match the number of beads in NN graph.")
  }

  do.norm <- TRUE
  if (weights == "none") {
    do.norm <- FALSE
  } else if (weights == "mean") {
    weight.mtx@x <- rep(1, length(weight.mtx@x))
  } else if (weights == "counts") {
    weight.mtx@x <- Matrix::colSums(dge.raw)[weight.mtx@i]
  } else if (weights == "dist") {
    weight.mtx@x <- exp(-weight.mtx@x)
  }

  # re-weight the weight mtx across each niche (column)
  if (do.norm) {
    weight.mtx <- normCounts.dgCMatrix(weight.mtx, scale.fac = Matrix::colSums(weight.mtx), method = "none")
    # weight.mtx@x = weight.mtx@x / rep.int(Matrix::colSums(weight.mtx), diff(weight.mtx@p))
  }

  dge.raw.imputed <- dge.raw %*% weight.mtx; gc()

  # weighted sum of scale factors as well
  scale.fac <- Matrix::Matrix(object@parameters$lib.size,
                              sparse = TRUE, nrow = 1,
                              dimnames = list(NULL, colnames(dge.raw)))
  scale.fac.imp <- scale.fac %*% weight.mtx

  scale.fac.imp <- as.numeric(scale.fac.imp)
  names(scale.fac.imp) <- names(object@parameters$lib.size)

  res.density <- sum(dge.raw.imputed != 0)/length(dge.raw.imputed)

  object@imputation[[imp.use]]@imp.data <- dge.raw.imputed
  object@imputation[[imp.use]]@nn.graph <- weight.mtx
  object@imputation[[imp.use]]@scale.fac <- scale.fac.imp
  object@imputation[[imp.use]]@log[["Density:"]] <- paste0(res.density*100, "%")

  return(object)
}


#' Normalize the data using the specified method
#' @rdname normCounts
#' @param object A CytoSignal object or a dgCMatrix object of imputation data.
#' @param method The method to use for normalization. \code{"default"}
#' normalizes each column by \code{scale.fac}, multiplies each by 1e4 and
#' finally takes log1p. \code{"cpm"} does the similar thing but multiples each
#' by 1e6. \code{"lib_size"} only normalize by \code{scale.fac}.
#' @param ... Arguments passed to other S3 methods.
#' @return For CytoSignal object, the original object with normalized data
#' updated. For dgCMatrix object, the normalized version of it.
#' @noRd
normCounts <- function(
    object,
    method = c("default", "cpm", "none", "scanpy"),
    ...
) {
  UseMethod(generic = 'normCounts', object = object)
}

#' @rdname normCounts
#' @param scale.fac Numeric vector of length equals to \code{ncol(object)}. For
#' CytoSignal method, this is pre-determined.
#' @method normCounts dgCMatrix
#' @noRd
normCounts.dgCMatrix <- function(
    object,
    scale.fac = NULL,
    method = c("default", "cpm", "none", "scanpy"),
    ...
){
  if (is.null(scale.fac)) {
    stop("scale.fac has to be specified!")
  }
  method <- match.arg(method)
  object@x <- object@x / rep.int(scale.fac, diff(object@p))
  if (method == "default"){
    object@x <- log1p(object@x * 1e4)
  } else if (method == "cpm"){
    object@x <- log1p(object@x * 1e6)
  } else if (method == "scanpy"){
    lib.median = stats::median(scale.fac)
    object@x <- object@x * lib.median
  } else if (method == "none"){
  } else {
    stop("Unknown method!")
  }

  return(object)
}

#' @rdname normCounts
#' @param imp.use The type of neighbors. The \code{imp.name} used when finding
#' nearest neighbors in the upstream step. If not customized, use
#' \code{"GauEps"} for NN found with \code{\link{findNNGauEB}}, or use
#' \code{"DT"} for NN found with \code{\link{findNNDT}}. Default use the most
#' lastly computed NN.
#' @param verbose Whether to print out the progress. Default \code{TRUE}.
#' @param method The method to use for normalization. default: seurat style \code{"default"}
#' @return A dgCMatrix of normalized data.
#' @method normCounts CytoSignal
#' @noRd
normCounts.CytoSignal <- function(
    object,
    method = c("default", "cpm", "none", "scanpy"),
    imp.use = NULL,
    verbose = TRUE,
    ...
){
  method <- match.arg(method)
  imp.use <- .checkImpUse(object, imp.use = imp.use)
  if (isTRUE(verbose))
  message(paste0("Normalizing method: ", method, " on Imputation slot: ", imp.use, "..."))

  mat <- object@imputation[[imp.use]]@imp.data
  scale.fac <- object@imputation[[imp.use]]@scale.fac

  if (method == "none"){
    return(mat)
  }

  mat <- normCounts(object = mat, scale.fac = scale.fac, method = method)

  # modified to output the normalized data instead of the object
  # object@imputation[[slot.use]]@imp.data <- mat
  # return(object)
  return(mat)

  # if (!all.equal(names(slot(object, slot.use)), names(object@imp.norm.data))){
  # 	warning("Names of imp.data and imp.norm.data not equal!")
  # }

}

#' Infer significance of LR scores
#' @description
#' After inferring the LR scores and conducting permutation tests, this function
#' computes p-values for each interaction on each bead and also produces
#' adjusted p-values.
#' @param object A CytoSignal object with \code{\link{inferIntrScore}}
#' performed, or with \code{\link{inferScoreLR}}, \code{\link{permuteLR}} and
#' \code{\link{inferNullScoreLR}} done in sequence.
#' @param lrscore.use Which LR score to use. Use the name specified with
#' \code{lrscore.name} when running \code{\link{inferScoreLR}}. Default
#' \code{NULL} apply specified significance and filtering criteria to all
#' available LR-score inferred.
#' @param fdr.method The false discovery rate method to use. Choose from
#' \code{"spatialFDR"} and \code{"fdr"}. Default \code{"spatialFDR"}.
#' @param p.thresh A numeric scalar. The adjusted p-value threshold to use for
#' filtering significant interactions. Default \code{0.05}.
#' @param reads.thresh A numeric scalar. The minimum number of reads for a
#' ligand-receptor interaction to be considered. Default \code{100}.
#' @param sig.thresh A numeric scalar. The minimum number of of beads for a
#' ligand-receptor interaction to be considered. Default \code{100}.
#' @return The input CytoSignal object with result updated in "lrscore" slot.
#' For each "lrScores" object inside, the "res.list" slot will be updated with:
#' \itemize{
#' \item{"result"}{A list, each element is named by a significant interaction,
#' the value of each element is a character vector of the beads ID where the
#' interaction is significant. This is derived after filtering with
#' \code{p.thresh}}
#' \item{"result.hq"}{A list of the same format as "result", but after further
#' quality control according to \code{reads.thresh} and \code{sig.thresh}}
#' }
#' @export
inferSignif <- function(
    object,
    lrscore.use = NULL,
    fdr.method = c("spatialFDR", "fdr"),
    p.thresh = 0.05,
    reads.thresh = 100,
    sig.thresh = 100
) {
  if (is.null(lrscore.use)) {
    # By default use all
    lrscore.use <- names(object@lrscore)
    lrscore.use <- lrscore.use[lrscore.use != "default"]
  } else {
    lrscore.use <- .checkSlotUse(object, slot.use = lrscore.use)
  }

  fdr.method <- match.arg(fdr.method)

  dge.raw <- object@counts
  gene_to_uniprot <- object@intr.valid[['gene_to_uniprot']]
  for (lr in lrscore.use) {
    imp.use <- object@lrscore[[lr]]@recep.slot
    imp.use <- .checkImpUse(object, imp.use = imp.use)

    message(sprintf(
      "Inferring significant beads on LR-score '%s', using imputation '%s'... ",
      lr, imp.use
    ))

    nn.graph <- object@imputation[[imp.use]]@nn.graph

    intr.db <- object@intr.valid[[object@lrscore[[lr]]@intr.slot]]

    lrscore.mtx <- object@lrscore[[lr]]@score
    null.lrscore.mtx <- object@lrscore[[lr]]@score.null

    pval.mtx <- getPvalues(lrscore.mtx, null.lrscore.mtx)

    if (fdr.method == "spatialFDR"){
      pval.spatial <- graphSpatialFDRNew(nn.graph, pval.mtx)
      dimnames(pval.spatial) <- dimnames(lrscore.mtx); gc()
    } else if (fdr.method == "fdr"){
      pval.spatial <- sapply(1:nrow(pval.mtx), function(i) {
        p.adjust(pval.mtx[i, ], method = "BH")
      })
      pval.spatial <- t(pval.spatial)
      dimnames(pval.spatial) <- dimnames(lrscore.mtx); gc()
    }

    res.list <- lapply(colnames(pval.spatial), function(cp){
      rownames(pval.spatial)[pval.spatial[, cp] < p.thresh]
    })

    names(res.list) <- colnames(pval.spatial)

    res.list <- res.list[lengths(res.list) != 0]
    message("- Number of interactions that have significant i-niche: ",
            length(res.list))

    res.list <- res.list[order(lengths(res.list), decreasing = TRUE)]

    res.list.hq = filterRes(dge.raw = dge.raw, res.list = res.list,
                            intr.db = intr.db, gene_to_uniprot = gene_to_uniprot,
                            reads.thresh = reads.thresh, sig.thresh = sig.thresh)

    res.list <- list(
      result = res.list,
      result.hq = res.list.hq
    )

    object@lrscore[[lr]]@res.list <- res.list

    object@lrscore[[lr]]@log[["Parameters for Significant beads"]] <- list(
      "p.value" = p.thresh,
      "reads.thresh" = reads.thresh,
      "sig.thresh" = sig.thresh
    )

    object@lrscore[[lr]]@log[["Significant Intrs"]] <- list(
      paste0("Number of interactions that have significant i-niche: ", length(res.list[["result"]])),
      paste0("Number of high quality intr: ", length(res.list[["result.hq"]]))
    )
  }

  return(object)
}


#' Identify spatially significant interactions using std-corrected pearson correlation

#' Normal Moran's I test is not applicable here since the total number of the cell is too large, causing
#' unnacceptable computation cost. Here we use a modified version of Moran's I test, which is to take only
#' the top KNNs to compute the Moran's I test.
#'
#' @param object A Cytosignal object
#' @param k The number of nearest neighbors to use
#' @param weight The weight of the nearest neighbors
#' @param score.slot Which LR score to use. Use the name specified with \code{tag}
#' when running \code{\link{inferLRScore}}.
#'
#' @return A Cytosignal object
#' @export
runPears.std <- function(
    object,
    k = 10,
    weight = 2,
    score.slot = NULL
) {
  if (is.null(score.slot)){
    score.slot = object@lrscore[["default"]]
  }

  if (!score.slot %in% names(object@lrscore)) {
    stop("LRscore obj not found. ")
  }

  message("Inferring spatial significant intrs on Score slot ", score.slot, "... ")

  # if (length(method) > 1 | !method %in% c("pear.var", "moranI")) {
  # 	stop("Method not supported. ")
  # }

  intr.hq <- names(object@lrscore[[score.slot]]@res.list[["result.hq"]])

  lr.mtx <- as.matrix(object@lrscore[[score.slot]]@score[, intr.hq])
  lt.mtx.imp <- dataImpKNN(lr.mtx, object@cells.loc, k = k, weight = weight)
  lt.mtx.imp <- as.matrix(Matrix::t(lt.mtx.imp))

  intr.order <- as.double(pearson_col_cpp(lt.mtx.imp, lr.mtx)) * as.double(stdMat_cpp(lt.mtx.imp)) / as.double(stdMat_cpp(lr.mtx))
  names(intr.order) <- colnames(lr.mtx)
  intr.order <- intr.order[order(intr.order, decreasing = T)]

  message("Reordering significant interactions...")

  res.list.pear <- object@lrscore[[score.slot]]@res.list[["result.hq"]][names(intr.order)]
  object@lrscore[[score.slot]]@res.list[["result.hq.pear"]] <- res.list.pear

  return(object)
}


#' Infer the correspondence between LR-scores and Significance
#'
#' Use DT neighbors to impute the p-value and LR-score for each cell, then compute the pearson
#' correlation between the imputed LR-score and p-value.
#'
#' @param object A Cytosignal object
#' @param correctBy Spatial FDR correction method. Choose from \code{"cell"} or
#' \code{"intr"}.
#' @param slot.use Which LR score to use. Use the name specified with \code{tag}
#' when running \code{\link{inferLRScore}}.
#' @return A Cytosignal object
#'
inferSpatialCorr <- function(
    object,
    correctBy = c("cell", "intr"),
    slot.use = NULL
) {
  if (is.null(slot.use)){
    signif.use = object@lrscore[["default"]]
  }

  lrscore.mtx = object@lrscore[[slot.use]]@score
  null.lrscore.mtx = object@lrscore[[slot.use]]@score.null
  nb.index = object@imputation[["DT"]]@nn.id

  # if Cell numbers in pvalue mtx and nb.fac do not match
  cellsNN = as.integer(levels(nb.index))
  if (!identical(nrow(lrscore.mtx), length(cellsNN))) {
    cellsNoNN <- setdiff(1:nrow(lrscore.mtx), cellsNN)
    lrscore.mtx <- lrscore.mtx[-cellsNoNN, ]
  }

  pval.mtx = getPvalues(lrscore.mtx, null.lrscore.mtx)

  if (correctBy == "cell"){
    pval.spatial = graphSpatialFDR(nb.fac, pval.mtx)
    dimnames(pval.spatial) = dimnames(lrscore.mtx); gc()
  } else if (correctBy == "intr"){
    pval.spatial = sapply(1:nrow(pval.mtx), function(i) {
      p.adjust(pval.mtx[i, ], method = "BH")
    })
    pval.spatial = t(pval.spatial)
    dimnames(pval.spatial) = dimnames(lrscore.mtx); gc()
  }

  pval.spatial[pval.spatial > 0.05] = 0
  pval.spatial[pval.spatial > 0] = 1

  nn.graph = object@imputation[["DT"]]@nn.graph
  nn.graph = nn.graph[-cellsNoNN, -cellsNoNN]
  nn.graph = to_mean(nn.graph)

  # impute pval and lrscore
  # pval.spatial = as(t(pval.spatial), "CsparseMatrix")
  # pval.spatial.imp = t(pval.spatial) %*% nn.graph
  # pval.spatial.imp = as.matrix(t(pval.spatial.imp))

  # lrscore.mtx.imp = t(lrscore.mtx) %*% nn.graph
  # lrscore.mtx.imp = as.matrix(t(lrscore.mtx.imp))

  pval.spatial.imp = pval.spatial
  lrscore.mtx.imp = as(lrscore.mtx, "matrix")

  # compute pearson correlation
  pearson.cor = as.double(pearson_col_cpp(pval.spatial.imp, lrscore.mtx.imp))
  names(pearson.cor) = colnames(lrscore.mtx)
  pearson.cor = pearson.cor[!is.na(pearson.cor)]
  pearson.cor = pearson.cor^2
  pearson.cor = pearson.cor[order(pearson.cor, decreasing = T)]

  res.list.corr = object@lrscore[[slot.use]]@res.list[["result.hq"]]
  intr.use = intersect(names(pearson.cor), names(res.list.corr))
  pearson.corr.use = pearson.cor[intr.use]
  pearson.corr.use = pearson.corr.use[order(pearson.corr.use, decreasing = T)]
  res.list.corr = res.list.corr[names(pearson.corr.use)]

  object@lrscore[[slot.use]]@res.list[["result.hq.corr"]] = res.list.corr

  return(object)
}



#' Infer the correspondence between LR-scores and Significance
#'
#' Min-max the LRscores, substract the significant ones, sum and average.
#'
#' @param object A Cytosignal object
#' @param correctBy Spatial FDR correction method. Choose from \code{"cell"} or
#' \code{"intr"}.
#' @param slot.use Which LR score to use. Use the name specified with \code{tag}
#' when running \code{\link{inferLRScore}}.
#' @return A Cytosignal object
#' @export
inferCorrScore <- function(
    object,
    correctBy = c("cell", "intr"),
    slot.use = NULL
) {
  correctBy <- match.arg(correctBy)
  if (is.null(slot.use)){
    slot.use = object@lrscore[["default"]]
  }

  lrscore.mtx = object@lrscore[[slot.use]]@score
  null.lrscore.mtx = object@lrscore[[slot.use]]@score.null
  nb.index = object@imputation[["DT"]]@nn.id

  # if Cell numbers in pvalue mtx and nb.fac do not match
  cellsNN = as.integer(levels(nb.index))
  if (!identical(nrow(lrscore.mtx), length(cellsNN))) {
    cellsNoNN <- setdiff(1:nrow(lrscore.mtx), cellsNN)
    lrscore.mtx <- lrscore.mtx[-cellsNoNN, ]
  }

  pval.mtx = getPvalues(lrscore.mtx, null.lrscore.mtx)

  if (correctBy == "cell"){
    nb.fac = list(
      id = object@imputation[["DT"]]@nn.id,
      dist = object@imputation[["DT"]]@nn.dist
    )
    pval.spatial = graphSpatialFDR(nb.fac, pval.mtx)
    dimnames(pval.spatial) = dimnames(lrscore.mtx); gc()
  } else if (correctBy == "intr"){
    pval.spatial = sapply(1:nrow(pval.mtx), function(i) {
      p.adjust(pval.mtx[i, ], method = "BH")
    })
    pval.spatial = t(pval.spatial)
    dimnames(pval.spatial) = dimnames(lrscore.mtx); gc()
  }

  pval.sig <- pval.spatial
  pval.sig[pval.sig >= 0.05] <- 1
  pval.sig[pval.sig < 0.05] <- 0

  # norm lrscore
  l.m.norm <- lrscore.mtx
  l.m.norm@x <- l.m.norm@x / rep.int(Matrix::colSums(l.m.norm), diff(l.m.norm@p))
  l.m.norm <- l.m.norm * pval.sig

  # sum and average
  l.m.norm@x <- l.m.norm@x / nrow(l.m.norm)
  corr.s <- Matrix::colSums(l.m.norm)

  res.list.corr = object@lrscore[[slot.use]]@res.list[["result.hq"]]
  intr.use = intersect(names(corr.s), names(res.list.corr))
  corr.s.hq = corr.s[intr.use]
  corr.s.hq = corr.s.hq[order(corr.s.hq)]
  res.list.corr = res.list.corr[names(corr.s.hq)]

  object@lrscore[[slot.use]]@res.list[["result.hq.corr"]] <- res.list.corr
  object@lrscore[[slot.use]]@res.list[["hq.corr.score"]] <- corr.s.hq

  return(object)
}


#' Rank the inferred high-quality interactions by their spatial variability
#' @description
#' This function utilizes SPARK package to calculate the spatial variability of
#' the high-quality interactions, using their LR scores. Please refer to
#' \href{https://doi.org/10.1186/s13059-021-02404-0}{Jiaqiang Zhu, et al., 2021, Genome Biology}
#' for more details of the method.
#' @param object A \linkS4class{CytoSignal} object, with
#' \code{\link{inferSignif}} already run.
#' @param slot.use Which LR score to use. Use the name specified with \code{tag}
#' when running \code{\link{inferLRScore}}.
#' @param numCores SPARK::sparkx parameter, an integer specifying the number of
#' threads.
#' @param verbose SPARK::sparkx parameter, a logical value indicating whether to
#' print details for debug purpose
#' @return The input \linkS4class{CytoSignal} object with the spatially variable
#' high-quality interaction list updated at
#' \code{object@lrscore[[slot.use]]@res.list$result.spx}
#' @export
#' @examples
#' \dontrun{
#' object <- findNN(object)
#' object <- imputeLR(object)
#' object <- inferScoreLR(object, lig.slot = "GauEps", recep.slot = "Raw",
#'                        intr.db.name = "diff_dep")
#' object <- permuteLR(object)
#' object <- inferNullScoreLR(object)
#' object <- inferSignif(object)
#' object <- rankIntrSpatialVar(object)
#' }
rankIntrSpatialVar <- function(
    object,
    slot.use = NULL,
    numCores = 1,
    verbose = FALSE
) {
  if (!requireNamespace("SPARK", quietly = TRUE)) {
    stop("Package SPARK is required for calculating the spatial variability. ",
         "Please install with:\ndevtools::install_github('xzhoulab/SPARK')")
  }
  if (is.null(slot.use)) {
    slots <- names(object@lrscore)
    slots <- slots[slots != "default"]
  } else {
    slots <- .checkSlotUse(object, slot.use = slot.use)
  }
  for (slot.use in slots) {
    if (!"result.hq" %in% names(object@lrscore[[slot.use]]@res.list)) {
      warning("Filtered high-quality interactions for ", slot.use, "not found. ",
              "Please run `inferSignif()` first.", immediate. = TRUE)
    }
    message("Ranking high-quality interactions by spatial variability using LR score at ", slot.use, "... ")

    spx.res <- SPARK::sparkx(count_in = t(object@lrscore[[slot.use]]@score),
                             locus_in = object@cells.loc,
                             numCores = numCores,
                             option = "mixture",
                             verbose = verbose)
    # head(spx.res$res_mtest)

    spx.pvals <- spx.res$res_mtest
    spx.pvals.hq <- spx.pvals[names(object@lrscore[[slot.use]]@res.list[["result.hq"]]), ]
    nVar <- sum(spx.pvals.hq$adjustedPval <= 0.05)
    message("- Number of spatially variable interactions: ", nVar, " out of ",
            nrow(spx.pvals.hq), " high-quality interactions.")
    spx.pvals.hq <- spx.pvals.hq[order(spx.pvals.hq$adjustedPval), ]

    object@lrscore[[slot.use]]@res.list[["result.spx"]] <- object@lrscore[[slot.use]]@res.list[["result.hq"]][rownames(spx.pvals.hq)]

  }
  return(object)
}



# #' Compute the LR velo for specific ligand-receptor imputation obj pairs
# #'
# #' @param object A Cytosignal object
# #' @param lig.slot The ligand slot to use
# #' @param recep.slot The receptor slot to use
# #' @param intr.db.name The intr database name to use
# #' @param nn.use The neighbor index as niche
# #'
# #' @return A Cytosignal object
# #' @export
# #'
# inferVeloLR <- function(
    #   object,
#   ...
# ) {
#   UseMethod(generic = 'inferVeloLR', object = object)
# }

# #' Sub function for inferVeloLR, input a sparse matrix
# #'
# #' @param dge.lig A sparse matrix for ligand
# #' @param dge.recep A sparse matrix for receptor
# #' @param nb.id.fac A factor of neighbor indices
# #' @param lig.fac A factor of ligand indices
# #' @param recep.fac A factor of receptor indices
# #'
# #' @return A sparse matrix
# #' @export
# #'
# inferVeloLR.matrix_like <- function(
    # 	dge.lig,
#     dge.recep,
# 	dge.lig.velo,
# 	dge.recep.velo,
# 	lig.fac,
# 	recep.fac
# ){

# 	### cavaet: remember to convert the Uniprot ids to indices!
# 	# convert nb fac
# 	# nb.id.fac = sort(nb.id.fac)
# 	# nb.index = facToIndex(nb.id.fac)

# 	if (max(as.integer(names(lig.fac))) > nrow(dge.lig)){
# 		stop("Intr index out of dge bounds.")
# 	}

# 	lig.index = facToIndex(lig.fac)
# 	recep.index = facToIndex(recep.fac)

# 	# compute velos
# 	res.mtx = inferVeloLR_cpp(
# 		unname(as.matrix(dge.lig)),
# 		unname(as.matrix(dge.recep)),
# 		unname(as.matrix(dge.lig.velo)),
# 		unname(as.matrix(dge.recep.velo)),
# 		lig.index[[1]], lig.index[[2]],
# 		recep.index[[1]], recep.index[[2]]
# 	)

# 	dimnames(res.mtx) = list(colnames(dge.lig), levels(lig.fac))
# 	# dimnames(res.mtx) = list(colnames(dge.lig), levels(lig.fac))
# 	# res.mtx = Matrix(res.mtx, sparse = T)

# 	return(res.mtx)
# }



# #' Sub function for inferVeloLR, input a CytoSignal object
# #'
# #' @param object A Cytosignal object
# #' @param lig.slot The ligand slot to use
# #' @param recep.slot The receptor slot to use
# #' @param intr.db.name The intr database name to use
# #' @param nn.use slot that the neighbor index should be taken from, by default is the same as
# #' 			the recep.slot. For example, if velo.obj = GauEps-DT, then nn.use = "DT".
# #' 			nn.use could also be a user-defind factor.
# #'
# #' @return A Cytosignal object
# #' @export
# inferVeloLR.CytoSignal <- function(
    # 	object,
# 	lig.slot,
# 	recep.slot,
# 	intr.db.name,
# 	tag = NULL
# ){
# 	if (!lig.slot %in% names(object@imputation)){
# 		stop("Ligand slot not found.")
# 	}

# 	if (!recep.slot %in% names(object@imputation)){
# 		stop("Receptor slot not found.")
# 	}

# 	message("Computing velos using ", intr.db.name, " database.")
# 	message("Ligand: ", lig.slot, ", Receptor: ", recep.slot, ".")

# 	if (!intr.db.name %in% c("diff_dep", "cont_dep")) {
# 		stop("intr.db.name must be either 'diff_dep' or 'cont_dep'.")
# 	}

# 	if (is.null(tag)) {
# 		tag <- paste0(lig.slot, "-", recep.slot)
# 	}

# 	if (tag %in% names(object@lrvelo)) {
# 		stop("Tag already exists.")
# 	}

# 	object@lrvelo[["default"]] <- tag

# 	dge.lig <- object@imputation[[lig.slot]]@imp.data
# 	dge.recep <- object@imputation[[recep.slot]]@imp.data
# 	dge.lig.velo <- object@imputation[[lig.slot]]@imp.velo
# 	dge.recep.velo <- object@imputation[[recep.slot]]@imp.velo

# 	# compare the dimnames of all four matrices
# 	if (!all.equal(dimnames(dge.lig), dimnames(dge.recep))){
# 		stop("dge.lig and dge.recep must have the same dimension names.")
# 	}

# 	if (!all.equal(dimnames(dge.lig.velo), dimnames(dge.recep.velo))){
# 		stop("dge.lig.velo and dge.recep.velo must have the same dimension names.")
# 	}

# 	use.genes <- intersect(rownames(dge.lig), rownames(dge.lig.velo))

# 	dge.lig <- dge.lig[use.genes, ]
# 	dge.recep <- dge.recep[use.genes, ]
# 	dge.lig.velo <- dge.lig.velo[use.genes, ]
# 	dge.recep.velo <- dge.recep.velo[use.genes, ]

# 	message("Number of velo genes: ", length(use.genes), " / ", nrow(object@imputation[[lig.slot]]@imp.data))

# 	intr.db.list <- checkIntr(use.genes, object@intr.valid[[intr.db.name]])

# 	res.mtx <- inferVeloLR.matrix_like(dge.lig, dge.recep,
# 				dge.lig.velo, dge.recep.velo,
# 				intr.db.list[["ligands"]], intr.db.list[["receptors"]])

# 	intr.hq <- names(object@lrscore[[tag]]@res.list[["result.hq.pear"]])
# 	intr.hq <- intr.hq[intr.hq %in% colnames(res.mtx)]

# 	lrvelo.obj <- new(
# 		"lrVelo",
# 		lig.slot = lig.slot,
# 		recep.slot = recep.slot,
# 		intr.slot = intr.db.name,
# 		intr.list = intr.db.list,
# 		intr.velo = res.mtx,
# 		velo.gene = list(
# 			"velo.genes" = use.genes,
# 			"velo.intr" = colnames(res.mtx),
# 			"velo.intr.hq" = intr.hq),
# 		# nn.id = velo.nn.list[["id"]],
# 		# nn.dist = velo.nn.list[["dist"]],
# 		log = list(
# 			"Used slot" = c(lig.slot, recep.slot)
# 		)
# 	)

# 	object@lrvelo[[tag]] <- lrvelo.obj

# 	return(object)
# }


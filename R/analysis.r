#' Infer the parameters of the Gaussian kernel
#'
#' @param object A Cytosignal object
#' @param scale.facor 1 µm equals to how many spatial coord units
#' @param r.eps.real The radius of the epsilon ball in tech resolution in um, default 200 um
#' @param thresh The total signal out of the epsilon ball
#'
#' @return A Cytosignal object
#' @export

inferEpsParams <- function(object, scale.factor = NULL, r.eps.real = 200, thresh = 0.001){
  r.eps.scale = r.eps.real/scale.factor # the radius of the epsilon ball in scaled resolution
  sigma.scale = r.eps.scale / sqrt(-2 * log(thresh)) # sigma in tech resolution

  object@parameters[["r.diffuse.scale"]] <- r.eps.scale
  object@parameters[["sigma.scale"]] <- sigma.scale

  return(object)
}


#' Find the neighbors of each cell in the Epsilon Ball
#' @rdname findNNGauEB
#' @param object A CytoSignal object or a matrix object for the cell spatial
#' location.
#' @param eps The radius of the epsilon ball. For CytoSignal object, by default
#' use parameter inferred with \code{\link{inferEpsParams}}.
#' @param sigma The sigma of the Gaussian kernel. Default \code{0.15} for matrix
#' method. For CytoSignal object, by default use parameter inferred with
#' \code{\link{inferEpsParams}}.
#' @param self.weight weight of the index cell. Use a number between 0-1, or
#' choose from \code{"auto"} or \code{"sum_1"}.
#' @param ... Arguments passed to other S3 methods.
#' @return For CytoSignal object, the orignal object with result updated. For
#' matrix object, a list of neighbors and their distances.
#' @export
findNNGauEB <- function(
    object,
    eps = NULL,
    sigma = NULL,
    self.weight = "auto",
    ...
) {
  UseMethod(generic = 'findNNGauEB', object = object)
}

#' @rdname findNNGauEB
#' @export
findNNGauEB.matrix <- function(
    object,
    eps,
    sigma = 0.15,
    self.weight = "auto",
    ...
){
  if (is.null(sigma)) stop("`sigma` has to be specified for matrix method.")
  # cat("Finding neighbors in epsilon circle...\n")
  nn <- dbscan::frNN(object, eps = eps, sort = F)

  if (is.numeric(self.weight)) {
    if (self.weight < 0 || self.weight > 1) {
      stop("Self weight should be within (0,1)!")
    }
    message("Using manual self weight: ", self.weight, "...")
  } else if (self.weight == "auto") {
    message("Determining self weight automatically...")
  } else if (self.weight == "sum_1") {
    message("Using self weight: all NB weights sum to 1...")
  } else {
    stop("Please set `self.weight` to a number within (0, 1) ",
         "or use \"auto\" or \"sum_1\".")
  }

  # get a sorted nn factor
  nn.fac <- factor(rep(seq_along(nn$id), lengths(nn$id)))
  names(nn.fac) = unlist(nn$id)
  # add index cell itself
  nn.fac <- addIndex(nn.fac)

  # get a sorted nn dist factor
  nn.dist.fac <- factor(rep(seq_along(nn$dist), lengths(nn$dist)))
  names(nn.dist.fac) = as.numeric(unlist(nn$dist))
  nn.dist.fac = addIndexOne(nn.dist.fac)

  # dist.list.gau = lapply(nn$dist, function(x){
  #     y = gauss_vec_cpp(c(x, 1e-8), sigma)[,1]
  #     return(y/sum(y))
  # })

  dist.list.gau = lapply(nn$dist, function(x){
    if (self.weight == "auto"){
      # times the index cell weight by 10
      y = gauss_vec_cpp(c(x, 1e-9), sigma)[,1]
      y[length(y)] <- y[length(y)] * 5
      return(y/sum(y))
    } else if (is.numeric(self.weight)) {
      # norm the sum except the index cell to 1
      y = gauss_vec_cpp(x, sigma)[,1]
      return(c(y/sum(y), self.weight))
    } else if (self.weight == "sum_1") {
      # norm the sum to 1
      y = gauss_vec_cpp(c(x, 1e-9), sigma)[,1]
      return(y/sum(y))
    }
  })

  # discard the cell that has no neighbors
  rm.index = which(lengths(dist.list.gau) == 1)
  if (length(rm.index) > 0) dist.list.gau = dist.list.gau[-rm.index]
  dist.list.gau = as.numeric(unlist(dist.list.gau))
  names(nn.dist.fac) = dist.list.gau

  # nn.dist.fac <- addIndexOne(nn.dist.fac)

  num.diff = nrow(object) - length(levels(nn.fac))
  if (num.diff > 0){
    message("A total of ", num.diff, " beads do not have NNs.")
  } else if (num.diff < 0) {
    stop("Result fac longer than original beads length!\n")
  }

  # message("Mean num of neighbors: ", ceiling(mean(table(nn.fac))))
  # message("Median num of neighbors: ", median(table(nn.fac)))

  return(list(
    "id" = nn.fac,
    "dist" = nn.dist.fac
  ))
}


#' @rdname findNNGauEB
#' @export
#' @param tag Name prefix of the analysis.
findNNGauEB.CytoSignal <- function(
    object,
    eps = NULL,
    sigma = NULL,
    self.weight = "auto",
    tag = NULL,
    ...
){
  cells.loc <- object@cells.loc

  if (is.null(tag)){tag <- "GauEps"}

  if (is.null(eps)){
    eps <- object@parameters$r.diffuse.scale
  } else{
    tag <- paste0(tag, "_eps-", eps)
  }

  if (is.null(sigma)){
    sigma <- object@parameters$sigma.scale
  } else{
    tag <- paste0(tag, "_sigma-", sigma)
  }

  message("Finding neighbors in epsilon circle with tag: ", tag, "...")

  nn <- findNNGauEB.matrix(cells.loc, eps, sigma, self.weight)

  nn.obj <- methods::new(
    "ImpData",
    method = tag,
    imp.data = new("dgCMatrix"),
    # imp.data.null = new("dgCMatrix"),
    # imp.data.null = list(),
    # intr.data = new("dgCMatrix"),
    imp.velo = new("matrix"),
    nn.graph = new("dgCMatrix"),
    nn.id = nn$id,
    nn.dist = nn$dist,
    scale.fac = new("numeric"),
    scale.fac.velo = new("numeric"),
    log = list(
      "Parameters" = paste0("eps: ", eps, ", sigma: ", sigma),
      "Num of neighbors" = paste0("Mean: ", mean(table(nn$id)), ", Median: ", median(table(nn$id)))
    )
  )

  object@imputation[[tag]] <- nn.obj
  object@imputation[["default"]] <- tag

  return(object)
}


#' Find the direct connected neighbor of each cell, using Delaunay triangulation
#' @rdname findNNDT
#' @param object A CytoSignal object or a matrix of cell location.
#' @param weight.sum The sum of the weights
#' @param max.r The maximum radius of the edges, by default is the r.diffuse.scale,
#' @param ... Arguments passed to other S3 methods
#' @return For CytoSignal object, the original object with result updated. For matrix
#' object, a list of neighbors.
#' @export
findNNDT <- function(
    object,
    ...
) {
  UseMethod(generic = 'findNNDT', object = object)
}

#' @rdname findNNDT
#' @export
#' @method findNNDT matrix
findNNDT.matrix <- function(
    object,
    weight.sum = 2,
    max.r = NULL,
    ...
){
  # cells.loc <- RTriangle::pslg(P=cells.loc)
  cells.dt <- RTriangle::triangulate(RTriangle::pslg(P=object))
  egdes <- cells.dt$E

  #### construct a double-corresponded factor
  nn.fac <- factor(c(egdes[, 1], egdes[, 2])) # neighbors.factor
  names(nn.fac) <- c(egdes[, 2], egdes[, 1])
  nn.fac <- sort(nn.fac)

  index.loc = object[as.integer(nn.fac),]
  nb.loc = object[as.integer(names(nn.fac)),]

  dist.euc.list <- as.double(euclidean_cpp(index.loc, nb.loc))
  nn.valid <- which(dist.euc.list <= max.r)
  nn.fac.old <- nn.fac
  nn.fac <- nn.fac[nn.valid, drop = T]

  # message("Filtering out ", length(dist.euc.list) - length(nn.valid), " edges out of range.")
  noNN <- as.numeric(setdiff(levels(nn.fac.old), levels(nn.fac)))

  if (length(noNN) > 0){
    message("A total of ", length(noNN), " beads do not have NNs!")
  }# else {
  #   message("All beads have NNs.")
  # }

  # set distance for each niche, based on the cell number of each niche
  nb.size <- table(nn.fac)
  dist.list <- lapply(levels(nn.fac), function(x){
    use.size <- unname(nb.size[x])
    if (weight.sum == 1){
      # norm the sum to 1
      y = rep(1/(use.size+1), use.size)
      return(y)
    } else if (weight.sum == 2) {
      # norm the sum except the index cell to 1
      y = rep(1/use.size, use.size)
      return(c(y, 1))
    }
  })
  dist.list <- unlist(dist.list)

  # add index bead it self to the neighbor list
  nn.fac <- addIndex(nn.fac)

  nn.dist.fac <- nn.fac
  names(nn.dist.fac) <- dist.list

  # cat("Mean num of neighbors: ", ceiling(mean(nb.size)), "\n")
  # cat("Median num of neighbors: ", median(nb.size), "\n")

  # return(nn.fac)
  return(list(
    "id" = nn.fac,
    "dist" = nn.dist.fac
  ))
}


#' @rdname findNNDT
#' @export
#' @method findNNDT CytoSignal
findNNDT.CytoSignal <- function(
    object,
    weight = 2,
    max.r = NULL,
    ...
) {
  tag <- "DT"

  # if (tag %in% names(object@imputation)){
  # 	stop("This imputation has been done before.")
  # }

  message("Finding neighbors using DT with tag: ", tag, "...")

  cells.loc <- object@cells.loc

  # filter out the edges that are over given radius
  if (is.null(max.r)) {
    max.r <- object@parameters$r.diffuse.scale
  }

  if (!is.numeric(max.r)){
    stop("max.r should be a numeric value.")
  }

  # message("Filtering out the edges over given radius: ", max.r, "...")

  nn <- findNNDT.matrix(cells.loc, weight.sum = weight, max.r = max.r)

  nn.obj <- methods::new(
    "ImpData",
    method = tag,
    imp.data = new("dgCMatrix"),
    # imp.data.null = new("dgCMatrix"),
    # imp.data.null = list(),
    # intr.data = new("dgCMatrix"),
    imp.velo = new("matrix"),
    nn.graph = new("dgCMatrix"),
    nn.id = nn$id,
    nn.dist = nn$dist,
    scale.fac = new("numeric"),
    scale.fac.velo = new("numeric"),
    log = list(
      "Parameters" = "Delauany Triangulation",
      "Num of neighbors" = paste0("Mean: ", mean(table(nn$id)), ", Median: ", median(table(nn$id)))
    )
  )

  object@imputation[[tag]] <- nn.obj
  object@imputation[["default"]] <- tag

  return(object)
}


#' Create a ImpData object using raw data without imputation
#'
#' @param object A Cytosignal object
#'
#' @return A Cytosignal object
#' @export
findNNRaw <- function(
    object
) {
  tag <- "Raw"

  message("Setting ImpData obj using NO imputation...")

  if ("DT" %in% names(object@imputation)){
    message("DT has been done before, taking the same neighbors.")
    nn <- new("list")
    nn$id <- object@imputation[["DT"]]@nn.id
    nn$dist <- object@imputation[["DT"]]@nn.dist
  } else {
    cells.loc <- object@cells.loc
    nn <- findNNDT.matrix(cells.loc)
  }

  nn.obj <- methods::new(
    "ImpData",
    method = tag,
    imp.data = object@counts,
    # imp.data.null = new("dgCMatrix"),
    # imp.norm.data = new("list"),
    imp.velo = new("matrix"),
    # intr.data = new("dgCMatrix"),
    nn.graph = Matrix::Diagonal(ncol(object@counts)),
    nn.id = nn$id,
    nn.dist = nn$dist,
    scale.fac = object@parameters[["lib.size"]],
    scale.fac.velo = new("numeric"),
    log = list(
      "Parameters" = "Raw data without imputation",
      "Num of neighbors" = paste0("Mean: ", mean(table(nn$id)), ", Median: ", median(table(nn$id)))
    )
  )

  object@imputation[[tag]] <- nn.obj
  object@imputation[["default"]] <- tag

  return(object)
}


#' Impute the data
#' @rdname imputeNiche
#' @param object A Cytosignal object or a dgCMatrix object of raw gene
#' expression matrix
#' @param weights The weight of the Delaunay triangulation. Choose from
#' \code{"mean"}, \code{"counts"}, \code{"dist"} or \code{"none"}.
#' @param ... Arguments passed to other S3 methods
#' @return For CytoSignal object, the original object with imputation result
#' updated. For dgCMatrix object, a sparse N x N graph presented as another
#' dgCMatrix object.
#' @export
imputeNiche <- function(
    object,
    weights = c("mean", "counts", "dist", "none"),
    ...
) {
  UseMethod(generic = 'imputeNiche', object = object)
}

#' @rdname imputeNiche
#' @param nb.id.fac A factor of neighbors
#' @param nb.dist.fac A factor of weights
#' @export
imputeNiche.dgCMatrix <- function(
    object,
    nb.id.fac,
    nb.dist.fac,
    weights = c("mean", "counts", "dist", "none"),
    ...
){
  weights <- match.arg(weights)
  i.list = as.integer(names(nb.id.fac)) # x coords, row number
  j.list = as.numeric(as.character(nb.id.fac)) # y coords, col number, hq beads

  if (weights == "avg") { # just do the mean
    x.list = rep(1, length(i.list))
    do.norm <- TRUE
  } else if (weights == "counts") { # use total counts in each bead as weights
    x.list = colSums(object)[i.list]
    do.norm <- TRUE
  } else if (weights == "dist") {
    x.list = exp(-as.numeric(names(nb.dist.fac))) # exp(-dist) as similarities
    do.norm <- TRUE
  } else if (weights == "none") {
    x.list = as.numeric(names(nb.dist.fac)) # gaussian distances has been calculated
    do.norm <- FALSE
  }

  weights.mtx = sparseMatrix(i = i.list, j = j.list, x = x.list,
                             dims = c(ncol(object), ncol(object)), dimnames = list(NULL, colnames(object)))

  # re-weight the weight mtx across each niche (column)
  if (do.norm) {
    weights.mtx@x = weights.mtx@x / rep.int(Matrix::colSums(weights.mtx), diff(weights.mtx@p))
  }

  dge.raw.imputed = object %*% weights.mtx

  # if (return.graph){
  # 	return(list(dge.raw.imputed, weights.mtx))
  # }

  # return(object.imputed)

  return(weights.mtx)
}

#' @rdname imputeNiche
#' @param nn.type The type of neighbors. The \code{tag} used when finding
#' nearest neighbors in the upstream step. If not customized, use
#' \code{"GauEps"} for NN found with \code{\link{findNNGauEB}}, or use
#' \code{"DT"} for NN found with \code{\link{findNNDT}}. Default use the most
#' lastly computed NN.
#' @export
imputeNiche.CytoSignal <- function(
    object,
    nn.type = NULL,
    weights = c("mean", "counts", "dist", "none"),
    ...
) {
  weights <- match.arg(weights)
  if (is.null(nn.type)) {
    nn.type <- object@imputation[["default"]]
  }

  message("Imputing using ", nn.type, "...")

  dge.raw <- object@counts

  if (nn.type %in% names(object@imputation)) {
    nb.id.fac <- object@imputation[[nn.type]]@nn.id
    nb.dist.fac <- object@imputation[[nn.type]]@nn.dist
  } else {
    stop("NN type not found. See options with `names(object@imputation)`.")
  }

  if (ncol(dge.raw) < length(levels(nb.id.fac))) {
    stop("Number of index beads larger than the number of beads in DGE.")
  }

  weights.mtx <- imputeNiche.dgCMatrix(
    dge.raw,
    nb.id.fac,
    nb.dist.fac,
    weights = weights
  )

  dge.raw.imputed <- dge.raw %*% weights.mtx; gc()

  # weighted sum of scale factors as well
  scale.fac <- Matrix::Matrix(object@parameters$lib.size,
                              sparse = TRUE, nrow = 1,
                              dimnames = list(NULL, colnames(dge.raw)))
  scale.fac.imp <- scale.fac %*% weights.mtx

  scale.fac.imp <- as.numeric(scale.fac.imp)
  names(scale.fac.imp) <- names(object@parameters$lib.size)

  res.density <- sum(dge.raw.imputed != 0)/length(dge.raw.imputed) # density 6.2%
  # message("Density after imputation: ", res.density*100, "%")

  object@imputation[[nn.type]]@imp.data <- dge.raw.imputed
  object@imputation[[nn.type]]@nn.graph <- weights.mtx
  object@imputation[[nn.type]]@scale.fac <- scale.fac.imp
  object@imputation[[nn.type]]@log[["Density:"]] <- paste0(res.density*100, "%")

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
#' @param slot.use The type of neighbors. The \code{tag} used when finding
#' nearest neighbors in the upstream step. If not customized, use
#' \code{"GauEps"} for NN found with \code{\link{findNNGauEB}}, or use
#' \code{"DT"} for NN found with \code{\link{findNNDT}}. Default use the most
#' lastly computed NN.
#' @param verbose Whether to print out the progress. Default \code{TRUE}.
#' @param method The method to use for normalization. default: seurat style \code{"default"}
#' @return A dgCMatrix of normalized data.
normCounts.CytoSignal <- function(
    object,
    method = c("default", "cpm", "none", "scanpy"),
    slot.use = NULL,
    verbose = TRUE,
    ...
){
  method <- match.arg(method)

  if (is.null(slot.use)) {
    slot.use <- object@imputation[["default"]]
  }

  if (!slot.use %in% names(object@imputation)) {
    stop("Data not found.")
  }
  if (isTRUE(verbose))
  message(paste0("Normalizing method: ", method, " on Imputation slot: ", slot.use, "..."))

  mat <- object@imputation[[slot.use]]@imp.data
  scale.fac <- object@imputation[[slot.use]]@scale.fac
  mat <- normCounts(object = mat, scale.fac = scale.fac, method = method)

  # modified to output the normalized data instead of the object
  # object@imputation[[slot.use]]@imp.data <- mat
  # return(object)
  return(mat)

  # if (!all.equal(names(slot(object, slot.use)), names(object@imp.norm.data))){
  # 	warning("Names of imp.data and imp.norm.data not equal!")
  # }

}


#' Sub function for inferSignif, input is a sparse matrix
#'
#' @param dge.raw A sparse matrix of raw counts
#' @param lrscore.mtx A matrix of LR scores
#' @param null.lrscore.mtx A matrix of NULL LR scores
#' @param nb.fac A factor of nearest neighbors index
#' @param intr.db A interaction database
#' @param gene_to_uniprot A dataframe of gene to uniprot mapping
#' @param p.thresh A numeric value of p-value threshold
#' @param reads.thresh A numeric value of reads threshold
#' @param sig.thresh A numeric value of significance threshold
#'
#' @return A list of indexes of significant cells
#' @export
#' @noRd
.inferSignif.matrix_like <- function(
    dge.raw,
    lrscore.mtx,
    null.lrscore.mtx,
    fdr.method = c("spatialFDR", "fdr"),
    nb.fac,
    intr.db,
    gene_to_uniprot,
    p.thresh = 0.05,
    reads.thresh = 100,
    sig.thresh = 100
){
  fdr.method <- match.arg(fdr.method)
  # if Cell numbers in pvalue mtx and nb.fac do not match
  cellsNN = as.integer(levels(nb.fac[["id"]]))
  if (!identical(nrow(lrscore.mtx), length(cellsNN))) {
    cellsNoNN <- setdiff(1:nrow(lrscore.mtx), cellsNN)
    lrscore.mtx <- lrscore.mtx[-cellsNoNN, ]
  }

  pval.mtx <- getPvalues(lrscore.mtx, null.lrscore.mtx)

  if (fdr.method == "spatialFDR"){
    pval.spatial <- graphSpatialFDR(nb.fac, pval.mtx)
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

  res.list <- res.list[order(lengths(res.list), decreasing = T)]

  res.list.hq = filterRes(dge.raw = dge.raw, res.list = res.list,
                          intr.db = intr.db, gene_to_uniprot = gene_to_uniprot,
                          reads.thresh = reads.thresh, sig.thresh = sig.thresh)

  # lrscore.mtx.hq = lrscore.mtx[, names(res.list.hq)]
  # moran.index = moranI(lrscore.mtx.hq, cells.loc)
  # res.list.hq.moran = res.list.hq[names(moran.index)]

  return(
    list(
      result = res.list,
      result.hq = res.list.hq
      # result.hq.moran = res.list.hq.moran
    )
  )

}

#' Infer significance of LR scores
#'
#' @param object A Cytosignal object
#' @param fdr.method Spatial FDR correction method. Choose from
#' \code{"spatialFDR"} or \code{"fdr"}.
#' @param p.value A numeric value of p-value threshold
#' @param reads.thresh The minimum number of reads for a ligand-receptor
#' interaction to be considered.
#' @param sig.thresh The minimum number of beads for a ligand-receptor
#' interaction to be considered.
#' @param slot.use Which LR score to use. Use the name specified with \code{tag}
#' when running \code{\link{inferLRScore}}.
#' @param nn.use Which imputation to use. Default the imputation used for
#' deriving the LRScore specified with \code{slot.use}. Use the name specified
#' with \code{tag} when running \code{\link{findNNGauEB}}; use \code{"DT"} for
#' imputation produced with \code{\link{findNNDT}}; or use \code{"Raw"} for
#' imputation produced with \code{\link{findNNRaw}}.
#'
#' @return A Cytosignal object
#' @export
inferSignif <- function(
    object,
    fdr.method = c("spatialFDR", "fdr"),
    p.value = 0.05,
    reads.thresh = 100,
    sig.thresh = 100,
    slot.use = NULL,
    nn.use = NULL
){
  if (is.null(slot.use)){
    slot.use = object@lrscore[["default"]]
  }

  if (!slot.use %in% names(object@lrscore)){
    stop("LRscores not found. ")
  }

  if (is.null(nn.use)) {
    nn.use <- object@lrscore[[slot.use]]@recep.slot
  }

  if (is.character(nn.use)) {
    if (!nn.use %in% names(object@imputation)){
      stop("Imputation slot not found.")
    }
    nb.id.fac <- object@imputation[[nn.use]]@nn.id
  } else {
    stop("`nn.use` must be either or a character.")
  }

  # else if (is.factor(nn.use)) {
  # 	if (length(nn.use) != ncol(object@imputation[[nn.use]]@intr.data))
  # 		stop("nn.use must have the same length as the number of cells.")
  # 	nb.id.fac <- nn.use
  # }

  message("Inferring significant beads on Score slot ", slot.use, "... ")

  # lrscore.mtx = object@lrscore[[slot.use]]@score
  # null.lrscore.mtx = object@lrscore[[slot.use]]@score.null

  nb.fac = list(
    id = object@imputation[[nn.use]]@nn.id,
    dist = object@imputation[[nn.use]]@nn.dist
  )

  use.intr.slot.name <- object@lrscore[[slot.use]]@intr.slot
  use.intr.db <- object@intr.valid[[use.intr.slot.name]]

  res.list = .inferSignif.matrix_like(object@counts, object@lrscore[[slot.use]]@score,
                                      object@lrscore[[slot.use]]@score.null, fdr.method, nb.fac,
                                      use.intr.db, object@intr.valid[["gene_to_uniprot"]], p.value,
                                      reads.thresh, sig.thresh)
  object@lrscore[[slot.use]]@res.list <- res.list

  object@lrscore[[slot.use]]@log[["Parameters for Significant beads"]] <- list(
    "p.value" = p.value,
    "reads.thresh" = reads.thresh,
    "sig.thresh" = sig.thresh
  )

  object@lrscore[[slot.use]]@log[["Significant Intrs"]] <- list(
    paste0("Number of interactions that have significant i-niche: ", length(res.list[["result"]])),
    paste0("Number of high quality intr: ", length(res.list[["result.hq"]]))
  )

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
  if (is.null(slot.use)){
    slot.use <- object@lrscore[["default"]]
  }
  if (!"result.hq" %in% names(object@lrscore[[slot.use]]@res.list)) {
    stop("Filtered high-quality interactions not found. Please run ",
         "`inferIntrScore()` or `inferSignif()` first.")
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


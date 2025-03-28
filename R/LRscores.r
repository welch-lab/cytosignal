#' Compute the LR score for specific ligand-receptor imputation pairs
#' @description
#' This function uses imputed ligand expression and receptor expression at each
#' spot to compute the ligand-receptor interaction scores (LR score). Usually,
#' for diffusion-dependent interactions, the ligand imputation is performed with
#' Gaussion-Epsilon (GauEps) neighbors, while the contact-dependent interactions
#' areinferred using Delaunay triangulation (DT) neighbors. Raw imputation can
#' be used for receptor, while DT imputation can also be used for receptor for
#' smoothing the imputation results.
#'
#' As explained above, the most common scenarios are: 1) Inferring
#' diffusion-dependent LR scores with \code{inferScoreLR(object, lig.imp =
#' "GauEps", recep.imp = "Raw", intr.db.name = "diff_dep")}, and 2) Inferring
#' contact-dependent LR scores with \code{inferScoreLR(object, lig.imp = "DT",
#' recep.imp = "Raw", intr.db.name = "cont_dep")}.
#' @param object A CytoSignal object, with \code{\link{imputeNiche}} performed
#' for desired analyses.
#' @param lig.imp The ligand imputation result to use. Specify the "imp.name"
#' used at \code{findNN*} functions.
#' @param recep.imp The receptor imputation result to use. Specify the
#' "imp.name" used at \code{findNN*} functions.
#' @param norm.method Method to normalize the data. Default is "default".
#' @param intr.db.use The interaction database name to use. Choose from
#' \code{"diff_dep"} for diffusion-dependent interactions, and \code{"cont_dep"}
#' for contact-dependent interactions.
#' @param lrscore.name Name of the result to be stored in object. Default a
#' combination of the ligand and receptor imputation result names.
#' @return Returns the input object with result updated within its "lrscore"
#' slot. A new "lrScores" object is created there and with LR-score stored in
#' "score" slot.
#' @export
inferScoreLR <- function(
    object,
    lig.imp,
    recep.imp,
    norm.method = c("default", "cpm", "none", "scanpy"),
    intr.db.use = c("diff_dep", "cont_dep"),
    lrscore.name = paste0(lig.imp, "-", recep.imp)
) {
  intr.db.use <- match.arg(intr.db.use)
  norm.method <- match.arg(norm.method)

  lig.imp <- .checkImpUse(object, lig.imp)
  recep.imp <- .checkImpUse(object, recep.imp)

  if (is.null(lrscore.name)) {
    stop("A `lrscore.name` name has to be specified for storing this result.")
  }

  message("Computing LR-scores using ", intr.db.use, " database.")
  message("- Ligand: ", lig.imp, ", Receptor: ", recep.imp, ".")

  # normalize using default method, normCount has been revised to internal function
  # dge.lig <- object@imputation[[lig.imp]]@imp.data
  # dge.recep <- object@imputation[[recep.imp]]@imp.data
  dge.lig <- normCounts(object, method = norm.method, imp.use = lig.imp,
                        verbose = FALSE)
  dge.recep <- normCounts(object, method = norm.method, imp.use = recep.imp,
                          verbose = FALSE)

  #----------- pre-computing the lrscores by averaging the DT scores, without norm -----------#
  if (norm.method != "none"){
    dt.avg.g <- object@imputation[[findImpByMethod(object, "DT")]]@nn.graph
    dt.avg.g <- to_mean(dt.avg.g)

    dge.lig <- dge.lig %*% dt.avg.g
    dge.recep <- dge.recep %*% dt.avg.g
  }
  #-------------------------------------------------------------------------------------------#

  if (!all.equal(dim(dge.lig), dim(dge.recep))){
    stop("Imputed ligand and receptor expression do not have the same dimension.")
  }

  intr.db.list <- checkIntr(unname(object@intr.valid[["symbols"]][["intr"]]),
                            object@intr.valid[[intr.db.use]])

  res.mtx <- .inferScoreLR.dgCMatrix(dge.lig, dge.recep,
                                     intr.db.list[["ligands"]], intr.db.list[["receptors"]])


  score.obj <- methods::new(
    "lrScores",
    lig.slot = lig.imp,
    recep.slot = recep.imp,
    intr.slot = intr.db.use,
    intr.list = intr.db.list,
    score = Matrix::Matrix(res.mtx, sparse = TRUE),
    log = list(
      Date = Sys.time()
    )
  )

  object@lrscore[["default"]] <- lrscore.name
  object@lrscore[[lrscore.name]] <- score.obj

  return(object)
}

.inferScoreLR.dgCMatrix <- function(
    dge.lig,
    dge.recep,
    lig.fac,
    recep.fac
){

  ### cavaet: remember to convert the Uniprot ids to indices!

  if (max(as.integer(names(lig.fac))) > nrow(dge.lig)){
    stop("Intr index out of dge bounds.")
  }

  lig.index = facToIndex(lig.fac)
  recep.index = facToIndex(recep.fac)

  # compute scores
  res.mtx = inferScoreLR_cpp(
    unname(as.matrix(dge.lig)),
    unname(as.matrix(dge.recep)),
    lig.index[[1]], lig.index[[2]],
    recep.index[[1]], recep.index[[2]]
  )

  dimnames(res.mtx) = list(colnames(dge.lig), levels(lig.fac))
  # res.mtx = Matrix(res.mtx, sparse = T)

  return(res.mtx)
}

#' Permute Imputation Results of specific imputation method
#' @description
#' This function permutes the imputation methods from which a given LR score is
#' calculated. Note that all rounds of permutation will use the same shuffle and
#' sample index.
#' @param object A CytoSignal object with \code{\link{inferLRScore}} run
#' beforehand.
#' @param slot.use Which LR score to use. Use the name specified with \code{tag}
#' when running \code{\link{inferLRScore}}.
#' @param perm.size Size of the permutation test. Defult \code{100000}.
#' @param norm.me Normalization method. See \code{\link{normCounts}} for detail.
#' @return A CytoSignal object
#' @export
permuteLR <- function(
    object,
    slot.use = NULL,
    norm.method = "default",
    perm.size = 100000
) {
  if (is.null(slot.use)){
    slot.use <- object@lrscore[["default"]]
  }

  if (!slot.use %in% names(object@lrscore)){
    stop("Cannot find corresponding LR score.")
  }

  # decide how many times to permute according to the size of the data
  if (!is.numeric(perm.size)){
    stop("perm.size must be a numeric value.")
  }

  dge.raw <- object@counts

  if (perm.size < ncol(dge.raw)){
    message("Permutation size too small, using ", ncol(dge.raw), " instead.")
    perm.size <- ncol(dge.raw)
  }

  #### get all data needed for permutation
  lig.slot <- object@lrscore[[slot.use]]@lig.slot
  recep.slot <- object@lrscore[[slot.use]]@recep.slot
  # cells.loc <- object@cells.loc
  intr.valid <- object@lrscore[[slot.use]]@intr.list

  times <- ceiling(perm.size / ncol(dge.raw))
  each.size <- ceiling(perm.size / times)

  message("Permuting whole dataset ", times, " times...")

  perm.idx.list <- lapply(1:times, function(x){ sample(ncol(dge.raw)) })
  sample.idx <- sample(ncol(dge.raw), each.size)

  message("- Permuting ligand Imp slot: ", lig.slot, "...")

  pb <- utils::txtProgressBar(min = 0, max = length(perm.idx.list), style = 3, file = stderr())
  null.lig.list <- list()
  for (i in seq_along(perm.idx.list)) {
    idx <- perm.idx.list[[i]]
    null.lig.list[[i]] <- permuteLR.sparse(object, nn.type = lig.slot, perm.index = idx, sample.idx, norm.method = norm.method)
    utils::setTxtProgressBar(pb, i)
  }
  cat('\n')
  # null.lig.list <- lapply(perm.idx.list, function(idx) {
  #   permuteLR.sparse(object, nn.type = lig.slot, perm.index = idx, sample.idx, norm.method = norm.method)
  # })
  null.lig.dge <- cbind_list(null.lig.list)

  message("- Permuting receptor Imp slot: ", recep.slot, "...")

  pb <- utils::txtProgressBar(min = 0, max = length(perm.idx.list), style = 3, file = stderr())
  null.recep.list <- list()
  for (i in seq_along(perm.idx.list)) {
    idx <- perm.idx.list[[i]]
    null.recep.list[[i]] <- permuteLR.sparse(object, nn.type = recep.slot, perm.index = idx, sample.idx, norm.method = norm.method)
    utils::setTxtProgressBar(pb, i)
  }
  cat('\n')

  null.recep.dge <- cbind_list(null.recep.list)

  message("- Calculating NULL scores...")

  object@lrscore[[slot.use]]@lig.null <- null.lig.dge
  object@lrscore[[slot.use]]@recep.null <- null.recep.dge
  object@lrscore[[slot.use]]@perm.idx <- list(
    perm = perm.idx.list,
    sample = sample.idx
  )

  rm(null.lig.list, null.recep.list, null.lig.dge, null.recep.dge); gc()

  return(object)

}




#' Permute imputation data on corresponded imputation methods, using the same shuffle index
#'
#' This function permutes the default or user-defined imputation method and stored the results in the
#' default or user-defined ImpData object. A index vector is needed to permute the data.
#'
#' @param object A Cytosignal object
#' @param nn.type The imputation method to use
#' @param perm.size Size of the permutation test
#'
#' @return A dgCMatrix object
#' @noRd
permuteLR.sparse <- function(
    object,
    nn.type = NULL,
    perm.index = NULL,
    sample.index = NULL,
    norm.method = "default"
){
  # check DT imputation first, this is for pre-computing lrscore.mtx
  DT.imp <- findImpByMethod(object, "DT")
  if (is.null(DT.imp)) {
    stop("Need to infer NN with DT method first.")
  }

  nn.type <- .checkImpUse(object, nn.type)

  if (!is.integer(perm.index)){
    stop("perm.index not valid.")
  } else if (length(perm.index) != ncol(object@counts)){
    stop("Length of perm.index not valid.")
  }

  if (!is.integer(sample.index)){
    stop("Sample index must be an integer vector.")
  }

  if (min(sample.index) < 1 || max(sample.index) > ncol(object@counts)){
    stop("Sample index out of range.")
  }

  # use.gau <- grepl("GauEps", nn.type)
  # use.dt <- grepl("DT", nn.type)
  # use.raw <- grepl("Raw", nn.type)

  dge.raw <- object@counts
  scale.fac <- object@parameters[["lib.size"]]
  scale.fac <- Matrix::Matrix(scale.fac, nrow = 1, byrow = T, sparse = T)

  if (nn.type == "Raw") {
    # The nn.graph stored in Raw imputation object is taken from DT, but we need
    # to use the Diagonal (i.e. Identity) matrix here
    nn.graph <- Matrix::Diagonal(ncol(dge.raw))
  } else {
    nn.graph <- object@imputation[[nn.type]]@nn.graph
  }
  null.graph <- nn.graph[perm.index, ]

  # #---------------------------------------------------------------------------------------------------#
  # #----------- if permuting the neighbors again, use either of the two following functions -----------#
  # null.graph <- shuffle_sp_mat_col(null.graph)
  # null.graph <- shuffleEdgeRandomNB(null.graph)
  # #---------------------------------------------------------------------------------------------------#

  # #---------------------------------------------------------------------------------------------------#
  # #-----------  if permuting dge.raw colnames instead of the graph, use the following code -----------#
  # null.dge.raw <- dge.raw[, perm.index]
  # scale.fac <- object@parameters[["lib.size"]][perm.index]
  # colnames(null.dge.raw) <- colnames(dge.raw)
  # names(scale.fac) <- colnames(dge.raw)
  # #---------------------------------------------------------------------------------------------------#

  # first: permute the dge.raw and scale.fac
  null.dge <- dge.raw %*% null.graph
  null.scale.fac <- scale.fac %*% null.graph

  null.dge <- normCounts(null.dge, scale.fac = as.numeric(null.scale.fac),
                         norm.method)

  # #----------- pre-computing the lrscores by averaging the DT scores, without norm -----------#
  dt.avg.g <- object@imputation[[DT.imp]]@nn.graph
  dt.avg.g <- to_mean(dt.avg.g)[perm.index, ]

  # sample the null graph to control the size of the permutation
  dt.avg.g <- dt.avg.g[, sample.index]
  # second: avg across null DT neighbors, without norm
  null.dge <- null.dge %*% dt.avg.g

  gc()

  return(null.dge)

}


#' Permute LR score for specific ligand-receptor imputation obj pairs
#'
#' This function is a follow-up function of inferScoreLR. It computes the NULL LR-scores
#' using the NULL imputation results and stores the results in the LR score object.
#' The null distribution of the LR scores can be used to test the significance of the LR scores.
#'
#' @param object A Cytosignal object
#' @param slot.use Which LR score to use. Use the name specified with \code{tag}
#' when running \code{\link{inferLRScore}}.
#'
#' @return A Cytosignal object
#' @export
#'
inferNullScoreLR <- function(
    object,
    slot.use = NULL
){
  slot.use <- .checkSlotUse(object, slot.use)

  message("Permuting scores on Score slot: ", slot.use, "...")

  # score.obj <- object@lrscore[[slot.use]]
  # lig.slot <- object@lrscore[[slot.use]]@lig.slot
  # recep.slot <- object@lrscore[[slot.use]]@recep.slot
  intr.valid <- object@lrscore[[slot.use]]@intr.list
  null.dge.lig <- object@lrscore[[slot.use]]@lig.null
  null.dge.recep <- object@lrscore[[slot.use]]@recep.null

  # if (ncol(null.dge.lig) != ncol(null.dge.recep)){
  # 	null.dge.recep <- null.dge.recep[, sample(ncol(null.dge.recep),
  # 						ncol(null.dge.lig), replace = T)]
  # }

  null.lrscore.mtx <- .inferScoreLR.dgCMatrix(null.dge.lig, null.dge.recep,
                                              intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()
  null.lrscore.mtx <- Matrix::Matrix(null.lrscore.mtx, sparse = T)

  colnames(null.lrscore.mtx) <- colnames(object@lrscore[[slot.use]]@score)
  rownames(null.lrscore.mtx) <- paste0("perm_", 1:nrow(null.lrscore.mtx))

  if (sum(colSums(null.lrscore.mtx) == 0) != 0){
    message("- A total of ", sum(colSums(null.lrscore.mtx) == 0), " interactions are empty in NULL scores.")
    null.lrscore.mtx = null.lrscore.mtx[, !colSums(null.lrscore.mtx) == 0]
  }

  intr.both <- intersect(colnames(null.lrscore.mtx), colnames(object@lrscore[[slot.use]]@score))

  if (length(intr.both) != ncol(null.lrscore.mtx)) {
    message("- Removing ", ncol(null.lrscore.mtx) - length(intr.both), " more intr from NULL scores.")
    null.lrscore.mtx = null.lrscore.mtx[, intr.both]
  }

  if (length(intr.both) != ncol(object@lrscore[[slot.use]]@score)) {
    message("- Removing ", ncol(object@lrscore[[slot.use]]@score) - length(intr.both), " corresponding intr from REAL scores.")
    object@lrscore[[slot.use]]@score = object@lrscore[[slot.use]]@score[, intr.both]
  }

  object@lrscore[[slot.use]]@score.null <- null.lrscore.mtx

  return(object)
}

#' Smooth existing LR-scores using spatial nearest neighbors
#' @description
#' This function takes an existing LR-score out and smooth it using the spatial
#' nearest neighbor graph inferred with Gaussian-Epsilon method. If the radius
#' queried (\code{eps}) is the same as any pre-calculated NN graph during
#' imputation stage, this function will directly use it. However, a different
#' \code{eps} can also be used. In this case, a new NN graph will be inferred
#' and stored in the object.
#' @param object A CytoSignal object with existing LR-scores available in
#' "lrscore" slot.
#' @param score.use Name of LR-scores to be smoothed. Default use the most
#' recently calculated LR-scores.
#' @param eps The radius of the epsilon ball. Default use parameter inferred
#' with \code{\link{inferEpsParams}}. See Description as well.
#' @param lrscore.name Name of the smoothed LR-scores to be stored in the object.
#' Default appending \code{"_smoothed"} after \code{score.use}.
#' @return Returns the original object with result updated in "lrscore" slot.
#' An "lrScores" object is created with its "score" slot filled a cell by
#' interaction sparse matrix indicating the smoothed LR-score. When a new GauEps
#' NN graph is inferred, it will be stored in the "imputation" slot. See
#' \code{\link{findNNGauEB}} for detail.
#' @export
smoothScoreLR <- function(
    object,
    score.use = NULL,
    eps = NULL,
    lrscore.name = NULL
) {
  # Find the score to be smoothed
  score.use <- .checkSlotUse(object, slot.use = score.use, velo = FALSE)
  # lrscore.name <- lrscore.name %||% paste0(score.use, "_smooth")
  # eps <- eps %||% object@parameters$r.diffuse.scale
  lrscore.name <- if (!is.null(lrscore.name)) lrscore.name else paste0(score.use, "_smooth")
  eps <- if (!is.null(eps)) eps else object@parameters$r.diffuse.scale
  # Check if there has been GauEps NN with the same `eps` pre-calculated
  nn.avail <- names(object@imputation)
  nn.avail <- nn.avail[nn.avail != "default"]
  nn.precal <- NULL
  for (nn.name in nn.avail) {
    if (object@imputation[[nn.name]]@method != "GauEps") next
    pre.eps <- object@imputation[[nn.name]]@log$parameters$eps
    if (is.null(pre.eps)) next
    if (!is.numeric(pre.eps)) next
    if (length(pre.eps) != 1) next
    if (pre.eps != eps) next
    nn.precal <- nn.name
    break
  }
  if (!is.null(nn.precal)) {
    message("GauEps NN with eps = ", eps, " already pre-calculated, using it.")
  } else {
    message("Capculating GauEps NN with eps = ", eps, "...")
    nn.precal <- paste0("GauEps_eps", round(eps, 2))
    object <- findNNGauEB(object, eps = eps, imp.name = nn.precal)
    message("New NN stored at ", nn.precal, ".")
  }
  message("Smoothing LR-scores from ", score.use, " using NN from ", nn.precal, "...")
  # nn.graph: N cell x N cell
  nn.graph <- object@imputation[[nn.precal]]@nn.graph
  # score: N cell x I intr
  score <- object@lrscore[[score.use]]@score
  score.smooth <- t(t(score) %*% nn.graph)

  # impute Null lrscore
  null.score <- t(object@lrscore[[score.use]]@score.null)
  perm.idx <- object@lrscore[[score.use]]@perm.idx[[1]]
  size.idx <- object@lrscore[[score.use]]@perm.idx[[2]]
  n_cell <- length(size.idx)

  message("Imputing Null LR-scores...")
  new.null.score <- lapply(1:length(perm.idx), function(i){
    start.idx <- 1+(i-1)*n_cell
    end.idx <- i*n_cell
    sub.null.score <- null.score[, start.idx:end.idx]

    null.graph <- nn.graph[sample(n_cell), sample(n_cell)]

    smooth.score <- sub.null.score %*% null.graph

    return(smooth.score)
  })

  smooth.null.score <- t(cbind_list(new.null.score))

  colnames(smooth.null.score) <- colnames(score)
  rownames(smooth.null.score) <- paste0("perm_", 1:nrow(smooth.null.score))
  # object@lrscore[[score.use]]@score.null <- t(smooth.null.score)

  lrscoreObj <- object@lrscore[[score.use]]
  lrscoreObj@score <- score.smooth
  lrscoreObj@score.null <- smooth.null.score
  lrscoreObj@res.list <- list()
  lrscoreObj@log$Smooth.NN <- nn.precal
  lrscoreObj@log$Date <- Sys.time()
  object@lrscore[["default"]] <- lrscore.name
  object@lrscore[[lrscore.name]] <- lrscoreObj

  return(object)
}

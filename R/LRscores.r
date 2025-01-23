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
    dt.avg.g <- object@imputation[["DT"]]@nn.graph
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
    lig.null = methods::new("dgCMatrix"),
    recep.null = methods::new("dgCMatrix"),
    intr.slot = intr.db.use,
    intr.list = intr.db.list,
    score = Matrix::Matrix(res.mtx, sparse = T),
    score.null = methods::new("matrix"),
    perm.idx = list(),
    res.list = list(),
    log = list(
      "Parmaeters:" = c(lig.imp, recep.imp),
      "Date:" = Sys.time()
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
  if (!"DT" %in% names(object@imputation)) {
    stop("Need to run DT imputation first.")
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
  dt.avg.g <- object@imputation[["DT"]]@nn.graph
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

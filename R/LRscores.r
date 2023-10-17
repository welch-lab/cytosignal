#' Compute the LR score for specific ligand-receptor imputation obj pairs
#' 
#' @param object A Cytosignal object
#' @param lig.slot The ligand slot to use
#' @param recep.slot The receptor slot to use
#' @param intr.db.name The intr database name to use
#' @param nn.use The neighbor index as niche
#' 
#' @return A Cytosignal object
#' @export
#' 
inferScoreLR <- function(
  object,
  ...
) {
  UseMethod(generic = 'inferScoreLR', object = object)
}

#' Sub function for inferScoreLR, input a sparse matrix
#' 
#' @param dge.lig A sparse matrix for ligand
#' @param dge.recep A sparse matrix for receptor
#' @param lig.fac A factor of ligand indices
#' @param recep.fac A factor of receptor indices
#' 
#' @return A sparse matrix
#' @export
#' 
inferScoreLR.dgCMatrix <- function(
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


#' Sub function for inferScoreLR, input a CytoSignal object
#' 
#' @param object A Cytosignal object
#' @param lig.slot The ligand slot to use
#' @param recep.slot The receptor slot to use
#' @param intr.db.name The intr database name to use
#' @param nn.use slot that the neighbor index should be taken from, by default is the same as
#' 			the recep.slot. For example, if score.obj = GauEps-DT, then nn.use = "DT".
#' 			nn.use could also be a user-defind factor.
#' 
#' @return A Cytosignal object
#' @export
inferScoreLR.CytoSignal <- function(
	object,
	lig.slot,
	recep.slot,
	intr.db.name,
	tag = NULL
){
	# check DT imputation first, this is for pre-computing lrscore.mtx
	if (!"DT" %in% names(object@imputation)) {
		stop("Need to run DT imputation first.")
	}

	if (!lig.slot %in% names(object@imputation)){
		stop("Ligand slot not found.")
	}

	if (!recep.slot %in% names(object@imputation)){
		stop("Receptor slot not found.")
	}

	if (!intr.db.name %in% c("diff_dep", "cont_dep")) {
		stop("intr.db.name must be either 'diff_dep' or 'cont_dep'.")
	}

	if (is.null(tag)) {
		tag <- paste0(lig.slot, "-", recep.slot)
	}

	message("Computing LR-scores using ", intr.db.name, " database.")
	message("Ligand: ", lig.slot, ", Receptor: ", recep.slot, ".")

	dge.lig <- object@imputation[[lig.slot]]@imp.data
	dge.recep <- object@imputation[[recep.slot]]@imp.data

	#----------- pre-computing the lrscores by averaging the DT scores, without norm -----------#
	message("Comfirming niche index...")

	dt.avg.g <- object@imputation[["DT"]]@nn.graph
	dt.avg.g <- to_mean(dt.avg.g)

	dge.lig <- dge.lig %*% dt.avg.g
	dge.recep <- dge.recep %*% dt.avg.g
	#-------------------------------------------------------------------------------------------#

	if (!all.equal(dim(dge.lig), dim(dge.recep))){
		stop("dge.lig and dge.recep must have the same dimension.")
	}

	intr.db.list <- checkIntr(unname(object@intr.valid[["symbols"]][["intr"]]), 
							object@intr.valid[[intr.db.name]])

	message("Calculating LR-scores...")

	res.mtx <- inferScoreLR.dgCMatrix(dge.lig, dge.recep,
				intr.db.list[["ligands"]], intr.db.list[["receptors"]])

	message("Done!")

	score.obj <- methods::new(
		"lrScores",
		lig.slot = lig.slot,
		recep.slot = recep.slot,
		lig.null = methods::new("dgCMatrix"),
		recep.null = methods::new("dgCMatrix"),
		intr.slot = intr.db.name,
		intr.list = intr.db.list,
		score = Matrix::Matrix(res.mtx, sparse = T),
		score.null = methods::new("matrix"),
		perm.idx = list(),
		res.list = list(),
		log = list(
			"Parmaeters:" = c(lig.slot, recep.slot),
			"Date:" = Sys.time()
		)
	)

	object@lrscore[["default"]] <- tag
	object@lrscore[[tag]] <- score.obj

	return(object)
}




#' Permute LR score for specific ligand-receptor imputation obj pairs
#' 
#' This function is a follow-up function of inferScoreLR. It computes the NULL LR-scores 
#' using the NULL imputation results and stores the results in the LR score object.
#' The null distribution of the LR scores can be used to test the significance of the LR scores.
#' 
#' @param object A Cytosignal object
#' @param nn.type The imputation method to use
#' 
#' @return A Cytosignal object
#' @export
#' 
inferNullScoreLR <- function(
	object,
	slot.use = NULL
){
	if (is.null(slot.use)){
		slot.use <- object@lrscore[["default"]]
	}

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

	null.lrscore.mtx <- inferScoreLR.dgCMatrix(null.dge.lig, null.dge.recep,
					intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()
    null.lrscore.mtx <- Matrix::Matrix(null.lrscore.mtx, sparse = T)

	colnames(null.lrscore.mtx) <- colnames(object@lrscore[[slot.use]]@score)
	rownames(null.lrscore.mtx) <- paste0("perm_", 1:nrow(null.lrscore.mtx))

    if (sum(Matrix::colSums(null.lrscore.mtx) == 0) != 0){
        message("A total of ", sum(Matrix::colSums(null.lrscore.mtx) == 0), " intr are empty in NULL scores.")
        null.lrscore.mtx = null.lrscore.mtx[, !Matrix::colSums(null.lrscore.mtx) == 0]
    }

	intr.both <- intersect(colnames(null.lrscore.mtx), colnames(object@lrscore[[slot.use]]@score))

	if (length(intr.both) != ncol(null.lrscore.mtx)) {
		message("Removing ", ncol(null.lrscore.mtx) - length(intr.both), " more intr from NULL scores.")
		null.lrscore.mtx = null.lrscore.mtx[, intr.both]
	}

	if (length(intr.both) != ncol(object@lrscore[[slot.use]]@score)) {
		message("Removing ", ncol(object@lrscore[[slot.use]]@score) - length(intr.both), " corresponding intr from REAL scores.")
		object@lrscore[[slot.use]]@score = object@lrscore[[slot.use]]@score[, intr.both]
	}

	object@lrscore[[slot.use]]@score.null <- null.lrscore.mtx

    return(object)
}




#' Permute Imputation Results of specific imputation method
#' 
#' This function is a conveninent wrapper of `permuteLR.sparse`, which permutes the imputation methods 
#' from which a given LR score is calculated. Note that all rounds of permutation will use the same
#' shuffle and sample index.
#'
#' 
#' @param object A Cytosignal object
#' @param slot.use The `lrscore` slot to use
#' @param perm.size Size of the permutation test
#' 
#' @return A Cytosignal object
#' @export
#' 
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

	message("Permuting ligand Imp slot: ", lig.slot, "...")

	null.lig.list <- lapply(perm.idx.list, function(idx){
		permuteLR.sparse(object, nn.type = lig.slot, perm.index = idx, sample.idx, norm.method = norm.method)
	})
	null.lig.dge <- cbind_list(null.lig.list)

	message("Permuting receptor Imp slot: ", recep.slot, "...")

	null.recep.list <- lapply(perm.idx.list, function(idx){
		permuteLR.sparse(object, nn.type = recep.slot, perm.index = idx, sample.idx, norm.method = norm.method)
	})
	null.recep.dge <- cbind_list(null.recep.list)
	
    message("Calculating NULL scores...")

    object@lrscore[[slot.use]]@lig.null <- null.lig.dge
    object@lrscore[[slot.use]]@recep.null <- null.recep.dge
    object@lrscore[[slot.use]]@perm.idx <- list(
        perm = perm.idx.list,
        sample = sample.idx
    )

    rm(null.lig.list, null.recep.list, null.lig.dge, null.recep.dge); gc()

    cat("Done!\n")

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
#' @export
#' 
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

	if (is.null(nn.type)){
		nn.type <- object@imputation[["default"]]
	}

	if (!nn.type %in% names(object@imputation)){
		stop("Cannot find corresponding imputation method.")
	}

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
    
	nn.graph <- object@imputation[[nn.type]]@nn.graph
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
	
	null.dge <- normCounts.list(list(mat=null.dge, scale.fac=as.numeric(null.scale.fac)), norm.method)

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


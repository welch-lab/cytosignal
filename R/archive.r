# #' Permute scores on corresponded imputation methods, using the same shuffle index
# #' 
# #' This function is a follow-up function of imputeNiche, which permutes the default or user-defined
# #' imputation method and stored the results in the ImpData object.
# #' 
# #' @param object A Cytosignal object
# #' @param nn.type The imputation method to use
# #' @param perm.size Size of the permutation test
# #' 
# #' @return A Cytosignal object
# #' @export
# permuteRandomNBScore <- function(
# 	object,
# 	slot.use = NULL,
# 	perm.size = 100000
# ){
# 	if (is.null(slot.use)){
# 		slot.use <- object@lrscore[["default"]]
# 	}

# 	message("Permuting scores on Score slot: ", slot.use, "...")

# 	# score.obj <- object@lrscore[[slot.use]]
# 	lig.slot <- object@lrscore[[slot.use]]@lig.slot
# 	recep.slot <- object@lrscore[[slot.use]]@recep.slot
# 	cells.loc <- object@cells.loc
# 	intr.valid <- object@lrscore[[slot.use]]@intr.list

# 	if (!nn.type %in% names(object@imputation)){
# 		stop("Cannot find corresponding imputation method.")
# 	}

# 	# call function permuteRandomScore separately, but add parameter fro shuffule index


# 	use.gau <- grepl("GauEps", nn.type)
# 	use.dt <- grepl("DT", nn.type)
# 	use.raw <- grepl("Raw", nn.type)

# 	# if (!use.gau & !use.dt){
# 	# 	stop("Cannot find corresponding imputation method.")
# 	# }

# 	message("Permuting imputation data on method ", nn.type, "...")

# 	dimnames.list <- dimnames(object@imputation[[nn.type]]@imp.data)
# 	cells.loc <- object@cells.loc
# 	dge.raw <- object@counts

# 	# decide how many times to permute according to the size of the data
# 	if (!is.numeric(perm.size)){
# 		stop("perm.size must be a numeric value.")
# 	}

# 	if (perm.size < ncol(dge.raw)){
# 		message("Permutation size too small, using ", ncol(dge.raw), " instead.")
# 		perm.size <- ncol(dge.raw)
# 	}
	
# 	# reason for computing by batch is to save costs
# 	times <- ceiling(perm.size / ncol(dge.raw))
# 	each.size <- ceiling(perm.size / times)

# 	message("Permuting whole dataset ", times, " times...")

# 	# permute cell locations and remember to change the rownames!!
# 	null.cells.loc = cells.loc[sample(nrow(cells.loc)), ]
# 	rownames(null.cells.loc) = rownames(cells.loc)

# 	# use the same imputation graph for permutation
# 	graph.nn <- object@imputation[[nn.type]]@nn.graph

# 	cat("Permuting NULL expressions...\nTimes No.")

# 	null.dge.list <- lapply(1:times, function(i){
# 		## shuffule dge, remember to change the colnames!!
# 		## could also shuffle the null graph, mathematically the same
# 		cat(i, ", ", sep = "")

# 		# # if shuffle dge.raw
# 		# perm.index <- sample(ncol(dge.raw))
# 		# null.dge.raw <- dge.raw[, perm.index]
# 		# scale.fac <- object@parameters[["lib.size"]][perm.index]
# 		# colnames(null.dge.raw) <- colnames(dge.raw)
# 		# names(scale.fac) <- colnames(dge.raw)

# 		scale.fac <- object@parameters[["lib.size"]]
# 		scale.fac <- Matrix::Matrix(scale.fac, nrow = 1, byrow = T, sparse = T) # computing scale factor
# 		# null.dge.raw <- changeUniprot.matrix_like(null.dge.raw.all, object@intr.valid[["gene_to_uniprot"]])[[1]]

# 		if (use.raw) {
# 			# sample the dge.raw to control the size of the permutation
# 			# not really necessary, but put here for consistency
# 			sample.index <- sample(ncol(dge.raw), each.size)
# 			null.dge.eps <- dge.raw[, sample.index]
# 			scale.fac.imp <- scale.fac[sample.index]
# 		} else {
# 			##### core permutation step #####
# 			# null.graph <- shuffle_sp_mat_col(graph.nn)
# 			null.graph <- shuffleEdgeRandomNB(graph.nn)
# 			##### core permutation step #####

# 			# sample the null graph to control the size of the permutation
# 			null.graph <- null.graph[, sample(ncol(null.graph), each.size)]

# 			null.dge.eps <- dge.raw %*% null.graph
# 			scale.fac.imp <- scale.fac %*% null.graph
# 		}
		
# 		null.dge.eps <- normCounts.list(list(mat=null.dge.eps, scale.fac=as.numeric(scale.fac.imp)), "default")
# 		# rm(null.dge.raw)

# 		return(null.dge.eps)
# 	})

# 	null.dge.imp <- cbind_list(null.dge.list)
# 	rm(null.dge.list); gc()

# 	cat("End.\nFinished!\n")

# 	object@imputation[[nn.type]]@imp.data.null <- null.dge.imp

# 	return(object)

# }







# #' Permute Imputation Results of specific imputation method
# #' 
# #' This function is a follow-up function of imputeNiche, which permutes the default or user-defined
# #' imputation method and stored the results in the ImpData object.
# #' 
# #' @param object A Cytosignal object
# #' @param nn.type The imputation method to use
# #' @param perm.size Size of the permutation test, by defualt NULL
# #' @param times Number of times to permute the whole dataset, by default 10
# #' 
# #' @return A Cytosignal object
# #' @export
# permuteImpLR_old <- function(
# 	object,
# 	nn.type = NULL,
# 	perm.size = NULL,
# 	times = 10
# ){
# 	if (is.null(nn.type)){
# 		nn.type <- object@imputation[["default"]]
# 	}

# 	if (!nn.type %in% names(object@imputation)){
# 		stop("Cannot find corresponding imputation method.")
# 	}

# 	use.gau <- grepl("GauEps", nn.type)
# 	use.dt <- grepl("DT", nn.type)
# 	use.raw <- grepl("Raw", nn.type)

# 	# if (!use.gau & !use.dt){
# 	# 	stop("Cannot find corresponding imputation method.")
# 	# }

# 	message("Permuting imputation data on method ", nn.type, "...")

# 	dimnames.list <- dimnames(object@imputation[[nn.type]]@imp.data)
# 	cells.loc <- object@cells.loc
# 	dge.raw <- object@counts

# 	if (!is.null(times)) {
# 		if (!is.numeric(times)){
# 			stop("times must be a numeric value.")
# 		}

# 		times <- ceiling(times)

# 		if (times < 1) {
# 			stop("times must be larger than 1.")
# 		}

# 	} else {
# 		# decide how many times to permute according to the size of the data
# 		if (is.null(perm.size)){
# 			message("Permutation size not specified, using 100000 instead.")
# 			perm.size <- 100000
# 		}

# 		if (!is.numeric(perm.size)){
# 			stop("perm.size must be a numeric value.")
# 		}

# 		if (perm.size < ncol(dge.raw)){
# 			message("Permutation size too small, using ", ncol(dge.raw), " instead.")
# 			perm.size <- ncol(dge.raw)
# 		}
		
# 		times <- ceiling(perm.size / ncol(dge.raw))
# 	}

# 	message("Permuting whole dataset ", times, " times...")

# 	# permute cell locations
# 	null.cells.loc.list = lapply(1:times, function(i){
# 		null.cells.loc = cells.loc[sample(nrow(cells.loc)), ]
# 		rownames(null.cells.loc) = rownames(cells.loc)
# 		return(null.cells.loc)
# 	})
	
# 	# Null imputations
# 	null.dge.list <- lapply(1:times, function(i){
# 		## MUST shuffule raw dge!!
# 		perm.index <- sample(ncol(dge.raw))
# 		null.dge.raw <- dge.raw[, perm.index]
# 		scale.fac <- object@parameters[["lib.size"]][perm.index]
		
# 		colnames(null.dge.raw) <- colnames(dge.raw)
# 		names(scale.fac) <- colnames(dge.raw)

# 		scale.fac <- Matrix::Matrix(scale.fac, nrow = 1, byrow = T, sparse = T) # computing scale factor
# 		# null.dge.raw <- changeUniprot.matrix_like(null.dge.raw.all, object@intr.valid[["gene_to_uniprot"]])[[1]]

# 		if (use.gau) {
# 			param.eps <- object@parameters$r.diffuse.scale
# 			param.sigma <- object@parameters$sigma.scale

# 			null.eps.nb <- findNNGauEB.matrix(null.cells.loc.list[[i]], eps = param.eps, 
# 							sigma = param.sigma, weight.sum = 2); gc()

# 			null.dge.eps <- imputeNiche.dgCMatrix(null.dge.raw, null.eps.nb$id, null.eps.nb$dist, weights = "none"); gc()
# 			scale.fac.imp <- imputeNiche.dgCMatrix(scale.fac, null.eps.nb$id, null.eps.nb$dist, weights = "none")
# 		} else if (use.dt) {
# 			null.nb.fac <- findNNDT.matrix(null.cells.loc.list[[i]])
#         	null.dge.eps <- imputeNiche.dgCMatrix(null.dge.raw, null.nb.fac$id, null.nb.fac$dist, weights = "none"); gc() # den: 31%
# 			scale.fac.imp <- imputeNiche.dgCMatrix(scale.fac, null.nb.fac$id, null.nb.fac$dist, weights = "none")
# 		} else if (use.raw) {
# 			null.dge.eps <- null.dge.raw
# 			scale.fac.imp <- scale.fac
# 		} else {
# 			stop("Cannot find corresponding imputation method.")
# 		}
		
# 		null.dge.eps = normCounts.list(list(mat=null.dge.eps, scale.fac=as.numeric(scale.fac.imp)), "default"); gc()
# 		rm(null.dge.raw); gc()

# 		return(null.dge.eps)
# 	})

# 	null.dge.eps.unpt = meanMat_cpp(null.dge.list, nrow(null.dge.list[[1]]), ncol(null.dge.list[[1]]))
# 	dimnames(null.dge.eps.unpt) <- dimnames.list
# 	rm(null.dge.list); gc()

#     object@imputation[[nn.type]]@imp.data.null <- null.dge.eps.unpt

# 	return(object)

# }



#' Compute the LR velo for specific ligand-receptor imputation obj pairs
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
veloNicheLR <- function(
  object,
  ...
) {
  UseMethod(generic = 'veloNicheLR', object = object)
}

#' Sub function for veloNicheLR, input a sparse matrix
#' 
#' @param dge.lig A sparse matrix for ligand
#' @param dge.recep A sparse matrix for receptor
#' @param nb.id.fac A factor of neighbor indices
#' @param lig.fac A factor of ligand indices
#' @param recep.fac A factor of receptor indices
#' 
#' @return A sparse matrix
#' @export
#' 
veloNicheLR.matrix_like <- function(
	dge.lig,
    dge.recep,
	velo.mtx,
    nb.id.fac,
	lig.fac,
	recep.fac
){

	### cavaet: remember to convert the Uniprot ids to indices!
	# convert nb fac
	nb.id.fac = sort(nb.id.fac)
	nb.index = facToIndex(nb.id.fac)

	if (max(as.integer(names(lig.fac))) > nrow(dge.lig)){
		stop("Intr index out of dge bounds.")
	}

	lig.index = facToIndex(lig.fac)
	recep.index = facToIndex(recep.fac)

	# compute velos
	res.mtx = VeloinferScoreLR_cpp(
		unname(as.matrix(dge.lig)),
		unname(as.matrix(dge.recep)),
		unname(as.matrix(velo.mtx)),
		lig.index[[1]], lig.index[[2]], 
		recep.index[[1]], recep.index[[2]],
		nb.index[[1]], nb.index[[2]]
	)

	dimnames(res.mtx) = list(colnames(dge.lig)[as.integer(levels(nb.id.fac))], levels(lig.fac))
	# dimnames(res.mtx) = list(colnames(dge.lig), levels(lig.fac))
	# res.mtx = Matrix(res.mtx, sparse = T)

	return(res.mtx)
}



#' Sub function for veloNicheLR, input a CytoSignal object
#' 
#' @param object A Cytosignal object
#' @param lig.slot The ligand slot to use
#' @param recep.slot The receptor slot to use
#' @param intr.db.name The intr database name to use
#' @param nn.use slot that the neighbor index should be taken from, by default is the same as
#' 			the recep.slot. For example, if velo.obj = GauEps-DT, then nn.use = "DT".
#' 			nn.use could also be a user-defind factor.
#' 
#' @return A Cytosignal object
#' @export
veloNicheLR.CytoSignal <- function(
	object,
	lig.slot,
	recep.slot,
	intr.db.name,
	tag = NULL
){
	if (!lig.slot %in% names(object@imputation)){
		stop("Ligand slot not found.")
	}

	if (!recep.slot %in% names(object@imputation)){
		stop("Receptor slot not found.")
	}

	message("Computing velos using ", intr.db.name, " database.")
	message("Ligand: ", lig.slot, ", Receptor: ", recep.slot, ".")

	if (!intr.db.name %in% c("diff_dep", "cont_dep")) {
		stop("intr.db.name must be either 'diff_dep' or 'cont_dep'.")
	}

	if (is.null(tag)) {
		tag <- paste0(lig.slot, "-", recep.slot)
	}

	if (tag %in% names(object@lrvelo)) {
		stop("Tag already exists.")
	}
	
	object@lrvelo[["default"]] <- tag

	# if (is.null(slot.use)) {
	# 	slot.use <- object@lrscore[["default"]]
	# }

	# if (!slot.use %in% names(object@lrscore)) {
	# 	stop("lrvelo slot not found.")
	# }

	# if (is.null(nn.use)) {
	# 	nn.use <- recep.slot
	# }
	
	# if (is.character(nn.use)) {
	# 	if (!nn.use %in% names(object@imputation)){
	# 		stop("Imputation slot not found.")
	# 	}
	# 	nb.id.fac <- object@imputation[[nn.use]]@nn.id
	# } else if (is.factor(nn.use)) {
	# 	if (length(nn.use) != ncol(object@imputation[[lig.slot]]@intr.data))
	# 		stop("nn.use must have the same length as the number of cells.")
	# 	nb.id.fac <- nn.use
	# } else {
	# 	stop("nn.use must be either a factor or a character.")
	# }

	dge.lig <- object@imputation[[lig.slot]]@intr.data
	dge.recep <- object@imputation[[recep.slot]]@intr.data

	if (!all.equal(dim(dge.lig), dim(dge.recep))){
		stop("dge.lig and dge.recep must have the same dimension.")
	}

	if (!all.equal(object@intr.valid[["symbols"]][[lig.slot]], 
					object@intr.valid[["symbols"]][[recep.slot]])){
		stop("Unpt symbols generated from imputations not equal!")
	}
	
	# velo and dge should have the same cells and genes
	use.cells <- intersect(colnames(dge.lig), colnames(object@velo[["velo.s"]]))
	use.genes <- intersect(rownames(dge.lig), rownames(object@velo[["velo.s"]]))

	# infer new neighbor index since the cells are filtered
	new.cells.loc <- object@cells.loc[use.cells, ]
	velo.nn.list <- findNNDT.matrix(new.cells.loc)

	# dge and velo should have the same cells and genes
	dge.lig <- dge.lig[use.genes, use.cells]
	dge.recep <- dge.recep[use.genes, use.cells]
	velo.s <- object@velo[["velo.s"]][use.genes, use.cells]
	velo.u <- object@velo[["velo.u"]][use.genes, use.cells]

	message("Number of velo cells: ", ncol(dge.lig), " / ", ncol(object@imputation[[lig.slot]]@intr.data))
	message("Number of velo genes: ", nrow(dge.lig), " / ", nrow(object@imputation[[lig.slot]]@intr.data))

	intr.db.list <- checkIntr(use.genes, object@intr.valid[[intr.db.name]])

	res.mtx <- veloNicheLR.matrix_like(dge.lig, dge.recep, 
				(velo.s + velo.u), velo.nn.list[["id"]], 
				intr.db.list[["ligands"]], intr.db.list[["receptors"]])

	lrvelo.obj <- new(
		"lrVelo",
		lig.slot = lig.slot,
		recep.slot = recep.slot,
		intr.slot = intr.db.name,
		intr.list = intr.db.list,
		dim.valid = list(
			"cells" = use.cells,
			"genes" = use.genes),
		intr.velo = res.mtx,
		nn.id = velo.nn.list[["id"]],
		nn.dist = velo.nn.list[["dist"]],
		log = list(
			"Used slot" = c(lig.slot, recep.slot)
		)
	)

	object@lrvelo[[tag]] <- lrvelo.obj

	return(object)
}






# #' Permute Imputation Results of specific imputation method
# #' 
# #' This function is a follow-up function of imputeNiche, which permutes the default or user-defined
# #' imputation method and stored the results in the ImpData object.
# #' 
# #' @param object A Cytosignal object
# #' @param nn.type The imputation method to use
# #' @param perm.size Size of the permutation test
# #' 
# #' @return A Cytosignal object
# #' @export
# permuteImpLRSort <- function(
# 	object,
# 	nn.type = NULL,
# 	perm.size = 20000
# ){
# 	if (is.null(nn.type)){
# 		nn.type <- object@imputation[["default"]]
# 	}

# 	if (!nn.type %in% names(object@imputation)){
# 		stop("Cannot find corresponding imputation method.")
# 	}

# 	use.gau <- grepl("GauEps", nn.type)
# 	use.dt <- grepl("DT", nn.type)
# 	use.raw <- grepl("Raw", nn.type)

# 	# if (!use.gau & !use.dt){
# 	# 	stop("Cannot find corresponding imputation method.")
# 	# }

# 	message("Permuting imputation data on method ", nn.type, "...")

# 	dimnames.list <- dimnames(object@imputation[[nn.type]]@imp.data)
# 	cells.loc <- object@cells.loc
# 	dge.raw <- object@raw.counts

# 	# decide how many times to permute according to the size of the data
# 	if (!is.numeric(perm.size)){
# 		stop("perm.size must be a numeric value.")
# 	}

# 	if (perm.size < ncol(dge.raw)){
# 		message("Permutation size too small, using ", ncol(dge.raw), " instead.")
# 		perm.size <- ncol(dge.raw)
# 	}
	
# 	times <- ceiling(perm.size / ncol(dge.raw))

# 	message("Permuting whole dataset ", times, " times...")

# 	# permute cell locations
# 	null.cells.loc.list = lapply(1:times, function(i){
# 		perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
# 		null.cells.loc = cells.loc[perm.index,]
# 		rownames(null.cells.loc) = rownames(cells.loc)
# 		return(null.cells.loc)
# 	})
	
# 	# Null imputations
# 	null.dge.list <- lapply(1:times, function(i){
# 		## MUST shuffule raw dge!!
# 		null.dge.raw.all <- dge.raw[, sample(seq_len(ncol(dge.raw)), ncol(dge.raw), replace = F)]
# 		scale.fac <- Matrix::Matrix(Matrix::colSums(null.dge.raw.all), nrow = 1, byrow = T, sparse = T) # computing scale factor
# 		null.dge.raw <- changeUniprot.matrix_like(null.dge.raw.all, object@intr.valid[["gene_to_uniprot"]])[[1]]

# 		if (use.gau) {
# 			param.eps <- object@parameters$r.diffuse.scale
# 			param.sigma <- object@parameters$sigma.scale

# 			null.eps.nb <- findNNGauEB.matrix(null.cells.loc.list[[i]], eps = param.eps, 
# 							sigma = param.sigma, weight.sum = 2); gc()

# 			null.dge.eps <- imputeNiche.dgCMatrix(null.dge.raw, null.eps.nb$id, null.eps.nb$dist, weights = "none"); gc()
# 			scale.fac.imp <- imputeNiche.dgCMatrix(scale.fac, null.eps.nb$id, null.eps.nb$dist, weights = "none")
# 		} else if (use.dt) {
# 			null.nb.fac <- findNNDT.matrix(null.cells.loc.list[[i]])
#         	null.dge.eps <- imputeNiche.dgCMatrix(null.dge.raw, null.nb.fac$id, null.nb.fac$dist, weights = "none"); gc() # den: 31%
# 			scale.fac.imp <- imputeNiche.dgCMatrix(scale.fac, null.nb.fac$id, null.nb.fac$dist, weights = "none")
# 		} else if (use.raw) {
# 			null.dge.eps <- null.dge.raw
# 			scale.fac.imp <- scale.fac
# 		} else {
# 			stop("Cannot find corresponding imputation method.")
# 		}
		
# 		null.dge.eps = normCounts.list(list(mat=null.dge.eps, scale.fac=as.numeric(scale.fac.imp)), "default"); gc()
# 		rm(null.dge.raw); gc()

# 		return(null.dge.eps)
# 	})

# 	# null.dge.eps.unpt = meanMat_cpp(null.dge.list, nrow(null.dge.list[[1]]), ncol(null.dge.list[[1]]))
#     null.dge.list <- lapply(null.dge.list, function(x) as.matrix(x))
#     null.dge.eps.unpt <- Matrix::Matrix(sample_rows_cpp(null.dge.list, 10000), sparse = T)
# 	# dimnames(null.dge.eps.unpt) <- dimnames.list
#     rownames(null.dge.eps.unpt) <- dimnames.list[[1]]
#     colnames(null.dge.eps.unpt) <- paste0("perm_", 1:10000)
# 	rm(null.dge.list); gc()

#     object@imputation[[nn.type]]@imp.data.null <- null.dge.eps.unpt

# 	return(object)

# }


#' Permute LR score for specific ligand-receptor imputation obj pairs
#' 
#' This function is a follow-up function of inferScoreLR. It computes the NUL LRscores 
#' using the NULL imputation results and stores the results in the LR score object.
#' The null distribution of the LR scores can be used to test the significance of the LR scores.
#' 
#' @param object A Cytosignal object
#' @param slot.use The imputation method to use
#' 
#' @return A Cytosignal object
#' @export
permuteNicheScoreLR <- function(
	object,
	slot.use = NULL
){
	if (is.null(slot.use)){
		slot.use <- object@lrscore[["default"]]
	}

	message("Permuting scores on Score slot: ", slot.use, "...")

	# score.obj <- object@lrscore[[slot.use]]
	lig.slot <- object@lrscore[[slot.use]]@lig.slot
	recep.slot <- object@lrscore[[slot.use]]@recep.slot
	cells.loc <- object@cells.loc
	intr.valid <- object@lrscore[[slot.use]]@intr.list

	null.dge.eps.unpt <- object@imputation[[lig.slot]]@imp.data.null
	null.dge.dt.unpt <- object@imputation[[recep.slot]]@imp.data.null

    perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
    null.cells.loc = cells.loc[perm.index,]
    rownames(null.cells.loc) = rownames(cells.loc)
    null.nb.fac = findNNDT.matrix(null.cells.loc); gc() # Mean num of neighbors: 44, median: 36

    null.lrscore.mtx = graphNicheLR.dgCMatrix(null.dge.eps.unpt, null.dge.dt.unpt, null.nb.fac[["id"]],
						intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()

	# null.lrscore.mtx <- null.dge.eps.unpt * null.dge.dt.unpt
 
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




# #' Permute Imputation Results of specific imputation method
# #' 
# #' This function is a follow-up function of imputeNiche, which permutes the default or user-defined
# #' imputation method and stored the results in the ImpData object.
# #' 
# #' @param object A Cytosignal object
# #' @param nn.type The imputation method to use
# #' @param perm.size Size of the permutation test
# #' 
# #' @return A Cytosignal object
# #' @export
# permuteListImpLR_ori <- function(
# 	object,
# 	nn.type = NULL,
# 	perm.size = 100000
# ){
# 	if (is.null(nn.type)){
# 		nn.type <- object@imputation[["default"]]
# 	}

# 	if (!nn.type %in% names(object@imputation)){
# 		stop("Cannot find corresponding imputation method.")
# 	}

# 	use.gau <- grepl("GauEps", nn.type)
# 	use.dt <- grepl("DT", nn.type)
# 	use.raw <- grepl("Raw", nn.type)

# 	# if (!use.gau & !use.dt){
# 	# 	stop("Cannot find corresponding imputation method.")
# 	# }

# 	message("Permuting imputation data on method ", nn.type, "...")

# 	dimnames.list <- dimnames(object@imputation[[nn.type]]@imp.data)
# 	cells.loc <- object@cells.loc
# 	dge.raw <- object@raw.counts

# 	# decide how many times to permute according to the size of the data
# 	if (!is.numeric(perm.size)){
# 		stop("perm.size must be a numeric value.")
# 	}

# 	if (perm.size < ncol(dge.raw)){
# 		message("Permutation size too small, using ", ncol(dge.raw), " instead.")
# 		perm.size <- ncol(dge.raw)
# 	}
	
# 	times <- ceiling(perm.size / ncol(dge.raw))
# 	each.size <- ceiling(perm.size / times)

# 	message("Permuting whole dataset ", times, " times...")

# 	# permute cell locations and remember to change the rownames!!
# 	null.cells.loc = cells.loc[sample(nrow(cells.loc)), ]
# 	rownames(null.cells.loc) = rownames(cells.loc)
	
# 	# get a random graph for permuting
# 	if (use.gau) {
# 		param.eps <- object@parameters$r.diffuse.scale
# 		param.sigma <- object@parameters$sigma.scale

# 		null.eps.nb <- findNNGauEB.matrix(null.cells.loc, eps = param.eps, 
# 						sigma = param.sigma, weight.sum = 2); gc()

# 		null.graph <- imputeNiche.dgCMatrix(dge.raw, null.eps.nb$id, null.eps.nb$dist, 
# 							weights = "none", return.graph = T); gc()
# 	} else if (use.dt) {
# 		null.nb.fac <- findNNDT.matrix(null.cells.loc)
# 		null.graph <- imputeNiche.dgCMatrix(dge.raw, null.nb.fac$id, null.nb.fac$dist, 
# 							weights = "none", return.graph = T); gc()
# 	} else if (use.raw) {
# 		cat("No need to generate a random graph for raw imputation.\n")
# 	} else {
# 		stop("Cannot find corresponding imputation method.")
# 	}

# 	cat("Permuting NULL expressions...\nTimes No.")

# 	null.dge.list <- lapply(1:times, function(i){
# 		## shuffule dge, remember to change the colnames!!
# 		## could also shuffle the null graph, mathematically the same
# 		cat(i, ", ", sep = "")
# 		perm.index <- sample(ncol(dge.raw))
# 		null.dge.raw <- dge.raw[, perm.index]
# 		scale.fac <- object@parameters[["lib.size"]][perm.index]
		
# 		colnames(null.dge.raw) <- colnames(dge.raw)
# 		names(scale.fac) <- colnames(dge.raw)

# 		scale.fac <- Matrix::Matrix(scale.fac, nrow = 1, byrow = T, sparse = T) # computing scale factor
# 		# null.dge.raw <- changeUniprot.matrix_like(null.dge.raw.all, object@intr.valid[["gene_to_uniprot"]])[[1]]

# 		if (use.raw) {
# 			# sample the dge.raw to control the size of the permutation
# 			sample.index <- sample(ncol(null.dge.raw), each.size)

# 			null.dge.eps <- null.dge.raw[, sample.index]
# 			scale.fac.imp <- scale.fac[sample.index]
# 		} else {
# 			# sample the null graph to control the size of the permutation
# 			null.graph <- null.graph[, sample(ncol(null.graph), each.size)]

# 			null.dge.eps <- null.dge.raw %*% null.graph; gc()
# 			scale.fac.imp <- scale.fac %*% null.graph; gc()
# 		}
		
# 		null.dge.eps <- normCounts.list(list(mat=null.dge.eps, scale.fac=as.numeric(scale.fac.imp)), "default"); gc()
# 		rm(null.dge.raw); gc()

# 		return(null.dge.eps)
# 	})

# 	cat("End.\nFinished!\n")

# 	# null.dge.eps.unpt = meanMat_cpp(null.dge.list, nrow(null.dge.list[[1]]), ncol(null.dge.list[[1]]))
# 	# null.dge.eps.unpt = combine_sparse_rows(null.dge.list)
# 	# dimnames(null.dge.eps.unpt) <- dimnames.list
#     # rownames(null.dge.eps.unpt) <- dimnames.list[[1]]
#     # colnames(null.dge.eps.unpt) <- paste0("perm_", 1:ncol(null.dge.eps.unpt))
# 	# rm(null.dge.list); gc()

#     # object@imputation[[nn.type]]@imp.data.null <- null.dge.eps.unpt
# 	object@imputation[[nn.type]]@imp.data.null <- null.dge.list

# 	return(object)

# }




#' Permute Imputation Results of specific imputation method
#' 
#' This function is a follow-up function of imputeNiche, which permutes the default or user-defined
#' imputation method and stored the results in the ImpData object.
#' 
#' @param object A Cytosignal object
#' @param nn.type The imputation method to use
#' @param perm.size Size of the permutation test
#' 
#' @return A Cytosignal object
#' @export
permuteRandomNB <- function(
	object,
	nn.type = NULL,
	perm.size = 100000
){
	if (is.null(nn.type)){
		nn.type <- object@imputation[["default"]]
	}

	if (!nn.type %in% names(object@imputation)){
		stop("Cannot find corresponding imputation method.")
	}

	use.gau <- grepl("GauEps", nn.type)
	use.dt <- grepl("DT", nn.type)
	use.raw <- grepl("Raw", nn.type)

	# if (!use.gau & !use.dt){
	# 	stop("Cannot find corresponding imputation method.")
	# }

	message("Permuting imputation data on method ", nn.type, "...")

	dimnames.list <- dimnames(object@imputation[[nn.type]]@imp.data)
	cells.loc <- object@cells.loc
	dge.raw <- object@counts

	# decide how many times to permute according to the size of the data
	if (!is.numeric(perm.size)){
		stop("perm.size must be a numeric value.")
	}

	if (perm.size < ncol(dge.raw)){
		message("Permutation size too small, using ", ncol(dge.raw), " instead.")
		perm.size <- ncol(dge.raw)
	}
	
	# reason for computing by batch is to save costs
	times <- ceiling(perm.size / ncol(dge.raw))
	each.size <- ceiling(perm.size / times)

	message("Permuting whole dataset ", times, " times...")

	# permute cell locations and remember to change the rownames!!
	null.cells.loc = cells.loc[sample(nrow(cells.loc)), ]
	rownames(null.cells.loc) = rownames(cells.loc)

	# use the same imputation graph for permutation
	graph.nn <- object@imputation[[nn.type]]@nn.graph

	cat("Permuting NULL expressions...\nTimes No.")

	null.dge.list <- lapply(1:times, function(i){
		## shuffule dge, remember to change the colnames!!
		## could also shuffle the null graph, mathematically the same
		cat(i, ", ", sep = "")

		# # if shuffle dge.raw
		# perm.index <- sample(ncol(dge.raw))
		# null.dge.raw <- dge.raw[, perm.index]
		# scale.fac <- object@parameters[["lib.size"]][perm.index]
		# colnames(null.dge.raw) <- colnames(dge.raw)
		# names(scale.fac) <- colnames(dge.raw)

		scale.fac <- object@parameters[["lib.size"]]
		scale.fac <- Matrix::Matrix(scale.fac, nrow = 1, byrow = T, sparse = T) # computing scale factor
		# null.dge.raw <- changeUniprot.matrix_like(null.dge.raw.all, object@intr.valid[["gene_to_uniprot"]])[[1]]

		if (use.raw) {
			# sample the dge.raw to control the size of the permutation
			# not really necessary, but put here for consistency
			sample.index <- sample(ncol(dge.raw), each.size)
			null.dge.eps <- dge.raw[, sample.index]
			scale.fac.imp <- scale.fac[sample.index]
		} else {
			##### core permutation step #####
			# null.graph <- shuffle_sp_mat_col(graph.nn)
			null.graph <- shuffleEdgeRandomNB(graph.nn)
			##### core permutation step #####

			# sample the null graph to control the size of the permutation
			null.graph <- null.graph[, sample(ncol(null.graph), each.size)]

			null.dge.eps <- dge.raw %*% null.graph
			scale.fac.imp <- scale.fac %*% null.graph
		}
		
		null.dge.eps <- normCounts.list(list(mat=null.dge.eps, scale.fac=as.numeric(scale.fac.imp)), "default")
		# rm(null.dge.raw)

		return(null.dge.eps)
	})

	null.dge.imp <- cbind_list(null.dge.list)
	rm(null.dge.list); gc()

	cat("End.\nFinished!\n")

	object@imputation[[nn.type]]@imp.data.null <- null.dge.imp

	return(object)

}




#' Permute LR score for specific ligand-receptor imputation obj pairs
#' 
#' This function is a follow-up function of inferScoreLR. It computes the NUL LRscores 
#' using the NULL imputation results and stores the results in the LR score object.
#' The null distribution of the LR scores can be used to test the significance of the LR scores.
#' 
#' @param object A Cytosignal object
#' @param slot.use The imputation method to use
#' 
#' @return A Cytosignal object
#' @export
permuteRandomScore <- function(
	object,
	slot.use = NULL
){
	if (is.null(slot.use)){
		slot.use <- object@lrscore[["default"]]
	}

	message("Permuting scores on Score slot: ", slot.use, "...")

	# score.obj <- object@lrscore[[slot.use]]
	lig.slot <- object@lrscore[[slot.use]]@lig.slot
	recep.slot <- object@lrscore[[slot.use]]@recep.slot
	cells.loc <- object@cells.loc
	intr.valid <- object@lrscore[[slot.use]]@intr.list

	null.dge.lig <- object@imputation[[lig.slot]]@imp.data.null
	null.dge.recep <- object@imputation[[recep.slot]]@imp.data

	if (ncol(null.dge.lig) != ncol(null.dge.recep)){
		null.dge.recep <- null.dge.recep[, sample(ncol(null.dge.recep), 
							ncol(null.dge.lig), replace = T)]
	}

	null.lrscore.mtx <- inferScoreLR.dgCMatrix(null.dge.lig, null.dge.recep,
					intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()

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
#' This function is a follow-up function of imputeNiche, which permutes the default or user-defined
#' imputation method and stored the results in the ImpData object.
#' 
#' @param object A Cytosignal object
#' @param nn.type The imputation method to use
#' @param perm.size Size of the permutation test
#' 
#' @return A Cytosignal object
#' @export
permuteListImpLR <- function(
	object,
	nn.type = NULL,
	perm.size = 100000
){
	if (is.null(nn.type)){
		nn.type <- object@imputation[["default"]]
	}

	if (!nn.type %in% names(object@imputation)){
		stop("Cannot find corresponding imputation method.")
	}

	use.gau <- grepl("GauEps", nn.type)
	use.dt <- grepl("DT", nn.type)
	use.raw <- grepl("Raw", nn.type)

	# if (!use.gau & !use.dt){
	# 	stop("Cannot find corresponding imputation method.")
	# }

	message("Permuting imputation data on method ", nn.type, "...")

	dimnames.list <- dimnames(object@imputation[[nn.type]]@imp.data)
	cells.loc <- object@cells.loc
	dge.raw <- object@counts

	# decide how many times to permute according to the size of the data
	if (!is.numeric(perm.size)){
		stop("perm.size must be a numeric value.")
	}

	if (perm.size < ncol(dge.raw)){
		message("Permutation size too small, using ", ncol(dge.raw), " instead.")
		perm.size <- ncol(dge.raw)
	}
	
	times <- ceiling(perm.size / ncol(dge.raw))
	each.size <- ceiling(perm.size / times)

	message("Permuting whole dataset ", times, " times...")

	# permute cell locations and remember to change the rownames!!
	null.cells.loc = cells.loc[sample(nrow(cells.loc)), ]
	rownames(null.cells.loc) = rownames(cells.loc)

	# use the same imputation graph for permutation
	null.graph <- object@imputation[[nn.type]]@nn.graph

	cat("Permuting NULL expressions...\nTimes No.")

	null.dge.list <- lapply(1:times, function(i){
		## shuffule dge, remember to change the colnames!!
		## could also shuffle the null graph, mathematically the same
		cat(i, ", ", sep = "")
		perm.index <- sample(ncol(dge.raw))
		null.dge.raw <- dge.raw[, perm.index]
		scale.fac <- object@parameters[["lib.size"]][perm.index]
		# scale.fac <- object@parameters[["lib.size"]]
		
		colnames(null.dge.raw) <- colnames(dge.raw)
		names(scale.fac) <- colnames(dge.raw)

		scale.fac <- Matrix::Matrix(scale.fac, nrow = 1, byrow = T, sparse = T) # computing scale factor
		# null.dge.raw <- changeUniprot.matrix_like(null.dge.raw.all, object@intr.valid[["gene_to_uniprot"]])[[1]]

		if (use.raw) {
			# sample the dge.raw to control the size of the permutation
			sample.index <- sample(ncol(null.dge.raw), each.size)

			null.dge.eps <- null.dge.raw[, sample.index]
			scale.fac.imp <- scale.fac[sample.index]
		} else {
			# try to permute the values of each column
			null.graph <- shuffle_sp_mat_col(null.graph)
			
			# sample the null graph to control the size of the permutation
			null.graph <- null.graph[, sample(ncol(null.graph), each.size)]

			null.dge.eps <- null.dge.raw %*% null.graph
			scale.fac.imp <- scale.fac %*% null.graph
		}
		
		null.dge.eps <- normCounts.list(list(mat=null.dge.eps, scale.fac=as.numeric(scale.fac.imp)), "default")
		rm(null.dge.raw)

		return(null.dge.eps)
	})

	null.dge.imp <- cbind_list(null.dge.list)
	rm(null.dge.list); gc()

	cat("End.\nFinished!\n")

	object@imputation[[nn.type]]@imp.data.null <- null.dge.imp

	return(object)

}



#' Permute LR score for specific ligand-receptor imputation obj pairs
#' 
#' This function is a follow-up function of inferScoreLR. It computes the NUL LRscores 
#' using the NULL imputation results and stores the results in the LR score object.
#' The null distribution of the LR scores can be used to test the significance of the LR scores.
#' 
#' @param object A Cytosignal object
#' @param slot.use The imputation method to use
#' 
#' @return A Cytosignal object
#' @export
permuteListScoreLR <- function(
	object,
	slot.use = NULL
){
	if (is.null(slot.use)){
		slot.use <- object@lrscore[["default"]]
	}

	message("Permuting scores on Score slot: ", slot.use, "...")

	# score.obj <- object@lrscore[[slot.use]]
	lig.slot <- object@lrscore[[slot.use]]@lig.slot
	recep.slot <- object@lrscore[[slot.use]]@recep.slot
	cells.loc <- object@cells.loc
	intr.valid <- object@lrscore[[slot.use]]@intr.list

	null.dge.eps <- object@imputation[[lig.slot]]@imp.data.null
	null.dge.dt <- object@imputation[[recep.slot]]@imp.data.null

	null.lrscore.mtx <- inferScoreLR.dgCMatrix(null.dge.eps, null.dge.dt,
					intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()

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





# #' Permute LR score for specific ligand-receptor imputation obj pairs
# #' 
# #' This function is a follow-up function of inferScoreLR. It computes the NUL LRscores 
# #' using the NULL imputation results and stores the results in the LR score object.
# #' The null distribution of the LR scores can be used to test the significance of the LR scores.
# #' 
# #' @param object A Cytosignal object
# #' @param slot.use The imputation method to use
# #' 
# #' @return A Cytosignal object
# #' @export
# permuteListScoreLR_ori_2 <- function(
# 	object,
# 	slot.use = NULL
# ){
# 	if (is.null(slot.use)){
# 		slot.use <- object@lrscore[["default"]]
# 	}

# 	message("Permuting scores on Score slot: ", slot.use, "...")

# 	# score.obj <- object@lrscore[[slot.use]]
# 	lig.slot <- object@lrscore[[slot.use]]@lig.slot
# 	recep.slot <- object@lrscore[[slot.use]]@recep.slot
# 	cells.loc <- object@cells.loc
# 	intr.valid <- object@lrscore[[slot.use]]@intr.list

# 	null.dge.eps.list <- object@imputation[[lig.slot]]@imp.data.null
# 	null.dge.dt.list <- object@imputation[[recep.slot]]@imp.data.null

#     # perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
#     # null.cells.loc = cells.loc[perm.index,]
#     # rownames(null.cells.loc) = rownames(cells.loc)
#     # null.nb.fac = findNNDT.matrix(null.cells.loc); gc() # Mean num of neighbors: 44, median: 36


# 	null.lrscore.mtx.list <- lapply(seq_along(null.dge.eps.list), function(i){
# 		null.dge.eps.unpt <- null.dge.eps.list[[i]]
# 		null.dge.dt.unpt <- null.dge.dt.list[[i]]

# 		# ## if permuting NULL mtx again
# 		# null.dge.eps.unpt = null.dge.eps.unpt[, sample(ncol(null.dge.eps.unpt))]
# 		# null.dge.dt.unpt = null.dge.dt.unpt[, sample(ncol(null.dge.dt.unpt))]

# 		null.lrscore.mtx <- inferScoreLR.dgCMatrix(null.dge.eps.unpt, null.dge.dt.unpt,
# 					intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()
		
# 		return( Matrix::Matrix(null.lrscore.mtx, sparse = T) )
# 	})

# 	null.lrscore.mtx <- combine_sparse_cols(null.lrscore.mtx.list)
# 	colnames(null.lrscore.mtx) <- colnames(object@lrscore[[slot.use]]@score)
# 	rownames(null.lrscore.mtx) <- paste0("perm_", 1:nrow(null.lrscore.mtx))
	
# 	rm(null.lrscore.mtx.list); gc()

#     if (sum(Matrix::colSums(null.lrscore.mtx) == 0) != 0){
#         message("A total of ", sum(Matrix::colSums(null.lrscore.mtx) == 0), " intr are empty in NULL scores.")
#         null.lrscore.mtx = null.lrscore.mtx[, !Matrix::colSums(null.lrscore.mtx) == 0]
#     }

# 	intr.both <- intersect(colnames(null.lrscore.mtx), colnames(object@lrscore[[slot.use]]@score))

# 	if (length(intr.both) != ncol(null.lrscore.mtx)) {
# 		message("Removing ", ncol(null.lrscore.mtx) - length(intr.both), " more intr from NULL scores.")
# 		null.lrscore.mtx = null.lrscore.mtx[, intr.both]
# 	}

# 	if (length(intr.both) != ncol(object@lrscore[[slot.use]]@score)) {
# 		message("Removing ", ncol(object@lrscore[[slot.use]]@score) - length(intr.both), " corresponding intr from REAL scores.")
# 		object@lrscore[[slot.use]]@score = object@lrscore[[slot.use]]@score[, intr.both]
# 	}

# 	object@lrscore[[slot.use]]@score.null <- null.lrscore.mtx

#     return(object)
# }



# nullDiffusibleLR <- function(
# 	object,
# 	times = 10,
# 	slot.use = NULL
# ){
# 	if (is.null(slot.use)){
# 		slot.use <- object@lrscore[["default"]]
# 	}

# 	# score.obj <- object@lrscore[[slot.use]]
# 	lig.slot <- object@lrscore[[slot.use]]@lig.slot
# 	recep.slot <- object@lrscore[[slot.use]]@recep.slot
# 	dimnames.list <- dimnames(object@imputation[[lig.slot]]@intr.data)
	
# 	cells.loc <- object@cells.loc
# 	dge.raw <- object@counts
# 	eps.params <- object@parameters

# 	null.cells.loc.list = lapply(1:times, function(i){
# 		perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
# 		null.cells.loc = cells.loc[perm.index,]
# 		rownames(null.cells.loc) = rownames(cells.loc)
# 		return(null.cells.loc)
# 	})

# 	null.dge.list <- lapply(1:times, function(i){
# 		# shuffule raw dge
# 		null.dge.raw = dge.raw[, sample(seq_len(ncol(dge.raw)), ncol(dge.raw), replace = F)]
# 		null.eps.nb = findNNGauEB.matrix(null.cells.loc.list[[i]], eps = eps.params[[1]], sigma = eps.params[[2]], weight.sum = 2)
# 		null.dge.eps = imputeNiche.dgCMatrix(null.dge.raw, null.eps.nb$id, null.eps.nb$dist, weights = "gaussian"); gc()

# 		null.dge.eps = normCounts.dgCMatrix(null.dge.eps, "default"); gc()
# 		null.dge.eps.unpt = changeUniprot.dgCMatrix(null.dge.eps)[[1]]
# 		rm(null.dge.eps); gc()

# 		return(null.dge.eps.unpt)
# 	})

# 	null.dge.eps.unpt = meanMat_cpp(null.dge.list, nrow(null.dge.list[[1]]), ncol(null.dge.list[[1]]))
# 	dimnames(null.dge.eps.unpt) <- dimnames.list
# 	rm(null.dge.list); gc()

#     object@imputation[[lig.slot]]@intr.data.null <- null.dge.eps.unpt

#     # permute DT as well
#     null.dge.dt.list <- lapply(1:times, function(i){
#         null.dge.raw = dge.raw[, sample(seq_len(ncol(dge.raw)), ncol(dge.raw), replace = F)]
#         null.nb.fac = findNNDT.matrix(null.cells.loc.list[[i]]); gc() # Mean num of neighbors: 44, median: 36
#         null.dge.raw.dt = imputeNiche.dgCMatrix(null.dge.raw, null.nb.fac$id, null.nb.fac$dist, weights = "dist"); gc() # den: 31%

#         null.dge.dt = normCounts.dgCMatrix(null.dge.raw.dt, "default"); gc()
#         null.dge.dt.unpt = changeUniprot.dgCMatrix(null.dge.dt)[[1]]
#         rm(null.dge.dt); gc()

#         return(null.dge.dt.unpt)
#     })

#     null.dge.dt.unpt = meanMat_cpp(null.dge.dt.list, nrow(null.dge.dt.list[[1]]), ncol(null.dge.dt.list[[1]]))
#     dimnames(null.dge.dt.unpt) <- dimnames.list
#     rm(null.dge.dt.list); gc()

#     object@imputation[[recep.slot]]@intr.data.null <- null.dge.dt.unpt

#     perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
#     null.cells.loc = cells.loc[perm.index,]
#     rownames(null.cells.loc) = rownames(cells.loc)
#     null.nb.fac = findNNDT.matrix(null.cells.loc); gc() # Mean num of neighbors: 44, median: 36

#     null.lrscore.mtx = graphNicheLR.dgCMatrix(null.dge.eps.unpt, null.dge.dt.unpt, null.nb.fac[["id"]],
# 						object@intr.valid[["ligands"]], object@intr.valid[["receptors"]]); gc()

#     if (sum(colSums(null.lrscore.mtx) == 0) != 0){
#         cat("Removing ", sum(colSums(null.lrscore.mtx) == 0), " NULL intr.\n")
#         null.lrscore.mtx = null.lrscore.mtx[, !colSums(null.lrscore.mtx) == 0]
#     }

# 	intr.both <- intersect(colnames(null.lrscore.mtx), colnames(object@lrscore[[slot.use]]@score))

# 	if (length(intr.both) != ncol(null.lrscore.mtx)) {
# 		cat("Removing ", ncol(null.lrscore.mtx) - length(intr.both), " NULL intr.\n")
# 		null.lrscore.mtx = null.lrscore.mtx[, intr.both]
# 	}

# 	if (length(intr.both) != ncol(object@lrscore[[slot.use]]@score)) {
# 		cat("Removing ", ncol(object@lrscore[[slot.use]]@score) - length(intr.both), " NULL intr.\n")
# 		object@lrscore[[slot.use]]@score = object@lrscore[[slot.use]]@score[, intr.both]
# 	}

# 	object@lrscore[[slot.use]]@score.null <- null.lrscore.mtx

#     return(object)

# }




# permuteImpLR_ori <- function(
# 	object,
# 	nn.type = NULL,
# 	times = 10
# ){
# 	if (is.null(nn.type)){
# 		nn.type <- object@imputation[["default"]]
# 	}

# 	use.gau <- grepl("GauEps", nn.type)
# 	use.dt <- grepl("DT", nn.type)

# 	if (!use.gau & !use.dt){
# 		stop("Cannot find corresponding imputation method.")
# 	}

# 	message("Permuting imputation data on method ", nn.type, "...")

# 	dimnames.list <- dimnames(object@imputation[[nn.type]]@intr.data)
# 	cells.loc <- object@cells.loc
# 	dge.raw <- object@counts
	
# 	# permute cell locations
# 	null.cells.loc.list = lapply(1:times, function(i){
# 		perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
# 		null.cells.loc = cells.loc[perm.index,]
# 		rownames(null.cells.loc) = rownames(cells.loc)
# 		return(null.cells.loc)
# 	})
	
# 	# Null imputations
# 	null.dge.list <- lapply(1:times, function(i){
# 		# shuffule raw dge
# 		null.dge.raw = dge.raw[, sample(seq_len(ncol(dge.raw)), ncol(dge.raw), replace = F)]

# 		if (use.gau){
# 			eps.params <- object@parameters
# 			null.eps.nb <- findNNGauEB.matrix(null.cells.loc.list[[i]], eps = eps.params[[1]], sigma = eps.params[[2]], weight.sum = 2)
# 			null.dge.eps <- imputeNiche.dgCMatrix(null.dge.raw, null.eps.nb$id, null.eps.nb$dist, weights = "gaussian"); gc()
# 		} else if (use.dt){
# 			null.nb.fac <- findNNDT.matrix(null.cells.loc.list[[i]])
#         	null.dge.eps <- imputeNiche.dgCMatrix(null.dge.raw, null.nb.fac$id, null.nb.fac$dist, weights = "dist"); gc() # den: 31%
# 		} else {
# 			stop("Cannot find corresponding imputation method.")
# 		}

# 		null.dge.eps = normCounts.dgCMatrix(null.dge.eps, "default"); gc()
# 		null.dge.eps.unpt = changeUniprot.dgCMatrix(null.dge.eps, object@intr.valid[["gene_to_uniprot"]])[[1]]
# 		rm(null.dge.eps); gc()

# 		return(null.dge.eps.unpt)
# 	})

# 	null.dge.eps.unpt = meanMat_cpp(null.dge.list, nrow(null.dge.list[[1]]), ncol(null.dge.list[[1]]))
# 	dimnames(null.dge.eps.unpt) <- dimnames.list
# 	rm(null.dge.list); gc()

#     object@imputation[[nn.type]]@intr.data.null <- null.dge.eps.unpt

# 	return(object)

# }



#' Sub function for graphNicheLR, input a CytoSignal object
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
graphNicheLR.CytoSignal <- function(
	object,
	lig.slot,
	recep.slot,
	intr.db.name,
	nn.use = NULL
){
	if (!lig.slot %in% names(object@imputation)){
		stop("Ligand slot not found.")
	}

	if (!recep.slot %in% names(object@imputation)){
		stop("Receptor slot not found.")
	}

	if (!intr.db.name %in% c("diff_dep", "cont_dep")) {
		stop("intr.db.name must be either 'diff_dep' or 'cont_dep'.")
	}

	if (is.null(nn.use)) {
		nn.use <- recep.slot
	}
	
	if (is.character(nn.use)) {
		if (!nn.use %in% names(object@imputation)){
			stop("Imputation slot not found.")
		}
		nb.id.fac <- object@imputation[[nn.use]]@nn.id
	} else if (is.factor(nn.use)) {
		if (length(nn.use) != ncol(object@imputation[[lig.slot]]@imp.data))
			stop("nn.use must have the same length as the number of cells.")
		nb.id.fac <- nn.use
	} else {
		stop("nn.use must be either a factor or a character.")
	}

	dge.lig <- object@imputation[[lig.slot]]@imp.data
	dge.recep <- object@imputation[[recep.slot]]@imp.data

	if (!all.equal(dim(dge.lig), dim(dge.recep))){
		stop("dge.lig and dge.recep must have the same dimension.")
	}

	if (ncol(dge.lig) < length(levels(nb.id.fac))){
		stop("Number of index beads larger than the number of beads in DGE.")
	}

	# if (!all.equal(object@intr.valid[["symbols"]][[lig.slot]], 
	# 				object@intr.valid[["symbols"]][[recep.slot]])){
	# 	stop("Unpt symbols generated from imputations not equal!")
	# }

	message("Computing scores using ", intr.db.name, " database.")
	message("Ligand: ", lig.slot, ", Receptor: ", recep.slot, ".")

	intr.db.list <- checkIntr(unname(object@intr.valid[["symbols"]][["intr"]]), 
							object@intr.valid[[intr.db.name]])

	res.mtx <- graphNicheLR.dgCMatrix(dge.lig, dge.recep, nb.id.fac, 
				intr.db.list[["ligands"]], intr.db.list[["receptors"]])

	score.obj <- methods::new(
		"lrScores",
		lig.slot = lig.slot,
		recep.slot = recep.slot,
		intr.slot = intr.db.name,
		intr.list = intr.db.list,
		score = res.mtx,
		score.null = methods::new("matrix"),
		res.list = list(),
		log = list(
			"Parmaeters:" = c(lig.slot, recep.slot),
			"Date:" = Sys.time()
		)
	)

	tag <- paste0(lig.slot, "-", recep.slot)
	object@lrscore[["default"]] <- tag
	object@lrscore[[tag]] <- score.obj

	return(object)
}


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
graphNicheLR <- function(
  object,
  ...
) {
  UseMethod(generic = 'graphNicheLR', object = object)
}

#' Sub function for graphNicheLR, input a sparse matrix
#' 
#' @param dge.lig A sparse matrix for ligand
#' @param dge.recep A sparse matrix for receptor
#' @param nb.id.fac A factor of neighbor indices
#' @param lig.fac A factor of ligand indices
#' @param recep.fac A factor of receptor indices
#' 
#' @return A sparse matrix
#' @export
#' 
graphNicheLR.dgCMatrix <- function(
	dge.lig,
    dge.recep,
    nb.id.fac,
	lig.fac,
	recep.fac
){

	### cavaet: remember to convert the Uniprot ids to indices!
	# convert nb fac
	nb.id.fac = sort(nb.id.fac)
	nb.index = facToIndex(nb.id.fac)

	# lig.fac = intr.valid[["ligands"]]
	# recep.fac = intr.valid[["receptors"]]

	if (max(as.integer(names(lig.fac))) > nrow(dge.lig)){
		stop("Intr index out of dge bounds.")
	}

	lig.index = facToIndex(lig.fac)
	recep.index = facToIndex(recep.fac)

	# compute scores
	res.mtx = graphNicheLR_cpp(
		unname(as.matrix(dge.lig)),
		unname(as.matrix(dge.recep)),
		lig.index[[1]], lig.index[[2]], 
		recep.index[[1]], recep.index[[2]],
		nb.index[[1]], nb.index[[2]]
	)

	dimnames(res.mtx) = list(colnames(dge.lig)[as.integer(levels(nb.id.fac))], levels(lig.fac))
	# dimnames(res.mtx) = list(colnames(dge.lig), levels(lig.fac))
	# res.mtx = Matrix(res.mtx, sparse = T)

	return(res.mtx)
}


#' Sub function for graphNicheLR, input a CytoSignal object
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
graphNicheLR.CytoSignal <- function(
	object,
	lig.slot,
	recep.slot,
	intr.db.name,
	nn.use = NULL
){
	if (!lig.slot %in% names(object@imputation)){
		stop("Ligand slot not found.")
	}

	if (!recep.slot %in% names(object@imputation)){
		stop("Receptor slot not found.")
	}

	if (!intr.db.name %in% c("diff_dep", "cont_dep")) {
		stop("intr.db.name must be either 'diff_dep' or 'cont_dep'.")
	}

	dge.lig <- object@imputation[[lig.slot]]@imp.data
	dge.recep <- object@imputation[[recep.slot]]@imp.data

	if (!all.equal(dim(dge.lig), dim(dge.recep))){
		stop("dge.lig and dge.recep must have the same dimension.")
	}

	message("Computing scores using ", intr.db.name, " database.")
	message("Ligand: ", lig.slot, ", Receptor: ", recep.slot, ".")

	intr.db.list <- checkIntr(unname(object@intr.valid[["symbols"]][["intr"]]), 
							object@intr.valid[[intr.db.name]])

	# res.mtx <- graphNicheLR.dgCMatrix(dge.lig, dge.recep, nb.id.fac, 
	# 			intr.db.list[["ligands"]], intr.db.list[["receptors"]])

	res.mtx <- dge.lig * dge.recep

	score.obj <- methods::new(
		"lrScores",
		lig.slot = lig.slot,
		recep.slot = recep.slot,
		intr.slot = intr.db.name,
		intr.list = intr.db.list,
		score = t(res.mtx),
		score.null = methods::new("matrix"),
		res.list = list(),
		log = list(
			"Parmaeters:" = c(lig.slot, recep.slot),
			"Date:" = Sys.time()
		)
	)

	tag <- paste0(lig.slot, "-", recep.slot)
	object@lrscore[["default"]] <- tag
	object@lrscore[[tag]] <- score.obj

	return(object)
}



# #' Sub function for changeUniprot, input a cytosignal object
# #' 
# #' @param object A Cytosignal object
# #' @param slot.use The slot to use for subsetting
# #' @param verbose Whether to print out the progress
# #' 
# #' @return A Cytosignal object
# #' @export
# changeUniprot.CytoSignal <- function(
# 	object,
# 	slot.use = NULL, 
# 	# mode = "unpt",
# 	verbose = T
# ){
# 	if (is.null(slot.use)){
# 		slot.use <- object@imputation[["default"]]
# 	}

# 	if (!slot.use %in% names(object@imputation)){
# 		stop("Data not found.")
# 	}

# 	gene_to_uniprot <- object@intr.valid[["gene_to_uniprot"]]
# 	dge.raw <- object@counts
# 	# check duplicate items exists in database

# 	unpt.list <- changeUniprot.matrix_like(dge.raw, gene_to_uniprot, verbose = verbose)
	
# 	object@imputation[[slot.use]]@intr.data <- unpt.list[[1]]
# 	object@intr.valid[["symbols"]][[slot.use]] <- unpt.list[[2]]

# 	scale.fac <- object@imputation[[slot.use]]@scale.fac$raw
# 	scale.fac <- scale.fac[colnames(unpt.list[[1]])]
# 	object@imputation[[slot.use]]@scale.fac$raw <- scale.fac

# 	# if (!all.equal(names(slot(object, slot.use)), names(object@intr.data))) {
# 	# 	warning("Names of imp.norm.data and intr.data not equal!")
# 	# }

# 	return(object)
# }
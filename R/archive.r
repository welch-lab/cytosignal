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
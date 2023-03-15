#' Impute the velocity mtx using the specified method
#' 
#' @param object A Cytosignal object
#' @param method The method to use for imputation
#' @param ... Other parameters
#' 
#' @return A Cytosignal object
#' @export
imputeNicheVelo <- function(
  object,
  ...
) {
  UseMethod(generic = 'imputeNicheVelo', object = object)
}



#' Sub function for imputeNicheVelo, input a Cytosignal object
#' 
#' @param object A Cytosignal object
#' @param nn.type The type of neighbors
#' @param weights The weight of the Delaunay triangulation
#' 
#' @return A Cytosignal object
#' @export
imputeNicheVelo.CytoSignal <- function(
	object, 
	nn.type = NULL, 
	weights = c("mean", "counts", "dist", "none")
){
	# dge.raw <- object@counts
	# Extract neighborhood factors

	if (is.null(nn.type)){
		nn.type <- object@imputation[["default"]]
	}

	if  (!weights %in% c("mean", "counts", "dist", "none")){
		stop("Incorret weighting method!\n")
	}

	message(paste0("Imputing using ", nn.type, "..."))

	dge.raw <- object@velo[["velo.s"]] + object@velo[["velo.u"]]

	if (nn.type %in% names(object@imputation)){
		nb.id.fac <- object@imputation[[nn.type]]@nn.id
		nb.dist.fac <- object@imputation[[nn.type]]@nn.dist
	} else {
		stop("NN type not found.")
	}

	if (nn.type == "Raw") {
		object@imputation[[nn.type]]@imp.velo <- dge.raw
		return(object)
	}

	if (ncol(dge.raw) < length(levels(nb.id.fac))) {
		stop("Number of index beads larger than the number of beads in DGE.")
	}

	nn.graph <- object@imputation[[nn.type]]@nn.graph
	dge.raw.imputed <- dge.raw %*% nn.graph

	# weighted sum of scale factors as well
	scale.fac <- Matrix::Matrix(object@parameters[["velo.lib.size"]], sparse = T, nrow = 1,
						dimnames = list(NULL, names(object@parameters[["velo.lib.size"]])))
	scale.fac.imp <- scale.fac %*% weights.mtx

	scale.fac.imp <- as.numeric(scale.fac.imp)
	names(scale.fac.imp) <- names(object@parameters[["velo.lib.size"]])

	res.density <- sum(dge.raw.imputed != 0)/length(dge.raw.imputed) # density 6.2%
	cat(paste0("Density after imputation: ", res.density*100, "%\n"))

	object@imputation[[nn.type]]@imp.velo <- dge.raw.imputed
	object@imputation[[nn.type]]@scale.fac.velo <- scale.fac.imp
	object@imputation[[nn.type]]@log[["Density:"]] <- paste0(res.density*100, "%")

	return(object)
}





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

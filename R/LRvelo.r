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
#' 
#' @return A Cytosignal object
#' @export
imputeNicheVelo.CytoSignal <- function(
	object, 
	nn.type = NULL
){
	if (is.null(nn.type)){
		nn.type <- object@imputation[["default"]]
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
	scale.fac.imp <- scale.fac %*% nn.graph

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
inferVeloLR <- function(
  object,
  ...
) {
  UseMethod(generic = 'inferVeloLR', object = object)
}

#' Sub function for inferVeloLR, input a sparse matrix
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
inferVeloLR.matrix_like <- function(
	dge.lig,
    dge.recep,
	dge.lig.velo,
	dge.recep.velo,
	lig.fac,
	recep.fac
){
	if (max(as.integer(names(lig.fac))) > nrow(dge.lig)){
		stop("Intr index out of dge bounds.")
	}

    if (max(as.integer(names(lig.fac))) > nrow(dge.lig.velo)){
		stop("Intr index out of dge bounds.")
	}

	lig.index = facToIndex(lig.fac)
	recep.index = facToIndex(recep.fac)

	# compute velos
	res.mtx = inferVeloLR_cpp(
		unname(as.matrix(dge.lig)),
		unname(as.matrix(dge.recep)),
		unname(as.matrix(dge.lig.velo)),
		unname(as.matrix(dge.recep.velo)),
		lig.index[[1]], lig.index[[2]], 
		recep.index[[1]], recep.index[[2]]
	)

	dimnames(res.mtx) = list(colnames(dge.lig), levels(lig.fac))
	res.mtx = Matrix(res.mtx, sparse = T)

	return(res.mtx)
}



#' Sub function for inferVeloLR, input a CytoSignal object
#' 
#' @param object A Cytosignal object
#' @param lig.slot The ligand slot to use
#' @param recep.slot The receptor slot to use
#' @param intr.db.name The intr database name to use
#' @param nn.use slot that the neighbor index should be taken from, by default is the same as
#' 			the recep.slot. For example, if velo.obj = GauEps-DT, then nn.use = "DT".
#' 			nn.use could also be a user-defind factor.
#' @param norm.method The normalization method to apply to the counts, need to be consistent with the
#' 			normalization method used for the RNA velocity method. Default is "scanpy".
#' 
#' @return A Cytosignal object
#' @export
inferVeloLR.CytoSignal <- function(
	object,
	lig.slot,
	recep.slot,
	intr.db.name,
	norm.method = "scanpy",
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

	message("Computing LR-velos using ", intr.db.name, " database.")
	message("Ligand: ", lig.slot, ", Receptor: ", recep.slot, ".")

	if (is.null(tag)) {
		tag <- paste0(lig.slot, "-", recep.slot)
	}

	# if (tag %in% names(object@lrvelo)) {
	# 	stop("Tag already exists.")
	# }

	# normalize using scanpy method, normCount has been revised to internal function
	# dge.lig <- object@imputation[[lig.slot]]@imp.data
	# dge.recep <- object@imputation[[recep.slot]]@imp.data
	dge.lig <- normCounts(object, method = norm.method, slot.use = lig.slot)
	dge.recep <- normCounts(object, method = norm.method, slot.use = recep.slot)

	dge.lig.velo <- object@imputation[[lig.slot]]@imp.velo
	dge.recep.velo <- object@imputation[[recep.slot]]@imp.velo

	# compare the dimnames of all four matrices
	if (!all.equal(dimnames(dge.lig), dimnames(dge.recep))){
		stop("dge.lig and dge.recep must have the same dimension names.")
	}

	if (!all.equal(dimnames(dge.lig.velo), dimnames(dge.recep.velo))){
		stop("dge.lig.velo and dge.recep.velo must have the same dimension names.")
	}

	#----------- pre-computing the lrscores by averaging the DT scores, without norm -----------#
	message("Comfirming niche index...")

	dt.avg.g <- object@imputation[["DT"]]@nn.graph
	dt.avg.g <- to_mean(dt.avg.g)

	dge.lig <- dge.lig %*% dt.avg.g
	dge.recep <- dge.recep %*% dt.avg.g

    dge.lig.velo <- dge.lig.velo %*% dt.avg.g
    dge.recep.velo <- dge.recep.velo %*% dt.avg.g
	#-------------------------------------------------------------------------------------------#

	intr.db.list <- checkIntr(unname(object@intr.valid[["symbols"]][["intr"]]), 
							object@intr.valid[[intr.db.name]])

	message("Calculating LR-velos...")

	res.mtx <- inferVeloLR.matrix_like(dge.lig, dge.recep, 
				dge.lig.velo, dge.recep.velo,
				intr.db.list[["ligands"]], intr.db.list[["receptors"]])

    message("Done!\n")

	lrvelo.obj <- new(
		"lrVelo",
		lig.slot = lig.slot,
		recep.slot = recep.slot,
		intr.slot = intr.db.name,
		intr.list = intr.db.list,
		intr.velo = res.mtx,
		# nn.id = velo.nn.list[["id"]],
		# nn.dist = velo.nn.list[["dist"]],
		log = list(
			"Used slot" = c(lig.slot, recep.slot)
		)
	)

    object@lrvelo[["default"]] <- tag
	object@lrvelo[[tag]] <- lrvelo.obj

	return(object)
}

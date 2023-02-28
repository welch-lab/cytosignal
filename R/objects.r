#sourceCpp("/home/alan/Documents/Welch_lab/pj19-cell_cell_int/cytosignal_v1_02142023/cytosignal/src/mat_exp.cpp")
#sourceCpp("/home/alan/Documents/Welch_lab/pj19-cell_cell_int/cytosignal_v1_02142023/cytosignal/src/utils_velo.cpp")

#' The CytoSignal Class
#'
#' The CytoSignal object is created from one ST dataset. To construct a
#' CytoSignal object, the user needs to provide both DGE and cell locations. The class provides functions for data
#' preprocessing, integrative analysis, and visualization.
#'
#' The key slots used in the CytoSignal object are described below.
#'
#' @slot counts raw dge, gene x cell
#' @slot cell.loc matrix of cell locations
#' @slot clusters cluster assignments
#' @slot imputation list of imputation objs
#' @slot lrscore list of lrScore objs
#' @slot velo list of velo objs
#' @slot intr.valid list of intr database
#' @slot cell.data data frame of cell data
#' @slot parameters List of parameters for gaussian imputation
#' @slot log Log of each main steps
#' @slot version Version of package used to create object
#'
#' @name CytoSignal-class
#' @rdname CytoSignal-class
#' @aliases CytoSignal-class
#' @exportClass CytoSignal


CytoSignal <- setClass(
  Class = "CytoSignal",
  slots = c(
	counts = "dgCMatrix",
	cells.loc = "matrix",
	clusters = "factor",
	imputation = "list",
	lrscore = "list",
	velo = "list",
	intr.valid = "list",
	parameters = "list",
	log = "list",
	version = "character"
  )
)

#' The ImpData Class
#'
#' The ImpData object is created from one ST dataset. User could choose a preferred imputation method
#' and the class stores the imputed data, the imputed normalized data, the intr database, the intr database.
#'
#' The key slots used in the ImpData object are described below.
#'
#' @slot method imputation method
#' @slot imp.data imputed data
#' @slot imp.norm.data imputed normalized data
#' @slot intr.data imputed normalized data subsetted by intr database
#' @slot intr.data.null permuted imputed normalized data subsetted by intr database
#' @slot nn.id nearest neighbor id
#' @slot nn.dist nearest neighbor distance
#' @slot log Log of each main steps
#' 
#' @name ImpData-class
#' @rdname ImpData-class
#' @aliases ImpData-class
#' @exportClass ImpData
 

ImpData <- setClass(
	Class = "ImpData",
	slots = c(
		method = "character",
		imp.data = "dgCMatrix",
		imp.norm.data = "dgCMatrix",
		intr.data = "dgCMatrix",
		intr.data.null = "dgCMatrix",
		nn.id = "factor",
		nn.dist = "factor",
		log = "list"
	)
)



#' The lrScores Class
#' 
#' The lrScores object is created from one ST dataset. User could choose two imputation methods to calculate
#' the ligand-receptor scores. The class stores the ligand, receptor, and interaction database, the ligand-receptor
#' scores, the ligand-receptor scores for permuted data, the ligand-receptor scores for permuted data, and the
#' log of each main steps.
#' 
#' The key slots used in the lrScores object are described below.
#' 
#' @slot lig.slot ligand database
#' @slot recep.slot receptor database
#' @slot intr.slot interaction database
#' @slot intr.list list of interaction database
#' @slot score ligand-receptor scores
#' @slot score.null permuted ligand-receptor scores
#' @slot res.list list of results
#' @slot log Log of each main steps
#' 
#' @name lrScores-class
#' @rdname lrScores-class
#' @aliases lrScores-class
#' @exportClass lrScores


lrScores <- setClass(
	Class = "lrScores",
	slots = c(
		lig.slot = "character",
		recep.slot = "character",
		intr.slot = "character",
		intr.list = "list",
		score = "matrix",
		score.null = "matrix",
		res.list = "list",
		log = "list"
	)
)


#' show method for CytoSignal
#'
#' @param object CytoSignal object
#' @export

setMethod(
  f = "show",
  signature = "CytoSignal",
  definition = function(object) {
	cat(
	  "An object of class",
	  class(object),
	  "\nwith raw data of dimension\n",
	  ncol(object@counts),
	  "cells and ",
	  nrow(object@counts),
	  "genes.\n",
	  "Head of the raw data:\n"
	)
	print(object@counts[1:5,1:5])
	invisible(x = NULL)
  }
)



#' show method for ImpData
#' 
#' @param object CytoSignal object
#' @param slot.use slot to use
#' @export 

showImp <- function(object, slot.use = NULL) {
	 if (is.null(slot.use)){
		 slot.use <- object@imputation[["default"]]
	 }

	 cat(
		 "Default: imputed data using: \n",
		 slot.use, "\n",
		 "Processing Log: \n"
	 )

	 print(object@imputation[[slot.use]]@log)

	 if (ncol(object@imputation[[slot.use]]@imp.data) > 0){
		 print(object@imputation[[slot.use]]@imp.data[1:5,1:5])
	 } else {
		 cat("No imputed data available.\n")
	 }
}

#' show intr.data in ImpData
#' 
#' @param object CytoSignal object
#' @param slot.use slot to use
#' @export
#' 
showUnpt <- function(object, slot.use = NULL){
	 if (is.null(slot.use)){
		 slot.use <- object@imputation[["default"]]
	 }

	 symbols = object@intr.valid[["symbols"]][[slot.use]]

	 cat(
		 "Default: Filtered data using: \n",
		 slot.use, "\n",
		 "Number of genes in the database: ", length(symbols), "\n"
	 )

	 if (ncol(object@imputation[[slot.use]]@intr.data) > 0){
		 print(object@imputation[[slot.use]]@intr.data[1:5,1:5])
	 } else {
		 cat("No imputed data available.\n")
	 }
}


#' show method for lrScores
#' 
#' @param object CytoSignal object
#' @param slot.use slot to use
#' 
#' @export
showScore <- function(object, slot.use = NULL){
	 if (is.null(slot.use)){
		 slot.use <- object@lrscore[["default"]]
	 }

	 cat(
		 "Default: Score data using: \n",
		 slot.use, "\n",
		 "Number of intrs in the mtx:", ncol(object@lrscore[[slot.use]]@score), "\n",
		 "Processing Log: \n"
	 )

	 print(object@lrscore[[slot.use]]@log)

	 if (ncol(object@lrscore[[slot.use]]@score) > 0){
		 print(object@lrscore[[slot.use]]@score[1:5,1:5])
	 } else {
		 cat("No score data available.\n")
	 }
}

#' show all current logs
#' 
#' @param object CytoSignal object
#' 
#' @export
showLog <- function(object){
	cat("Pre-processing Log: \n")
	print(object@log)

	cat("Ligand imputation Log: \n")
	print(object@imputation[[1]]@log)

	cat("Receptor imputation Log: \n")
	print(object@imputation[[3]]@log)

	cat("Score Log: \n")
	print(object@lrscore[[2]]@log)
}

#' show method for cytosignal obj
#' 
setMethod(
  f = "show",
  signature = "ImpData",
  definition = function(object) {
	cat(
		"Imputed data using: \n",
		object@method, "\n",
		"Parameters: \n",
		object@log[[1]],
		"\nNeighbors: \n",
		object@log[[2]],
		"\n"
	)
	print(object@imp.data[1:5,1:5])
	invisible(x = NULL)
  }
)

#' Create a CytoSignal object
#' 
#' @param raw.data raw data matrix
#' @param cells.loc cells location matrix
#' @param clusters cluster information
#' @param name name of the dataset
#' @param parameters parameters used
#' @param log log of the processing
#' @param version version of the package
#' 
#' @return a CytoSignal object
#' 
#' @export
#' 
createCytoSignal <- function(
	raw.data,
	cells.loc,
	clusters = NULL,
	name = NULL,
	parameters = NULL,
	log = NULL,
	version = NULL) {

	if (is.null(clusters)) {
		clusters <- factor(rep(1, ncol(raw.data)))
	}
	if (is.null(parameters)) {
		parameters <- list()
	}
	if (is.null(log)) {
		log <- list()
		log[["Dataset"]] <- name
	}
	if (is.null(version)) {
		version <- "v1.0.0"
	}

	# check the class of raw.data
	if (!inherits(raw.data, "CsparseMatrix")) {
		raw.data <- as(raw.data, "CsparseMatrix")
	}

	cells.loc <- as.matrix(cells.loc)[colnames(raw.data), ]
	if (!all.equal(colnames(raw.data), rownames(cells.loc))) {
		stop("The cell names in raw.data and cells.loc are not the same.")
	}

	clusters <- clusters[colnames(raw.data), drop = T]

	object <- methods::new(
		"CytoSignal",
		counts = raw.data,
		cells.loc = cells.loc,
		clusters = clusters,
		imputation = list(),
		lrscore = list(),
		velo = list(),
		intr.valid = list(),
		parameters = list(),
		log = list(),
		version = version
	)

	return(object)
}

#' Add interaction database to CytoSignal object
#' 
#' @param object CytoSignal object
#' @param gene_to_uniprot df, gene to uniprot symbols mapping
#' @param intr.db.diff_dep intr.db for diffusable ligands + non-diffusable receptors
#' @param intr.db.cont_dep intr.db for contact dependent ligands + receptors
#' @param inter.index df, collection of intraction information
#' 
#' @return a CytoSignal object
#' 
#' @export
addIntrDB <- function(
	object,
	gene_to_uniprot,
	intr.db.diff_dep,
	intr.db.cont_dep,
	inter.index
) {
	object@intr.valid <- list(
		"gene_to_uniprot" = gene_to_uniprot,
		"diff_dep" = intr.db.diff_dep,
		"cont_dep" = intr.db.cont_dep,
		"intr.index" = inter.index
	)

	return(object)
}

#' Remove low quality cells and genes from raw counts
#' 
#' @param object CytoSignal object
#' @param counts.thresh threshold for cell counts
#' @param gene.thresh threshold for gene counts
#' 
#' @return a CytoSignal object
#' 
#' @export
removeLowQuality <- function(object, counts.thresh = 300, gene.thresh = 50) {
	dge.raw.filter <- object@counts
	cells.loc <- object@cells.loc

	del.genes = filterGene(dge.raw.filter, object@intr.valid[["gene_to_uniprot"]], gene.thresh)
	del.cells = which(colSums(dge.raw.filter) < counts.thresh)

	# Removed  13017 / 33611  low quality cells.
	if (length(del.cells) == 0 & length(del.genes) != 0){
		text1 = "No cells removed."
		text2 = paste0("Removed ", length(del.genes), "/", nrow(dge.raw.filter), " low quality genes.")
		dge.raw.filter = dge.raw.filter[-del.genes, ]

	} else if (length(del.genes) == 0 & length(del.cells) != 0){
		text1 = "No genes removed."
		# cat("No genes removed.\n")
		text2 = paste0("Removed ", length(del.cells), "/", ncol(dge.raw.filter), " low quality cells.")
		# cat("Removed ", length(del.cells), "/", ncol(dge.raw.filter), " low quality cells.\n")
		dge.raw.filter = dge.raw.filter[, -del.cells]
		cells.loc = cells.loc[-del.cells,]

	} else if (length(del.genes) != 0 & length(del.cells) != 0){
		text1 = paste0("Removed ", length(del.genes), "/", nrow(dge.raw.filter), " low quality genes.")
		text2 = paste0("Removed ", length(del.cells), "/", ncol(dge.raw.filter), " low quality cells.")
		# cat("Removed ", length(del.genes), "/", nrow(dge.raw.filter), " low quality genes.\n")
		# cat("Removed ", length(del.cells), "/", ncol(dge.raw.filter), " low quality cells.\n")
		dge.raw.filter = dge.raw.filter[-del.genes, -del.cells]
		cells.loc = cells.loc[-del.cells,]
	} else {
		text1 = "No cells removed."
		text2 = "No genes removed."
		# cat("No genes or cells removed.\n")
	}

	if (!all.equal(rownames(cells.loc), colnames(dge.raw.filter))){
		stop("The cell names in raw data and cells locations are not the same.")
	}

	object@clusters <- object@clusters[colnames(dge.raw.filter), drop = T]
	object@counts <- dge.raw.filter
	object@cells.loc <- cells.loc

	log_text = c(
		paste0("Gene count threshold: ", gene.thresh),
		text1,
		paste0("Cell count threshold: ", counts.thresh),
		text2
	)

	object@log[["Step_1"]] <- log_text

	cat(log_text, sep = "\n")

	return(object)
}

#' Remove imputed data and normalized imputed data from CytoSignal object to save disk space
#' 
#' @param object CytoSignal object
#' 
#' @return a CytoSignal object
#' 
#' @export
purgeBeforeSave <- function(object) {
	# get slot names
	imp.names = setdiff(names(object@imputation), "default")

	for (i in 1:length(imp.names)) {
		object@imputation[[imp.names[i]]]@imp.data <- new("dgCMatrix")
		object@imputation[[imp.names[i]]]@imp.norm.data <- new("dgCMatrix")
	}

	return(object)
}

#' Suggest cell intervels for scaling
#' 
#' This function is used to estimate the actual cell intervel in cells.loc after scaling
#' returns the top 5 least intervels, users can choose to which to use.
#' 
#' @param object CytoSignal object
#' 
#' @return a vector of intervels
#' 
#' @export
suggestInterval <- function(object) {
	nn = RANN::nn2(object@cells.loc, object@cells.loc, k = 6, searchtype = "priority")
	# nn.idx <- t(nn[["nn.idx"]]) # k X N
	nn.dist <- t(nn[["nn.dists"]]) # k X N
	nn.dist <- nn.dist[-1, ] # remove the first row (self)

	# find the minimum within each row
	nn.dist.min <- apply(nn.dist, 1, min)
	nn.dist.min <- nn.dist.min[order(nn.dist.min)]

	return(nn.dist.min)
}

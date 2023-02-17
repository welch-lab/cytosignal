#' The CytoSignal Class
#'
#' The CytoSignal object is created from one ST dataset. To construct a
#' CytoSignal object, the user needs to provide both DGE and cell locations. The class provides functions for data
#' preprocessing, integrative analysis, and visualization.
#'
#' The key slots used in the CytoSignal object are described below.
#'
#' @slot raw.data raw dge, gene x cell
#' @slot imp.data normalized dge
#' @slot intr.data normed dge of only intr genes
#' @slot cell.loc matrix of cell locations
#' @slot clusters cluster assignments
#' @slot cell.data data frame of cell data
#' @slot parameters List of parameters for gaussian imputation
#' @slot log Log of each main steps
#' @slot version Version of package used to create object
#'
#' @name CytoSignal-class
#' @rdname CytoSignal-class
#' @aliases CytoSignal-class
#' @exportClass CytoSignal
#' @importFrom Rcpp evalCpp
#' @useDynLib rCytoSignal
#' 
#' 

source("/home/alan/Documents/Welch_lab/pj19-cell_cell_int/cytosignal_v1_02142023/master_pipelines/sub_utils.r")

# CytoSignal <- setClass(
#   Class = "CytoSignal",
#   slots = c(
# 	raw.data = "CsparseMatrix",
# 	imp.data = "list",
# 	imp.norm.data = "list",
# 	intr.data = "list",
# 	lrscore = "list",
# 	lrscore.null = "list",
# 	res.list = "list",
# 	cells.loc = "matrix",
# 	velo.data = "list"
# 	clusters = "factor",
# 	intr.db.valid = "list",
# 	nn = "list",
# 	parameters = "list",
# 	log = "list",
# 	version = "character"
#   )
# )


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


ImpData <- setClass(
	Class = "ImpData",
	slots = c(
		method = "character",
		imp.data = "dgCMatrix",
		imp.norm.data = "dgCMatrix",
		intr.data = "dgCMatrix",
		nn.id = "factor",
		nn.dist = "factor",
		log = "list"
	)
)


lrScores <- setClass(
	Class = "lrScores",
	slots = c(
		lig.slot = "character",
		recep.slot = "character",
		score = "matrix",
		score.null = "matrix",
		res.list = "list",
		log = "list"
	)
)


#' show method for CytoSignal
#'
#' @param object CytoSignal object
#' @name show
#' @aliases show,CytoSignal-method
#' @docType methods
#' @rdname show-methods

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


removeLowQuality <- function(object, intr.df, counts.thresh = 300, gene.thresh = 50) {
	dge.raw.filter <- object@counts
	cells.loc <- object@cells.loc

	del.genes = filterGene(dge.raw.filter, intr.df, gene.thresh)
	del.cells = which(colSums(dge.raw.filter) < counts.thresh)

	# Removed  13017 / 33611  low quality cells.
	if (length(del.cells) == 0 & length(del.genes) != 0){
		test1 = "No cells removed."
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
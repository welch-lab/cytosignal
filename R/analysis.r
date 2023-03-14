#' Infer the parameters of the Gaussian kernel
#' 
#' @param object A Cytosignal object
#' @param tech.res The resolution of the technology
#' @param bin_size The size of the bin
#' @param loc.d The scaled distance between two cells
#' @param r.eps.real The radius of the epsilon ball in tech resolution in um, default 200 um
#' @param thresh The total signal out of the epsilon ball
#' 
#' @return A Cytosignal object
#' @export

inferEpsParams <- function(object, tech.res, bin_size, loc.d, r.eps.real = 200, thresh = 0.001){
	# tech.res: distance between two cells, ori res * number of bins
	# loc.d: SCALED distance between two cells
	# r.eps.real: the radius of the epsilon ball in tech resolution in um, default 200 um
	# thresh: total signal out of the epsilon ball
	tech.res = tech.res * bin_size
	# loc.d = object@cells.loc[1,1] - cells.loc[2,1]
	r.eps.scale = (loc.d*r.eps.real)/tech.res # the radius of the epsilon ball in scaled resolution
	sigma.real = r.eps.real / sqrt(-2 * log(thresh)) # sigma in tech resolution
	sigma.scale = (loc.d*sigma.real)/tech.res # sigma in scaled resolution
	# object@parameters <- list(r.diffuse.scale = r.eps.scale, sigma.scale = sigma.scale)
	object@parameters[["r.diffuse.scale"]] <- r.eps.scale
	object@parameters[["sigma.scale"]] <- sigma.scale

	return(object)
}


#' Find the neighbors of each cell in the epsilon ball
#' 
#' @param object A Cytosignal object
#' @param eps The radius of the epsilon ball
#' @param sigma The sigma of the Gaussian kernel
#' @param self.weight weight of the index cell
#' 
#' @return A list of neighbors and their distances
#' @export
#' 
findNNGauEB <- function(
  object,
  ...
) {
  UseMethod(generic = 'findNNGauEB', object = object)
}


#' Sub function for findNNGauEB, input a sparse matrix
#' 
#' @param object A Cytosignal object
#' @param eps The radius of the epsilon ball
#' @param sigma The sigma of the Gaussian kernel
#' @param self.weight weight of the index cell
#' 
#' @return A Cytosignal object
#' @export
findNNGauEB.matrix <- function(
	cells.loc,
	eps,
	sigma = 0.15,
	self.weight = "auto"
){
	# cat("Finding neighbors in epsilon circle...\n")
	nn <- dbscan::frNN(cells.loc, eps = eps, sort = F)

	if (is.numeric(self.weight)) {
		if (self.weight < 0 || self.weight > 1) {
			stop("Self weight should be within (0,1)!\n")
		}
		message("Using manual self weight: ", self.weight, "...")
	} else if (self.weight == "auto") {
		message("Determining self weight automatically...")
	} else if (self.weight == "sum_1") {
		message("Using self weight: all NB weights sum to 1...")
	} else {
		message("Unknown self weight, determining automatically...")
		self.weight = "auto"
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
		else {
			stop("Unknown self weight!\n")
		}
	})

	# discard the cell that has no neighbors
	rm.index = which(lengths(dist.list.gau) == 1)
	if (length(rm.index) > 0) dist.list.gau = dist.list.gau[-rm.index]
	dist.list.gau = as.numeric(unlist(dist.list.gau))
	names(nn.dist.fac) = dist.list.gau

	# nn.dist.fac <- addIndexOne(nn.dist.fac)

	num.diff = nrow(cells.loc) - length(levels(nn.fac))
	if (num.diff > 0){
		cat("A total of ", num.diff, " beads do not have NNs.\n")
	} else if (num.diff < 0) {
	   stop("Result fac longer than original beads length!\n")
	}

	cat(paste0("Mean num of neighbors: ", ceiling(mean(table(nn.fac))), "\n"))
	cat(paste0("Median num of neighbors: ", median(table(nn.fac)), "\n"))

	return(list(
		"id" = nn.fac,
		"dist" = nn.dist.fac
	))
}


#' Sub function for findNNGauEB, input a Cytosignal object
#' 
#' @param object A Cytosignal object
#' @param eps The radius of the epsilon ball
#' @param sigma The sigma of the Gaussian kernel
#' @param self.weight weight of the index cell
#' 
#' @return A Cytosignal object
#' @export
findNNGauEB.CytoSignal <- function(
	object,
	eps = NULL,
	sigma = NULL,
	self.weight = "auto",
	tag = NULL
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
#' 
#' @param object A Cytosignal object
#' @param weight.sum The sum of the weights
#' 
#' @return A Cytosignal object
#' @export
findNNDT <- function(
  object,
  ...
) {
  UseMethod(generic = 'findNNDT', object = object)
}

#' Sub function for findNNDT, input a matrix
#' 
#' @param cells.loc A matrix of cells location
#' @param weight.sum The sum of the weights
#' @importFrom RTriangle triangulate pslg
#' @return A list of neighbors
#' @export
findNNDT.matrix <- function(
	cells.loc,
	weight.sum = 2
){
	cat("Finding DT neighbors...\n")
	# cells.loc <- pslg(P=cells.loc)
	cells.dt <- triangulate(pslg(P=cells.loc))
	egdes <- cells.dt$E

	#### construct a double-corresponded factor
	nn.fac <- factor(c(egdes[, 1], egdes[, 2])) # neighbors.factor
	names(nn.fac) <- c(egdes[, 2], egdes[, 1])
	nn.fac <- sort(nn.fac)
	
	index.loc = cells.loc[as.integer(nn.fac),]
	nb.loc = cells.loc[as.integer(names(nn.fac)),]
	
	# set distance for each niche, based on the cell number of each niche
	nb.size <- table(nn.fac)
	dist.list = lapply(levels(nn.fac), function(x){
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

	cat(paste0("Mean num of neighbors: ", ceiling(mean(nb.size)), "\n"))
	cat(paste0("Median num of neighbors: ", median(nb.size), "\n"))

	# return(nn.fac)
	return(list(
		"id" = nn.fac,
		"dist" = nn.dist.fac
	))
}


#' Sub function for findNNDT, input a Cytosignal object
#' 
#' @param object A Cytosignal object
#' @param weight The weight of the Delaunay triangulation
#' 
#' @return A Cytosignal object
#' @export
findNNDT.CytoSignal <- function(
	object,
	weight = 2
){
	tag <- "DT"

	# if (tag %in% names(object@imputation)){
	# 	stop("This imputation has been done before.")
	# }

	message("Finding neighbors using DT with tag: ", tag, "...")

	cells.loc <- object@cells.loc
	nn <- findNNDT.matrix(cells.loc, weight.sum = weight)

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

	message("Setting Imp obj using NO imputation...")

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


#' Impute the data using the specified method
#' 
#' @param object A Cytosignal object
#' @param method The method to use for imputation
#' @param ... Other parameters
#' 
#' @return A Cytosignal object
#' @export
imputeNiche <- function(
  object,
  ...
) {
  UseMethod(generic = 'imputeNiche', object = object)
}


#' Sub function for imputeNiche, input a sparse matrix
#' 
#' @param dge.raw A sparse matrix
#' @param nb.id.fac A factor of neighbors
#' @param nb.dist.fac A factor of weights
#' @param weights The weight of the Delaunay triangulation
#' 
#' @return A sparse matrix
#' @export
imputeNiche.dgCMatrix <- function(
	dge.raw, 
	nb.id.fac,
	nb.dist.fac,
	weights = c("mean", "counts", "dist", "none"),
	return.graph = F
){
	cat("Imputing...\n")

	i.list = as.integer(names(nb.id.fac)) # x coords, row number
	j.list = as.numeric(as.character(nb.id.fac)) # y coords, col number, hq beads

	if (weights == "avg") { # just do the mean
		x.list = rep(1, length(i.list))
		do.norm <- T
	} else if (weights == "counts") { # use total counts in each bead as weights
		x.list = Matrix::colSums(dge.raw)[i.list]
		do.norm <- T
	} else if (weights == "dist") {
		x.list = exp(-as.numeric(names(nb.dist.fac))) # exp(-dist) as similarities
		do.norm <- T
	} else if (weights == "none") {
		x.list = as.numeric(names(nb.dist.fac)) # gaussian distances has been calculated
		do.norm <- F
	} else{
		stop("Incorret weighting method!\n")
	}

	weights.mtx = sparseMatrix(i = i.list, j = j.list, x = x.list, 
		dims = c(ncol(dge.raw), ncol(dge.raw)), dimnames = list(NULL, colnames(dge.raw)))

	# re-weight the weight mtx across each niche (column)
	if (do.norm) {
		weights.mtx@x = weights.mtx@x / rep.int(Matrix::colSums(weights.mtx), diff(weights.mtx@p))
	}

	dge.raw.imputed = dge.raw %*% weights.mtx

	if (return.graph){
		return(list(dge.raw.imputed, weights.mtx))
	}

	return(dge.raw.imputed)
}



#' Sub function for imputeNiche, input a Cytosignal object
#' 
#' @param object A Cytosignal object
#' @param nn.type The type of neighbors
#' @param weights The weight of the Delaunay triangulation
#' 
#' @return A Cytosignal object
#' @export
imputeNiche.CytoSignal <- function(object, 
	nn.type = NULL, 
	weights = "none"
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

	dge.raw <- object@counts

	if (nn.type %in% names(object@imputation)){
		nb.id.fac <- object@imputation[[nn.type]]@nn.id
		nb.dist.fac <- object@imputation[[nn.type]]@nn.dist
	} else {
		stop("NN type not found.")
	}

	if (ncol(dge.raw) < length(levels(nb.id.fac))) {
		stop("Number of index beads larger than the number of beads in DGE.")
	}

	imp.list <- imputeNiche.dgCMatrix(
		dge.raw, 
		nb.id.fac,
		nb.dist.fac,
		weights = weights,
		return.graph = T
	)

	dge.raw.imputed <- imp.list[[1]]
	weights.mtx <- imp.list[[2]]
	rm(imp.list); gc()

	# weighted sum of scale factors as well
	scale.fac <- Matrix(object@parameters$lib.size, sparse = T, nrow = 1,
					dimnames = list(NULL, colnames(dge.raw)))
	scale.fac.imp <- scale.fac %*% weights.mtx

	scale.fac.imp <- as.numeric(scale.fac.imp)
	names(scale.fac.imp) <- names(object@parameters$lib.size)
	
	res.density <- sum(dge.raw.imputed != 0)/length(dge.raw.imputed) # density 6.2%
	cat(paste0("Density after imputation: ", res.density*100, "%\n"))

	object@imputation[[nn.type]]@imp.data <- dge.raw.imputed
	object@imputation[[nn.type]]@nn.graph <- weights.mtx
	object@imputation[[nn.type]]@scale.fac <- scale.fac.imp
	object@imputation[[nn.type]]@log[["Density:"]] <- paste0(res.density*100, "%")

	return(object)
}



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



#' Normalize the data using the specified method
#' 
#' @param object A Cytosignal object
#' @param method The method to use for normalization
#' @param ... Other parameters
#' 
#' @return A Cytosignal object
#' @export
#' 
normCounts <- function(
  object,
  ...
) {
  UseMethod(generic = 'normCounts', object = object)
}


#' Sub function for normCounts, input a sparse matrix
#' 
#' @param mat A sparse matrix
#' @param method The method to use for normalization
#' 
#' @return A sparse matrix
#' @export
normCounts.dgCMatrix <- function(
	mat, 
	method = c("default", "cpm")
){
	if (method == "default"){
		mat@x <- mat@x / rep.int(Matrix::colSums(mat), diff(mat@p))
		mat@x <- log1p(mat@x * 1e4)
	} else if (method == "cpm"){
		mat@x <- mat@x / rep.int(Matrix::colSums(mat), diff(mat@p))
		mat@x <- log1p(mat@x * 1e6)
	} else {
		stop("Method not found.")
	}

	return(mat)
}


#' Sub function for normCounts, input a sparse matrix
#' 
#' @param mat.list list of (mat, scale.fac)
#' @param method The method to use for normalization
#' 
#' @return A sparse matrix
#' @export
normCounts.list <- function(
	mat.list, 
	method = c("default", "cpm")
){
	mat <- mat.list[["mat"]]
	scale.fac <- mat.list[["scale.fac"]]

	if (method == "default"){
		mat@x <- mat@x / rep.int(scale.fac, diff(mat@p))
		mat@x <- log1p(mat@x * 1e4)
	} else if (method == "cpm"){
		mat@x <- mat@x / rep.int(scale.fac, diff(mat@p))
		mat@x <- log1p(mat@x * 1e6)
	} else {
		stop("Method not found.")
	}

	return(mat)
}


#' Sub function for normCounts, input a Cytosignal object
#' 
#' @param object A Cytosignal object
#' @param method The method to use for normalization
#' @param slot.use The slot to use for normalization
#' @param velo Whether to normalize the velocity data instead
#' 
#' @return A Cytosignal object
#' @export
normCounts.CytoSignal <- function(
	object, 
	method = c("default", "cpm"),
	slot.use = NULL,
	velo = FALSE
){
	if (is.null(slot.use)){
		slot.use <- object@imputation[["default"]]
	}

	if (!slot.use %in% names(object@imputation)) {
		stop("Data not found.")
	}

	message(paste0("Normalizing on Imp slot: ", slot.use, "..."))

	if (velo){
		mat <- object@imputation[[slot.use]]@imp.velo
		scale.fac <- object@imputation[[slot.use]]@scale.fac.velo
	} else {
		mat <- object@imputation[[slot.use]]@imp.data
		scale.fac <- object@imputation[[slot.use]]@scale.fac
	}

	mat <- normCounts.list(
		list(mat = mat, scale.fac = scale.fac), 
		method = method)

	if (velo){
		object@imputation[[slot.use]]@imp.velo <- mat
	} else {
		object@imputation[[slot.use]]@imp.data <- mat
	}

	# if (!all.equal(names(slot(object, slot.use)), names(object@imp.norm.data))){
	# 	warning("Names of imp.data and imp.norm.data not equal!")
	# }
	
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
	intr.db.name
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

	res.mtx <- inferScoreLR.dgCMatrix(dge.lig, dge.recep,
				intr.db.list[["ligands"]], intr.db.list[["receptors"]])

	# res.mtx <- dge.lig * dge.recep

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

	tag <- paste0(lig.slot, "-", recep.slot)
	object@lrscore[["default"]] <- tag
	object@lrscore[[tag]] <- score.obj

	return(object)
}



#' Permute Imputation Results of specific imputation method
#' 
#' This function is a follow-up function of imputeNiche, which permutes the default or user-defined
#' imputation method and stored the results in the ImpData object.
#' 
#' @param object A Cytosignal object
#' @param nn.type The imputation method to use
#' @param perm.size Size of the permutation test, by defualt NULL
#' @param times Number of times to permute the whole dataset, by default 10
#' 
#' @return A Cytosignal object
#' @export
permuteImpLR <- function(
	object,
	nn.type = NULL,
	perm.size = NULL,
	times = 10
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

	message("Permuting imputation data on method ", nn.type, "...")

	dimnames.list <- dimnames(object@imputation[[nn.type]]@imp.data)
	cells.loc <- object@cells.loc
	dge.raw <- object@counts

	if (!is.null(times)) {
		if (!is.numeric(times)){
			stop("times must be a numeric value.")
		}

		times <- ceiling(times)

		if (times < 1) {
			stop("times must be larger than 1.")
		}

	} else {
		# decide how many times to permute according to the size of the data
		if (is.null(perm.size)){
			message("Permutation size not specified, using 100000 instead.")
			perm.size <- 100000
		}

		if (!is.numeric(perm.size)){
			stop("perm.size must be a numeric value.")
		}

		if (perm.size < ncol(dge.raw)){
			message("Permutation size too small, using ", ncol(dge.raw), " instead.")
			perm.size <- ncol(dge.raw)
		}
		
		times <- ceiling(perm.size / ncol(dge.raw))
	}

	message("Permuting whole dataset ", times, " times...")

	null.cells.loc = cells.loc[sample(nrow(cells.loc)), ]
	rownames(null.cells.loc) = rownames(cells.loc)

	# get a random graph for permuting
	if (use.gau) {
		param.eps <- object@parameters$r.diffuse.scale
		param.sigma <- object@parameters$sigma.scale

		null.eps.nb <- findNNGauEB.matrix(null.cells.loc, eps = param.eps, 
						sigma = param.sigma, weight.sum = 2); gc()

		null.graph <- imputeNiche.dgCMatrix(dge.raw, null.eps.nb$id, null.eps.nb$dist, 
							weights = "none", return.graph = T); gc()
	} else if (use.dt) {
		null.nb.fac <- findNNDT.matrix(null.cells.loc)
		null.graph <- imputeNiche.dgCMatrix(dge.raw, null.nb.fac$id, null.nb.fac$dist, 
							weights = "none", return.graph = T); gc()
	} else if (use.raw) {
		cat("No need to generate a random graph for raw imputation.\n")
	} else {
		stop("Cannot find corresponding imputation method.")
	}

	cat("Permuting NULL expressions...\nTimes No.")

	null.dge.list <- lapply(1:times, function(i){
		## MUST shuffule raw dge!!
		cat(i, ", ", sep = "")
		perm.index <- sample(ncol(dge.raw))
		null.dge.raw <- dge.raw[, perm.index]
		scale.fac <- object@parameters[["lib.size"]][perm.index]
		
		colnames(null.dge.raw) <- colnames(dge.raw)
		names(scale.fac) <- colnames(dge.raw)

		scale.fac <- Matrix::Matrix(scale.fac, nrow = 1, byrow = T, sparse = T) # computing scale factor
		# null.dge.raw <- changeUniprot.matrix_like(null.dge.raw.all, object@intr.valid[["gene_to_uniprot"]])[[1]]

		if (use.raw) {
			null.dge.eps <- null.dge.raw
			scale.fac.imp <- scale.fac
		} else {
			null.dge.eps <- null.dge.raw %*% null.graph; gc()
			scale.fac.imp <- scale.fac %*% null.graph; gc()
		}
		
		null.dge.eps <- normCounts.list(list(mat=null.dge.eps, scale.fac=as.numeric(scale.fac.imp)), "default"); gc()
		rm(null.dge.raw); gc()

		return(null.dge.eps)
	})

	cat("End.\nFinished!\n")


	null.dge.eps.unpt = meanMat_cpp(null.dge.list, nrow(null.dge.list[[1]]), ncol(null.dge.list[[1]]))
	dimnames(null.dge.eps.unpt) <- dimnames.list
	rm(null.dge.list); gc()

    object@imputation[[nn.type]]@imp.data.null <- null.dge.eps.unpt

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
permuteScoreLR <- function(
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

    # perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
    # null.cells.loc = cells.loc[perm.index,]
    # rownames(null.cells.loc) = rownames(cells.loc)
    # null.nb.fac = findNNDT.matrix(null.cells.loc); gc() # Mean num of neighbors: 44, median: 36

	### test permuting NULL mtx again
	# null.dge.eps.unpt = null.dge.eps.unpt[, sample(ncol(null.dge.eps.unpt))]
	# null.dge.dt.unpt = null.dge.dt.unpt[, sample(ncol(null.dge.dt.unpt))]

    null.lrscore.mtx = inferScoreLR.dgCMatrix(null.dge.eps.unpt, null.dge.dt.unpt,
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



#' Infer significance of LR scores
#' 
#' @param object A Cytosignal object
#' @param lrscore.mtx A matrix of LR scores
#' @param null.lrscore.mtx A matrix of NULL LR scores
#' @param nb.fac A factor of nearest neighbors index
#' @param intr.db A interaction database
#' @param gene_to_uniprot A dataframe of gene to uniprot mapping
#' @param p.thresh A numeric value of p-value threshold
#' @param reads.thresh A numeric value of reads threshold
#' @param sig.thresh A numeric value of significance threshold
#' @param ... Other arguments
#' 
#' @return A list of indexes of significant cells
#' @export
#' 
inferSignif <- function(
    object,
    ...
) {
    UseMethod(generic = 'inferSignif', object = object)
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
#' 
inferSignif.matrix_like <- function(
    dge.raw,
    lrscore.mtx,
    null.lrscore.mtx,
    nb.fac,
	intr.db,
	gene_to_uniprot,
	p.thresh = 0.05,
	reads.thresh = 100,
	sig.thresh = 100
){
    pval.mtx = getPvalues(lrscore.mtx, null.lrscore.mtx)
    pval.spatial = graphSpatialFDR(nb.fac, pval.mtx)
    # pval.spatial = graphSpatialFDR(nb.fac, pval.mtx, null.lrscore.mtx, reciprocal = F)
    dimnames(pval.spatial) = dimnames(lrscore.mtx); gc()

    res.list = lapply(colnames(pval.spatial), function(cp){
        bead.index = names(which(pval.spatial[, cp] < p.thresh))
    })

    names(res.list) = colnames(pval.spatial)

    num.sig.intr = sum(lengths(res.list) != 0) # 534 / 744
    cat("Number of interactions that have significant i-niche: ", num.sig.intr, "\n")

    res.list = res.list[lengths(res.list) != 0]
    res.list = res.list[order(lengths(res.list), decreasing = T)]

    res.list.hq = filterRes(dge.raw, res.list, intr.db, gene_to_uniprot,
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


#' Sub function for inferSignif, input is a Cytosignal object
#' 
#' @param object A Cytosignal object
#' @param p.value A numeric value of p-value threshold
#' @param reads.thresh A numeric value of reads threshold
#' @param sig.thresh A numeric value of significance threshold
#' @param slot.use A string of slot name to use
#' @param nn.use The lrscore obj to infer from
#' 
#' @return A Cytosignal object
#' @export
inferSignif.CytoSignal <- function(
    object,
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
	} 
	# else if (is.factor(nn.use)) {
	# 	if (length(nn.use) != ncol(object@imputation[[nn.use]]@intr.data))
	# 		stop("nn.use must have the same length as the number of cells.")
	# 	nb.id.fac <- nn.use
	# } 
	else {
		stop("nn.use must be either a factor or a character.")
	}

	message("Inferring significant beads on Score slot ", slot.use, "... ")

    # lrscore.mtx = object@lrscore[[slot.use]]@score
    # null.lrscore.mtx = object@lrscore[[slot.use]]@score.null

    nb.fac = list(
        id = object@imputation[[nn.use]]@nn.id,
        dist = object@imputation[[nn.use]]@nn.dist
    )

	use.intr.slot.name <- object@lrscore[[slot.use]]@intr.slot
	use.intr.db <- object@intr.valid[[use.intr.slot.name]]

    res.list = inferSignif.matrix_like(object@counts, object@lrscore[[slot.use]]@score, 
				object@lrscore[[slot.use]]@score.null, nb.fac,
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
#' @param score.slot The slot name of the lrscore obj to use
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
	lt.mtx.imp <- dataImpKNN(as.matrix(lr.mtx), object@cells.loc, k = k, weight = weight)
	lt.mtx.imp <- as.matrix(Matrix::t(lt.mtx.imp))

	intr.order <- as.double(pearson_col_cpp(lt.mtx.imp, lr.mtx)) * as.double(stdMat_cpp(lt.mtx.imp)) / as.double(stdMat_cpp(lr.mtx))
    names(intr.order) <- colnames(lr.mtx)
	intr.order <- intr.order[order(intr.order, decreasing = T)]

	cat("Reordering significant interactions... \n")
		
	res.list.pear <- object@lrscore[[score.slot]]@res.list[["result.hq"]][names(intr.order)]
	object@lrscore[[score.slot]]@res.list[["result.hq.pear"]] <- res.list.pear

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

	### cavaet: remember to convert the Uniprot ids to indices!
	# convert nb fac
	# nb.id.fac = sort(nb.id.fac)
	# nb.index = facToIndex(nb.id.fac)

	if (max(as.integer(names(lig.fac))) > nrow(dge.lig)){
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
	# dimnames(res.mtx) = list(colnames(dge.lig), levels(lig.fac))
	# res.mtx = Matrix(res.mtx, sparse = T)

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
#' 
#' @return A Cytosignal object
#' @export
inferVeloLR.CytoSignal <- function(
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

	dge.lig <- object@imputation[[lig.slot]]@imp.data
	dge.recep <- object@imputation[[recep.slot]]@imp.data
	dge.lig.velo <- object@imputation[[lig.slot]]@imp.velo
	dge.recep.velo <- object@imputation[[recep.slot]]@imp.velo

	# compare the dimnames of all four matrices
	if (!all.equal(dimnames(dge.lig), dimnames(dge.recep))){
		stop("dge.lig and dge.recep must have the same dimension names.")
	}

	if (!all.equal(dimnames(dge.lig.velo), dimnames(dge.recep.velo))){
		stop("dge.lig.velo and dge.recep.velo must have the same dimension names.")
	}

	use.genes <- intersect(rownames(dge.lig), rownames(dge.lig.velo))

	dge.lig <- dge.lig[use.genes, ]
	dge.recep <- dge.recep[use.genes, ]
	dge.lig.velo <- dge.lig.velo[use.genes, ]
	dge.recep.velo <- dge.recep.velo[use.genes, ]

	message("Number of velo genes: ", length(use.genes), " / ", nrow(object@imputation[[lig.slot]]@imp.data))
	
	intr.db.list <- checkIntr(use.genes, object@intr.valid[[intr.db.name]])

	res.mtx <- inferVeloLR.matrix_like(dge.lig, dge.recep, 
				dge.lig.velo, dge.recep.velo,
				intr.db.list[["ligands"]], intr.db.list[["receptors"]])

	intr.hq <- names(object@lrscore[[tag]]@res.list[["result.hq.pear"]])
	intr.hq <- intr.hq[intr.hq %in% colnames(res.mtx)]

	lrvelo.obj <- new(
		"lrVelo",
		lig.slot = lig.slot,
		recep.slot = recep.slot,
		intr.slot = intr.db.name,
		intr.list = intr.db.list,
		intr.velo = res.mtx,
		velo.gene = list(
			"velo.genes" = use.genes,
			"velo.intr" = colnames(res.mtx),
			"velo.intr.hq" = intr.hq),
		# nn.id = velo.nn.list[["id"]],
		# nn.dist = velo.nn.list[["dist"]],
		log = list(
			"Used slot" = c(lig.slot, recep.slot)
		)
	)

	object@lrvelo[[tag]] <- lrvelo.obj

	return(object)
}


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
	object@parameters <- list(r.diffuse.scale = r.eps.scale, sigma.scale = sigma.scale)
	return(object)
}


findNNGauEB <- function(
  object,
  ...
) {
  UseMethod(generic = 'findNNGauEB', object = object)
}


### Epsilon Circle
findNNGauEB.matrix <- function(
	cells.loc,
	eps,
	sigma = 0.15,
	weight.sum = 1
){
	cat("Finding neighbors in epsilon circle...\n")
	nn <- dbscan::frNN(cells.loc, eps = eps, sort = F)

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
		if (weight.sum == 1){
			# norm the sum to 1
			y = gauss_vec_cpp(c(x, 1e-8), sigma)[,1]
			return(y/sum(y))
		} else if (weight.sum == 2) {
			# norm the sum except the index cell to 1
			y = gauss_vec_cpp(x, sigma)[,1]
			return(c(y/sum(y), 1))
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

	cat(paste0("Mean num of neighbors: ", mean(table(nn.fac)), "\n"))
	cat(paste0("Median num of neighbors: ", median(table(nn.fac)), "\n"))

	return(list(
		"id" = nn.fac,
		"dist" = nn.dist.fac
	))
}



### Epsilon Circle
findNNGauEB.CytoSignal <- function(
	object,
	eps = NULL,
	sigma = NULL,
	weight.sum = 2
){
	cells.loc <- object@cells.loc

	tag <- "GauEps"

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

	nn <- findNNGauEB.matrix(cells.loc, eps, sigma, weight.sum)

	nn.obj <- methods::new(
		"ImpData",
		method = tag,
		imp.data = new("dgCMatrix"),
		imp.norm.data = new("dgCMatrix"),
		intr.data = new("dgCMatrix"),
		nn.id = nn$id,
		nn.dist = nn$dist,
		log = list(
			"Parameters" = paste0("eps: ", eps, ", sigma: ", sigma),
			"Num of neighbors" = paste0("Mean: ", mean(table(nn$id)), ", Median: ", median(table(nn$id)))
		)
	)

	object@imputation[[tag]] <- nn.obj
	object@imputation[["default"]] <- tag

	return(object)
}


findNNDT <- function(
  object,
  ...
) {
  UseMethod(generic = 'findNNDT', object = object)
}


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

	cat(paste0("Mean num of neighbors: ", mean(nb.size), "\n"))
	cat(paste0("Median num of neighbors: ", median(nb.size), "\n"))

	# return(nn.fac)
	return(list(
		"id" = nn.fac,
		"dist" = nn.dist.fac
	))
}



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
		imp.norm.data = new("dgCMatrix"),
		intr.data = new("dgCMatrix"),
		nn.id = nn$id,
		nn.dist = nn$dist,
		log = list(
			"Parameters" = "Delauany Triangulation",
			"Num of neighbors" = paste0("Mean: ", mean(table(nn$id)), ", Median: ", median(table(nn$id)))
		)
	)

	object@imputation[[tag]] <- nn.obj
	object@imputation[["default"]] <- tag

	return(object)
}


# This function is used to norm and filter by UNPT only
# Corresponding to the situation that user does not want to do imputation
findNNRaw <- function(
	object
) {
	tag <- "Raw"

	message("Setting Imp obj using NO imputation...")

	cells.loc <- object@cells.loc
	nn <- findNNDT.matrix(cells.loc)

	nn.obj <- methods::new(
		"ImpData",
		method = tag,
		imp.data = object@counts,
		imp.norm.data = new("dgCMatrix"),
		intr.data = new("dgCMatrix"),
		nn.id = nn$id,
		nn.dist = nn$dist,
		log = list(
			"Parameters" = "Raw data without imputation",
			"Num of neighbors" = paste0("Mean: ", mean(table(nn$id)), ", Median: ", median(table(nn$id)))
		)
	)

	object@imputation[[tag]] <- nn.obj
	object@imputation[["default"]] <- tag

	return(object)
}



imputeNiche <- function(
  object,
  ...
) {
  UseMethod(generic = 'imputeNiche', object = object)
}


imputeNiche.dgCMatrix <- function(
	dge.raw, 
	nb.id.fac,
	nb.dist.fac,
	weights = c("mean", "counts", "dist", "none")
){
	cat("Imputing...\n")

	i.list = as.integer(names(nb.id.fac)) # x coords, row number
	j.list = as.numeric(as.character(nb.id.fac)) # y coords, col number, hq beads

	if (weights == "avg") { # just do the mean
		x.list = rep(1, length(i.list))
		do.norm <- T
	} else if (weights == "counts") { # use total counts in each bead as weights
		x.list = colSums(dge.raw)[i.list]
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
		weights.mtx@x = weights.mtx@x / rep.int(colSums(weights.mtx), diff(weights.mtx@p))
	}

	dge.raw.imputed = dge.raw %*% weights.mtx

	return(dge.raw.imputed)
}


imputeNiche.CytoSignal <- function(object, 
	nn.type = NULL, 
	weights = c("mean", "counts", "dist", "gaussian"),
	save.obj = TRUE
){
	dge.raw <- object@counts
	# Extract neighborhood factors

	if (is.null(nn.type)){
		nn.type <- object@imputation[["default"]]
	}

	message(paste0("Imputing using ", nn.type, "..."))

	if (nn.type %in% names(object@imputation)){
		nb.id.fac <- object@imputation[[nn.type]]@nn.id
		nb.dist.fac <- object@imputation[[nn.type]]@nn.dist
	} else {
		stop("NN type not found.")
	}

	if (ncol(dge.raw) < length(levels(nb.id.fac))) {
		stop("Number of index beads larger than the number of beads in DGE.")
	}

	dge.raw.imputed <- imputeNiche.dgCMatrix(
		dge.raw, 
		nb.id.fac,
		nb.dist.fac,
		weights = weights
	)

	# object@imp.data[[nn.type]] <- dge.raw.imputed

	res.density <- sum(dge.raw.imputed != 0)/length(dge.raw.imputed) # density 6.2%
	cat(paste0("Density after imputation: ", res.density*100, "%\n"))

	if (save.obj){
		object@imputation[[nn.type]]@imp.data <- dge.raw.imputed
	}

	object@imputation[[nn.type]]@log[["Density:"]] <- paste0(res.density*100, "%")

	return(object)
}


normCounts <- function(
  object,
  ...
) {
  UseMethod(generic = 'normCounts', object = object)
}



normCounts.dgCMatrix <- function(
	mat, 
	method = c("default", "cpm")
){
	if (method == "default"){
		mat@x <- mat@x / rep.int(colSums(mat), diff(mat@p))
		mat@x <- log1p(mat@x * 1e4)
	} else if (method == "cpm"){
		mat@x <- mat@x / rep.int(colSums(mat), diff(mat@p))
		mat@x <- log1p(mat@x * 1e6)
	} else {
		stop("Method not found.")
	}

	return(mat)
}



normCounts.CytoSignal <- function(
	object, 
	method = c("default", "cpm"),
	slot.use = NULL,
	save.obj = TRUE
){
	if (is.null(slot.use)){
		slot.use <- object@imputation[["default"]]
	}

	if (!slot.use %in% names(object@imputation)) {
		stop("Data not found.")
	}

	message(paste0("Normalizing using ", slot.use, "..."))

	mat <- object@imputation[[slot.use]]@imp.data
	mat <- normCounts.dgCMatrix(mat, method = method)

	if (save.obj){
		object@imputation[[slot.use]]@imp.norm.data <- mat
	}

	# if (!all.equal(names(slot(object, slot.use)), names(object@imp.norm.data))){
	# 	warning("Names of imp.data and imp.norm.data not equal!")
	# }
	
	return(object)

}


# normCounts_list <- function(object, 
# 	slot.use = NULL, 
# 	method = c("default", "cpm")
# ){
# 	if (is.null(slot.use)){
# 		slot.use <- object@imputation[["default"]]
# 	}

# 	if (!slot.use %in% slotNames(object)) {
# 		stop("Data not found.")
# 	}

# 	mat_list_norm <- lapply(object@imputation[[slot.use]], function(mat) {
# 		if (method == "default") {
# 			mat@x = mat@x / rep.int(colSums(mat), diff(mat@p))
# 			mat@x = lop1p(mat@x * 1e4)
# 		} else if (method == "cpm") {
# 			mat@x = mat@x / rep.int(colSums(mat), diff(mat@p))
# 			mat@x = lop1p(mat@x * 1e6)
# 		} else {
# 			stop("Method not found.")
# 		}
		
# 		return(mat)
# 	})
	
# 	object@imp.norm.data <- mat_list_norm
	
# 	if (!all.equal(names(slot(object, slot.use)), names(object@imp.norm.data))) {
# 		warning("Names of imp.data and imp.norm.data not equal!")
# 	}
	
# 	return(object)
# }


changeUniprot <- function(
  object,
  ...
) {
  UseMethod(generic = 'changeUniprot', object = object)
}


changeUniprot.dgCMatrix <- function(
    dge.raw, 
    gene_to_uniprot,
    # mode = "unpt",
    verbose = T
){
    # uppercase the gene to match database
    rownames(dge.raw) = toupper(rownames(dge.raw))

	unique.check = sapply(rownames(dge.raw), function(x){
		length(unique(gene_to_uniprot$uniprot[which(gene_to_uniprot$gene_name == x)]))
	})

	if (max(unique.check) > 1){ # remove duplicated genes if nany
		dup.num = sum(unique.check > 1)
		stop(paste0("Must remove duplicated Uniprot items: ", dup.num, " genes."))
	}
    
    uniprot.genes = rownames(dge.raw)[rownames(dge.raw) %in% gene_to_uniprot$gene_name] # 738/22683

    if (verbose){
        cat(paste0("Number of genes in the database: ", length(uniprot.genes), "\n"))
    }

    uniprot.names = sapply(uniprot.genes, function(x){
        unique(gene_to_uniprot$uniprot[which(gene_to_uniprot$gene_name == x)])
    })
    
    dge.raw = dge.raw[uniprot.genes, ]

    # rownames(dge.raw) = unname(lowwords(rownames(dge.raw)))
	rownames(dge.raw) = unname(uniprot.names)

    # if (mode == "unpt"){
    #     rownames(dge.raw) = unname(uniprot.names)
    # }

    if (anyDuplicated(rownames(dge.raw)) != 0){
		del.num <- sum(duplicated(rownames(dge.raw)))
		message(paste0("Removing ", del.num, " duplicated unpt genes."))

        if (del.num > 10){
            stop("More than 10 duplicated Uniprot items detected.")
        }

        dup.index = which(duplicated(rownames(dge.raw)) == T)
        dge.raw = dge.raw[-dup.index, ]
		uniprot.names = uniprot.names[-dup.index]
		uniprot.genes = uniprot.genes[-dup.index]
    }
    
	names(uniprot.names) <- lowwords(uniprot.genes)

    return(list(
		dge.raw,
		uniprot.names
	))
}


changeUniprot.CytoSignal <- function(
	object,
	slot.use = NULL, 
	# mode = "unpt",
	verbose = T
){
	if (is.null(slot.use)){
		slot.use <- object@imputation[["default"]]
	}

	if (!slot.use %in% names(object@imputation)){
		stop("Data not found.")
	}

	gene_to_uniprot <- object@intr.valid[["gene_to_uniprot"]]
	dge.raw <- object@imputation[[slot.use]]@imp.norm.data
	# check duplicate items exists in database

	unpt.list <- changeUniprot.dgCMatrix(dge.raw, gene_to_uniprot, verbose = verbose)
	
	object@imputation[[slot.use]]@intr.data <- unpt.list[[1]]
	object@intr.valid[["symbols"]][[slot.use]] <- unpt.list[[2]]

	# if (!all.equal(names(slot(object, slot.use)), names(object@intr.data))) {
	# 	warning("Names of imp.norm.data and intr.data not equal!")
	# }

	return(object)
}



graphNicheLR <- function(
  object,
  ...
) {
  UseMethod(generic = 'graphNicheLR', object = object)
}


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

# nn.use: slot that the neighbor index should be taken from, by default is the same as
# the recep.slot. For example, if score.obj = GauEps-DT, then nn.use = "DT".
# nn.use could also be a user-defind factor.
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
		if (length(nn.use) != ncol(object@imputation[[lig.slot]]@intr.data))
			stop("nn.use must have the same length as the number of cells.")
		nb.id.fac <- nn.use
	} else {
		stop("nn.use must be either a factor or a character.")
	}

	dge.lig <- object@imputation[[lig.slot]]@intr.data
	dge.recep <- object@imputation[[recep.slot]]@intr.data

	if (!all.equal(dim(dge.lig), dim(dge.recep))){
		stop("dge.lig and dge.recep must have the same dimension.")
	}

	if (ncol(dge.lig) < length(levels(nb.id.fac))){
		stop("Number of index beads larger than the number of beads in DGE.")
	}

	if (!all.equal(object@intr.valid[["symbols"]][[lig.slot]], 
					object@intr.valid[["symbols"]][[recep.slot]])){
		stop("Unpt symbols generated from imputations not equal!")
	}

	message("Computing scores using ", intr.db.name, " database.")
	message("Ligand: ", lig.slot, ", Receptor: ", recep.slot, ".")

	intr.db.list <- checkIntr(unname(object@intr.valid[["symbols"]][[lig.slot]]), 
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



permuteImpLR <- function(
	object,
	nn.type = NULL,
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

	# if (!use.gau & !use.dt){
	# 	stop("Cannot find corresponding imputation method.")
	# }

	message("Permuting imputation data on method ", nn.type, "...")

	dimnames.list <- dimnames(object@imputation[[nn.type]]@intr.data)
	cells.loc <- object@cells.loc
	dge.raw <- object@counts

	# permute cell locations
	null.cells.loc.list = lapply(1:times, function(i){
		perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
		null.cells.loc = cells.loc[perm.index,]
		rownames(null.cells.loc) = rownames(cells.loc)
		return(null.cells.loc)
	})
	
	# Null imputations
	null.dge.list <- lapply(1:times, function(i){
		## MUST shuffule raw dge!!
		null.dge.raw <- dge.raw[, sample(seq_len(ncol(dge.raw)), ncol(dge.raw), replace = F)]

		if (use.gau) {
			eps.params <- object@parameters
			null.eps.nb <- findNNGauEB.matrix(null.cells.loc.list[[i]], eps = eps.params[[1]], sigma = eps.params[[2]], weight.sum = 2)
			null.dge.eps <- imputeNiche.dgCMatrix(null.dge.raw, null.eps.nb$id, null.eps.nb$dist, weights = "none"); gc()
		} else if (use.dt) {
			null.nb.fac <- findNNDT.matrix(null.cells.loc.list[[i]])
        	null.dge.eps <- imputeNiche.dgCMatrix(null.dge.raw, null.nb.fac$id, null.nb.fac$dist, weights = "none"); gc() # den: 31%
		} else if (use.raw) {
			null.dge.eps <- null.dge.raw
		} else {
			stop("Cannot find corresponding imputation method.")
		}

		null.dge.eps = normCounts.dgCMatrix(null.dge.eps, "default"); gc()
		null.dge.eps.unpt = changeUniprot.dgCMatrix(null.dge.eps, object@intr.valid[["gene_to_uniprot"]])[[1]]
		rm(null.dge.eps, null.dge.raw); gc()

		return(null.dge.eps.unpt)
	})

	null.dge.eps.unpt = meanMat_cpp(null.dge.list, nrow(null.dge.list[[1]]), ncol(null.dge.list[[1]]))
	dimnames(null.dge.eps.unpt) <- dimnames.list
	rm(null.dge.list); gc()

    object@imputation[[nn.type]]@intr.data.null <- null.dge.eps.unpt

	return(object)

}


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

	null.dge.eps.unpt <- object@imputation[[lig.slot]]@intr.data.null
	null.dge.dt.unpt <- object@imputation[[recep.slot]]@intr.data.null

    perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
    null.cells.loc = cells.loc[perm.index,]
    rownames(null.cells.loc) = rownames(cells.loc)
    null.nb.fac = findNNDT.matrix(null.cells.loc); gc() # Mean num of neighbors: 44, median: 36

    null.lrscore.mtx = graphNicheLR.dgCMatrix(null.dge.eps.unpt, null.dge.dt.unpt, null.nb.fac[["id"]],
						intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()

    if (sum(colSums(null.lrscore.mtx) == 0) != 0){
        message("A total of ", sum(colSums(null.lrscore.mtx) == 0), " intr are empty in NULL scores.")
        null.lrscore.mtx = null.lrscore.mtx[, !colSums(null.lrscore.mtx) == 0]
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


inferSignif <- function(
    object,
    ...
) {
    UseMethod(generic = 'inferSignif', object = object)
}

inferSignif.dgCMatrix <- function(
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

inferSignif.CytoSignal <- function(
    object,
	p.value = 0.05, 
	reads.thresh = 100,
	sig.thresh = 100,
    slot.use = NULL
){
    if (is.null(slot.use)){
        slot.use = object@lrscore[["default"]]
    }

    if (!slot.use %in% names(object@lrscore)){
        stop("LRscores not found. ")
    }

	message("Inferring significant beads on Score slot ", slot.use, "... ")

    # lrscore.mtx = object@lrscore[[slot.use]]@score
    # null.lrscore.mtx = object@lrscore[[slot.use]]@score.null

    nb.fac = list(
        id = object@imputation[["DT"]]@nn.id,
        dist = object@imputation[["DT"]]@nn.dist
    )

	use.intr.slot.name <- object@lrscore[[slot.use]]@intr.slot
	use.intr.db <- object@intr.valid[[use.intr.slot.name]]

    res.list = inferSignif(object@counts, object@lrscore[[slot.use]]@score, 
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



# Normal Moran's I test is not applicable here since the total number of the cell is too large, causing 
# unnacceptable computation cost. Here we use a modified version of Moran's I test, which is to take only 
# the top KNNs to compute the Moran's I test.
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

    lr.mtx <- object@lrscore[[score.slot]]@score[, intr.hq]
	lt.mtx.imp <- dataImpKNN(lr.mtx, object@cells.loc, k = k, weight = weight)
	lt.mtx.imp <- as.matrix(t(lt.mtx.imp))

	intr.order <- as.double(pearson_col_cpp(lt.mtx.imp, lr.mtx)) * as.double(stdMat_cpp(lt.mtx.imp)) / as.double(stdMat_cpp(lr.mtx))
    names(intr.order) <- colnames(lr.mtx)
	intr.order <- intr.order[order(intr.order, decreasing = T)]

	cat("Reordering significant interactions... \n")
		
	res.list.pear <- object@lrscore[[score.slot]]@res.list[["result.hq"]][names(intr.order)]
	object@lrscore[[score.slot]]@res.list[["result.hq.pear"]] <- res.list.pear

	return(object)
}


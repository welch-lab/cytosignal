inferEpsParams <- function(object, tech.res, bin_size, loc.d, r.eps.real = 200, thresh = 0.01){
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


### Epsilon Circle
findNNGauEB <- function(
	object,
	eps = NULL,
	sigma = NULL,
	weight.sum = 2
){
	cells.loc <- object@cells.loc

	tag <- "GauEps"

	if (tag %in% names(object@imputation)){
		stop("This imputation has been done before.")
	}

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

	# object@nn[["epsilon"]] <- list("id" = nn.fac, "dist" = nn.dist.fac)

	nn.obj <- methods::new(
		"ImpData",
		method = tag,
		imp.data = new("dgCMatrix"),
		imp.norm.data = new("dgCMatrix"),
		intr.data = new("dgCMatrix"),
		nn.id = nn.fac,
		nn.dist = nn.dist.fac,
		log = list(
			"Parameters" = paste0("eps: ", eps, ", sigma: ", sigma),
			"Num of neighbors" = paste0("Mean: ", mean(table(nn.fac)), ", Median: ", median(table(nn.fac)))
		)
	)

	object@imputation[[tag]] <- nn.obj
	object@imputation[["default"]] <- tag

	return(object)
}



findNNDT <- function(
    object
){
	tag <- "DT"

	if (tag %in% names(object@imputation)){
		stop("This imputation has been done before.")
	}

	cells.loc <- object@cells.loc
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

    dist.list <- as.double(euclidean_cpp(index.loc, nb.loc))

    nn.dist.fac = nn.fac
    names(nn.dist.fac) = dist.list

    nn.fac = addIndex(nn.fac) # add index bead it self to the neighbor list
    nn.dist.fac = addIndexZero(nn.dist.fac)

    cat(paste0("Mean num of neighbors: ", mean(table(nn.fac)), "\n"))
    cat(paste0("Median num of neighbors: ", median(table(nn.dist.fac)), "\n"))

	nn.obj <- methods::new(
		"ImpData",
		method = tag,
		imp.data = new("dgCMatrix"),
		imp.norm.data = new("dgCMatrix"),
		intr.data = new("dgCMatrix"),
		nn.id = nn.fac,
		nn.dist = nn.dist.fac,
		log = list(
			"Parameters" = "None",
			"Num of neighbors" = paste0("Mean: ", mean(table(nn.fac)), ", Median: ", median(table(nn.fac)))
		)
	)

	object@imputation[[tag]] <- nn.obj
	object@imputation[["default"]] <- tag

	return(object)
}



imputeNiche <- function(object, 
	nn.type = NULL, 
	weights = c("mean", "counts", "dist", "gaussian")
){
	dge.raw <- object@counts
	# Extract neighborhood factors

	if (is.null(nn.type)){
		nn.type <- object@imputation[["default"]]
	}

	cat(paste0("Imputing using ", nn.type, "...\n"))

	if (nn.type %in% names(object@imputation)){
		nb.id.fac <- object@imputation[[nn.type]]@nn.id
		nb.dist.fac <- object@imputation[[nn.type]]@nn.dist
	} else {
		stop("NN type not found.")
	}

	if (ncol(dge.raw) < length(levels(nb.id.fac))) {
		stop("Number of index beads larger than the number of beads in DGE.")
	}

	cat("Imputing...\n")

	i.list = as.integer(names(nb.id.fac)) # x coords, row number
	j.list = as.numeric(as.character(nb.id.fac)) # y coords, col number, hq beads

	if (weights == "mean") { # just do the mean
		x.list = rep(1, length(i.list))
	} else if (weights == "counts") { # use total counts in each bead as weights
		x.list = colSums(dge.raw)[i.list]
	} else if (weights == "dist") {
		x.list = exp(-as.numeric(names(nb.dist.fac))) # exp(-dist) as similarities
	} else if (weights == "gaussian") {
		x.list = as.numeric(names(nb.dist.fac)) # gaussian distances has been calculated
	} else {
		stop("Incorret weighting method!\n")
	}

	weights.mtx = sparseMatrix(i = i.list, j = j.list, x = x.list, 
		dims = c(ncol(dge.raw), ncol(dge.raw)), dimnames = list(NULL, colnames(dge.raw)))

	# re-weight the weight mtx across each niche (column)
	weights.mtx@x = weights.mtx@x / rep.int(colSums(weights.mtx), diff(weights.mtx@p))
	dge.raw.imputed = dge.raw %*% weights.mtx

	res.density = sum(dge.raw.imputed != 0) / length(dge.raw.imputed) # density 6.2%
	cat(paste0("Density after imputation: ", res.density*100, "%\n"))

	# object@imp.data[[nn.type]] <- dge.raw.imputed

	object@imputation[[nn.type]]@imp.data <- dge.raw.imputed
	object@imputation[[nn.type]]@log[["Density:"]] <- paste0(res.density*100, "%")

	return(object)
}



normCounts <- function(
	object, 
	method = c("default", "cpm"),
	slot.use = NULL
){
	if (is.null(slot.use)){
		slot.use <- object@imputation[["default"]]
	}

	if (!slot.use %in% names(object@imputation)) {
		stop("Data not found.")
	}

	mat <- object@imputation[[slot.use]]@imp.data

	if (method == "default"){
		mat@x = mat@x / rep.int(colSums(mat), diff(mat@p))
		mat@x = log1p(mat@x * 1e4)
	} else if (method == "cpm"){
		mat@x = mat@x / rep.int(colSums(mat), diff(mat@p))
		mat@x = log1p(mat@x * 1e6)
	} else {
		stop("Method not found.")
	}

	object@imputation[[slot.use]]@imp.norm.data <- mat

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
	slot.use = NULL, 
    mode = "unpt",
    verbose = T
){
	if (is.null(slot.use)){
		slot.use <- object@imputation[["default"]]
	}

	if (!slot.use %in% names(object@imputation)){
		stop("Data not found.")
	}

	dge.raw <- object@imputation[[slot.use]]@imp.norm.data
    # check duplicate items exists in database
    unique.check = sapply(rownames(dge.raw), function(x){
        length(unique(gene_to_uniprot$uniprot[which(gene_to_uniprot$gene_name == x)]))
    })

    if (max(unique.check) > 1){
        # dup.genes = which(unique.check > 1)
        stop("Duplicated Uniprot items detected.")
    }

    # uppercase the gene to match database
    rownames(dge.raw) = toupper(rownames(dge.raw))
    
    uniprot.genes = rownames(dge.raw)[rownames(dge.raw) %in% gene_to_uniprot$gene_name] # 738/22683

    if (verbose){
        cat(paste0("Number of genes in the database: ", length(uniprot.genes), "\n"))
    }

    uniprot.names = sapply(uniprot.genes, function(x){
        unique(gene_to_uniprot$uniprot[which(gene_to_uniprot$gene_name == x)])
    })
    
    dge.raw = dge.raw[uniprot.genes, ]

    rownames(dge.raw) = unname(lowwords(rownames(dge.raw)))
	# uniprot.genes = unname(lowwords(uniprot.genes))

    if (mode == "unpt"){
        rownames(dge.raw) = unname(uniprot.names)
    }

    if (anyDuplicated(rownames(dge.raw)) != 0){
        cat("Duplicated Uniprot items detected.")

        if (sum(duplicated(rownames(dge.raw))) > 10){
            stop("More than 10 duplicated Uniprot items detected.")
        }

        dup.index = which(duplicated(rownames(dge.raw)) == T)
        dge.raw = dge.raw[-dup.index, ]
    }
    
	object@imputation[[slot.use]]@intr.data <- dge.raw
	names(uniprot.names) <- lowwords(uniprot.genes)
	object@intr.valid[["symbols"]][[slot.use]] <- uniprot.names

	# if (!all.equal(names(slot(object, slot.use)), names(object@intr.data))) {
	# 	warning("Names of imp.norm.data and intr.data not equal!")
	# }

    return(object)
}


# subset intr.fac, lig.fac, recep.fac with unpts in DGE
checkIntr <- function(object, intr.db){
	if (!all.equal(object@intr.valid[["symbols"]][[1]], object@intr.valid[["symbols"]][[2]])){
		stop("Symbols not equal!")
	}

	symbols.list = unname(object@intr.valid[["symbols"]][[1]])
    
	intr.fac = intr.db[['combined']]
    lig.fac = intr.db[['ligands']]
    recep.fac = intr.db[['receptors']]

    # filter all empty intrs
    if.all.in = sapply(levels(intr.fac), function(intr) {
        intr.cmp = names(intr.fac[intr.fac == intr, drop = T])
        return(sum(intr.cmp %in% symbols.list) == length(intr.cmp))
    })
    cat(paste0("Number of valid intrs: ", sum(if.all.in), " / ", length(if.all.in), "\n"))
    levels.valid = levels(intr.fac)[if.all.in]
    # se.cmp.fac = cmp.fac[cmp.fac %in% levels(cmp.fac)[if.pair.in], drop = T] # 294 / 1396 levels

    intr.fac = intr.fac[intr.fac %in% levels.valid, drop = T]
    lig.fac = lig.fac[lig.fac %in% levels.valid, drop = T]
    recep.fac = recep.fac[recep.fac %in% levels.valid, drop = T]

    # convert lig and recep
    names(lig.fac) = unname(sapply(names(lig.fac), function(gene){which(symbols.list == gene)}))
    names(recep.fac) = unname(sapply(names(recep.fac), function(gene){which(symbols.list == gene)}))

	object@intr.valid[["combined"]] <- intr.fac
	object@intr.valid[["ligands"]] <- lig.fac
	object@intr.valid[["receptors"]] <- recep.fac

	return(object)
}


graphNicheLR <- function(
  object,
  ...
) {
  UseMethod(generic = 'optimizeALS', object = object)
}


graphNicheLR.list <- function(
	object
){

	dge.lig <- object[[1]]
	dge.recep <- object[[2]]
	nb.id.fac <- object[[3]]

    ### cavaet: remember to convert the Uniprot ids to indices!
    # convert nb fac
    nb.id.fac = sort(nb.id.fac)
    nb.index = facToIndex(nb.id.fac)

    lig.fac = object@intr.valid[["ligands"]]
    recep.fac = object@intr.valid[["receptors"]]

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


graphNicheLR.CytoSignal <- function(
	object,
	lig.slot,
	recep.slot
){

	dge.lig <- object@imputation[[lig.slot]]@intr.data
	dge.recep <- object@imputation[[recep.slot]]@intr.data
	nb.id.fac <- object@imputation$DT@nn.id

    if (!all.equal(dim(dge.lig), dim(dge.recep))){
        stop("dge.lig and dge.recep must have the same dimension.")
    }

    if (ncol(dge.lig) < length(levels(nb.id.fac))){
        stop("Number of index beads larger than the number of beads in DGE.")
    }

    # if (!all(unique(names(lig.fac)) %in% rownames(dge))){
    #     stop("Not all Uniprot names in intr.fac are in input dge.")
    # }

	res.mtx <- graphNicheLR.list(
		list(dge.lig, dge.recep, nb.id.fac)
	)

	score.obj <- methods::new(
		"lrScores",
		lig.slot = lig.slot,
		recep.slot = recep.slot,
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




nullNicheLR <- function(
	object,
	slot.use = NULL
){
	if (is.null(slot.use)){
		slot.use <- object@lrscore[["default"]]
	}

	score.obj <- object@lrscore[[slot.use]]
	
}
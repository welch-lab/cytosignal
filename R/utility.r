filterGene <- function(dge.raw, gene_to_uniprot, thresh = 50){
    # find intr gene index
    rownames(dge.raw) = toupper(rownames(dge.raw))
    # uniprot.genes = rownames(dge.raw)[rownames(dge.raw) %in% gene_to_uniprot$gene_name] # 738/22683
    intr.genes.index = which(rownames(dge.raw) %in% gene_to_uniprot$gene_name)
    cat(paste0("Number of genes in the database: ", length(intr.genes.index), "\n"))

    # find mito gene index
    mt_idx <- grep("MT-",rownames(dge.raw))

    # find low quality gene index
    low_quality_idx <- which(Matrix::rowSums(dge.raw) < thresh)

    cat(paste0("Number of low-quality intr genes: ", sum(intr.genes.index %in% low_quality_idx), "\n"))

    # return non-intr filtered gene index
    fin.index = setdiff(union(low_quality_idx, mt_idx), intr.genes.index)
    cat(paste0("Number of genes to be removed: ", length(fin.index), " / ", nrow(dge.raw), "\n"))

    return(fin.index)
}

addIndex <- function(fac.use){
    fac.use = sort(fac.use)

    new.fac = as.integer(c(as.character(fac.use), levels(fac.use)))
    new.fac = factor(new.fac)

    names(new.fac) = as.integer(c(names(fac.use), levels(fac.use)))
    new.fac = sort(new.fac)

    return(new.fac)
}

addIndexZero <- function(fac.use){
    fac.use = sort(fac.use)

    new.fac = as.integer(c(as.character(fac.use), levels(fac.use)))
    new.fac = factor(new.fac)

    names(new.fac) = as.numeric(c(names(fac.use), rep(0, length(levels(fac.use)))))
    new.fac = sort(new.fac)

    return(new.fac)
}

addIndexOne <- function(fac.use){
    fac.use = sort(fac.use)

    new.fac = as.integer(c(as.character(fac.use), levels(fac.use)))
    new.fac = factor(new.fac)

    names(new.fac) = as.numeric(c(names(fac.use), rep(1, length(levels(fac.use)))))
    new.fac = sort(new.fac)

    return(new.fac)
}

lowwords <- function(w){
    paste0(
        substr(w, 1, 1),
        tolower(substr(w, 2, 999))
    )
}



facToIndex <- function(fac.use){
    # Note that i_index and nb_index start with 0 --> add 0 to the front
    nb.index = c(0, cumsum(as.integer(table(fac.use))))
    nb.list = as.integer(names(fac.use))-1
    return(list(
        index = nb.index,
        nb = nb.list
    ))
}


# subset intr.fac, lig.fac, recep.fac with unpts in DGE
checkIntr <- function(symbols.list, intr.db) {
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

	# convert lig and recep into index in the matrix
	names(lig.fac) = unname(sapply(names(lig.fac), function(gene){which(symbols.list == gene)}))
	names(recep.fac) = unname(sapply(names(recep.fac), function(gene){which(symbols.list == gene)}))

	return(list(
		"combined" = intr.fac,
		"ligands" = lig.fac,
		"receptors" = recep.fac
	))
}




getPvalues <- function(
    i.mtx,
    null.i.mtx,
    adjust = F
){
    if (ncol(i.mtx) < ncol(null.i.mtx)){
        cat("Real lr-score mtx has fewer features than NULL lr-score mtx.\n")
        cat("Removing extra", ncol(null.i.mtx) - ncol(i.mtx), "interactions from NULL lr-score mtx.\n")
        null.i.mtx <- null.i.mtx[, colnames(i.mtx)]
    }

    if (ncol(i.mtx) > ncol(null.i.mtx)){
        cat("NULL lr-score mtx has fewer features than Real lr-score mtx.\n")
        cat("Removing extra", ncol(i.mtx) - ncol(null.i.mtx), "interactions from Real lr-score mtx.\n")
        i.mtx <- i.mtx[, colnames(null.i.mtx)]
    }

    p.val.mtx = sapply(seq_len(ncol(i.mtx)), function(n){
        p.values = 1 - findInterval(i.mtx[, n], sort(null.i.mtx[, n]), left.open = T)/nrow(null.i.mtx)
        if (adjust){p.values = p.adjust(p.values, "BH")}
        # return(p.adjust(p.values, "BH"))
        return(p.values)
    })
    dimnames(p.val.mtx) = dimnames(i.mtx)

    return(p.val.mtx)
}



graphSpatialFDR <- function(nb.fac, pval.mtx, spatial.coords=NULL, weighting='DT',
                            k=NULL, reciprocal=T){

    # # Discarding NA pvalues.
    # haspval <- !is.na(pvalues)

    # if (!all(haspval)) {
    #     coords <- coords[haspval, , drop=FALSE]
    #     pvalues <- pvalues[haspval]
    # }

    # if(weighting[1] == "none"){
    #     return(rep(NA_real_, length(pvalues)))
    # }
    # # if the weighting vector length > 1 then just take the first element as this is the default
    # # is this a bit hacky?
    # if(length(weighting) > 1){
    #     weighting <- weighting[1]
    # }

    nn.index = nb.fac[["id"]]
    nn.dist = nb.fac[["dist"]]

    # spatial.coords = cells.loc
    # pvalues = pval.mtx

    if (length(levels(nn.index)) != nrow(pval.mtx)){
        warning("Cell numbers in pvalue mtx and neighbor factor do not match!!\n")
    }

    if (weighting == "DT"){
        # for DT, we use the max distance within the nb as the weight
        cat("Using DT weighting.\n")
        # t.connect <- sapply(levels(nn.dist), function(x){
        #     nb.dist = as.numeric(names(nn.dist[nn.dist == x]))
        #     return(max(nb.dist))
        # })

        t.connect <- sapply(split(as.double(names(nn.dist)), nn.dist), function(x){
            return(max(x))
        })

    } else if (weighting == "EB"){

    } else if (weighting == "KNN"){

    } else{
        stop("Weighting must be one of 'DT', 'EB', 'KNN'.")
    }

    # use 1/connectivity as the weighting for the weighted BH adjustment from Cydar
    if (reciprocal){
        w <- 1/t.connect
    } else{
        w <- t.connect
    }
    # w <- 1/t.connect
    w[is.infinite(w)] <- 1
    
    # Compute a adjusted pvalues for every interaction
    padj.mtx = sapply(1:ncol(pval.mtx), function(intr){
        pvalues = unname(pval.mtx[ ,intr])
        haspval <- !is.na(pvalues)
        o <- order(pvalues)
        pvalues <- pvalues[o]
        w.sub <- w[o]

        adjp <- numeric(length(o))
        adjp[o] <- rev(cummin(rev(sum(w.sub)*pvalues/cumsum(w.sub))))
        adjp <- pmin(adjp, 1)

        if (!all(haspval)) {
            refp <- rep(NA_real_, length(haspval))
            refp[haspval] <- adjp
            adjp <- refp
        }
        return(adjp)
    })

    # # Computing a density-weighted q-value.
    # o <- order(pvalues)
    # pvalues <- pvalues[o]
    # w <- w[o]

    # adjp <- numeric(length(o))
    # adjp[o] <- rev(cummin(rev(sum(w)*pvalues/cumsum(w))))
    # adjp <- pmin(adjp, 1)

    # if (!all(haspval)) {
    #     refp <- rep(NA_real_, length(haspval))
    #     refp[haspval] <- adjp
    #     adjp <- refp
    # }
    # return(adjp)

    return(padj.mtx)
}


# This is the second layer of filtering after the initial filtering of the total counts
# within each cell. This function filters out interactions that have: 1) low total counts
# in the cells, 2) low number of significant interactions.
filterRes <- function(dge.raw, res.list, intr.db, gene_to_uniprot, reads.thresh = 100, sig.thresh = 100){
    low_genes = diff(Matrix::t(dge.raw)@p)
    
    # # if dge.raw is genes X cell
    # low_genes = toupper(rownames(dge.raw)[which(low_genes < reads.thresh)])
    # low_genes = intersect(low_genes, gene_to_uniprot$gene_name)
    # low_unpt = gene_to_uniprot$uniprot[gene_to_uniprot$gene_name %in% low_genes]

    # if dge is UNPT X cell
    low_unpt = rownames(dge.raw)[which(low_genes < reads.thresh)]

    low_unpt = low_unpt[!duplicated(low_unpt)]
    all.intrs = intr.db[["combined"]]
    low.intrs = levels(all.intrs[names(all.intrs) %in% low_unpt, drop = T])

    intr.hq = names(res.list)[!names(res.list) %in% low.intrs]

    res.list.hq = res.list[intr.hq]
    res.list.hq = res.list.hq[lengths(res.list.hq) != 0]
    res.list.hq = res.list.hq[order(lengths(res.list.hq), decreasing = T)]

    res.list.hq = res.list.hq[lengths(res.list.hq) >= sig.thresh]
    cat("Number of high quality intr: ", length(res.list.hq), "/ ", length(res.list), "\n")

    return(res.list.hq)
}


# impute data using KNN
dataImpKNN <- function(
    data, 
    cells.loc, 
    # density = 0.05, 
    k = NULL, 
    weight = 2
) {
    # standard version 1: if alpha% of cells < 500, use 500; otherwise use alpha = 5%
    # standard version 2: if alpha% of cells < 500, use 500; if alpha% of cells > 2000, use 2000; otherwise use alpha = 5%

    # den.thresh = round(nrow(cells.loc)*density)
    # if (is.null(k)){
    #     if (den.thresh < 500){
    #         k = 500
    #     } else if (den.thresh > 2000){
    #         k = 2000
    #     } else {
    #         k = den.thresh
    #     }
    # } else {k = round(k)}

    message(paste0("Using ", k, " neighbors for imputation..."))
    # attention that the cells.loc is already scaled but not centered!!!!!!
    nn <- RANN::nn2(cells.loc, cells.loc, k = k, searchtype = "priority")
    # cat("Finished!\n")

    # cat("Gaussian transforming...")
    nn.idx <- t(nn[["nn.idx"]]) # k X N
    # nn.dist <- t(nn[["nn.dists"]]) # k X N

    if (weight == 2){
        nn.dist <- matrix(1/(k-1), nrow = k, ncol = ncol(nn.idx)) # k X N
        nn.dist[1,] <- 1
    } else if (weight == 1) {
        nn.dist <- matrix(1/(k), nrow = k, ncol = ncol(nn.idx)) # k X N
    }
    
    rm(nn); gc()
    # cat("Finished!\n")
    
    # cat("Creating weight matrix...")
    # construct the weight matrix direcly
    i <- c(nn.idx) # join each column and convert to int vector
    j <- as.integer(rep_each_cpp(ncol(nn.idx), k)) # rep (1,2,3,...,n), each k times
    x <- c(nn.dist)
    dist.graph <- sparseMatrix(i, j, x = x) # every COLUMN represents a cell's distances to its k nearest neighbors
    # all.equal(which(dist.graph[,1] != 0), sort(nn.idx[,1]))

    # dist.graph <- quickNorm(dist.graph)

    #### this matrix is NOT symmetric!!!!!
    # cat("Finished!\n")
    rm(i, j, x, nn.idx, nn.dist); gc()

    cat("Imputing...")

    # dge.intr.imp <- dge.intr %*% dist.graph
    # dge.mix.imp <- rbind(dge.intr.imp, dge.other)

    data.imp <- Matrix::crossprod(data, dist.graph)
    data.imp <- as(data.imp, "CsparseMatrix")

    cat("Finished!\n")

    dimnames(data.imp) <- list(
        colnames(data), 
        rownames(data)
    )

    # sum(dge.raw!=0)/length(dge.raw) # density 98%

    return(data.imp)
}


getIntrNames <- function(res.intr.list, inter.index){
    intr.names = sapply(seq_along(res.intr.list), function(i){
        pair.index = which(inter.index$id_cp_interaction == res.intr.list[i])

        ligand.name = inter.index[pair.index, 4]
        if (ligand.name == ""){ # if complex
            ligand.name = inter.index[pair.index, 2]
        }
        ligand.name = gsub("_HUMAN", "", ligand.name)

        receptor.name = inter.index[pair.index, 5]
        if (receptor.name == ""){
            receptor.name = inter.index[pair.index, 3]
        }
        receptor.name = gsub("_HUMAN", "", receptor.name)

        # generate interaction names
        # intr.name = paste0(gsub("_HUMAN", "", ligand.name), "-", gsub("_HUMAN", "", receptor.name))
        intr.name = paste0(ligand.name, "-", receptor.name)
    })

    return(intr.names)
}


# x is the new column names of the sparse matrix
# y is the sparse matrix
increase_columns_sparse <- function(x, y) {
	# Get the column names of the sparse matrix
    colnames_y <- colnames(y)
	
	# Determine which columns are missing
	missing_cols <- setdiff(x, colnames_y)
	
	# Add missing columns filled with 0s
	if (length(missing_cols) > 0) {
		n_missing_cols <- length(missing_cols)
		n_rows <- nrow(y)
		y_new <- cbind(y, sparseMatrix(i = integer(0), j = integer(0), x = double(0), 
										dims = c(n_rows, n_missing_cols), 
										dimnames = list(NULL, missing_cols)))
	} else {
		y_new <- y
	}
	
	# Reorder columns to match x
	y_new[, x]
}



# x is the new column names of the matrix
# y is the matrix
increase_columns <- function(x, y) {
	# Get the column names of the matrix
	colnames_y <- colnames(y)
	
	# Determine which columns are missing
	missing_cols <- setdiff(x, colnames_y)
	
	# Add missing columns filled with 0s
	if (length(missing_cols) > 0) {
		n_missing_cols <- length(missing_cols)
		n_rows <- nrow(y)
		y_new <- cbind(y, matrix(0, nrow = n_rows, ncol = n_missing_cols,
								 dimnames = list(NULL, missing_cols)))
	} else {
		y_new <- y
	}
	
	# Reorder columns to match x
	y_new[, x]
}


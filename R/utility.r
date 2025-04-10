#' Process Ligand–Receptor Data Without Assigning New Random CPI-CC IDs
#'
#' This function processes user-supplied ligand–receptor data frames (ligands, receptors, and interaction types) and a gene-to-UniProt mapping data frame. 
#' Rather than generating new CPI-CC IDs, it retains the existing IDs (assumed to be in the first column of each input data frame). The function also 
#' removes duplicated interactions that have the exact same sets of ligands and receptors, as well as the same interaction type.
#'
#' @param interaction_type_df A \code{data.frame} with columns \code{c("cpi_cc_id", "interaction_type")}.
#' @param ligands_df A \code{data.frame} where the first column is \code{cpi_cc_id} and subsequent columns are ligand genes.
#' @param receptors_df A \code{data.frame} where the first column is \code{cpi_cc_id} and subsequent columns are receptor genes.
#' @param g_to_u A \code{data.frame} with columns \code{c("gene_name", "uniprot")}.
#'
#' @details
#' \enumerate{
#'   \item Merges the provided data frames by \code{cpi_cc_id} into a single data frame.
#'   \item Builds a "signature" for each interaction by sorting the ligand genes, sorting the receptor genes, 
#'         and combining them with the interaction type.
#'   \item Removes duplicate rows with matching signatures, retaining only the first occurrence.
#'   \item Preserves the user-provided \code{cpi_cc_id} rather than assigning new random IDs.
#' }
#'
#' @return A list with three elements:
#' \describe{
#'   \item{db.diff}{A \code{data.frame} containing diffusion-dependent interactions 
#'                  (ligands, receptors, and factor mappings from UniProt IDs to CPI-IDs).}
#'   \item{db.cont}{A \code{data.frame} containing contact-dependent interactions, similarly structured.}
#'   \item{intr.index}{A \code{data.frame} with columns \code{"id_cp_interaction"}, 
#'                    \code{"partner_a"}, and \code{"partner_b"}.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- process_LR_data(
#'   interaction_type_df,  # data.frame with columns [cpi_cc_id, interaction_type]
#'   ligands_df,           # data.frame with first column = cpi_cc_id, subsequent columns = ligand genes
#'   receptors_df,         # data.frame with first column = cpi_cc_id, subsequent columns = receptor genes
#'   g_to_u                # data.frame with columns [gene_name, uniprot]
#' )
#' }
#'
#' @export
formatLRDB <- function(interaction_type_df, ligands_df, receptors_df, g_to_u) {
  ############################
  # 1) Basic Merging (base R)
  ############################
  
  # For consistency, rename columns so that the first column is cpi_cc_id
  # and the second column in the interaction_type_df is "interaction_type".
  colnames(interaction_type_df)[1:2] <- c("cpi_cc_id", "interaction_type")
  colnames(ligands_df)[1] <- "cpi_cc_id"
  colnames(receptors_df)[1] <- "cpi_cc_id"
  
  # Manually merge data frames by cpi_cc_id using base R's merge.
  merged_df <- merge(interaction_type_df, ligands_df, by = "cpi_cc_id", all.x = TRUE)
  merged_df <- merge(merged_df, receptors_df, by = "cpi_cc_id", all.x = TRUE)
  
  # Identify columns for ligands and receptors
  ligand_cols <- colnames(ligands_df)[-1]  # all except first col
  receptor_cols <- colnames(receptors_df)[-1]
  
  #####################################
  # 2) Gather ligand/receptor genes row by row,
  #    build a signature, then remove duplicates.
  #####################################
  
  n_rows <- nrow(merged_df)
  ligand_genes_raw_list <- vector("list", n_rows)
  receptor_genes_raw_list <- vector("list", n_rows)
  signature_vec <- character(n_rows)
  
  gather_non_na <- function(row_data, cols) {
    vals <- row_data[cols]
    vals <- vals[!is.na(vals)]
    as.character(vals)
  }
  
  for (i in seq_len(n_rows)) {
    row_data <- merged_df[i, ]
    
    # gather ligand genes
    lig_vec <- gather_non_na(row_data, ligand_cols)
    # gather receptor genes
    rec_vec <- gather_non_na(row_data, receptor_cols)
    
    ligand_genes_raw_list[[i]] <- lig_vec
    receptor_genes_raw_list[[i]] <- rec_vec
    
    # sort them for signature
    lig_sorted <- paste(sort(lig_vec), collapse = "|")
    rec_sorted <- paste(sort(rec_vec), collapse = "|")
    type_val <- as.character(row_data$interaction_type)
    
    signature_vec[i] <- paste(lig_sorted, rec_sorted, type_val, sep = "##")
  }
  
  # Mark duplicates (keep the first occurrence)
  dup_mask <- duplicated(signature_vec)
  keep_idx <- !dup_mask
  
  merged_df <- merged_df[keep_idx, ]
  ligand_genes_raw_list <- ligand_genes_raw_list[keep_idx]
  receptor_genes_raw_list <- receptor_genes_raw_list[keep_idx]
  
  ##########################
  # 3) Build output structures
  ##########################
  
  # We'll define a helper function to map gene -> uniprot from g_to_u
  gene_to_uniprot <- function(gene) {
    match_idx <- match(gene, g_to_u$gene_name)
    if (!is.na(match_idx)) {
      g_to_u$uniprot[match_idx]
    } else {
      NA
    }
  }
  
  # browser()
  # create empty structures for diffusion- and contact-dependent data
  db.diff <- list(ligands = factor(), receptors = factor(), combined = factor())
  db.cont <- list(ligands = factor(), receptors = factor(), combined = factor())
  
  intr.index <- data.frame(
    id_cp_interaction = character(0),
    partner_a = character(0),
    partner_b = character(0),
    stringsAsFactors = FALSE
  )
  
  # helper function to append factors
  assign_factors <- function(uniprot_vec, cpi_id, existing_factor) {
    new_factor <- factor(
      rep(cpi_id, length(uniprot_vec)),
      levels = unique(c(levels(existing_factor), cpi_id))
    )
    names(new_factor) <- uniprot_vec
    c(existing_factor, new_factor)
  }
  
  # iterate over each deduplicated row
  n_dedup <- nrow(merged_df)
  for (i in seq_len(n_dedup)) {
    cpi_id_val <- as.character(merged_df$cpi_cc_id[i])
    type_val <- as.character(merged_df$interaction_type[i])
    
    ligand_genes <- ligand_genes_raw_list[[i]]
    receptor_genes <- receptor_genes_raw_list[[i]]
    
    # map genes to uniprot
    ligand_uniprot <- sapply(ligand_genes, gene_to_uniprot)
    receptor_uniprot <- sapply(receptor_genes, gene_to_uniprot)
    
    # Build partner_a, partner_b
    # If there's only 1 ligand gene, omit "-ligand" suffix
    if (length(ligand_genes) == 1) {
      partner_a_val <- paste(ligand_genes, collapse = "_")
    } else {
      partner_a_val <- paste(paste(ligand_genes, collapse = "_"), "-ligand", sep = "")
    }
    
    # If there's only 1 receptor gene, omit "-receptor" suffix
    if (length(receptor_genes) == 1) {
      partner_b_val <- paste(receptor_genes, collapse = "_")
    } else {
      partner_b_val <- paste(paste(receptor_genes, collapse = "_"), "-receptor", sep = "")
    }
    
    intr.index <- rbind(
      intr.index,
      data.frame(
        id_cp_interaction = cpi_id_val,
        partner_a = partner_a_val,
        partner_b = partner_b_val,
        stringsAsFactors = FALSE
      )
    )
    
    # assign to db.diff or db.cont based on interaction type
    if (length(ligand_uniprot) > 0) {
      if (type_val == "diffusion-dependent") {
        db.diff$ligands <- assign_factors(ligand_uniprot, cpi_id_val, db.diff$ligands)
      } else {
        db.cont$ligands <- assign_factors(ligand_uniprot, cpi_id_val, db.cont$ligands)
      }
    }
    
    if (length(receptor_uniprot) > 0) {
      if (type_val == "diffusion-dependent") {
        db.diff$receptors <- assign_factors(receptor_uniprot, cpi_id_val, db.diff$receptors)
      } else {
        db.cont$receptors <- assign_factors(receptor_uniprot, cpi_id_val, db.cont$receptors)
      }
    }
    
    # combined factor has both ligand & receptor uniprots
    combined_uniprots <- unique(c(ligand_uniprot, receptor_uniprot))
    if (length(combined_uniprots) > 0) {
      if (type_val == "diffusion-dependent") {
        db.diff$combined <- assign_factors(combined_uniprots, cpi_id_val, db.diff$combined)
      } else {
        db.cont$combined <- assign_factors(combined_uniprots, cpi_id_val, db.cont$combined)
      }
    }
  }
  
  # return final structures
  list(
    db.diff = db.diff,
    db.cont = db.cont,
    intr.index = intr.index
  )
}

filterGene <- function(dge.raw, gene_to_uniprot, thresh = 50){
  # find intr gene index
  rownames(dge.raw) <- toupper(rownames(dge.raw))
  # uniprot.genes = rownames(dge.raw)[rownames(dge.raw) %in% gene_to_uniprot$gene_name] # 738/22683
  # intr.genes.index = which(rownames(dge.raw) %in% gene_to_uniprot$gene_name)
  intr.genes.index <- rownames(dge.raw) %in% gene_to_uniprot$gene_name
  # message("Number of genes in the database: ", length(intr.genes.index))
  message("Number of genes in the database: ", sum(intr.genes.index))

  # find mito gene index
  # mt_idx <- grep("MT-",rownames(dge.raw))
  mt_idx <- startsWith(rownames(dge.raw), "MT-")

  # find low quality gene index
  # low_quality_idx <- which(rowSums(dge.raw) < thresh)
  low_quality_idx <- rowSums(dge.raw) < thresh

  # message("Number of low-quality intr genes: ", sum(intr.genes.index %in% low_quality_idx))
  message("Number of low-quality intr genes: ", sum(intr.genes.index & low_quality_idx))

  # return non-intr filtered gene index
  # fin.index = setdiff(union(low_quality_idx, mt_idx), intr.genes.index)
  gene.del <- (low_quality_idx | mt_idx) & (!intr.genes.index)
  # return(fin.index)
  return(!gene.del)
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
  message("- Number of valid intrs: ", sum(if.all.in), " out of ", length(if.all.in), " from database.")
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
    warning("Real lr-score mtx has fewer features than NULL lr-score mtx.")
    warning("Removing extra", ncol(null.i.mtx) - ncol(i.mtx), "interactions from NULL lr-score mtx.")
    null.i.mtx <- null.i.mtx[, colnames(i.mtx)]
  }

  if (ncol(i.mtx) > ncol(null.i.mtx)){
    warning("NULL lr-score mtx has fewer features than Real lr-score mtx.")
    warning("Removing extra", ncol(i.mtx) - ncol(null.i.mtx), "interactions from Real lr-score mtx.")
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

graphSpatialFDRNew <- function(
    nn.graph,
    pval.mtx,
    spatial.coords = NULL,
    weighting = "DT",
    k = NULL,
    reciprocal = TRUE
) {
  if (ncol(nn.graph) != nrow(pval.mtx)){
    stop("Cell numbers in pvalue mtx and neighbor graph do not match!")
  }

  if (weighting == "DT"){
    # for DT, we use the max distance within the nb as the weight
    t.connect <- rep(NA, ncol(nn.graph))
    for (j in seq(ncol(nn.graph))) {
      weights <- nn.graph[, j]
      # Do not include cell without neighbors
      if (sum(weights) == 0) next
      # Remove the self-weight
      weights <- weights[-j]
      t.connect[j] <- max(weights)
    }
    t.connect <- t.connect[!is.na(t.connect)]
  } else if (weighting == "EB"){
    stop("EB weighting not developed yet.")
  } else if (weighting == "KNN"){
    stop("KNN weighting not developed yet.")
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

  return(padj.mtx)
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
    stop("Cell numbers in pvalue mtx and neighbor factor do not match!")
  }

  if (weighting == "DT"){
    # for DT, we use the max distance within the nb as the weight
    # message("Using DT weighting.")
    # t.connect <- sapply(levels(nn.dist), function(x){
    #     nb.dist = as.numeric(names(nn.dist[nn.dist == x]))
    #     return(max(nb.dist))
    # })

    t.connect <- sapply(
      split(as.numeric(names(nn.dist)), nn.dist),
      function(x){
        # return(max(x))
        return(max(x[-which(x == max(x))[1]]))
      }
    )

  } else if (weighting == "EB"){
    stop("EB weighting not developed yet.")
  } else if (weighting == "KNN"){
    stop("KNN weighting not developed yet.")
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
  message("- Number of high quality interactions: ", length(res.list.hq))

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

  message("Using ", k, " neighbors for imputation...")
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
  # construct the weight matrix direcly
  i <- c(nn.idx) # join each column and convert to int vector
  j <- as.integer(rep_each_cpp(ncol(nn.idx), k)) # rep (1,2,3,...,n), each k times
  x <- c(nn.dist)
  dist.graph <- sparseMatrix(i, j, x = x) # every COLUMN represents a cell's distances to its k nearest neighbors
  # all.equal(which(dist.graph[,1] != 0), sort(nn.idx[,1]))

  # dist.graph <- quickNorm(dist.graph)

  #### this matrix is NOT symmetric!!!!!
  rm(i, j, x, nn.idx, nn.dist); gc()

  message("Imputing...")

  # dge.intr.imp <- dge.intr %*% dist.graph
  # dge.mix.imp <- rbind(dge.intr.imp, dge.other)
  data.imp <- Matrix::crossprod(data, dist.graph)
  data.imp <- as(data.imp, "CsparseMatrix")

  dimnames(data.imp) <- list(
    colnames(data),
    rownames(data)
  )

  # sum(dge.raw!=0)/length(dge.raw) # density 98%

  return(data.imp)
}


getIntrNames <- function(object, intr){
  ligandNames <- getLigandNames(object, intr)
  receptorNames <- getReceptorNames(object, intr)
  intrName <- paste0(ligandNames, "-", receptorNames)
  names(intrName) <- intr
  return(intrName)
}

getLigandNames <- function(object, intr){
  inter.index <- object@intr.valid$intr.index
  sapply(intr, function(x) {
    pair.index <- which(inter.index$id_cp_interaction == x)

    ligand.name <- inter.index[pair.index, 4]
    if (ligand.name == "") { # if complex
      ligand.name <- inter.index[pair.index, 2]
    }
    gsub("_HUMAN", "", ligand.name)
  })
}

getReceptorNames <- function(object, intr){
  inter.index <- object@intr.valid$intr.index
  sapply(intr, function(x) {
    pair.index <- which(inter.index$id_cp_interaction == x)
    receptor.name <- inter.index[pair.index, 5]
    if (receptor.name == "") {
      receptor.name <- inter.index[pair.index, 3]
    }
    gsub("_HUMAN", "", receptor.name)
  })
}


getIntrRanks <- function(object, intr.list, slot.use=NULL, signif.use=NULL){
  slot.use <- .checkSlotUse(object, slot.use)
  signif.use <- .checkSignifUse(object, signif.use, slot.use)
  intr.list <- .checkIntrAvail(object, intr.list, slot.use, signif.use)

  rank.list <- lapply(intr.list, function(intr) {
    match(intr, names(object@lrscore[[slot.use]]@res.list[[signif.use]]))
  })
  names(rank.list) <- intr.list

  return(rank.list)
}

getCPIs <- function(object, intr.idx, slot.use=NULL, signif.use=NULL){
  if (!is.numeric(intr.idx)) {
    stop("intr.idx must be numeric")
  }

  slot.use <- .checkSlotUse(object, slot.use)
  signif.use <- .checkSignifUse(object, signif.use, slot.use)
  intr.list <- names(object@lrscore[[slot.use]]@res.list[[signif.use]])

  if (max(intr.idx) > length(intr.list)) {
    stop("intr.idx contains index out of range")
  }

  intr.list <- intr.list[intr.idx]

  return(intr.list)
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

increase_rows <- function(x, y) {
  # Get the row names of the matrix
  rownames_y <- rownames(y)

  # Determine which rows are missing
  missing_rows <- setdiff(x, rownames_y)

  # Add missing rows filled with 0s
  if (length(missing_rows) > 0) {
    n_missing_rows <- length(missing_rows)
    n_cols <- ncol(y)
    y_new <- rbind(y, matrix(0, nrow = n_missing_rows, ncol = n_cols,
                             dimnames = list(missing_rows, NULL)))
  } else {
    y_new <- y
  }

  # Reorder rows to match x
  y_new[x, ]
}


# mat_a is a matrix, mat_b is another matrix
increase_dims <- function(mat_a, mat_b) {
  mat_a <- increase_rows(rownames(mat_b), mat_a)
  mat_a <- increase_columns(colnames(mat_b), mat_a)
  return(mat_a)
}



# Shuffling strategy 1: find a random set of cells to be the new neighbors; same weights, but different cells.
# shuffle the values in very column of a sparse matrix, keep the column index the same
shuffle_sp_mat_col <- function(mat) {
  new.i.len <- diff(mat@p) # num of non-zero elements in each column
  nr <- nrow(mat)

  new.i <- lapply(new.i.len, function(n){
    return(sample(0:(nr-1), n))
  })

  i.list <- unlist(new.i, use.names = FALSE)

  mat@i <- i.list

  return(mat)

}


# Shuffling strategy 2: keep the neighbors and corresponding weights the same, but infer the signals from a
# random set of cells.
shuffleEdgeRandomNB <- function(mat) {
  size <- ncol(mat)

  new.list <- lapply(seq_len(size), function(j){
    # get the @i indices for column j
    start <- mat@p[j] + 1
    end   <- mat@p[j+1]
    indx  <- start:end
    len <- length(indx)

    # sample a new set of cells to be the new sender cells
    new.i <- sample(size, len)
    keep.idx <- which(mat@i[indx] %in% new.i)

    # if no new sender cells are found, return NULL
    if (length(keep.idx) == 0){
      return(list(x = NULL, i = NULL, j = NULL))
    }

    # if new sender cells are found, return the new @x, @i, @j
    keep.x <- mat@x[keep.idx]
    keep.i <- mat@i[keep.idx]
    keep.j <- rep(j, length(keep.idx))

    return(list(x = keep.x, i = keep.i, j = keep.j))
  })

  new.x <- unlist(lapply(new.list, function(x) x$x ))
  new.i <- unlist(lapply(new.list, function(x) x$i ))
  new.j <- unlist(lapply(new.list, function(x) x$j )) - 1 # 0-based index

  new.mat <- sparseMatrix(i = new.i, x = new.x, j = new.j, dims = c(size, size),
                          dimnames = list(NULL, colnames(mat)), index1=F) # 0-based index

  return(new.mat)
}


to_mean <- function(mat) {
  new.x <- rep(1/(diff(mat@p)), diff(mat@p))
  mat@x <- new.x
  return(mat)
}


# # helper function for plotting
# sample_null_dge <- function(object, imp.slot){
#     n.cells <- ncol(object@counts)
#     n.null <- ncol(object@imputation[[imp.slot]]@imp.data.null[[1]])
#     need.times <- ceiling(n.cells/n.null)
#     each.size <- ceiling(n.null/need.times)

#     sel.index <- sample(n.null, each.size)
#     null.data <- lapply(object@imputation[[imp.slot]]@imp.data.null, function(x){
#         return( x[ ,sel.index] )
#     })

#     null.data <- combine_sparse_rows(null.data); gc()
#     null.data <- null.data[, 1:n.cells]

#     dimnames(null.data) <- dimnames(object@counts)

#     return(null.data)
# }
.checkImpUse <- function(object, imp.use = NULL) {
  avail <- names(object@imputation)
  avail <- avail[avail != "default"]
  if (is.null(imp.use)) {
    imp.use <- object@imputation[["default"]]
  }
  if (!imp.use %in% avail) {
    stop("Specified imputation (`imp.use = '", imp.use, "'`) is not available. ",
         "Available ones are: ",
         paste(paste0('"', avail, '"'), collapse = ", "))
  }
  return(imp.use)
}

.checkSlotUse <- function(object, slot.use = NULL, velo = FALSE) {
  avail <- if (velo) names(object@lrvelo) else names(object@lrscore)
  avail <- avail[avail != "default"]
  if (is.null(slot.use)) {
    slot.use <- if (velo) object@lrvelo[["default"]] else object@lrscore[["default"]]
  }
  infoType <- ifelse(velo, "VeloLRScore", "LRScore")
  if (length(avail) == 0 || is.null(slot.use)) {
    funName <- ifelse(velo, "inferIntrVelo()", "inferIntrScore()")
    stop("No ", infoType, " available or specified. ",
         "Please run `", funName, "` first.")
  }
  if (!slot.use %in% avail) {
    stop("Specified ", infoType, " (`slot.use = '", slot.use, "'`) is not available. ",
         "Please choose from: ",
         paste(paste0('"', avail, '"'), collapse = ", "))
  }
  return(slot.use)
}

.checkSignifUse <- function(object, signif.use, slot.use) {
  avail <- names(object@lrscore[[slot.use]]@res.list)
  if (is.null(signif.use)) {
    signif.use <- avail[length(avail)]
    # message("`signif.use` not specified, using the latest available option: ",
    #         signif.use)
  }

  if (!signif.use %in% avail) {
    stop("Invalid `signif.use`. Available options: ",
         paste(avail, collapse = ", "))
  }
  return(signif.use)
}

.checkIntrAvail <- function(object, intr, slot.use, signif.use) {
  intrList <- showIntr(object, slot.use, signif.use)
  if (!all(intr %in% intrList)) {
    nf <- intr[!intr %in% intrList]
    stop("Specified interaction not found, see available options with ",
         "`showIntr(object)`. \nUnavailable ones: ",
         paste(nf, collapse = ", "))
  }
  return(intr)
}

.checkArgLen <- function(arg, len) {
  argname <- deparse(substitute(arg))
  if (!is.null(arg)) {
    if (!is.character(arg)) {
      stop("`", argname, "` has to be character vector/scalar.")
    }
    if (length(arg) != len) {
      stop("`", argname, "` has to be a vector of length ", len)
    }
  }
  return(arg)
}

findImpByMethod <- function(object, method) {
  nn.avail <- names(object@imputation)
  nn.avail <- nn.avail[nn.avail != "default"]
  nn.target <- NULL
  for (nn in nn.avail) {
    if (object@imputation[[nn]]@method == method) {
      nn.target <- nn
      break
    }
  }
  return(nn.target)
}

filterGene <- function(dge.raw, gene_to_uniprot, thresh = 50){
    # find intr gene index
    rownames(dge.raw) = toupper(rownames(dge.raw))
    # uniprot.genes = rownames(dge.raw)[rownames(dge.raw) %in% gene_to_uniprot$gene_name] # 738/22683
    intr.genes.index = which(rownames(dge.raw) %in% gene_to_uniprot$gene_name)
    cat(paste0("Number of genes in the database: ", length(intr.genes.index), "\n"))

    # find mito gene index
    mt_idx <- grep("MT-",rownames(dge.raw))

    # find low quality gene index
    low_quality_idx <- which(rowSums(dge.raw) < thresh)

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


checkIntr <- function(dge, intr.db){
    intr.fac = intr.db[['combined']]
    lig.fac = intr.db[['ligands']]
    recep.fac = intr.db[['receptors']]

    # filter all empty intrs
    if.all.in = sapply(levels(intr.fac), function(intr) {
        intr.cmp = names(intr.fac[intr.fac == intr, drop = T])
        return(sum(intr.cmp %in% rownames(dge)) == length(intr.cmp))
    })
    cat(paste0("Number of valid intrs: ", sum(if.all.in), " / ", length(if.all.in), "\n"))
    levels.valid = levels(intr.fac)[if.all.in]
    # se.cmp.fac = cmp.fac[cmp.fac %in% levels(cmp.fac)[if.pair.in], drop = T] # 294 / 1396 levels

    intr.fac = intr.fac[intr.fac %in% levels.valid, drop = T]
    lig.fac = lig.fac[lig.fac %in% levels.valid, drop = T]
    recep.fac = recep.fac[recep.fac %in% levels.valid, drop = T]

    # convert lig and recep
    names(lig.fac) = unname(sapply(names(lig.fac), function(gene){which(rownames(dge) == gene)}))
    names(recep.fac) = unname(sapply(names(recep.fac), function(gene){which(rownames(dge) == gene)}))

    return(list(
        "combined" = intr.fac,
        "ligands" = lig.fac,
        "receptors" = recep.fac
    ))
}
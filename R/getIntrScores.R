getIntrValue <- function(
        object,
        intr,
        type = c("ligand", "ligand_ori", "ligand_null",
                 "receptor", "receptor_ori", "receptor_null",
                 "score", "score_null"),
        slot.use = NULL,
        signif.use = NULL,
        ...
) {
    slot.use <- .checkSlotUse(object, slot.use)
    signif.use <- .checkSignifUse(object, signif.use, slot.use)
    intr <- .checkIntrAvail(object, intr, slot.use, signif.use)

    score.obj <- object@lrscore[[slot.use]]
    intr.db <- object@intr.valid[[score.obj@intr.slot]]

    # Identify which interactions are queried
    # if (is.null(intr)) {
    #     index.len <- length(score.obj@res.list[[signif.use]])
    #     # By default use the first 20 and last 10 if more than 30, otherwise all
    #     intr <- intersect(union(seq(20), seq(index.len - 9, index.len)),
    #                           seq(index.len))
    # }

    # Initialize data.frame with location
    cells.loc <- as.data.frame(object@cells.loc)
    nCells <- nrow(cells.loc)

    # Prepare data to retrieve from
    res.list <- score.obj@res.list[[signif.use]]
    res.intr.list <- names(res.list)

    lig.slot <- score.obj@lig.slot
    # use normalized ligand data for visualization
    # dge.lig <- object@imputation[[lig.slot]]@imp.data
    dge.lig <- normCounts(object, imp.use = lig.slot, verbose = FALSE)
    recep.slot <- score.obj@recep.slot
    # use normalized receptor data for visualization
    # dge.recep <- object@imputation[[recep.slot]]@imp.data
    dge.recep <- normCounts(object, imp.use = recep.slot, verbose = FALSE)
    dge.raw <- object@counts
    score.mtx <- score.obj@score

    sample.index <- sample(ncol(score.obj@lig.null), nCells)
    null.dge.gau <- score.obj@lig.null[, sample.index]
    null.dge.dt <- score.obj@recep.null[, sample.index]
    rownames(null.dge.dt) <- rownames(null.dge.gau) <- rownames(dge.lig)
    colnames(null.dge.dt) <- colnames(null.dge.gau) <- colnames(dge.lig)
    null.lrscore.mtx <- score.obj@score.null[sample.index, ]
    rownames(null.lrscore.mtx) <- rownames(score.mtx)

    # For each interaction, update the data.frame with queried information
    intrList <- list()
    for (intrx in intr) {
        scoreDF <- cells.loc
        for (info in type) {
            if (!startsWith(info, "score")) {
                if (startsWith(info, "ligand")) {
                    idx <- names(intr.db[[2]][intr.db[[2]] == intrx])
                    if (info == "ligand") value <- dge.lig[idx, ]
                    if (info == "ligand_ori") value <- dge.raw[idx, ]
                    if (info == "ligand_null") value <- null.dge.gau[idx, ]
                } else if (startsWith(info, "recep")) {
                    idx <- names(intr.db[[3]][intr.db[[3]] == intrx])
                    if (info == "receptor") value <- dge.recep[idx, ]
                    if (info == "receptor_ori") value <- dge.raw[idx, ]
                    if (info == "receptor_null") value <- null.dge.dt[idx, ]
                }
                if (!is.null(nrow(value))) value = Matrix::colSums(value)
                scoreDF[[info]] <- value
            } else {
                scoreDF[[info]] <- 0
                if (info == "score") {
                    scoreDF[rownames(score.mtx), info] <-
                        score.mtx[, intrx]
                }
                if (info == "score_null") {
                    scoreDF[rownames(null.lrscore.mtx), info] <-
                        null.lrscore.mtx[, intrx]
                }
            }
        }
        # Identify ligand-recepter names
        pair.index <- which(object@intr.valid$intr.index$id_cp_interaction ==
                                intrx)

        ligand.name <- object@intr.valid$intr.index[pair.index, 4]
        if (ligand.name == "") {
            ligand.name <- object@intr.valid$intr.index[pair.index, 2]
        }
        ligand.name <- gsub("_HUMAN", "", ligand.name)

        receptor.name <- object@intr.valid$intr.index[pair.index, 5]
        if (receptor.name == "") {
            receptor.name <- object@intr.valid$intr.index[pair.index, 3]
        }
        receptor.name <- gsub("_HUMAN", "", receptor.name)
        intr.name <- paste0(ligand.name, "-", receptor.name)
        intrList[[intr.name]] <- scoreDF
    }
    return(intrList)
}

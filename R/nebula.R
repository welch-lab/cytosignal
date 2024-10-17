#' mergedCytoSignal Class for cross-dataset analysis of spatial cell-cell
#' communication
#' @slot lrscore A sparse matrix (dgCMatrix class) with LR scores from all
#' dataset merged into one.
#' @slot metadata A data.frame with location level metadata covariates.
setClass(
    Class = "mergedCytoSignal",
    representation = representation(
        metadata = "data.frame",
        diff.lrscore = "dgCMatrix",
        cont.lrscore = "dgCMatrix"
    ),
    prototype = prototype(
        metadata = NULL,
        diff.lrscore = NULL,
        cont.lrscore = NULL
    )
)

#' Show method for mergedCytoSignal object
#' @param object A mergedCytoSignal object
#' @return Show a message with object summary, no return value
#' @export
setMethod(
    f = "show",
    signature = "mergedCytoSignal",
    definition = function(object) {
        cat("An object of class mergedCytoSignal\n")
        cat("  - Number of spots: ", nrow(object@metadata), "\n")
        cat("  - Number of datasets:", nlevels(object@metadata$dataset), "\n")
        ncovar <- ncol(object@metadata)
        if ("x" %in% colnames(object@metadata)) ncovar <- ncovar - 1
        if ("y" %in% colnames(object@metadata)) ncovar <- ncovar - 1
        cat("  - Number of covariates: ", ncovar, "\n")
        cat("  - Number of diffusible ligand-receptor pairs: ", nrow(object@diff.lrscore), "\n")
        cat("  - Number of contact-dependent ligand-receptor pairs: ", nrow(object@cont.lrscore), "\n")
    }
)

#' dim Method for mergedCytoSignal object
#' @description Only returns nrow(x) which indicates the number of spots/cells
#' @param x A mergedCytoSignal object
#' @return A numeric vector with the number of spots/cells and NA
#' @export
#' @method dim mergedCytoSignal
dim.mergedCytoSignal <- function(x) {
    c(nrow(x@metadata), NA)
}

#' Dataset names of a mergedCytoSignal object
#' @description
#' Allows getting and setting dataset names of a mergedCytoSignal object.
#' @param x A mergedCytoSignal object
#' @param value A character vector of up to the same length as x to rename
#' datasets. or NULL to reset to default names.
#' @return
#' For \code{names}, a character vector of dataset names.
#'
#' For \code{names<-}, a mergedCytoSignal object with dataset names updated
#' in both metadata and the rownames of lrscore matrices
#' @rdname names.mergedCytoSignal
#' @export
#' @method names mergedCytoSignal
names.mergedCytoSignal <- function(x) {
    levels(x@metadata$dataset)
}

#' @rdname names.mergedCytoSignal
#' @export
#' @method names<- mergedCytoSignal
`names<-.mergedCytoSignal` <- function(x, value) {
    if (!is.null(value)) {
        if (length(value) != nlevels(x@metadata$dataset)) {
            stop("Length of `value` must match the number of datasets.")
        }
    } else {
        value <- paste0("dataset", seq(nlevels(x@metadata$dataset)))
    }
    levels(x@metadata$dataset) <- value
    newID <- paste0(x@metadata$dataset, "_", x@metadata$originalID)
    rownames(x@metadata) <- newID
    colnames(x@diff.lrscore) <- newID
    colnames(x@cont.lrscore) <- newID
    x
}

#' Extract metadata variable from mergedCytoSignal object
#' @description
#' Fast expression for getting and setting metadata variables in a
#' \linkS4class{mergedCytoSignal} object.
#' @rdname Extract.mergedCytoSignal
#' @export
#' @method $ mergedCytoSignal
`$.mergedCytoSignal` <- function(x, name) {
    stats::setNames(x@metadata[[name]], rownames(x@metadata))
}

#' @rdname Extract.mergedCytoSignal
#' @export
#' @method $<- mergedCytoSignal
`$<-.mergedCytoSignal` <- function(x, name, value) {
    x@metadata[[name]] <- value
    x
}

#' @export
#' @method .DollarNames mergedCytoSignal
.DollarNames.mergedCytoSignal <- function(x, pattern = "") {
    grep(pattern, colnames(x@metadata), value = TRUE)
}

.valid.mergedCytoSignal <- function(object) {
    # Check if spot IDs are identical from difflrscore, contlrscore and metadata
    metarownames <- rownames(object@metadata)
    diffcolnames <- colnames(object@diff.lrscore)
    contcolnames <- colnames(object@cont.lrscore)
    if (!identical(metarownames, diffcolnames) ||
        !identical(metarownames, contcolnames)) {
        return("Spot/cell IDs in metadata, diff.lrscore and cont.lrscore do not match.")
    }

    # Check if the spot IDs equal to "{dataset}_{originalID}"
    if (!identical(metarownames, paste0(object$dataset, "_", object$originalID))) {
        return("Spot/cell IDs in metadata do not follow the format '{dataset}_{originalID}'.")
    }
}

setValidity("mergedCytoSignal", .valid.mergedCytoSignal)


#' Create mergedCytoSignal object
#' @description
#' This function accepts a list of CytoSignal objects together with a metadata
#' \code{data.frame} that provides dataset level metadata covariates, which will
#' be required for the cross-dataset analysis.
#'
#' The function will preprocess each CytoSignal object with cross-dataset
#' analysis specific steps, including imputing diffusible LR values using
#' Gaussion Epsilon model with weight 1.
#' @param objList A list of CytoSignal object.
#' @param metadata A data.frame object. Default \code{NULL}. If provided, must
#' have as many rows as the length of \code{objList} with rows in matching
#' order to \code{objList}.
#' @param name.by A column name in \code{metadata} indicating where to find the
#' dataset names. Default \code{NULL} names datasets by "dataset1", "dataset2"
#' and etc, and will create column "dataset" in the metadata.
#' @param counts.thresh A number indicating the threshold for removing low-quality
#' spots. Default 0.
#' @param scale.factor A number indicating 1 spatial coord unit equals to how 
#' many µm. 
#' @return A \linkS4class{mergedCytoSignal} object preprocessed.
#' @export
mergeCytoSignal <- function(
        objList,
        metadata = NULL,
        name.by = NULL,
        counts.thresh = 0,
        scale.factor = NULL
) {
    # Sanity checks on objList
    for (i in seq_along(objList)) {
        obj <- objList[[i]]
        if (!inherits(obj, "CytoSignal")) {
            stop(sprintf("Dataset %d (%s) is not a CytoSignal object", i, names(objList)[i]))
        }
        if (length(obj@clusters) == 0) {
            stop(sprintf("No clusters added in dataset %d (%s)", i, names(objList)[i]))
        }
        if (any(dim(obj@cells.loc) == 0)) {
            stop(sprintf("No cell locations added in dataset %d (%s)", i, names(objList)[i]))
        }
    }

    # Check metadata
    if (is.null(metadata)) {
        sid <- paste0("dataset", seq_along(objList))
        metadata <- data.frame(
            dataset = factor(sid, levels = sid),
        )
    } else {
        if (nrow(metadata) != length(objList)) {
            stop("`metadata` contains different number of rows then number of ",
                 "datasets in `objList`.")
        }
        if (is.null(name.by)) {
            sid <- paste0("dataset", seq_along(objList))
            if ("dataset" %in% colnames(metadata)) {
                stop('Column name "dataset" is preserved by downstream ',
                     'analysis, please rename/remove it, or set `name.by = ',
                     '"dataset"`.')
            }
        } else {
            if (!(name.by %in% colnames(metadata))) {
                stop(sprintf('Column name "%s" not found in metadata.', name.by))
            }
            if (name.by != "dataset" &&
                "dataset" %in% colnames(metadata)) {
                stop('Column name "dataset" is preserved by downstream ',
                     'analysis, please rename/remove it, or set `name.by = ',
                     '"dataset"`.')
            }
            sid <- as.character(metadata[[name.by]])
        }
        metadata$dataset <- factor(sid, levels = sid)

        preserved <- c("originalID", "clusters", "x", "y")
        for (p in preserved) {
            if (p %in% colnames(metadata)) {
                stop(sprintf('Column name "%s" is preserved by downstream analysis, please rename or remove it.', p))
            }
        }
    }
    
    if (is.null(scale.factor)) {
        stop("Please provide `scale.factor`.")
    }

    # broadcast dataset-level metadata to location level
    spotIDs <- lapply(objList, function(x) colnames(x@raw.counts))
    sizes <- lengths(spotIDs)
    meta.full <- metadata[rep(seq_along(objList), sizes), , drop = FALSE]
    meta.full$originalID <- unlist(spotIDs)
    rownames(meta.full) <- paste0(meta.full$dataset, "_", meta.full$originalID)

    clusterUnion <- Reduce(union, lapply(objList, function(x) levels(x@clusters)))
    meta.full$clusters <- factor(NA, levels = clusterUnion)
    meta.full$x <- NA
    meta.full$y <- NA

    diffList <- list()
    contList <- list()
    for (i in seq_along(objList)) {
        obj <- objList[[i]]
        ids <- colnames(obj@raw.counts)
        newids <- paste0(sid[i], "_", ids)

        # Insert single-loc level metadata
        meta.full[newids, "clusters"] <- obj@clusters
        meta.full[newids, "x"] <- obj@cells.loc[,1]
        meta.full[newids, "y"] <- obj@cells.loc[,2]

        # Preprocessing
        message(Sys.time(), " - Preprocessing dataset: ", i)
        # Hard coded parameters, pay attention to these
        obj <- addIntrDB(obj, g_to_u, db.diff, db.cont, inter.index)
        obj <- removeLowQuality(obj, counts.thresh = counts.thresh)
        obj <- changeUniprot(obj)
        obj <- inferEpsParams(obj, scale.factor = scale.factor)
        obj <- findNN(obj, diff.weight = 1)
        obj <- imputeLR(obj)

        obj <- inferScoreLR(obj, lig.slot = "GauEps", recep.slot = "Raw",
                            norm.method = "none", intr.db.name = "diff_dep")
        diffscore <- obj@lrscore[["GauEps-Raw"]]@score
        diff.dimnames <- list(
            newids,
            unname(getIntrNames(obj, colnames(diffscore)))
        )
        diffscore <- cleanLRscore_sparse_cpp(
            diffscore@i, diffscore@p, diffscore@x, nrow(diffscore), ncol(diffscore)
        )
        dimnames(diffscore) <- diff.dimnames
        diffList[[i]] <- diffscore

        obj <- inferScoreLR(obj, lig.slot = "DT", recep.slot = "Raw",
                            norm.method = "none", intr.db.name = "cont_dep")
        contscore <- obj@lrscore[["DT-Raw"]]@score
        cont.dimnames <- list(
            newids,
            unname(getIntrNames(obj, colnames(contscore)))
        )
        contscore <- cleanLRscore_sparse_cpp(
            contscore@i, contscore@p, contscore@x, nrow(contscore), ncol(contscore)
        )
        dimnames(contscore) <- cont.dimnames
        contList[[i]] <- contscore

        message()
    }

    message(Sys.time(), " - Merging LR scores")

    # Merge diff lrscore
    diffIntrIsec <- Reduce(intersect, lapply(diffList, colnames))
    diffList <- lapply(diffList, `[`, i = , j = diffIntrIsec, drop = FALSE)
    diffLR <- Reduce(rbind, diffList)

    # Merge cont lrscore
    contIntrIsec <- Reduce(intersect, lapply(contList, colnames))
    contList <- lapply(contList, `[`, i = , j = contIntrIsec, drop = FALSE)
    contLR <- Reduce(rbind, contList)

    methods::new(
        Class = "mergedCytoSignal",
        metadata = meta.full,
        diff.lrscore = t(diffLR),
        cont.lrscore = t(contLR)
    )
}

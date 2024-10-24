#' mergedCytoSignal Class for cross-dataset analysis of spatial cell-cell
#' communication
#' @description
#' An object that holds necessary information from multiple CytoSignal objects
#' that is required for performing cross-dataset analysis. Results from the
#' cross-dataset analysis will also be held in this object.
#'
#' Please use \code{\link{mergeCytoSignal}()} to create a
#' \code{mergedCytoSignal} object
#' @slot metadata A data.frame with spot/cell level metadata covariates. Column
#' name "dataset", "originalID", "clusters", "x" and "y" are reserved for
#' downstream analysis and object validity.
#' @slot diff.lrscore A sparse matrix (dgCMatrix) of LR scores for
#' diffusion-dependent interactions, across all datasets.
#' @slot diff.results A list of data.frame objects, each containing the results
#' of NEBULA analysis for a specific covariate as specified when running the
#' analysis.
#' @slot cont.lrscore A sparse matrix (dgCMatrix) of LR scores for
#' contact-dependent interactions, across all datasets.
#' @slot cont.results A list of data.frame objects, each containing the results
#' of NEBULA analysis for a specific covariate as specified when running the
#' analysis.
#' @slot cov.use A list that organizes the covariates being used for the
#' modeling. Each element is named by a variable name selected from metadata,
#' and is a list with the following structure
#' \itemize{
#'   \item{"type": "factor" or "numeric" etc. Original class of the variable}
#'   \item{"levels" (if a factor): levels of the factor that'll be found in the
#'   model.matrix outcome}
#'   \item{"range" (if numeric): range of the variable}
#' }
#' @seealso [mergeCytoSignal()] [runNEBULA()]
setClass(
    Class = "mergedCytoSignal",
    representation = representation(
        metadata = "data.frame",
        diff.lrscore = "dgCMatrix",
        diff.results = "list",
        cont.lrscore = "dgCMatrix",
        cont.results = "list",
        cov.use = "list"
    ),
    prototype = prototype(
        metadata = NULL,
        diff.lrscore = NULL,
        cont.lrscore = NULL,
        diff.results = list(),
        cont.results = list(),
        cov.use = list()
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
        cat(sprintf(
            "  - Number of datasets (%d): %s\n",
            length(unique(object@metadata$dataset)),
            paste0(paste0('"', unique(object@metadata$dataset), '"'), collapse = ", ")
        ))
        cat(sprintf(
            "  - Metadata variables (%d): %s\n",
            ncol(object@metadata),
            paste0(paste0('"', colnames(object@metadata), '"'), collapse = ", ")
        ))
        cat("  - Number of diffusible ligand-receptor pairs: ", nrow(object@diff.lrscore), "\n")
        cat("  - Number of contact-dependent ligand-receptor pairs: ", nrow(object@cont.lrscore), "\n")
        if (length(object@cov.use)) {
            all_repr <- lapply(seq_along(object@cov.use), function(i) {
                struct <- object@cov.use[[i]]
                name <- names(object@cov.use)[i]
                if (struct$type != "numeric") {
                    text <- sprintf('%s (%d levels)', name, length(struct$levels))
                } else {
                    text <- sprintf('%s (%3f - %3f)', name, struct$range[1], struct$range[2])
                }
                return(text)
            })
            cat(sprintf("  - Analyzed with covariates (%d): %s\n",
                        length(all_repr),
                        paste0(all_repr, collapse = ", ")))
        }
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

    # Check if "gene names" in result data.frame all came from lrscore rownames
    if (length(object@diff.results) > 0) {
        for (i in seq_along(object@diff.results)) {
            if (!all(object@diff.results[[i]]$interaction %in% rownames(object@diff.lrscore))) {
                return(sprintf("Unindentified interaction found from NEBULA result for diffusion-dependent interactions, in covariate %s", names(object@diff.results)[i]))
            }
        }
    }
    if (length(object@cont.results) > 0) {
        for (i in seq_along(object@cont.results)) {
            if (!all(object@cont.results[[i]]$interaction %in% rownames(object@cont.lrscore))) {
                return(sprintf("Unindentified interaction found from NEBULA result for contact-dependent interactions, in covariate %s", names(object@cont.results)[i]))
            }
        }
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
#' @param objList A list of CytoSignal objects or paths to CytoSignal objects
#' stored as RDS format. Can be mixed of two types.
#' @param metadata A data.frame object of dataset level metadata, containing
#' information like gender, age, treatment and etc. Default \code{NULL}. If
#' provided, must have as many rows as the length of \code{objList} with rows
#' in matching order to \code{objList}.
#' @param name.by A column name in \code{metadata} indicating where to find the
#' dataset names. Default \code{NULL} names datasets by "dataset1", "dataset2"
#' and etc, and will create column "dataset" in the metadata.
#' @param counts.thresh A number indicating the threshold for removing low-quality
#' spots. Default 0.
#' @param scale.factor A number indicating 1 spatial coord unit equals to how 
#' many µm. 
#' @return A \linkS4class{mergedCytoSignal} object preprocessed.
#' @export
#' @seealso [runNEBULA()]
mergeCytoSignal <- function(
        objList,
        metadata = NULL,
        name.by = NULL
        # counts.thresh = 0,
        # scale.factor = NULL
) {
    # Check metadata
    if (is.null(metadata)) {
        sid <- paste0("dataset", seq_along(objList))
        metadata <- data.frame(
            dataset = factor(sid, levels = sid),
        )
    } else {
        if (nrow(metadata) != length(objList)) {
            stop("`metadata` contains different number of rows than number of ",
                 "datasets in `objList`.")
        }
        if (is.null(name.by)) {
            sid <- paste0("dataset", seq_along(objList))
            if ("dataset" %in% colnames(metadata)) {
                stop('Column name "dataset" is reserved by downstream ',
                     'analysis, please rename/remove it, or set `name.by = ',
                     '"dataset"`.')
            }
        } else {
            if (!(name.by %in% colnames(metadata))) {
                stop(sprintf('Column name "%s" not found in metadata.', name.by))
            }
            if (name.by != "dataset" &&
                "dataset" %in% colnames(metadata)) {
                stop('Column name "dataset" is reserved by downstream ',
                     'analysis, please rename/remove it, or set `name.by = ',
                     '"dataset"`.')
            }
            sid <- as.character(metadata[[name.by]])
        }
        metadata$dataset <- factor(sid, levels = sid)

        reserved <- c("originalID", "clusters", "x", "y")
        for (p in reserved) {
            if (p %in% colnames(metadata)) {
                stop(sprintf('Column name "%s" is reserved by downstream analysis, please rename or remove it.', p))
            }
        }
    }
    
    # if (is.null(scale.factor)) {
    #     stop("Please provide `scale.factor`.")
    # }

    # Initialize lists that hold information needed
    diffList <- list()
    contList <- list()
    originalIDs <- list()
    clusters <- list()
    locations <- list()
    # Whether an object is read from RDS. If TRUE, need to remove from memory
    # after acquiring necessary information
    read <- FALSE

    for (i in seq_along(objList)) {
        obj <- objList[[i]]

        # Sanity check first
        if (inherits(obj, "CytoSignal")) {
            read <- FALSE
        } else if (is.character(obj)) {
            if (!file.exists(obj)) {
                stop(sprintf("File not found: %s", obj))
            }
            read <- TRUE
        } else {
            stop(sprintf("Invalid object type in objList at index %d", i))
        }

        if (read) {
            message("Loading object from file: ", obj)
            obj <- readRDS(obj)
            if (!inherits(obj, "CytoSignal")) {
                stop("Loaded object is not of CytoSignal class.")
            }
            # objList[[i]] <- obj
        }

        if (length(obj@clusters) == 0) {
            stop(sprintf("No clusters added in dataset %d (%s)", i, names(objList)[i]))
        }
        if (any(dim(obj@cells.loc) == 0)) {
            stop(sprintf("No cell locations added in dataset %d (%s)", i, names(objList)[i]))
        }


        # Preprocessing

        ids <- colnames(obj@raw.counts)
        originalIDs[[i]] <- ids
        newids <- paste0(sid[i], "_", ids)
        clusters[[i]] <- obj@clusters
        locations[[i]] <- obj@cells.loc

        # Preprocessing
        message(Sys.time(), " - Preprocessing dataset: ", i)
        
        # Hard coded parameters, pay attention to these
        #------------------- steps already run, skip for now ---------------#
        # obj <- addIntrDB(obj, g_to_u, db.diff, db.cont, inter.index)
        # obj <- removeLowQuality(obj, counts.thresh = counts.thresh)
        # obj <- changeUniprot(obj)
        # obj <- inferEpsParams(obj, scale.factor = scale.factor)
        # obj <- findNN(obj, diff.weight = 1)
        
        #----------------- Running Diff only, Hard coded params ------------#
        obj <- findNNGauEB(obj, self.weight = 1)
        # obj <- imputeLR(obj)
        obj <- imputeNiche(obj, nn.type = "GauEps", weights = "none")

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

        if (read) rm(obj)
        gc()
        message()
    }

    message(Sys.time(), " - Merging")
    # broadcast dataset-level metadata to location level
    sizes <- lengths(originalIDs)
    meta.full <- metadata[rep(seq_along(originalIDs), sizes), , drop = FALSE]
    meta.full$originalID <- unlist(originalIDs, use.names = FALSE)
    rownames(meta.full) <- paste0(meta.full$dataset, "_", meta.full$originalID)

    clusterUnion <- Reduce(union, lapply(clusters, levels))
    meta.full$clusters <- factor(NA, levels = clusterUnion)
    meta.full$x <- NA
    meta.full$y <- NA
    for (i in seq_along(objList)) {
        cellidx <- paste0(sid[i], '_', originalIDs[[i]])
        meta.full[cellidx, "clusters"] <- clusters[[i]]
        meta.full[cellidx, "x"] <- locations[[i]][, 1]
        meta.full[cellidx, "y"] <- locations[[i]][, 2]
    }

    # Merge diff lrscore
    diffIntrIsec <- Reduce(intersect, lapply(diffList, colnames))
    message(" - Found ", length(diffIntrIsec), " common interactions in diffusion-dependent LR scores.")
    diffList <- lapply(diffList, `[`, i = , j = diffIntrIsec, drop = FALSE)
    diffLR <- Reduce(rbind, diffList)

    # Merge cont lrscore
    contIntrIsec <- Reduce(intersect, lapply(contList, colnames))
    message(" - Found ", length(contIntrIsec), " common interactions in contact-dependent LR scores.")
    contList <- lapply(contList, `[`, i = , j = contIntrIsec, drop = FALSE)
    contLR <- Reduce(rbind, contList)

    methods::new(
        Class = "mergedCytoSignal",
        metadata = meta.full,
        diff.lrscore = t(diffLR),
        cont.lrscore = t(contLR)
    )
}

checkDependence <- function(X){
  # Perform QR decomposition
  qrX <- qr(X)
  
  # Get the rank
  matrix_rank <- qrX$rank
  
  # Get the pivot indices
  pivot_indices <- qrX$pivot
  
  # identify linear dependent columns
  # Get the R matrix
  R_matrix <- qr.R(qrX)
  
  # Extract the absolute values of the diagonal elements
  diag_R <- abs(diag(R_matrix))
  
  # Set tolerance level
  tolerance <- max(dim(X)) * .Machine$double.eps * max(diag_R)
  
  # Identify dependent diagonals
  dependent_diagonals <- diag_R < tolerance
  
  # Get the indices of dependent columns
  dependent_indices <- pivot_indices[which(dependent_diagonals)]
  
  # Get the names of dependent columns
  dependent_columns <- colnames(X)[dependent_indices]
  
  # # Output dependent columns
  # if (length(dependent_columns) > 0) {
  #   cat("Linearly dependent columns detected:\n")
  #   print(dependent_columns)
  #   return(dependent_columns)
  # } else {
  #   cat("No linearly dependent columns detected.\n")
  #   return()
  # }

  return(dependent_columns)
}


# Setup the model matrix and organizes the object internal information for
# downstream use.
# Returns list of (1) the modal matrix (2) the internal structure to save
.setup.model <- function(object, covariates) {
    if (!inherits(object, "mergedCytoSignal")) {
        stop("`object` must be a mergedCytoSignal object.")
    }

    if (!all(covariates %in% colnames(object@metadata))) {
        stop("Not all covariates specified are available in the metadata.")
    }


    model <- stats::model.matrix(
        stats::formula(paste0("~ ", paste0(c(-1, covariates), collapse = " + "))),
        data = object@metadata
    )
    m.names <- colnames(model)
    
    model <- cbind(1, model)
    colnames(model) <- c("(Intercept)", m.names)
    
    # check for linear dependence
    dep.cols <- checkDependence(model)
    
    if (length(dep.cols) > 0) {
      cat("Removing linearly dependent columns from model matrix.\n")
      print(dep.cols)
      remove.idx <- which(colnames(model) %in% dep.cols)
      model <- model[, -remove.idx]
    }
    
    # Organize cov.use for future use
    cov.use <- lapply(covariates, function(cov) {
        var <- object@metadata[[cov]]
        if (is.factor(var)) {
            struct <- list(
                type = "factor",
                levels = levels(var)
            )
        } else if (is.character(var)) {
            struct <- list(
                type = "character",
                levels = unique(var)
            )
        } else if (is.numeric(var)) {
            struct <- list(
                type = class(var)[1],
                range = range(var)
            )
        } else if (is.logical(var)) {
            struct <- list(
                type = "logical",
                levels = c("TRUE", "FALSE")
            )
        } else {
            struct <- list(
                type = class(var)[1]
            )
        }
    })
    cov.use <- lapply(seq_along(cov.use), function(i) {
        struct <- cov.use[[i]]
        type <- struct$type
        if (type != "numeric") {
            modelvar <- paste0(covariates[i], struct$levels)
            keep <- modelvar %in% colnames(model)
            struct$levels <- struct$levels[keep]
        }
        return(struct)
    })
    names(cov.use) <- covariates

    # Last check
    lapply(seq_along(cov.use), function(i) {
        struct <- cov.use[[i]]
        if (struct$type != "numeric") {
            paste0(covariates[i], struct$levels)
        } else {
            covariates[i]
        }
    }) |>
    unlist(use.names = FALSE) -> all_from_covuse
    modelcnCheck <- colnames(model)
    modelcnCheck <- modelcnCheck[modelcnCheck != "(Intercept)"]
    if (!all(all_from_covuse %in% modelcnCheck) ||
        !all(modelcnCheck %in% all_from_covuse)) {
        warning("Failed organizing internal record for covariate used")
    }

    return(list(model, cov.use))
}

.covlevelmap <- function(object) {
    covs <- character()
    suffix <- character()
    for (i in seq_along(object@cov.use)) {
        varname <- names(object@cov.use)[i]
        struct <- object@cov.use[[i]]
        if (struct$type != "numeric") {
            levels <- struct$levels
            covs <- c(covs, rep(varname, length(levels)))
            suffix <- c(suffix, levels)
        } else {
            covs <- c(covs, varname)
            suffix <- c(suffix, "")
        }
    }
    data.frame(covs = covs, suffix = suffix, stringsAsFactors = FALSE)
}

#' Run NEBULA's association analysis to compare datasets within mergedCytoSignal object
#' @description
#' After merging multiple CytoSignal objects and adding metadata covariates as
#' needed for modeling the comparison, this function runs NEBULA's association
#' analysis to identify significant interactions between cells from different
#' datasets and of different covariate levels.
#' @param object A \linkS4class{mergedCytoSignal} object
#' @param covariates A character vector of metadata covariates to include in the
#' model.
#' @param ncore Integer, number of cores to use for parallel computation.
#' Default \code{1L}.
#' @param verbose Logical, print progress messages. Default \code{TRUE}.
#' @export
#' @return Input \code{object} with slot \code{diff.results} and
#' \code{cont.results} updated with the results from NEBULA analysis.
#' @seealso [mergeCytoSignal()]
runNEBULA <- function(
        object,
        covariates,
        ncore = 1L,
        verbose = TRUE
) {
    if ("dataset" %in% covariates) {
        warning("'dataset' will be ignored from covariates.")
        covariates <- covariates[covariates != "dataset"]
    }
    modelsetup <- .setup.model(object, covariates)
    model <- modelsetup[[1]]
    object@cov.use <- modelsetup[[2]]

    if (isTRUE(verbose)) {
        message(Sys.time(), " - Running on Diffusion-dependent interactions...")
    }
    diff.res <- nebula::nebula(object@diff.lrscore, object$dataset, pred = model, ncore=ncore, verbose = verbose)

    if (isTRUE(verbose)) {
        message(Sys.time(), " - Running on Contact-dependent interactions...")
    }
    cont.res <- nebula::nebula(object@cont.lrscore, object$dataset, pred = model, ncore=ncore, verbose = verbose)

    diff.res.summary <- diff.res$summary
    cont.res.summary <- cont.res$summary

    # Since `model.matrix()` is stupid, we want to modify the column names it
    # produces in a robust way. A named character vector helps by directly
    # converting the model colnames to what I want
    covlevelmap <- .covlevelmap(object)
    modelcnmap <- stats::setNames(
        object = paste0(covlevelmap[[1]], '_', covlevelmap[[2]]),
        nm = paste0(covlevelmap[[1]], covlevelmap[[2]])
    )

    # split the summary tables by covariate levels, remove meta columms first,
    # then trim the prefixes and take uniq suffixes.
    summary.cn <- colnames(diff.res.summary)
    # Ignore the Intercept term
    summary.cn <- summary.cn[!grepl("_\\(Intercept\\)$", summary.cn)]
    # Ignore the "gene_id" and "gene" columns
    summary.cn <- summary.cn[!startsWith(summary.cn, "gene")]
    # Get the colnames trimmed removing the stats prefixes
    covnames <- unique(gsub("^(logFC|se|p)_", "", summary.cn))
    summary.cn <- lapply(covnames, function(cov) {
        paste0(c("logFC", "se", "p"), "_", cov)
    })

    diff.res.list <- lapply(seq_along(summary.cn), function(i) {
        sub <- diff.res.summary[, summary.cn[[i]]]
        colnames(sub) <- c("logFC", "se", "p")
        sub$padj <- stats::p.adjust(sub$p, method = "fdr")
        sub$interaction <- diff.res.summary$gene
        sub <- sub[, c(5, 1, 2, 3, 4)]
        rownames(sub) <- NULL
        return(sub)
    })
    cont.res.list <- lapply(seq_along(summary.cn), function(i) {
        sub <- cont.res.summary[, summary.cn[[i]]]
        colnames(sub) <- gsub(sprintf("_%s$", covnames[i]), "", colnames(sub))
        sub$padj <- stats::p.adjust(sub$p, method = "fdr")
        sub$interaction <- cont.res.summary$gene
        sub <- sub[, c(5, 1, 2, 3, 4)]
        rownames(sub) <- NULL
        return(sub)
    })

    covnames <- modelcnmap[covnames]
    names(diff.res.list) <- covnames
    names(cont.res.list) <- covnames

    object@diff.results <- diff.res.list
    object@cont.results <- cont.res.list

    methods::validObject(object)
    return(object)
}

#' Show available covariate levels used in NEBULA analysis
#' @param object A \linkS4class{mergedCytoSignal} object
#' @param intr.type A character string specifying the type of interaction to
#' show covariate levels for. Choose from \code{"diff"} for diffusion-dependent
#' or \code{"cont"} for contact-dependent interactions. Default \code{NULL}
#' returns the covariate levels used for both if they are the same, which is
#' true most of the time.
#' @return A character vector
#' @export
showCov <- function(object, intr.type = NULL) {
    diff_covs <- names(object@diff.results)
    cont_covs <- names(object@cont.results)
    if (is.null(intr.type)) {
        if (!identical(diff_covs, cont_covs)) {
            warning("Convariate levels found from diff-dep and cont-dep analyses differ. Showing covariate levels used for analysis of diff-dep intrs for now. Choose with `intr.type`.")
        }
        return(diff_covs)
    } else {
        switch(
            EXPR = intr.type,
            diff = diff_covs,
            cont = cont_covs,
            stop("Invalid `intr.type` argument. Choose from 'diff' or 'cont'.")
        )
    }
}

#' Make volcano plot for NEBULA analysis result
#' @description
#' This function allows drawing a volcano plot for the Nebula analysis on a
#' single covariate level, for either diffusion-dependent or contact-dependent
#' interactions. Either static ggplot figure or interactive plotly figure can
#' be generated. Interactions of interests can be labeled with: 1. choosing
#' specific interactions by name or index. 2. choosing top N interactions ranked
#' by significance (FDR and abs(logFC)). A non-NULL \code{highlight} argument is
#' considered before \code{topN}. When both are not specified, highlight all
#' significant interactions.
#' @param object A \linkS4class{mergedCytoSignal} object, with
#' \code{\link{runNEBULA}()} run in advance.
#' @param covariate A character string of the covariate to plot. Must be
#' available in \code{showCov(object)}.
#' @param intr.type A character string specifying the type of interaction to
#' plot. Choose from \code{"diff"} (default) for diffusion-dependent or
#' \code{"cont"} for contact-dependent interactions.
#' @param highlight A character vector of interaction names, a numeric or
#' logical vector of indices to select specific interactions to highlight; or
#' a single \code{TRUE} to highlight all significant interactions. See
#' \code{showIntr(object, intr.type)} for available options. Default \code{NULL}.
#' @param topN An integer to select top N interactions to highlight ranked by
#' significance. Rank by FDR first and abs(logFC) then. Default \code{NULL}.
#' @param fdrThresh A numeric value to set the FDR threshold for significance.
#' FDR value higher than this threshold will be considered as not significant.
#' Default \code{0.05}.
#' @param logfcThresh A numeric value to set the threshold on the absolute value
#' if logFC. Absoluate logFC value lower than this threshold will be considered
#' as not significant. Default \code{1}.
#' @param interactive Logical, whether to generate an interactive plotly figure.
#' Default \code{TRUE}. Generate a static ggplot figure if \code{FALSE}.
#' @param dotSize A numeric value to set the size of the dots in the plot.
#' Default \code{1.3}.
#' @param dotAlpha A numeric value to set the transparency of the dots in the
#' plot. Default \code{0.8}.
#' @param xTextSize A numeric value to set the size of the text on x-axis.
#' Default \code{10}.
#' @param xTitleSize A numeric value to set the size of the title on x-axis.
#' Default \code{12}.
#' @param yTextSize A numeric value to set the size of the text on y-axis.
#' Default \code{10}.
#' @param yTitleSize A numeric value to set the size of the title on y-axis.
#' Default \code{12}.
#' @param legendTextSize A numeric value to set the size of the text in the
#' legend. Default \code{12}.
#' @param legendTitleSize A numeric value to set the size of the title in the
#' legend. Default \code{12}.
#' @param labelTextSize A numeric value to set the size of the text on the
#' highlight labels. Default \code{4}.
#' @param titleSize A numeric value to set the size of the main title of the
#' plot. Default \code{16}.
#' @return
#' \itemize{
#'  \item A plotly figure (\code{plotly}, \code{htmlwidget} object) if \code{interactive = TRUE}.
#'  \item A ggplot figure (\code{gg}, \code{ggplot} object) if \code{interactive = FALSE}.
#' }
#' @export
plotNebulaVolcano <- function(
        object,
        covariate,
        intr.type = c("diff", "cont"),
        highlight = NULL,
        topN = NULL,
        fdrThresh = 0.05,
        logfcThresh = 1,
        interactive = TRUE,
        dotSize = 1.3,
        dotAlpha = 0.8,
        xTextSize = 10,
        xTitleSize = 12,
        yTextSize = 10,
        yTitleSize = 12,
        legendTextSize = 12,
        legendTitleSize = 12,
        labelTextSize = 4,
        titleSize = 16
) {
    intr.type <- match.arg(intr.type)
    res.list <- switch(
        EXPR = intr.type,
        diff = object@diff.results,
        cont = object@cont.results
    )
    if (length(covariate) != 1) stop("Only one covariate is allowed at a time.")
    if (!covariate %in% names(res.list)) {
        stop("Covariate given not found in the results. Available options are: ",
             paste0(paste0('"', names(res.list), '"'), collapse = ", "))
    }
    df <- res.list[[covariate]]
    if (any(df$padj == 0, na.rm = TRUE)) {
        warning("Adjusted p-value of 0 will be replaced with  0.1 * the minimum non-zero value.")
    }
    minNonZeroFDR <- min(df$padj[df$padj > 0], na.rm = TRUE)
    maxNeglog10FDR <- -log10(minNonZeroFDR) + 1
    df <- df |>
        dplyr::rename(FDR = .data[["padj"]]) |>
        dplyr::mutate(
            Significance = factor(dplyr::case_when(
                .data[['FDR']] < fdrThresh & abs(.data[['logFC']]) > logfcThresh ~ "FDR & LogFC",
                .data[['FDR']] < fdrThresh ~ "FDR",
                abs(.data[['logFC']]) > logfcThresh ~ "LogFC",
                .default = "Not significant"
            ), levels = c("Not significant", "FDR", "LogFC", "FDR & LogFC")),
            `-log10(FDR)` = dplyr::case_when(
                .data[['FDR']] == 0 ~ maxNeglog10FDR,
                .default = -log10(.data[['FDR']])
            )
        )
    if (!is.null(highlight)) {
        if (is.numeric(highlight)) {
            if (max(highlight) > nrow(df)) {
                stop("`highlight` can be numeric vector of numbers within 1 - ", nrow(df))
            }
            labelDF <- df[highlight, ,drop = FALSE]
        } else if (is.logical(highlight)) {
            if (length(highlight) == 1) {
                labelDF <- df |>
                    dplyr::filter(.data[["Significance"]] == "FDR & LogFC")
            } else if (length(highlight) == nrow(df)) {
                labelDF <- df[highlight, ,drop = FALSE]
            } else {
                stop("`highlight` can be a logical vector of length 1 or equal to the number of available interactions. See `showIntr(object, intr.type)` for availability.")
            }
        } else if (is.character(highlight)) {
            if (any(!highlight %in% df$interaction)) {
                stop("Some interactions in `highlight` are not available. See `showIntr(object, intr.type)` for available options.")
            }
            labelDF <- df |>
                dplyr::filter(.data[['interaction']] %in% highlight)
        } else {
            stop("`highlight` must be a numeric vector, a logical vector or a character vector.")
        }
    } else {
        if (is.null(topN)) topN <- sum(df$Significance == "FDR & LogFC")
        labelDF <- df |>
            dplyr::filter(.data[['Significance']] == "FDR & LogFC") |>
            dplyr::arrange(.data[["FDR"]], -dplyr::desc(.data[['logFC']])) |>
            dplyr::slice_head(n = topN)
    }

    colors <- c("grey20", "#4360DD", "#EA2116", "#2F872B")

    covlevelmap <- .covlevelmap(object)
    rownames(covlevelmap) <- paste0(covlevelmap$covs, "_", covlevelmap$suffix)
    varname <- covlevelmap[covariate, "covs"]
    level <- covlevelmap[covariate, "suffix"]
    titlelast <- if (nchar(level) > 0) paste0(varname, ": ", level) else varname
    titleText <- sprintf(
        "%s interaction, by %s",
        ifelse(intr.type == "diff", "Diffusion-dependent", "Contact-dependent"),
        titlelast
    )


    if (isFALSE(interactive)) {
        # Static ggplot figure
        p <- ggplot2::ggplot(
            data = df,
            mapping = ggplot2::aes(
                x = .data[['logFC']],
                y = .data[['-log10(FDR)']],
                color = .data[['Significance']]
            )
        ) +
            ggplot2::geom_point(alpha = 0.6, size = dotSize) +
            ggrepel::geom_text_repel(
                mapping = ggplot2::aes(label = .data[['interaction']],
                                       x = .data[['logFC']],
                                       y = .data[['-log10(FDR)']]),
                data = labelDF,
                inherit.aes = FALSE,
                force = 1, min.segment.length = 0, nudge_y = 1, # This line sets the small tick from point to label
                size = labelTextSize
            ) +
            ggplot2::guides(color = ggplot2::guide_legend(title = "Significance", override.aes = list(size = 4))) +
            ggplot2::geom_hline(yintercept = -log10(fdrThresh), linetype = "dashed") +
            ggplot2::geom_vline(xintercept = c(-logfcThresh, logfcThresh), linetype = "dashed") +
            ggplot2::scale_color_manual(values = c("grey20", "#4360DD", "#EA2116", "#2F872B"), drop = FALSE) +
            ggplot2::xlim(c(-max(abs(df$logFC)), max(abs(df$logFC)))) +
            ggplot2::theme_classic() +
            ggplot2::labs(title = titleText) +
            ggplot2::theme(
                legend.position = "bottom",
                plot.title = ggplot2::element_text(size = titleSize, hjust = 0.5),
                axis.text.x = ggplot2::element_text(size = xTextSize),
                axis.title.x = ggplot2::element_text(size = xTitleSize),
                axis.text.y = ggplot2::element_text(size = yTextSize),
                axis.title.y = ggplot2::element_text(size = yTitleSize),
                legend.text = ggplot2::element_text(size = legendTextSize),
                legend.title = ggplot2::element_text(size = legendTitleSize, face = "bold")
            )
    } else {
        # Interactive plotly figure
        p <- plotly::plot_ly(
            data = df,
            x = df[['logFC']],
            y = df[['-log10(FDR)']],
            type = "scatter",
            mode = "markers",
            color = df[['Significance']],
            colors = colors,
            text = paste0(
                "<b>", df[['interaction']], '</b><br>',
                '  FDR: ', signif(df[['FDR']], digits = 3), '<br>',
                'logFC: ', signif(df[['logFC']], digits = 3), '<br>'
            ),
            marker = list(
                size = dotSize*5,
                opacity = dotAlpha
            ),
            hovertemplate = "%{text}"
        ) |>
            plotly::layout(
                shapes = list(
                    list(type = "line",
                         x0 = -logfcThresh, y0 = min(df$`-log10(FDR)`),
                         x1 = -logfcThresh, y1 = max(df$`-log10(FDR)`),
                         line = list(color = "grey2", width = 1, dash = "dash")),
                    list(type = "line",
                         x0 = logfcThresh, y0 = min(df$`-log10(FDR)`),
                         x1 = logfcThresh, y1 = max(df$`-log10(FDR)`),
                         line = list(color = "grey2", width = 1, dash = "dash")),
                    list(type = "line",
                         x0 = -max(abs(df$logFC)), y0 = -log10(fdrThresh),
                         x1 = max(abs(df$logFC)), y1 = -log10(fdrThresh),
                         line = list(color = "grey2", width = 1, dash = "dash"))
                ),
                legend = list(
                    title = list(text = "<b>Significance</b>",
                                 font = list(size = legendTitleSize*1.2)),
                    font = list(size = legendTextSize*1.2),
                    orientation = "h",  # Horizontal orientation for the legend
                    x = 0.5,  # Center the legend horizontally
                    xanchor = "center",  # Anchor the legend at the center
                    y = -0.2  # Position the legend below the plot
                ),
                title = list(
                    text = titleText,
                    font = list(size = titleSize*1.5, weight = "bold")
                ),
                margin = list(t = 90, b = 90, l = 90, r = 90),
                annotations = lapply(seq_len(nrow(labelDF)), function(i) {
                    list(
                        x = labelDF$logFC[i],
                        y = labelDF$`-log10(FDR)`[i],
                        text = labelDF$interaction[i],
                        font = list(size = labelTextSize*4, color = "black"),
                        showarrow = TRUE,
                        arrowhead =0.5,
                        ax = 15,
                        ay = -30
                    )
                }),
                xaxis = list(
                    title = "LogFC",  # Add x-axis title
                    titlefont = list(size = xTitleSize*1.4)  # Customize font size of the title
                ),
                yaxis = list(
                    title = plotly::TeX("-log_{10}FDR"),  # Add y-axis title
                    titlefont = list(size = yTitleSize*1.4)  # Customize font size of the title
                )
            ) |>
            plotly::config(mathjax = "cdn")
    }
    return(p)
}


#' Visualize NEBULA analysis result with dot plot that involves all covariates
#' @description
#' This visualization method uses the result for either diffusion-dependent or
#' contact-dependent interactions from NEBULA analysis to generate a dot plot.
#' Each dot represents an interaction. The x axis is grouped by covariate levels
#' and jittered to avoid overlap. Therefore, each interaction appears in all the
#' covariate groups. The y axis is the logFC value of the interaction. The size
#' of the dots represents the significance (\eqn{-log_{10}FDR}) of an
#' interaction evaluated on a covariate level. Transparency of the dots shows
#' overall significance combining the FDR threshold and logFC threshold.
#' @param object A \linkS4class{mergedCytoSignal} object, with
#' \code{\link{runNEBULA}()} run in advance.
#' @param intr.type A character string specifying the type of interaction to
#' plot. Choose from \code{"diff"} (default) for diffusion-dependent or
#' \code{"cont"} for contact-dependent interactions.
#' @param fdrThresh A numeric value to set the FDR threshold for significance.
#' FDR value higher than this threshold will be considered as not significant.
#' Default \code{0.05}.
#' @param logfcThresh A numeric value to set the threshold on the absolute value
#' of logFC. Absoluate logFC value lower than this threshold will be considered
#' as not significant. Default \code{1}.
#' @param title A character string to set the title of the plot. Default
#' \code{NULL} does not add a title.
#' @param titleSize A numeric value to set the size of the title of the plot.
#' Default \code{16}.
#' @param xTextSize A numeric value to set the size of the text on x-axis.
#' Default \code{10}.
#' @param xTextAngle A numeric value to set the angle (in degree) to rotate the
#' x-axis text labels. Default \code{40}.
#' @param yTextSize A numeric value to set the size of the text on y-axis.
#' Default \code{10}.
#' @param yTitleSize A numeric value to set the size of the title on y-axis.
#' Default \code{12}.
#' @param legendTextSize A numeric value to set the size of the text in the
#' legend. Default \code{12}.
#' @param legendTitleSize A numeric value to set the size of the title in the
#' legend. Default \code{12}.
#' @param stripTitleSize A numeric value to set the size of the title in the
#' strip on top of each grid for the covariate groups. Default \code{12}.
#' @param ylim Numeric vector with 2 numbers to set the limits of y-axis.
#' Default \code{NULL} shows the full figure.
#' @return A ggplot figure (\code{gg}, \code{ggplot} object).
#' @export
plotNebulaAll <- function(
        object,
        intr.type = c("diff", "cont"),
        fdrThresh = 0.05,
        logfcThresh = 1,
        title = NULL,
        titleSize = 16,
        xTextSize = 10,
        xTextAngle = 40,
        yTextSize = 10,
        yTitleSize = 12,
        legendTextSize = 12,
        legendTitleSize = 12,
        stripTitleSize = 12,
        ylim = NULL
) {
    # object <- mergedcs
    intr.type <- match.arg(intr.type)
    all_list <- switch(
        EXPR = intr.type,
        diff = object@diff.results,
        cont = object@cont.results
    )

    all_list <- lapply(seq_along(all_list), function(i) {
        all_list[[i]][['covlevel']] <- names(all_list)[i]
        all_list[[i]]
    })
    full_df <- Reduce(rbind, all_list)

    covlevelmap <- .covlevelmap(object)
    rownames(covlevelmap) <- paste0(covlevelmap[[1]], '_', covlevelmap[[2]])
    covlevelmap$suffix[nchar(covlevelmap$suffix) == 0] <- covlevelmap$covs[nchar(covlevelmap$suffix) == 0]
    full_df$covariate_grouping <- covlevelmap[full_df$covlevel, "covs"]
    full_df$levels <- covlevelmap[full_df$covlevel, "suffix"]

    full_df$`-log10(FDR)` <- -log10(full_df$padj)
    nonzerononinfmax <- max(full_df$`-log10(FDR)`[is.finite(full_df$`-log10(FDR)`)], na.rm = TRUE)
    full_df$`-log10(FDR)`[is.infinite(full_df$`-log10(FDR)`)] <- nonzerononinfmax + 1
    # Scale the -log10FDR values to a range of 0.6 to 0.8 for alpha values
    # and set non-significant points to alpha 0.1
    old_min <- min(full_df$`-log10(FDR)`, na.rm = TRUE)
    old_max <- max(full_df$`-log10(FDR)`, na.rm = TRUE)
    new_min <- 0.6
    new_max <- 0.8
    alphas <- (full_df$`-log10(FDR)` - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
    alphas[full_df$padj > fdrThresh | abs(full_df$logFC) < logfcThresh] <- 0.2

    ylim <- if (is.null(ylim)) range(full_df$logFC, na.rm = TRUE) else ylim
    ggplot2::ggplot(
        full_df,
        ggplot2::aes(
            x = .data[['levels']],
            y = .data[['logFC']],
            size = .data[['-log10(FDR)']],
            color = .data[['levels']]
        )
    ) +
        ggplot2::geom_point(
            position = ggplot2::position_jitter(width = 0.4, height = 0),
            alpha = alphas,
            stroke = 0.2
        ) +
        # Do not show color legend
        ggplot2::guides(color = "none") +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(yintercept = logfcThresh, linetype = "dashed") +
        ggplot2::geom_hline(yintercept = -logfcThresh, linetype = "dashed") +
        ggplot2::facet_grid(
            stats::formula("~covariate_grouping"),
            scales = 'free_x',
            space = 'free_x'
        ) +
        ggplot2::ylim(ylim) +
        ggplot2::labs(title = title) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = titleSize, hjust = 0.5),
            axis.text.x = ggplot2::element_text(
                angle = xTextAngle,
                hjust = ifelse(xTextAngle > 0, 1, 0),
                size = xTextSize
            ),
            axis.title.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = yTextSize),
            axis.title.y = ggplot2::element_text(size = yTitleSize),
            axis.line = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(fill = NA, color = 'black'),
            panel.grid = ggplot2::element_blank(),
            strip.text.x = ggplot2::element_text(size = stripTitleSize),
            legend.title = ggplot2::element_text(size = legendTitleSize),
            legend.text = ggplot2::element_text(size = legendTextSize)
        )
}

methods::setOldClass('package_version')
methods::setClassUnion(
    name = 'sparse_or_null',
    members = c('dgCMatrix', 'NULL')
)
methods::setClassUnion(
    name = 'array_or_null',
    members = c('array', 'NULL')
)


#' @name CytoSignal2-class
#' @title CytoSignal2-class
#' @description
#' Container object for CytoSignal analysis that holds input data, parameters,
#' and results.
#' @slot rawData A sparse matrix of gene expression data with genes as rows and spots as columns.
#' @slot spatial A matrix of spatial coordinates with rows as spots and two columns (x and y).
#' @slot metadata A tibble (data.frame) containing metadata for each spot.
#' @slot intrDB An interaction database object of class `csdb`. See \link{CytoSignal-database} for how to construct.
#' @slot parameters A list of parameters used in the analysis.
#' @slot neighborDiff A sparse matrix representing the diffusion-dependent neighborhood model.
#' @slot neighborCont A sparse matrix representing the contact-dependent neighborhood model.
#' @slot LRScore A sparse matrix of ligand-receptor interaction scores.
#' @slot nullECDF A sparse matrix of empirical cumulative distribution functions
#' for interaction scores. A permutation p-value can be derived by subtracting
#' these values from 1.
#' @slot significance A list containing significance results for the interaction scores.
#' @slot version A package_version object indicating the version of CytoSignal
#' used to create the object.
#' @export
#' @rdname CytoSignal2-class
methods::setClass(
    Class = 'cytosignal2',
    representation = methods::representation(
        rawData = 'dgCMatrix',
        spatial = 'matrix',
        metadata = 'tbl_df',
        intrDB = 'ANY',
        parameters = 'list',
        neighborDiff = 'sparse_or_null',
        neighborCont = 'sparse_or_null',
        LRScore = 'sparse_or_null',
        clusterLRScore = 'array_or_null',
        significance = 'list',
        version = 'package_version'
    ),
    prototype = methods::prototype(
        intrDB = NULL,
        neighborDiff = NULL,
        neighborCont = NULL,
        LRScore = NULL,
        clusterLRScore = NULL,
        significance = list(),
        version = packageVersion('cytosignal')
    )
)

.valid.cytosignal2 <- function(object) {
    bc <- colnames(object@rawData)
    if (!identical(rownames(object@spatial), bc)) {
        return('Row names of spatial must be identical to column names of rawData.')
    }
    if (ncol(object@spatial) != 2) {
        return('Spatial must have exactly 2 columns representing x and y coordinates.')
    }
    if (!'barcodes' %in% colnames(object@metadata)) {
        return('Metadata must contain a column named "barcodes".')
    }
    if (!identical(object@metadata$barcodes, bc)) {
        return('Metadata barcodes must be identical to column names of rawData.')
    }
    if (!is.null(object@neighborDiff)) {
        if (ncol(object@neighborDiff) != length(bc) ||
            nrow(object@neighborDiff) != length(bc)) {
            return('Diffusion-dependent neighbor graph must be a square matrix with number of rows/columns equal to number of spots in rawData.')
        }
    }
    if (!is.null(object@neighborCont)) {
        if (ncol(object@neighborCont) != length(bc) ||
            nrow(object@neighborCont) != length(bc)) {
            return('Contact-dependent neighbor graph must be a square matrix with number of rows/columns equal to number of spots in rawData.')
        }
    }
    if (!is.null(object@intrDB)) {
        dbgenes <- .uniqGeneInDB(object@intrDB)
        if (length(intersect(rownames(object@rawData), dbgenes)) == 0) {
            return('No genes in rawData are found in the interaction database intrDB.')
        }
    }

    if (!is.null(object@LRScore)) {
        if (!identical(rownames(object@LRScore), colnames(object@rawData)) ||
            !identical(colnames(object@LRScore), object@intrDB$interactors)) {
            return('LRScore matrix do not match to spots in rawData and/or interactions in intrDB.')
        }
    }

    # TODO: whether or not to check clusterLRScore dimensions?
    # User might be switching the default cluster variable after running one
    # so things won't match naturally.
    return(TRUE)
}

methods::setValidity(
    Class = 'cytosignal2',
    method = .valid.cytosignal2
)

#' @rdname CytoSignal2-class
#' @param rawData A sparse matrix of gene expression data with genes as rows and
#' spots as columns.
#' @param spatial The spatial coordinates coercible to a matrix, with rows as
#' spots and two columns for x and y.
#' @param cluster A vector/factor of cluster assignments for each spot. The
#' length must equal to the number of spots (columns) in \code{rawData}.
#' @param ... Additional arguments directly passed to object slots.
#' @export
#' @examples
#' spatial <- matrix(runif(100), ncol = 2)
#' rownames(spatial) <- paste0('Spot', 1:50)
#' cs <- createCytoSignal2(
#'    rawData = Matrix::rsparsematrix(
#'        nrow = 100, ncol = 50, density = 0.1,
#'        dimnames = list(
#'            paste0('Gene', 1:100),
#'            paste0('Spot', 1:50)
#'        )
#'    ),
#'    spatial = spatial
#' )
#' cs$total_counts
#' cs$cluster <- sample(letters[1:3], 50, replace = TRUE)
createCytoSignal2 <- function(
        rawData,
        spatial,
        cluster = NULL,
        ...
) {
    if (!inherits(rawData, 'dgCMatrix')) {
        rawData <- methods::as(rawData, 'CsparseMatrix')
    }
    barcodes <- colnames(rawData)
    genes <- rownames(rawData)

    spatial <- as.matrix(spatial)
    if (nrow(spatial) != length(barcodes)) {
        cli::cli_abort(
            'Number of rows of {.code spatial} must equal to number of spots ({length(barcodes)}), while {nrow(spatial)} was provided'
        )
    }
    if (ncol(spatial) != 2) {
        cli::cli_abort(
            '{.code spatial} must have exactly 2 columns representing x and y coordinates, while {ncol(spatial)} were provided.'
        )
    } else {
        colnames(spatial) <- c('x', 'y')
    }
    if (!is.null(rownames(spatial))) {
        if (!identical(rownames(spatial), barcodes)) {
            cli::cli_abort(
                'Row names of {.code spatial} must be identical to column names of {.code rawData}.'
            )
        }
    } else {
        cli::cli_warn(
            'Row names of {.code spatial} do not exist. Setting them to column names of {.code rawData}.'
        )
        rownames(spatial) <- barcodes
    }

    metadata <- tibble::tibble(
        barcodes = barcodes,
        total_counts = Matrix::colSums(rawData),
        gene_detected = Matrix::colSums(rawData > 0)
    )

    object <- methods::new(
        Class = 'cytosignal2',
        rawData = rawData,
        spatial = spatial,
        metadata = metadata,
        ...
    )
    if (!is.null(cluster)) {
        object <- setCluster(object, cluster)
    }
    return(object)
}

#' @rdname CytoSignal2-class
#' @export
#' @param object,x The `cytosignal2` object.
methods::setMethod(
    f = 'show',
    signature = 'cytosignal2',
    definition = function(object) {
        cat('An object of class "cytosignal2"\n')
        cat(cli::format_inline(
            '- Raw data: {nrow(object@rawData)} genes x {ncol(object@rawData)} spots\n'
        ))
        cat(cli::format_inline(
            '- Metadata columns ({ncol(object@metadata)}): {.val {colnames(object@metadata)}}\n'
        ))
        cat('- Interaction database: ')
        if (is.null(object@intrDB)) {
            cat('(do `setIntrDB()`)\n')
        } else {
            cat(cli::format_inline('{nrow(object@intrDB)} interaction{?s}, {length(.uniqGeneInDB(object@intrDB))} gene{?s}.\n'))
        }

        cat('- Parameters: ')
        if (any(!c('micronPerUnit', 'ballRadiusMicron') %in% names(object@parameters))) {
            cat('(do `setParams()`)\n')
        } else {
            cat('\n  - Microns per coordinate unit: ')
            if (is.null(object@parameters$micronPerUnit)) {
                cat('\n')
            } else {
                cat(object@parameters$micronPerUnit, '\n')
            }
            cat('  - Epsilon ball radius in micron: ')
            if (is.null(object@parameters$ballRadiusMicron)) {
                cat('\n')
            } else {
                cat(object@parameters$ballRadiusMicron, '\n')
            }
        }


        cat('- Neighborhood detection\n')
        cat('  - Diffusion-dependent model:')
        if (is.null(object@neighborDiff)) {
            cat('\n')
        } else {
            cat('Done\n')
        }
        cat('  - Contact-dependent model:')
        if (is.null(object@neighborCont)) {
            cat('\n')
        } else {
            cat('Done\n')
        }
    }
)

#' @rdname CytoSignal2-class
#' @export
#' @param pattern A character string containing a regular expression to match column names of metadata.
.DollarNames.cytosignal2 <- function(x, pattern = "") {
    grep(pattern, colnames(x@metadata), value = TRUE)
}


#' @rdname CytoSignal2-class
#' @export
#' @param name The name of the metadata column to access.
`$.cytosignal2` <- function(x, name) {
    if (!name %in% colnames(x@metadata)) {
        NULL
    } else {
        value <- x@metadata[[name]]
        names(value) <- colnames(x@rawData)
        value
    }
}

#' @rdname CytoSignal2-class
#' @export
#' @param value The value to assign.
`$<-.cytosignal2` <- function(x, name, value) {
    if (!is.null(names(value))) {
        if (!identical(names(value), colnames(x@rawData))) {
            cli::cli_warn(
                "Assigning named variable while the names do not match to spot barcodes."
            )
        }
    }
    x@metadata[[name]] <- value
    x
}

#' Set default cluster variable
#' @description
#' Set an existing metadata variable as the default cluster annotation, or add
#' a new vector/factor as cluster annotation. This variable will be used by
#' default in cluster-wise inference and visualization functions.
#' @param object A \linkS4class{cytosignal2} object.
#' @param cluster A single string for selecting a metadata column, or a vector
#' or factor for cluster assignment.
#' @return The input \linkS4class{cytosignal2} object with updated parameters
#' at \code{object@parameters$cluster}.
#' @export
setCluster <- function(object, cluster) {
    if (length(cluster) == 1) {
        cluster <- rlang::arg_match(
            arg = cluster,
            values = colnames(object@metadata),
            multiple = FALSE
        )
        object@parameters[['cluster']] <- cluster
    } else {
        if (length(cluster) != nrow(object@metadata)) {
            cli::cli_abort(
                "Vector/factor cluster variable must have {nrow(object@metadata)} elements. Given has {length(cluster)}."
            )
        }
        if (!is.null(names(cluster))) {
            if (!identical(names(cluster), colnames(object@rawData))) {
                cli::cli_warn(
                    "Assigning named cluster variable while the names do not match to spot barcodes."
                )
            }
        }
        object@metadata[['cluster']] <- cluster
        object@parameters[['cluster']] <- 'cluster'
    }
    return(object)
}


#' Subsetting cytosignal2 object by spots
#' @description
#' Subset a \linkS4class{cytosignal2} object by spots.
#' @param x A \linkS4class{cytosignal2} object.
#' @param subset Logical, numeric, or character vector specifying the spots to
#' keep. Usually passed with expression using metadata variables.
#' @param drop Logical, whether to drop unused factor levels in metadata after
#' subsetting. Default \code{FALSE}.
#' @param keepNeighbor Logical, whether to keep the subset neighbor graphs
#' regardless of whether the subset spots still have the same neighbor structure.
#' Default \code{FALSE} remove graphs from the subset object.
#' @param ... Additional arguments (not used).
#' @return A \linkS4class{cytosignal2} object subsetted by spots.
#' @export
#' @method subset cytosignal2
subset.cytosignal2 <- function(
        x,
        subset,
        drop = FALSE,
        keepNeighbor = FALSE,
        ...
) {
    rawData <- x@rawData[, subset, drop = FALSE]
    spatial <- x@spatial[subset, , drop = FALSE]
    metadata <- x@metadata[subset, , drop = drop]
    x@rawData <- rawData
    x@spatial <- spatial
    x@metadata <- metadata
    if (isTRUE(keepNeighbor)) {
        if (!is.null(x@neighborDiff)) {
            x@neighborDiff <- x@neighborDiff[subset, subset, drop = FALSE]
        }
        if (!is.null(x@neighborCont)) {
            x@neighborCont <- x@neighborCont[subset, subset, drop = FALSE]
        }
    } else {
        cli::cli_warn('Neighbor graphs are removed after subsetting. Please rerun {.fn findNeighbor} or set {.code keepNeighbor = TRUE}.')
        x@neighborDiff <- NULL
        x@neighborCont <- NULL
    }

    return(x)
}


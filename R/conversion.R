#' Conversion between CytoSignal object and Seurat object
#' @description
#' Converts a CytoSignal object to a Seurat object or from a Seurat object to a
#' CytoSignal object.
#' @rdname csToSeurat
#' @export
#' @param object A CytoSignal object or a Seurat object.
#' @param assay Name of a Seurat assay to extract count matrix from or put into.
#' Default \code{"Spatial"}
#' @param layer Name of Seurat V5 layer name to extract the raw count matrix
#' from the assay. Default \code{"counts"}.
#' @param image Name of the spatial image object in the Seurat object where the
#' spatial coordinates are stored. Default \code{NULL} uses the first image
#' available to the specified assay.
#' @param cluster Name of the metadata variable in the Seurat object that
#' contains the cluster identities of each spot. Default \code{NULL} uses the
#' activate identity in the Seurat object.
#' @return
#' \itemize{
#' \item {
#' \code{SeuratToCS} returns a CytoSignal object that includes only the
#' raw counts, coordinates, and cluster identities of the spatial assay.
#' }
#' \item {
#' \code{csToSeurat} returns a Seurat object with "Spatial" assay that contains
#' the raw counts. When LR-score and signifcance inferrence results are present,
#' another assay called "CytoSignal" is created to store the LR-score matrices
#' and sparse binary matrices indicating whether an interaction is significant
#' at a spot. The spatial coordinates are coerced to a "SlideSeq" object
#' regardless of the exact technology being imported from. A interaction
#' metadata data frame is stored at the feature metadata slot of the assay. It
#' contains interpretable interaction names, and the ranking of the interaction
#' at different significance inferrence levels (NA indicating not significant).
#' }
#' }
csToSeurat <- function(object) {
    if (!inherits(object, "CytoSignal")) stop("Input object must be a CytoSignal object.")
    .checkDep(c("Seurat", "SeuratObject"), c("5.0.0", "5.0.0"))
    loc <- as.data.frame(object@cells.loc)
    cluster <- object@clusters

    raw.counts <- object@raw.counts
    seu <- SeuratObject::CreateSeuratObject(
        raw.counts,
        assay = "Spatial",
        project = "CytoSignal"
    )
    seu[['cluster']] <- cluster
    SeuratObject::Idents(seu) <- cluster

    coordAsasy <- "Spatial"
    lrscores <- names(object@lrscore)
    lrscores <- lrscores[lrscores != "default"]
    if (length(lrscores) > 0) {
        scoreMats <- lapply(lrscores, function(scoreName) {
            t(object@lrscore[[scoreName]]@score)
        })
        names(scoreMats) <- lrscores

        signifMats <- lapply(lrscores, function(scoreName) {
            signifLists <- object@lrscore[[scoreName]]@res.list
            resUse <- names(signifLists)
            resUse <- resUse[startsWith(resUse, "result.")]
            signifLists <- signifLists[resUse]
            names(signifLists) <- paste0("isSignif.", scoreName, ".", gsub("result.", "", resUse))
            signifLists
        })
        signifMats <- do.call(c, signifMats)
        signifMats <- lapply(signifMats, makeSignifMat, spotNames = rownames(loc))
        csAssay <- SeuratObject::CreateAssay5Object(data = c(scoreMats, signifMats))

        # Managing the interaction metadata
        intrNames <- getIntrNames(object, rownames(csAssay))
        csAssay[["name"]] <- intrNames
        for (scoreName in lrscores) {
            signifLists <- object@lrscore[[scoreName]]@res.list
            resUse <- names(signifLists)
            resUse <- resUse[startsWith(resUse, "result.")]
            for (resName in resUse) {
                metaName <- paste0("rank_", scoreName, "_", resName)
                resList <- signifLists[[resName]]
                resCol <- setNames(rep(NA, length(intrNames)), intrNames)
                resCol[names(resList)] <- seq_along(resList)
                csAssay[[metaName]] <- resCol
            }
        }

        # Insert the assay to the seurat object
        Seurat::Key(csAssay) <- "cytosignal_"
        seu@assays[["CytoSignal"]] <- csAssay
        coordAsasy <- "CytoSignal"
        SeuratObject::DefaultAssay(seu) <- "CytoSignal"
    }

    spatialImageObj <- new(
        Class = 'SlideSeq',
        assay = coordAsasy,
        coordinates = loc,
        key = 'slice_'
    )
    seu@images[["csimage"]] <- spatialImageObj

    seu
}

#' @rdname csToSeurat
#' @export
SeuratToCS <- function(
        object,
        layer = "counts",
        assay = "Spatial",
        image = NULL,
        cluster = NULL
) {
    if (!inherits(object, "Seurat")) stop("Input object must be a Seurat object.")
    .checkDep(c("Seurat", "SeuratObject"), c("5.0.0", "5.0.0"))
    imageAvail <- SeuratObject::Images(object, assay = assay)
    if (length(imageAvail) == 0) {
        stop("No image data found in Seurat object for assay '", assay, "'.")
    }
    if (is.null(image)) {
        if (length(imageAvail) > 1) {
            warning("Multiple images found in Seurat object for assay '", assay,
                    "'. Using the first one.")
        }
        image <- imageAvail[1]
    } else {
        if (!(image %in% imageAvail)) {
            stop("Image '", image, "' not found in Seurat object for assay '",
                 assay, "'.")
        }
    }
    imageObj <- object@images[[image]]
    loc <-SeuratObject::GetTissueCoordinates(imageObj)[, c('x', 'y')]
    rownames(loc) <- SeuratObject::Cells(imageObj)
    loc <- as.matrix(loc)

    counts <- SeuratObject::LayerData(object, layer = layer, assay = assay)
    if (is.null(cluster)) cluster <- SeuratObject::Idents(object)
    else cluster <- object[[cluster, drop = TRUE]]
    spotUse <- intersect(colnames(counts), rownames(loc))
    if (!all(colnames(counts) %in% spotUse)) {
        warning("Not all spots have location information in the image. ",
                "Dropping spots without location information.")
    }
    if (!all(rownames(loc) %in% spotUse)) {
        warning("Not all spots have count matrix information. ",
                "Dropping spots without count matrix information.")
    }
    loc <- loc[spotUse, , drop = TRUE]
    counts <- counts[, spotUse, drop = TRUE]
    cluster <- cluster[spotUse]
    createCytoSignal(
        raw.data = counts,
        cells.loc = loc,
        clusters = cluster,
        name = assay
    )
}

#' Conversion between CytoSignal object and AnnData object in H5AD file
#' @description
#' Extract data from H5AD file to a CytoSignal object or write a CytoSignal
#' object to a H5AD file.
#' @rdname csToH5AD
#' @export
#' @param object A CytoSignal object.
#' @param filename Name of the H5AD file.
#' @param overwrite Whether to overwrite the existing H5AD file. Default
#' \code{FALSE}.
#' @return
#' \itemize{
#' \item{
#' \code{H5ADToCS} returns a CytoSignal object, within only the raw counts
#' expression matrix, spatial coordinates, and optional cluster identities
#' extracted.
#' }
#' \item{
#' \code{csToH5AD} writes an H5AD file to disk and returns \code{NULL}. with the
#' following information written to corresponding location:
#' \itemize{
#'   \item{Raw counts - \code{adata.X}}
#'   \item{Clusters - \code{adata.obs['cluster']}}
#'   \item{Spatial coordinates - \code{adata.obsm['spatial']}}
#'   \item{LR-score for analysis type, e.g. 'GauEps-Raw' - \code{adata.obsm['GauEps-Raw_score']}}
#'   \item{Significant spots for analysis type, e.g. 'GauEps-Raw', inferred with
#'   level, e.g. 'result.spx' - \code{adata.obsm['GauEps-Raw_signif_result.spx']}.
#'   These are stored in form of sparse matrices, where value 1 means the
#'   interaction is tested to be significant at the spot.}
#'   \item{Interaction metadata for analysis type, e.g. 'GauEps-Raw' - \code{adata.uns['CytoSignal']['GauEps-Raw']}}
#'   This is stored as a DataFrame containing
#'   \itemize{
#'     \item{id - Interaction ID}
#'     \item{name - Interpretable name of the interaction}
#'     \item{rank_result.spx - Rank of the interaction in the significance
#'     inferrence level 'result.spx'}
#'   }
#' }
#' }
#' }
csToH5AD <- function(object, filename, overwrite = FALSE) {
    if (!inherits(object, "CytoSignal")) stop("Input object must be a CytoSignal object.")
    .checkDep("hdf5r")
    if (file.exists(filename)) {
        if (isTRUE(overwrite)) {
            file.remove(filename)
        } else {
            stop("H5AD file exists at: ", normalizePath(filename),
                 "\nPlease set `overwrite = TRUE` to overwrite the file.")
        }
    }

    finished <- FALSE
    adata <- hdf5r::H5File$new(filename, mode = "w")
    on.exit({
        adata$close_all()
        if (!finished) warning("File closed before writing completed. File may be corrupted.")
    })
    message("Writing .X")
    Xmat <- object@raw.counts
    .writeMatrixToH5AD(x = Xmat, dfile = adata, dname = "X")

    message("Writing .var")
    varDF <- data.frame(name = rownames(Xmat), row.names = rownames(Xmat))
    .writeDataFrameToH5AD(df = varDF, dfile = adata, dname = "var")

    message("Writing .obs")
    obsDF <- data.frame(cluster = object@clusters, row.names = colnames(Xmat))
    .writeDataFrameToH5AD(df = obsDF, dfile = adata, dname = "obs")

    message("Writing .obsm['spatial']")
    .createMapping(adata, 'obsm')
    adata[['obsm']]$create_dataset(
        name = 'spatial',
        robj = t(object@cells.loc),
        dtype = .H5AD.guessDType(object@cells.loc),
        chunk_dims = NULL
    )

    .createMapping(adata, 'uns')
    .createMapping(adata, 'uns/CytoSignal')

    scoreAvail <- names(object@lrscore)
    scoreAvail <- scoreAvail[scoreAvail != "default"]

    for (scoreName in scoreAvail) {
        obsmScoreKey <- paste0(scoreName, "_score")
        message("Writing LR-score from ", scoreName, " to .obsm['", obsmScoreKey, "']")
        .writeMatrixToH5AD(
            x = t(object@lrscore[[scoreName]]@score),
            dfile = adata,
            dname = paste0('obsm/', obsmScoreKey)
        )
        allIntrs <- colnames(object@lrscore[[scoreName]]@score)

        intrDF <- data.frame(
            id = allIntrs,
            name = getIntrNames(object, allIntrs),
            row.names = allIntrs
        )

        signifLists <- object@lrscore[[scoreName]]@res.list
        resUse <- names(signifLists)
        resUse <- resUse[startsWith(resUse, "result")]
        signifLists <- signifLists[resUse]
        names(signifLists) <- sprintf("%s_signif_%s", scoreName, resUse)
        signifMats <- lapply(
            X = signifLists, FUN = makeSignifMat,
            spotNames = colnames(Xmat),
            intrNames = allIntrs
        )
        for (i in seq_along(resUse)) {
            metaName <- paste0("rank_", resUse[i])
            intrDF[[metaName]] <- -1L
            resList <- signifLists[[i]]
            intrDF[names(resList), metaName] <- seq_along(resList)
        }
        for (i in seq_along(signifMats)) {
            obsmSignifKey <- names(signifMats)[i]
            message("Writing significant spot indication to .obsm['", obsmSignifKey, "']")
            .writeMatrixToH5AD(
                x = signifMats[[i]],
                dfile = adata,
                dname = paste0('obsm/', obsmSignifKey)
            )
        }

        message("Writing interaction metadata to .uns['CytoSignal']['", scoreName, "']")
        .writeDataFrameToH5AD(df = intrDF, dfile = adata, dname = paste0('uns/CytoSignal/', scoreName))
    }

    finished <- TRUE
    return(invisible(NULL))
}

#' @rdname csToH5AD
#' @export
#' @param clusterKey Variable name in \code{adata.obs} that provides categorical
#' clustering.
#' @param layer Data from the H5AD file to extract. Default \code{"X"} extracts
#' \code{adata.X}. Use \code{"raw/X"} to extract \code{adata.raw.X}. Use
#' \code{"layers/key"} to extract layer data at \code{adata.layers['key']}.
#' @param spatialKey Key to extract spatial coordinate matrix from
#' \code{adata.obsm}. Default \code{"spatial"}.
H5ADToCS <- function(filename, clusterKey, layer = "X", spatialKey = "spatial") {
    .checkDep("hdf5r")
    if (!file.exists(filename)) stop("File '", filename, "' not found.")
    adata <- hdf5r::H5File$new(filename, mode = "r")
    # Checking if specified fields exist
    clusterKey <- match.arg(clusterKey, choices = names(adata[['obs']]))
    spatialKey <- match.arg(spatialKey, choices = names(adata[['obsm']]))
    matOptions <- "X"
    if (adata$exists('raw/X')) matOptions <- c(matOptions, "raw/X")
    if (adata$exists('layers'))
        matOptions <- c(matOptions, paste0('layers/', names(adata[['layers']])))
    layer <- match.arg(layer, choices = matOptions)

    # Load expression matrix
    matObj <- adata[[layer]]
    if (!inherits(matObj, "H5Group") ||
        !'data' %in% names(matObj) ||
        !'indices' %in% names(matObj) ||
        !'indptr' %in% names(matObj)) {
        warning("Selected layer '", layer, "' does not contain a sparse matrix.")
    }
    x <- matObj[['data']][]
    if (!rlang::is_integerish(x)) {
        otherAvail <- matOptions[matOptions != layer]
        otherAvailStr <- paste(paste0("'", otherAvail, "'"), collapse = ", ")
        warning(
            "CytoSignal recommends using raw counts as input, where as the ",
            "expression values in layer '", layer, "' are not integerish. ",
            "Other available options include: ", otherAvailStr
        )
    }
    i <- matObj[['indices']][]
    p <- matObj[['indptr']][]
    if ('shape' %in% names(matObj)) {
        shape <- rev(matObj[['shape']][])
    } else if ('shape' %in% hdf5r::h5attr_names(matObj)) {
        shape <- rev(hdf5r::h5attr(matObj, 'shape'))
    } else {
        stop("Unsupported matrix encoding")
    }
    mat <- Matrix::sparseMatrix(
        i = i + 1, p = p, x = x,
        dims = shape, repr = "C"
    )
    obsNamesKey <- hdf5r::h5attr(adata[['obs']], '_index')
    obsNames <- adata[[paste0('obs/', obsNamesKey)]][]
    varUse <- ifelse(layer == 'raw/X', 'raw/var', 'var')
    varNamesKey <- hdf5r::h5attr(adata[[varUse]], '_index')
    varNames <- adata[[paste0(varUse, '/', varNamesKey)]][]
    dimnames(mat) <- list(varNames, obsNames)

    # Load cluster key
    clusterObj <- adata[[paste0('obs/', clusterKey)]]
    if (inherits(clusterObj, "H5Group") &&
        'categories' %in% names(clusterObj) &&
        'codes' %in% names(clusterObj)) {
        # Encoding-version 0.2.0, 'categories' and 'code' in a group for the categorical variable
        cluster <- clusterObj[['categories']][clusterObj[['codes']][] + 1]
        cluster <- factor(cluster, levels = clusterObj[['categories']][])
    } else {
        encodingVer <- hdf5r::h5attr(clusterObj, 'encoding-version')
        stop("adata.obs categorical encoding (version ", encodingVer,
             ") is not supported. Please submit an issue.")
    }
    names(cluster) <- obsNames

    # Load location
    loc <- adata[[paste0('obsm/', spatialKey)]][,]
    loc <- t(loc)
    dimnames(loc) <- list(obsNames, c('x', 'y'))

    cs <- createCytoSignal(raw.data = mat, cells.loc = loc, clusters = cluster)
    return(cs)
}

.checkDep <- function(pkgs, version = NULL) {
    for (i in seq_along(pkgs)) {
        pkg <- pkgs[i]
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Package ", pkg, " is required but not installed.")
        }
        if (!is.null(version[i])) {
            hasVer <- utils::packageVersion(pkg)
            requiredVer <- package_version(version[i])
            if (hasVer < requiredVer) {
                stop("Package ", pkg, " version ", version[i], " or higher is required.")
            }
        }
    }
}

makeSignifMat <- function(res, spotNames, intrNames = NULL) {
        # Here, res is a list, where each element is a character vector of
        # spot IDs that are significant for a given interaction.
        allIntr <- names(res)
        # intrNames <- intrNames %||% allIntr
        intrNames <- if (!is.null(intrNames)) intrNames else allIntr
        sparseIJ <- list()
        for (intr in allIntr) {
            IJ <- data.frame(
                i = rep(match(intr, intrNames), length(res[[intr]])),
                j = match(res[[intr]], spotNames)
            )
            sparseIJ[[intr]] <- IJ
        }
        sparseIJ <- do.call(rbind, sparseIJ)
        Matrix::sparseMatrix(
            i = sparseIJ$i, j = sparseIJ$j, x = 1,
            dims = c(length(intrNames), length(spotNames)),
            dimnames = list(intrNames, spotNames), repr = "C"
        )
}

.writeMatrixToH5AD <- function(
        x,
        dfile,
        dname
) {
    dnamePaths <- unlist(strsplit(dname, '/'))
    # Recursively create groups
    for (i in seq_along(dnamePaths)) {
        search <- paste0(dnamePaths[1:i], collapse = '/')
        if (!dfile$exists(name = search)) {
            dfile$create_group(name = search)
        }
    }
    dfile[[dname]]$create_dataset(
        name = 'data',
        robj = x@x,
        dtype = .H5AD.guessDType(x@x)
    )
    dfile[[dname]]$create_dataset(
        name = 'indices',
        robj = x@i,
        dtype = .H5AD.guessDType(x@i)
    )
    dfile[[dname]]$create_dataset(
        name = 'indptr',
        robj = x@p,
        dtype = .H5AD.guessDType(x@p)
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-type',
        robj = 'csr_matrix',
        dtype = .H5AD.guessDType('csr_matrix'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-version',
        robj = '0.1.0',
        dtype = .H5AD.guessDType('0.1.0'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'shape',
        robj = rev(dim(x)),
        dtype = .H5AD.guessDType(dim(x))
    )

    return(invisible(NULL))
}

.writeDataFrameToH5AD <- function(
        df,
        dfile,
        dname
) {
    dnamePaths <- unlist(strsplit(dname, '/'))
    # Recursively create groups
    for (i in seq_along(dnamePaths)) {
        search <- paste0(dnamePaths[1:i], collapse = '/')
        if (!dfile$exists(name = search)) {
            dfile$create_group(name = search)
        }
    }
    # Add index
    dfile[[dname]]$create_dataset(
        name = '_index',
        robj = rownames(df),
        dtype = .H5AD.guessDType(rownames(df))
    )
    dfile[[dname]]$create_attr(
        attr_name = '_index',
        robj = '_index',
        dtype = .H5AD.guessDType('index'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    # Add columns
    for (i in colnames(df)) {
        path <- paste0(dname, '/', i)
        if (is.factor(df[[i]])) {
            # Writing categorical
            .writeCategroricalToH5AD(df[[i]], dfile, path)
        } else {
            dfile[[dname]]$create_dataset(
                name = i,
                robj = df[[i]],
                dtype = .H5AD.guessDType(df[[i]])
            )
        }
    }
    dfile[[dname]]$create_attr(
        attr_name = 'column-order',
        robj = colnames(df),
        dtype = .H5AD.guessDType(colnames(df))
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-type',
        robj = 'dataframe',
        dtype = .H5AD.guessDType('dataframe'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-version',
        robj = '0.2.0',
        dtype = .H5AD.guessDType('0.2.0'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    return(invisible(NULL))
}

.writeCategroricalToH5AD <- function(x, dfile, dname) {
    dfile$create_group(name = dname)
    dfile[[dname]]$create_dataset(
        name = 'categories',
        robj = levels(x),
        dtype = .H5AD.guessDType(levels(x)),
        chunk_dims = NULL
    )
    dfile[[dname]]$create_dataset(
        name = "codes",
        robj = as.integer(x) - 1L,
        dtype = .H5AD.guessDType(as.integer(x) - 1L)
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-type',
        robj = 'categorical',
        dtype = .H5AD.guessDType('categorical'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-version',
        robj = '0.2.0',
        dtype = .H5AD.guessDType('0.2.0'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'ordered',
        robj = FALSE,
        dtype = .H5AD.guessDType(FALSE),
        space = hdf5r::H5S$new(type = "scalar")
    )
    return(invisible(NULL))
}

.createMapping <- function(dfile, dname) {
    dfile$create_group(name = dname)
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-type',
        robj = 'dict',
        dtype = .H5AD.guessDType('dict'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-version',
        robj = '0.1.0',
        dtype = .H5AD.guessDType('0.1.0'),
        space = hdf5r::H5S$new(type = "scalar")
    )
}

.H5AD.guessDType <- function(x, stype = "utf8", ...) {
    dtype <- hdf5r::guess_dtype(x = x, ...)
    if (inherits(dtype, "H5T_STRING")) {
        dtype <- .H5AD.stringType(stype = stype)
    }
    else if (inherits(dtype, "H5T_COMPOUND")) {
        cpd.dtypes <- dtype$get_cpd_types()
        for (i in seq_along(cpd.dtypes)) {
            if (inherits(cpd.dtypes[[i]], "H5T_STRING")) {
                cpd.dtypes[[i]] <- .H5AD.stringType(stype = stype)
            }
        }
        dtype <- hdf5r::H5T_COMPOUND$new(
            labels = dtype$get_cpd_labels(),
            dtypes = cpd.dtypes,
            size = dtype$get_size()
        )
    }
    else if (inherits(dtype, "H5T_LOGICAL")) {
        dtype <- hdf5r::guess_dtype(x = .boolToInt(x = x), ...)
    }
    return(dtype)
}

.boolToInt <- function(x) {
    x <- as.integer(x)
    x[is.na(x)] <- 2L
    return(x)
}

.H5AD.stringType <- function(stype = c("utf8", "ascii7"))
{
    stype <- match.arg(arg = stype)
    switch(
        EXPR = stype,
        utf8 = hdf5r::H5T_STRING$new(size = Inf)$set_cset(
            cset = hdf5r::h5const$H5T_CSET_UTF8
        ),
        ascii7 = hdf5r::H5T_STRING$new(size = 7L)
    )
}

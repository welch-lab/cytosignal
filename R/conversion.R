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
#' and binary matrices indicating whether an interaction is significant at a
#' spot. The spatial coordinates are coerced to a "SlideSeq" object regardless
#' of the exact technology being imported from. Interpretable interaction names,
#' rather than the CPI IDs of each interaction, are available as feature
#' metadata of the CytoSignal assay.
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
        signifMats <- lapply(signifMats, function(res) {
            # Here, res is a list, where each element is a character vector of
            # spot IDs that are significant for a given interaction.
            allIntr <- names(res)
            sparseIJ <- list()
            for (i in seq_along(res)) {
                IJ <- data.frame(
                    i = rep(i, length(res[[i]])),
                    j = match(res[[i]], rownames(loc))
                )
                sparseIJ[[i]] <- IJ
            }
            sparseIJ <- do.call(rbind, sparseIJ)
            methods::as(Matrix::sparseMatrix(
                i = sparseIJ$i, j = sparseIJ$j, x = 1,
                dims = c(length(allIntr), nrow(loc)),
                dimnames = list(allIntr, rownames(loc))
            ), "CsparseMatrix")
        })
        csAssay <- SeuratObject::CreateAssay5Object(data = c(scoreMats, signifMats))
        intrNames <- getIntrNames(object, rownames(csAssay))
        csAssay[["name"]] <- intrNames
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
#' @param layer Data from the H5AD file to extract. Default \code{"X"} extracts
#' \code{adata.X}. Use \code{"raw/X"} to extract \code{adata.raw.X}. Use
#' \code{"layers/key"} to extract layer data with key \code{"key"}.
#' @return \code{csToH5AD} writes an H5AD file to disk and returns \code{NULL}.
#' \code{H5ADToCS} returns a CytoSignal object.
csToH5AD <- function(object, filename, overwrite = FALSE) {

}

#' @rdname csToH5AD
#' @export
H5ADToCS <- function(filename, layer = "X") {

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

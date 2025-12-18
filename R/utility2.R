
# index - integer, character, or logical selector vector
# allAvail - vector of all elements to be selected from
# NSelect - NULL any number, or an exact number to check.
.checkValid.Index <- function(index, allAvail, NSelect = NULL) {
    if (length(index) > 0) {
        if (is.numeric(index)) {
            if (index != as.integer(index)) {
                cli::cli_abort('Numeric index should be an integer.')
            }
            if (any(index < 1) || any(index > length(allAvail))) {
                cli::cli_abort('Numeric index out of bounds.')
            }
        } else if (is.character(index)) {
            if (any(!index %in% allAvail)) {
                cli::cli_abort('Character ID not found.')
            }
            index <- match(index, allAvail)
        } else if (is.logical(index)) {
            if (length(index) != allAvail) {
                cli::cli_abort('Logical index length mismatch.')
            }
            index <- which(index)
        } else {
            cli::cli_abort('Index should be integer, character, or logical.')
        }
    }
    if (!is.null(NSelect)) {
        if (length(index) != NSelect) {
            cli::cli_abort('Index should select exactly {NSelect} element{?s}.')
        }
    }
    return(index)
}

.sparsity <- function(x, digits = 2) {
    if (inherits(x, 'dgCMatrix')) {
        s <- 1 - length(x@x) / (nrow(x) * ncol(x))
    } else {
        s <- 1 - sum(x != 0) / (length(x))
    }
    round(s, digits = digits)
}

validGenes <- function(object) {
    rawData <- object@rawData
    dbGenes <- .uniqGeneInDB(object@intrDB)
    geneUse <- intersect(rownames(rawData), dbGenes)
    sort(geneUse)
}

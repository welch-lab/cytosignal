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

# Sanity check for interaction selection
# Only use in an user exposed function!
# csdb - the csdb from object
# interaction - all possible types of interaction selection
# signifOnly - Whether to only consider significant interactions when `interaction = NULL`
# fdrThresh, minExp, minSpot - parameters for filtering significant interactions
# error - TRUE, abort when invalid selection is made; FALSE, warn and ignore invalid selection
.checkIntrSelection <- function(
        object,
        interaction,
        signifOnly = FALSE,
        fdrThresh = 0.05,
        minExp = 100,
        minSpot = 100,
        error = TRUE
) {
    csdb <- intrDB(object)
    if (is.null(csdb)) {
        cli::cli_abort("No interaction database available.")
    }

    if (is.null(interaction)) {
        if (signifOnly) {
            return(getTopIntr(
                object,
                fdrThresh = fdrThresh,
                minExp = minExp,
                minSpot = minSpot
            ))
        } else {
            return(csdb)
        }
    }

    if (error) messager <- cli::cli_abort
    else messager <- function(texts) {
        for (i in seq_along(texts)) {
            txt <- texts[i]
            type <- names(texts)[i]
            if (is.null(type)) cli::cli_alert_warning(txt)
            else if (type == 'x') cli::cli_alert_danger(txt)
            else if (type == 'i') cli::cli_alert_info(txt)
            else if (type == 'v') cli::cli_alert_success(txt)
            else cli::cli_alert_warning(txt)
        }
    }

    if (inherits(interaction, 'csdb')) {
        if (nrow(csdb) == 0) {
            messager("Empty interaction database provided.")
        }
        interaction <- interaction$interactors
        if (is.null(interaction)) {
            cli::cli_abort("Field {.field interactors} not found from database specification")
        }
    }

    if (is.character(interaction)) {
        if (any(!interaction %in% csdb$interactors)) {
            notfound <- interaction[!interaction %in% csdb$interactors]
            messager(c(
                x = "{length(notfound)} (out of {length(interaction)}) interaction{?s} {?is/are} not available. ",
                i = "Unavailable ones: {.val {notfound}}"
            ))
        }
        interaction <- interaction[interaction %in% csdb$interactors]
        csdb <- csdb[match(interaction, csdb$interactors), , drop = FALSE]
    } else if (is.numeric(interaction)) {
        if (any(interaction > nrow(csdb) | interaction < 1)) {
            cli::cli_abort("Interaction index out of bounds.")
        }
        csdb <- csdb[interaction, , drop = FALSE]
    } else if (is.logical(interaction)) {
        if (length(interaction) != nrow(csdb)) {
            cli::cli_abort("Logical vector length does not match number of interactions in database.")
        }
        csdb <- csdb[interaction, , drop = FALSE]
    } else {
        # Get the name of the function that user called
        caller <- as.character(sys.call(-1)[1])
        cli::cli_abort(c(
            x = "Argument {.field interaction} does not accept a {.cls {class(interaction)[1]}} object.",
            i = "Please check {.code ?{caller}}"
        ))
    }
    return(csdb)
}


# This function returns a logical sparse matrix, indicating whether a spot
# is considered as a valid sender for each interaction, based on ligand
# component (potentially multiple genes) expression and interaction type.
# object - cytosignal2 object
# Return - sparse logical matrix of dimension nSpot x nIntr
senderLigandMask <- function(object) {
    csdb <- object@intrDB
    raw <- object@rawData
    # Convert to binary to mark whether a gene is expressed
    raw <- raw[validGenes(object), , drop = FALSE] > 0
    Lmap <- intrGeneMap(object, component = 'ligands')
    # senderNGene: nSender x nIntr, for a sum of ligand components expressed
    senderNGene <- t(raw) %*% Lmap
    # How many components are involved for the ligand of each interaction
    ligNGene <- colSums(Lmap)
    # ligNExpect: length nIntr, expected nunmber of ligand components
    # For contact dependent, all components in a complex must be found in a sender
    # For diffusion, at least one component is required, assuming that other
    # components may come from other senders in a neighborhood.
    ligNExpect <- ifelse(csdb$type == 'contact', ligNGene, 1)

    sweep(x = senderNGene, MARGIN = 2, STATS = ligNExpect, FUN = ">=")
}

# object - cytosignal2 object
# csdb - Returned valid-only interaction selection from .checkIntrSelection
# Return - data.frame of fields:
# sender - 1-based integer index, in the range of all spots
# receiver - 1-based integer index, in the range of all spots
# interaction - 1-based integer index, in the range of all interactions
# type - 1-based integer index, 1 for contact 2 for diffusion
getAllEdges <- function(
        object,
        csdb,
        fdrThresh = 0.05
) {
    raw <- object@rawData
    genes <- validGenes(object)
    raw <- raw[genes, , drop = FALSE]
    fdr <- object@significance$spatialFDR
    if (is.null(fdr)) {
        cli::cli_abort('Spatial FDR not found. Please run {.fn inferLRScore} first.')
    }
    # signifReceiverMtx - nSpot x nIntr logical matrix, whether a spot receives
    # significant signal for each interaction
    signifReceiverMtx <- fdr[, csdb$interactors, drop = FALSE] < fdrThresh
    signifReceiverMtx[is.na(signifReceiverMtx)] <- FALSE

    nSignifSpots <- colSums(signifReceiverMtx)
    if (any(nSignifSpots == 0)) {
        cli::cli_alert_warning('{sum(nSignifSpots == 0)} interaction{?s} do{?es/} not show significance in any spot. {?It/They} will be ignored.')
    }
    csdb <- csdb[nSignifSpots > 0, , drop = FALSE]
    # senderMask - nSpot x nIntr logical matrix, whether a spot is a valid
    # sender for each interaction (i.e. expresses ligand components)
    senderMask <- senderLigandMask(object)
    graphCont <- object@neighborCont
    graphDiff <- object@neighborDiff
    dimnames(graphCont) <- list(colnames(raw), colnames(raw))
    dimnames(graphDiff) <- list(colnames(raw), colnames(raw))
    edgeDFs <- list()

    cli::cli_progress_bar(sprintf('Finding all edge for %d significant interactions', nrow(csdb)), total = nrow(csdb))
    for (i in seq_len(nrow(csdb))) {
        intr <- csdb$interactors[i]
        intrType <- csdb$type[i]
        receiverIdxBool <- signifReceiverMtx[, intr]

        if (intrType == 'contact') graph <- graphCont
        else if (intrType == 'diffusion') graph <- graphDiff
        else {
            cli::cli_warn("Unknown interaction type {.val {intrType}} for {.val {intr}}. Ask maintainer to fix it.")
            next
        }

        # senderIdxBool - eventually a logical masking the spots that are in the
        # valid neighborhood and express ligand components
        senderIdx <- graph[, receiverIdxBool, drop = FALSE]@i + 1
        if (length(senderIdx) == 0) {
            cli::cli_alert_warning('No ligand-expressing senders found connecting to significant receivers for interaction {.val {intr}}. Skipping.')
            next
        }
        senderIdxBool <- seq_len(ncol(raw)) %in% senderIdx
        senderIdxBool <- senderIdxBool & senderMask[, intr]

        graphSub <- graph[senderIdxBool, receiverIdxBool, drop = FALSE]
        edgeDFs[[i]] <- data.frame(
            sender = which(senderIdxBool)[graphSub@i + 1],
            receiver = rep(which(receiverIdxBool), diff(graphSub@p)),
            intr = factor(intr),
            type = intrType
        )
        cli::cli_progress_update()
    }
    cli::cli_progress_done()
    dplyr::bind_rows(edgeDFs) %>% dplyr::as_tibble()
}

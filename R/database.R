#' CytoSignal Database class
#' @rdname csdb-class
#' @name csdb
#' @description
#' A tibble data frame with the following columns:
#' \itemize{
#'   \item{\code{interactors}: Interaction ID in the format "partner_a-partner_b",
#'   where \code{partner_a} and \code{partner_b} are gene names of ligand and receptor,
#'   respectively. For complexes, gene names are joined by "+".}
#'   \item{\code{ligands}: Gene names of ligand(s), separated by semicolon if multiple.}
#'   \item{\code{receptors}: Gene names of receptor(s), separated by semicolon if multiple.}
#'   \item{\code{type}: Type of interaction, either \code{'diffusion'} or \code{'contact'}.}
#'   \item{\code{alt_name}: Alternative name for the interaction, using complex names
#'   if applicable.}
#' }
NULL


#' Download and make DB from public CellPhoneDB data
#' @export
#' @param path URL or local path to CellPhoneDB-data folder.
#' @param saveLoc Local path to save downloaded CellPhoneDB data if downloading
#' from URL. Defaults to a temporary directory.
#' @param interactionInput Relative path to interaction input CSV file within
#' the CellPhoneDB data folder.
#' @param proteinInput Relative path to protein input CSV file within
#' the CellPhoneDB data folder.
#' @param complexInput Relative path to complex input CSV file within
#' the CellPhoneDB data folder.
#' @param geneInput Relative path to gene input CSV file within
#' the CellPhoneDB data folder.
#' @return
#' A \code{\link{csdb-class}} object.
getCellPhoneDB <- function(
        path = 'https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v5.0.0.zip',
        saveLoc = paste0(tempdir(), '_cellphonedb'),
        interactionInput = 'data/interaction_input.csv',
        proteinInput = 'data/protein_input.csv',
        complexInput = 'data/complex_input.csv',
        geneInput = 'data/gene_input.csv'
) {
    # Create the save location directory if it doesn't exist
    if (!dir.exists(saveLoc)) {
        dir.create(saveLoc, recursive = TRUE)
    }

    # Is path an online url or local path
    if (grepl("^https?://", path)) {
        # Download the ZIP or tarball
        download.file(path, file.path(saveLoc, 'cellphonedb_data.zip'))
        # Unzip the file
        content_paths <- unzip(file.path(saveLoc, 'cellphonedb_data.zip'), exdir = saveLoc)
        path_split <- strsplit(content_paths, split = .Platform$file.sep)
        shortest <- min(lengths(path_split))
        common_path_split <- path_split[[1]][1:(shortest - 1)]
        path <- do.call(file.path, as.list(common_path_split))
    } else if (dir.exists(path)) {
        path <- normalizePath(path)
    } else {
        cli::cli_abort("The provided path is neither a valid URL nor an existing local directory.")
    }

    cli::cli_alert_success("Using data at folder: {.file {path}}")
    interactionInput <- file.path(path, interactionInput)
    if (file.exists(interactionInput)) {
        cli::cli_alert_success("Found interaction input file: {.file {interactionInput}}")
    } else {
        cli::cli_abort("Interaction input file not found at: {.file {interactionInput}}")
    }
    proteinInput <- file.path(path, proteinInput)
    if (file.exists(proteinInput)) {
        cli::cli_alert_success("Found protein input file: {.file {proteinInput}}")
    } else {
        cli::cli_abort("Protein input file not found at: {.file {proteinInput}}")
    }
    complexInput <- file.path(path, complexInput)
    if (file.exists(complexInput)) {
        cli::cli_alert_success("Found complex input file: {.file {complexInput}}")
    } else {
        cli::cli_abort("Complex input file not found at: {.file {complexInput}}")
    }
    geneInput <- file.path(path, geneInput)
    if (file.exists(geneInput)) {
        cli::cli_alert_success("Found gene input file: {.file {geneInput}}")
    } else {
        cli::cli_abort("Gene input file not found at: {.file {geneInput}}")
    }

    interactionInput <- read.csv(interactionInput)
    proteinInput <- read.csv(proteinInput)
    complexInput <- read.csv(complexInput)
    geneInput <- read.csv(geneInput)

    interactionInput <- interactionInput %>%
        filter(
            # .data[['directionality']] == 'Ligand-Receptor',
            .data[['is_ppi']] == 'True'
        ) %>%
        mutate(
            a_is_complex = .data[['partner_a']] %in% complexInput[,1],
            b_is_complex = .data[['partner_b']] %in% complexInput[,1]
        ) %>%
        mutate(
            a_tm_or_sec = case_when(
                .data[['a_is_complex']] & complexInput[match(.data[['partner_a']], complexInput[,1]), 'transmembrane'] ~ 'Transmembrane',
                .data[['a_is_complex']] & complexInput[match(.data[['partner_a']], complexInput[,1]), 'secreted'] ~ 'Secreted',
                !.data[['a_is_complex']] & proteinInput[match(.data[['partner_a']], proteinInput[,1]), 'transmembrane'] ~ 'Transmembrane',
                !.data[['a_is_complex']] & proteinInput[match(.data[['partner_a']], proteinInput[,1]), 'secreted'] ~ 'Secreted',
                TRUE ~ 'Other'
            ),
            b_tm = case_when(
                .data[['b_is_complex']] & complexInput[match(.data[['partner_b']], complexInput[,1]), 'transmembrane'] ~ 'Transmembrane',
                !.data[['b_is_complex']] & proteinInput[match(.data[['partner_b']], proteinInput[,1]), 'transmembrane'] ~ 'Transmembrane',
                TRUE ~ 'Other'
            ),
            type = case_when(
                .data[['a_tm_or_sec']] == 'Secreted' & .data[['b_tm']] == 'Transmembrane' ~ 'diffusion',
                .data[['a_tm_or_sec']] == 'Transmembrane' & .data[['b_tm']] == 'Transmembrane' ~ 'contact',
                TRUE ~ 'other'
            )
        ) %>%
        filter(.data[['type']] %in% c('diffusion', 'contact')) %>%
        mutate(
            ligands = NA, receptors = NA, alt_name = '',
            type = factor(type)
        )
    for (i in seq_len(nrow(interactionInput))) {
        partner_a <- interactionInput$partner_a[i]
        do_alt <- FALSE
        if (interactionInput$a_is_complex[i]) {
            a_table <- complexInput
            a_alt <- partner_a
            do_alt <- TRUE
        } else {
            a_table <- proteinInput
            a_alt <- geneInput[match(partner_a, geneInput$uniprot), 'gene_name']
        }
        a_table_uniprot_col <- grep('uniprot', colnames(a_table), value = TRUE)
        a_interactors <- a_table[match(partner_a, a_table[,1]), a_table_uniprot_col] %>%
            unlist() %>%
            unname()
        a_interactors <- a_interactors[a_interactors != '']
        a_genes <- geneInput[match(a_interactors, geneInput$uniprot), 'gene_name']
        interactionInput$ligands[i] <- paste(a_genes, collapse = ';')

        partner_b <- interactionInput$partner_b[i]
        if (interactionInput$b_is_complex[i]) {
            b_table <- complexInput
            b_alt <- partner_b
            do_alt <- TRUE
        } else {
            b_table <- proteinInput
            b_alt <- geneInput[match(partner_b, geneInput$uniprot), 'gene_name']
        }
        b_table_uniprot_col <- grep('uniprot', colnames(b_table), value = TRUE)
        b_interactors <- b_table[match(partner_b, b_table[,1]), b_table_uniprot_col] %>%
            unlist() %>%
            unname()
        b_interactors <- b_interactors[b_interactors != '']
        b_genes <- geneInput[match(b_interactors, geneInput$uniprot), 'gene_name']
        interactionInput$receptors[i] <- paste(b_genes, collapse = ';')

        if (do_alt) {
            interactionInput$alt_name[i] <- paste0(a_alt, '-', b_alt)
        }
    }
    csdb <- interactionInput %>%
        select(
            .data[['interactors']],
            .data[['ligands']],
            .data[['receptors']],
            .data[['type']],
            .data[['alt_name']]
        ) %>%
        filter(!duplicated(.data[['interactors']])) %>%
        tibble::as_tibble()
    class(csdb) <- c('csdb', class(csdb))
    return(csdb)
}


#' Acess interaction database of a CytoSignal2 object
#' @rdname intrDB
#' @description
#' The setter method \code{intrDB<-} allows user to set the interaction
#' database for a \code{\linkS4class{cytosignal2}} object. "Lowword" case
#' conversion is allowed when using different species. The getter method
#' \code{intrDB} retrieves the interaction database in a
#' \code{\linkS4class{cytosignal2}} object, with optional subsetting by
#' ligand/receptor genes or interaction type.
#' @param object A \code{\linkS4class{cytosignal2}} object.
#' @param toLowwords Logical, whether to convert uppercase gene names from
#' database to "Lowword" format (e.g., "TP53" to "Tp53"). Default \code{FALSE}.
#' @param value A `csdb` object.
#' @return
#' For \code{intrDB<-}, the updated \code{\linkS4class{cytosignal2}} object with
#' the interaction database set.
#'
#' For \code{intrDB}, a \code{csdb} object, possibly subsetted by specified
#' ligand/receptor genes or interaction type.
#' @export
`intrDB<-` <- function(object, toLowwords = FALSE, value) {
    if (!inherits(object, 'cytosignal2')) {
        cli::cli_abort('{.field object} must be of class {.cls cytosignal2}.')
    }
    if (!inherits(value, 'csdb')) {
        cli::cli_abort('Assigned value must be of class {.cls csdb}.')
    }

    if (isTRUE(toLowwords)) {
        value$ligands <- strsplit(value$ligands, split = ';') %>%
            lapply(lowwords) %>%
            sapply(paste, collapse = ';')
        value$receptors <- strsplit(value$receptors, split = ';') %>%
            lapply(lowwords) %>%
            sapply(paste, collapse = ';')
    }

    dbgenes <- .uniqGeneInDB(value)
    common <- intersect(
        dbgenes,
        rownames(object@rawData)
    )
    if (length(common) < 100) {
        cli::cli_warn(
            c(`!` = 'Only {length(common)} gene{?s} {?overlaps/overlap} with the expression data. See argument {.field toLowwords} ',
              i = 'First three genes in database are: {.val {dbgenes[1:3]}}',
              i = 'First three genes in expression data are: {.val {rownames(object@rawData)[1:3]}}'
             )
        )
    }

    # Filter DB to include only intrs with all genes present in data
    maskL <- strsplit(value$ligands, split = ';') %>%
        sapply(function(genes) all(genes %in% rownames(object@rawData)))
    maskR <- strsplit(value$receptors, split = ';') %>%
        sapply(function(genes) all(genes %in% rownames(object@rawData)))
    value <- value[maskL & maskR, , drop = FALSE]
    if (nrow(value) == 0) {
        cli::cli_warn('No interaction has all genes present in expression data.')
    } else {
        cli::cli_alert_success('{nrow(value)} valid interaction{?s} loaded from database.')
    }
    methods::slot(object, 'intrDB') <- value
    return(object)
}

#' @rdname intrDB
#' @export
#' @inheritParams subsetCSDB
intrDB <- function(
        object,
        ligandGenes = NULL,
        receptorGenes = NULL,
        type = NULL,
        LorR = FALSE
) {
    if (!inherits(object, 'cytosignal2')) {
        cli::cli_abort('{.field object} must be of class {.cls cytosignal2}.')
    }
    csdb <- methods::slot(object, 'intrDB')
    if (is.null(csdb)) {
        cli::cli_abort('Interaction database not set in the object. Please set one using {.fn intrDB<-}.')
    }
    csdb <- subsetCSDB(
        csdb = csdb,
        ligandGenes = ligandGenes,
        receptorGenes = receptorGenes,
        type = type,
        LorR = LorR
    )
    return(csdb)
}

.uniqGeneInDB <- function(csdb) {
    if (!inherits(csdb, 'csdb')) {
        cli::cli_abort('Input csdb must be of class {.cls csdb}.')
    }
    c(csdb$ligands, csdb$receptors) %>%
        strsplit(split = ';') %>%
        unlist() %>%
        unname() %>%
        unique()
}


#' Subset cytosignal database object by genes
#' @description
#' This function allows user to limit the interaction database to only those
#' involving specified ligand and/or receptor genes, or type
#' (diffusion-/contact-dependent).
#' @param csdb A \code{\link{csdb-class}} object.
#' @param ligandGenes A character vector of regular expression to match ligand
#' genes to keep. Complexes are kept if any of their member genes match. Default
#' \code{NULL} skip the filtering by ligand genes.
#' @param receptorGenes A character vector of regular expression to match
#' receptor genes to keep. Complexes are kept if any of their member genes
#' match. Default \code{NULL} skip the filtering by receptor genes.
#' @param type Type of interaction to keep. Either \code{'diffusion'} or
#' \code{'contact'}. Default \code{NULL} skip the filtering by type.
#' @param LorR Logical, whether to keep interactions that match the ligand gene
#' pattern or the receptor gene pattern. Default \code{FALSE} keeps interactions
#' that match both patterns at the same time.
#' @return
#' The subset \code{csdb}.
#' @export
subsetCSDB <- function(
        csdb,
        ligandGenes = NULL,
        receptorGenes = NULL,
        type = NULL,
        LorR = FALSE
) {
    if (!inherits(csdb, 'csdb')) {
        cli::cli_abort('Input object must be of class {.cls csdb}.')
    }
    if (!is.null(type)) {
        type <- match.arg(arg = type, choices = c('diffusion', 'contact'))
        csdb <- csdb[csdb$type == type, ]
    }
    if (!is.null(ligandGenes)) {
        maskL <- sapply(ligandGenes, function(pattern) {
            grepl(pattern, csdb$ligands)
        }) %>%
            rowSums() > 0
    } else {
        maskL <- rep(TRUE, nrow(csdb))
    }
    if (!is.null(receptorGenes)) {
        maskR <- sapply(receptorGenes, function(pattern) {
            grepl(pattern, csdb$receptors)
        }) %>%
            rowSums() > 0
    } else {
        maskR <- rep(TRUE, nrow(csdb))
    }

    if (isTRUE(LorR)) {
        if (!is.null(ligandGenes) && !is.null(receptorGenes)) {
            mask <- maskL | maskR
        } else if (!is.null(ligandGenes)) {
            mask <- maskL
        } else if (!is.null(receptorGenes)) {
            mask <- maskR
        } else {
            mask <- rep(TRUE, nrow(csdb))
        }
    } else {
        mask <- maskL & maskR
    }
    csdb <- csdb[mask, , drop = FALSE]
    if (nrow(csdb) == 0) {
        cli::cli_warn('No interaction left after subsetting the database.')
    }
    return(csdb)
}

intrGeneMap <- function(
        object,
        type = c('diffusion', 'contact'),
        component = c('ligands', 'receptors')
) {
    typeUse <- match.arg(type, several.ok = TRUE) # rename variable to avoid dplyr confusion
    component <- match.arg(component)
    allGenes <- validGenes(object)
    db <- object@intrDB %>% filter(.data[['type']] %in% typeUse)
    intrGeneIdxList <- db %>%
        pull(.data[[component]]) %>%
        strsplit(split = ';') %>%
        lapply(match, table = allGenes)
    i <- unname(unlist(intrGeneIdxList))
    p <- c(0, cumsum(lengths(intrGeneIdxList)))
    x <- rep(1, length(i))
    Matrix::sparseMatrix(
        i = i,
        p = p,
        x = x,
        dims = c(length(allGenes), nrow(db)),
        dimnames = list(allGenes, db$interactors)
    )
}

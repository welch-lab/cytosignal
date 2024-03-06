#' Submit query to REVIGO and get result to R
#' @description
#' This function gets an input of vector of GO terms and submit query to REVIGO
#' via HTTP, and retrieves the result from the job.
#' The code is based on example provided by REVIGO, originally located at
#' http://revigo.irb.hr/FAQ
#'
#' Please visit http://revigo.irb.hr/ for the web-based interactive tool, and
#' please see https://doi.org/10.1371/journal.pone.0021800 for detailed
#' introduction to the functionality.
#' @param GOterms Character vector of GO term IDs.
#' @param cutoff Size of the result, between 0.4 to 0.9. Larger value returns
#' more result. Default \code{0.7}.
#' @param value A vector of additional information for each GO term. Default
#' \code{NULL} does not include additional information. Type of the value has
#' to be specified with \code{valueType}.
#' @param valueType Type of \code{value}. Default \code{"PValue"}. Optionally,
#' \code{"Higher"}, \code{"Lower"}, \code{"HigherAbsolute"},
#' \code{"HigherAbsLog2"}.
#' @param speciesTaxon Taxon ID of the species. Default \code{"0"} queries to
#' all available species. \code{"10090"} is for mouse and \code{"9606"} is for
#' human.
#' @param measure The semantic similarity measure to use. Default
#' \code{"SIMREL"}. Optionally, \code{"LIN"}, \code{"RESNIK"}, \code{"JIANG"}.
#' @param removeObsolete Logical. Whether to remove obsolete GO terms. Default
#' \code{TRUE}.
#' @return A data.frame with many information
#' @export
#' @examples
#' rvg <- revigo(c("GO:0003002", "GO:0007389", "GO:0007423", "GO:0009653"),
#'               speciesTaxon = 10090)
revigo <- function(
    GOterms,
    cutoff = 0.7,
    value = NULL,
    valueType =  c("PValue", "Higher", "Lower", "HigherAbsolute", "HigherAbsLog2"),
    speciesTaxon = "0",
    measure = c("SIMREL", "LIN", "RESNIK", "JIANG"),
    removeObsolete = TRUE
) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for this function to work. ",
         "Please install with command: `install.packages('httr')`")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for this function to work. ",
         "Please install with command: `install.packages('jsonlite')`")
  }
  if (!requireNamespace("stringi", quietly = TRUE)) {
    stop("Package 'stringi' is required for this function to work. ",
         "Please install with command: `install.packages('stringi')`")
  }
  if (!is.character(GOterms)) {
    stop("`GOterms` should be a character vector.")
  }
  if (cutoff < 0.4 || cutoff > 0.9) {
    warning("cutoff should be in [0.4, 0.9]. Taking default 0.7.")
    cutoff <- 0.7
  }
  cutoff <- as.character(cutoff)
  valueType <- match.arg(valueType)
  speciesTaxon <- as.character(speciesTaxon)
  measure <- match.arg(measure)
  removeObsolete <- tolower(as.character(removeObsolete))

  if (is.null(value)) {
    userData <- paste0(paste(GOterms, collapse = "\n"), "\n")
  } else {
    value <- as.character(value)
    rows <- paste0(GOterms, "\t", value)
    userData <- paste0(paste(rows, collapse = "\n"), "\n")
  }

  # Submit job to Revigo
  httr::POST(
    url = "http://revigo.irb.hr/StartJob",
    body = list(
      goList = userData,
      cutoff = cutoff,
      valueType = valueType,
      speciesTaxon = speciesTaxon,
      measure = measure,
      removeObsolete = removeObsolete
    ),
    # application/x-www-form-urlencoded
    encode = "form"
  ) -> res


  dat <- httr::content(res, encoding = "UTF-8")
  # fix when httr::content automatically converts json to list
  if (typeof(dat) != "list") {
    jobid <- jsonlite::fromJSON(dat, bigint_as_char = TRUE)$jobid
  } else {
    jobid <- dat$jobid
  }

  # Check job status
  running <- "1"
  while (running != "0" ) {
    httr::GET(
      url = "http://revigo.irb.hr/QueryJob",
      query = list( jobid = jobid, type = "jstatus" )
    ) -> res2
    dat2 <- httr::content(res2, encoding = "UTF-8")
    # fix when httr::content automatically converts json to list
    if (typeof(dat2) != "list") {
      running <- jsonlite::fromJSON(dat2)$running
    } else {
      running <- dat2$running
    }
    Sys.sleep(1)
  }

  # Fetch results
  httr::GET(
    url = "http://revigo.irb.hr/QueryJob",
    query = list(
      jobid = jobid,
      namespace = "1",
      type = "Scatterplot"
    )
  ) -> res3

  dat3 <- httr::content(res3, encoding = "UTF-8")

  # Write results to a file
  dat3 <- stringi::stri_replace_all_fixed(dat3, "\r", "")

  tmpres <- tempfile(fileext = ".tsv")
  dat3 <- gsub("null", "NA", dat3)
  cat(dat3, file = tmpres, fill = FALSE)

  final <- utils::read.table(tmpres, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)
  class(final) <- c("revigo", class(final))
  attributes(final)$valueType <- valueType
  attributes(final)$cutoff <- cutoff
  attributes(final)$speciesTaxon <- speciesTaxon
  attributes(final)$measure <- measure
  attributes(final)$removeObsolete <- removeObsolete
  return(final)
}

#' Create scatter plot from REVIGO result
#' @param rvg A data.frame from \code{\link{revigo}}. Must contain columns:
#' \code{PC_0}, \code{PC_1}, \code{LogSize}, \code{Frequency}, and \code{Name}.
#' @param labelBy From which field of the revigo result to retrieve label text.
#' Choose from \code{"Name"} or \code{"TermID"}. Default \code{"Name"}.
#' @param labelIdx Index matching to rows of \code{rvg} to indicate whether a
#' label has to be shown for the bubble. Default \code{NULL} automatically
#' shows tidy view with ggrepel.
#' @param labelSize Size of text label for marked dots. Default \code{4}.
#' @param max.overlaps \code{?ggrepel::geom_text_repel} parameter. Default
#' \code{5}.
#' @return A ggplot object
#' @export
#' @examples
#' rvg <- revigo(c("GO:0003002", "GO:0007389", "GO:0007423", "GO:0009653"),
#'               speciesTaxon = 10090)
#' plotREVIGO(rvg)
plotREVIGO <- function(
    rvg,
    labelBy = c("Name", "TermID"),
    labelIdx = NULL,
    labelSize = 4,
    max.overlaps = 5)
{
  labelBy <- match.arg(labelBy)
  original <- rvg[[labelBy]]
  rvg$LABEL <- NA
  if (!is.null(labelIdx)) {
    rvg$LABEL[labelIdx] <- original[labelIdx]
  } else {
    rvg$LABEL <- original
  }
  p <- ggplot2::ggplot(rvg, ggplot2::aes(x = .data[["PC_0"]],
                                         y = .data[["PC_1"]],
                                         size = .data[["LogSize"]],
                                         color = .data[["Value"]],
                                         label = .data[["LABEL"]])) +
    ggplot2::geom_point(alpha = 0.7)
  if (is.null(labelIdx)) {
    p <- p + ggrepel::geom_text_repel(color = "black", size = labelSize,
                                      max.overlaps = max.overlaps,
                                      min.segment.length = 0)
  } else {
    p <- p + ggrepel::geom_text_repel(color = "black", size = labelSize,
                                      max.overlaps = max.overlaps,
                                      min.segment.length = 0)
  }
  if (!is.null(attributes(rvg)$valueType)) {
    p <- p + ggplot2::scale_color_viridis_c(name = attributes(rvg)$valueType)
  } else {
    p <- p + ggplot2::scale_color_viridis_c()
  }
  p +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          legend.position = "bottom",
          axis.line = ggplot2::element_blank())
}

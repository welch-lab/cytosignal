#' Identify nearest neighbors for each location using different strategies for
#' different types of interactions
#' @description
#' This is a wrapper function of \code{\link{findNNGauEB}},
#' \code{\link{findNNDT}} and \code{\link{findNNRaw}}. We use a Gaussian Epsilon
#' ball to identify the nearest neighbors that can contribute to the diffusible
#' ligands. And then uses a Delaunay triangulation to identify the nearest
#' neighbors that can contribute to the contact-dependent ligands. The
#' identified nearest neighbors will then be used for imputing the \eqn{L} or
#' \eqn{R} value for each location.
#' @param object A \linkS4class{CytoSignal} object.
#' @param eps Gaussian Epsilon Ball method parameter. A numeric scalar. The
#' radius of the Gaussian Epsilon ball. Default \code{NULL} uses parameter
#' inferred with \code{\link{inferEpsParams}}.
#' @param sigma Gaussian Epsilon Ball method parameter. The \eqn{\sigma} of the
#' Gaussian kernel. Default \code{NULL} uses parameter inferred with
#' \code{\link{inferEpsParams}}.
#' @param self.weight Gaussian Epsilon Ball method parameter. Weight of the
#' index cell. Use a number between 0-1 or a string \code{"auto"} or
#' \code{"sum_1"}. Default \code{"auto"}.
#' @param weight Delaunay Triangulation method parameter. A numeric scalar for
#' the sum of the weights of the edges of the Delaunay Triangulation. Default
#' \code{2}.
#' @param max.r Delaunay Triangulation method parameter. A numeric scalar for
#' the maximum radius of the edges. Default \code{NULL} uses parameter inferred
#' with \code{\link{inferEpsParams}}.
#' @return A \linkS4class{CytoSignal} object updated. \code{object@imputation}
#' slot will be updated with three new entries: \code{object@imputation$GauEps},
#' \code{object@imputation$DT} and \code{object@imputation$Raw}.
#' @export
#' @examples
#' \dontrun{
#' object <- findNN(object)
#' }
findNN <- function(
    object,
    # GauEB params
    eps = NULL,
    sigma = NULL,
    diff.weight = "auto",
    # DT params
    dt.weight = 2,
    max.r = NULL
) {
    object <- findNNGauEB(object, eps = eps, sigma = sigma, self.weight = diff.weight)
    object <- findNNDT(object, weight = dt.weight, max.r = max.r)
    object <- findNNRaw(object)
    return(object)
}

#' Impute the \eqn{L} or \eqn{R} value from the nearest neighbors of each
#' location
#' @description
#' After running \code{\link{findNN}}, we can impute the \eqn{L} or \eqn{R}
#' value from the nearest neighbors of each location basing on the types of
#' nearest neighbors.
#' @param object A \linkS4class{CytoSignal} object, with \code{\link{findNN}}
#' already run.
#' @param weights The method to transform distances to weights. Choose from \code{"none"}, \code{"mean"},
#' \code{"counts"} and \code{"dist"}. Default \code{"none"} since weights are already calculated by Gaussian
#' kernel and DT mean.
#' @return A \linkS4class{CytoSignal} object updated. Entries in
#' \code{object@imputation} slot will be updated with the imputation values.
#' @export
#' @examples
#' \dontrun{
#' object <- findNN(object)
#' object <- imputeLR(object)
#' }
imputeLR <- function(
    object,
    weights = c("none", "mean", "counts", "dist")
) {
  weights = match.arg(weights)
  object <- imputeNiche(object, nn.type = "GauEps", weights = weights)
  object <- imputeNiche(object, nn.type = "DT", weights = weights)

  return(object)
}

#' Calculate LRScore from the imputed L and R values
#' @description
#' After running \code{\link{imputeLR}}, we can calculate the LR score for each
#' location. The LR score is calculated as the product of the imputed \eqn{L}
#' and \eqn{R} values. With the LR score inferred, we subsequently perform
#' permutation tests to construct the null distribution for testing the
#' significance of the interactions.
#' @param object A \linkS4class{CytoSignal} object, with \code{\link{imputeLR}}
#' already run.
#' @param recep.smooth A logical scalar. Whether to use the smoothed \eqn{R}
#' values which is imputed with DT method. Default \code{FALSE}.
#' @param intr.type A character vector. The type of interactions to calculate
#' the LR score. Choose from one or both of \code{"diff"} and \code{"cont"} for
#' diffusion-dependent and contact-dependent interactions, respectively. Default
#' uses both.
#' @param perm.size A numeric scalar. The number of permutations to perform.
#' Default \code{1e5} times.
#' @param norm.method The normalization method to apply to the imputed L and R
#' @param numCores SPARK::sparkx parameter. The number of cores to use. Default
#' \code{1}.
#' @return A \linkS4class{CytoSignal} object updated. Entries in
#' \code{object@lrscore} slot will be updated with the LR scores and the
#' significance inferrence. When \code{recep.smooth} is by default \code{FALSE},
#' the \code{object@lrscore} slot will be updated with
#' \code{object@lrscore$`GauEps-Raw`} and \code{object@lrscore$`DT-Raw`}. When
#' \code{recep.smooth} is \code{TRUE}, the \code{object@lrscore} slot will be
#' updated with \code{object@lrscore$`GauEps-DT`} and
#' \code{object@lrscore$`DT-DT`}.
#' @export
#' @examples
#' \dontrun{
#' object <- findNN(object)
#' object <- imputeLR(object)
#' object <- inferIntrScore(object)
#' }
inferIntrScore <- function(
    object,
    recep.smooth = FALSE,
    intr.type = c("diff", "cont"),
    perm.size = 1e5,
    # fdr.method = c("spatialFDR", "fdr"),
    # p.value = 0.05,
    # reads.thresh = 100,
    # sig.thresh = 100,
    norm.method = c("default", "cpm", "none", "scanpy"),
    numCores = 1
) {
  # fdr.method <- match.arg(fdr.method)
  norm.method <- match.arg(norm.method)
  if (isTRUE(recep.smooth)) {
    recep.slot <- "DT"
  } else {
    recep.slot <- "Raw"
  }
  if (any(!intr.type %in% c("diff", "cont")) ||
      is.null(intr.type)) {
    stop("Please specify `intr.type` with one or two of the following values: 'diff', 'cont'")
  }
  if ("diff" %in% intr.type) {
    message("== Calculating diffusible ligand-receptor scores ==")
    # Calculate the ligand-receptor scores for diffusible ligands and raw receptors
    object <- inferScoreLR(object, lig.slot = "GauEps", recep.slot = recep.slot, 
                           norm.method = norm.method, intr.db.name = "diff_dep")
    # Permutation test to calculate the null distribution of the imputed ligands and receptors
    object <- permuteLR(object, perm.size = perm.size)
    # Calculate the null distribution of the ligand-receptor scores
    object <- inferNullScoreLR(object)
    # Infer the significant ligand-receptor interactions by comparing real scores with the null distribution
    # object <- inferSignif(object, fdr.method = fdr.method, p.value = p.value, reads.thresh = reads.thresh, sig.thresh = sig.thresh)
    # object <- rankIntrSpatialVar(object, numCores = numCores, verbose = FALSE)
  }
  if ("cont" %in% intr.type) {
    message("== Calculating contact-dependent ligand-receptor scores ==")
    object <- inferScoreLR(object, lig.slot = "DT", recep.slot = recep.slot, intr.db.name = "cont_dep",
                           norm.method = norm.method)
    object <- permuteLR(object, perm.size = perm.size)
    object <- inferNullScoreLR(object)
    if (!isTRUE(recep.smooth)) {
      object <- inferScoreLR(object, lig.slot = "Raw", recep.slot = "Raw", intr.db.name = "cont_dep",
                             norm.method = norm.method)
      object <- permuteLR(object, perm.size = perm.size)
      object <- inferNullScoreLR(object)
    }
    # object <- inferSignif(object, fdr.method = fdr.method, p.value = p.value, reads.thresh = reads.thresh, sig.thresh = sig.thresh)
    # object <- rankIntrSpatialVar(object, numCores = numCores, verbose = FALSE)
  }
  return(object)
}

#' Impute time derivative of \eqn{L} or \eqn{R} from the nearest neighbors of
#' each location
#' @description
#' After running \code{\link{findNN}}, we can impute the temporal \eqn{L} or
#' \eqn{R} change from the nearest neighbors of each location basing on the
#' types of nearest neighbors.
#' @param object A \linkS4class{CytoSignal} object, with \code{\link{findNN}}
#' already run and velocity information added with \code{\link{addVelo}}.
#' @return A \linkS4class{CytoSignal} object updated. Entries in
#' \code{object@imputation} slot will be updated with the imputation values.
#' @export
#' @examples
#' \dontrun{
#' object <- findNN(object)
#' object <- imputeVeloLR(object)
#' }
imputeVeloLR <- function(
    object
) {
  object <- imputeNicheVelo(object, nn.type = "GauEps")
  object <- imputeNicheVelo(object, nn.type = "DT")
  object <- imputeNicheVelo(object, nn.type = "Raw")
  return(object)
}

#' Calculate the interaction velocity from the imputed time derivative of
#' \eqn{L} or \eqn{R} values
#' @description
#' After running \code{\link{imputeVeloLR}}, we can calculate the interaction
#' velocity for each location. Please refer to the manualscript for detail of
#' the calculation.
#' @param object A \linkS4class{CytoSignal} object, with
#' \code{\link{imputeVeloLR}} already run.
#' @param recep.smooth A logical scalar. Whether to use the smoothed \eqn{R}
#' values which is imputed with DT method. Default \code{FALSE}.
#' @param norm.method The normalization method to apply to the velocity data,
#' need to be consistent with the normalization method used when generating the
#' input RNA velocity. Please consult the tool used for it, e.g. scVelo,
#' VeloVAE. Default is \code{"scanpy"}. Can choose from \code{"scanpy"},
#' \code{"cpm"}, \code{"default"} or \code{"none"}.
#' @return A \linkS4class{CytoSignal} object updated. Entries in
#' \code{object@lrvelo} slot will be updated with the velocity scores. When
#' \code{recep.smooth} is by default \code{FALSE}, the \code{object@lrvelo}
#' slot will be updated with \code{object@lrvelo$`GauEps-Raw`} and
#' \code{object@lrvelo$`DT-Raw`}. When \code{recep.smooth} is \code{TRUE}, the
#' \code{object@lrvelo} slot will be updated with
#' \code{object@lrvelo$`GauEps-DT`} and \code{object@lrvelo$`DT-DT`}.
#' @export
#' @examples
#' \dontrun{
#' object <- addVelo(object, velo.s, velo.u)
#' object <- findNN(object)
#' object <- imputeVeloLR(object)
#' object <- inferIntrVelo(object)
#' }
inferIntrVelo <- function(
    object,
    recep.smooth = FALSE,
    norm.method = "scanpy" # depending on velocity method
) {
  if (isTRUE(recep.smooth)) {
    recep.slot <- "DT"
  } else {
    recep.slot <- "Raw"
  }
  message("== Calculating diffusion-dependent interaction velocity ==")
  object <- inferVeloLR(object, norm.method = norm.method, lig.slot = "GauEps", recep.slot = recep.slot, intr.db.name = "diff_dep")
  message("== Calculating contact-dependent interaction velocity ==")
  object <- inferVeloLR(object, norm.method = norm.method, lig.slot = "DT", recep.slot = recep.slot, intr.db.name = "cont_dep")
  return(object)
}



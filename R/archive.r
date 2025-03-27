# #' Plot significant interactions ranked by the user-specified metric
# #'
# #' @param object A cytosignal object
# #' @param num.plot Number of interactions to plot
# #' @param res_dir Directory to save the plots
# #' @param plot.details Whether to plot NULL imputed values and scores
# #' @param slot.use The LRscore slot to use for plotting
# #' @param signif.use The metric used to rank the interactions, by default "result.hq.pear"
# #' @param plot.clusters Whether to plot the clusters
# #' @param plot.velo Whether to plot the velocity
# #' @param colors.list A list of colors to use for plotting
# #' @param pt.size Size of the points
# #' @param pt.stroke Stroke of the points
# #' @param u_width Width of the plot
# #' @param u_hgt Height of the plot
# #' @param set.res Resolution of the plot
# #' @param return.plot Whether to return the plot
# #'
# #' @return A plot if return.plot is TRUE. Otherwise, plots are saved to the specified directory.
# #'
# #' @export
# plotSignif <- function(object, num.plot = NULL, res_dir, plot.details = T, slot.use = NULL, signif.use = NULL, plot.clusters = T,
#                     plot.velo = F, colors.list = NULL, pt.size=0.1, pt.stroke = 0.2, u_width = 6, u_hgt = 5, set.res = 200,
#                     return.plot = F
# ){
#     if (length(object@lrscore) == 0){
#         stop("No score matrix found!\n")
#     }

#     if (is.null(slot.use)){
#         slot.use = object@lrscore[["default"]]
#     }

#     score.obj = object@lrscore[[slot.use]]
#     intr.db <- object@intr.valid[[score.obj@intr.slot]]

#     if (is.null(signif.use)){
#         signif.use = "result.hq.pear"
#     }

#     if (!signif.use %in% names(score.obj@res.list)){
#         stop("No such significance level: ", signif.use, " found!\n")
#     }

#     if (is.null(num.plot)) {
#         index.len <- length(object@lrscore[[slot.use]]@res.list[[signif.use]])
#         num.plot <- intersect(union(1:20, (index.len-10):index.len), 1:index.len)
#     }

#     dge.raw <- object@counts
#     clusters <- object@clusters
#     lig.slot <- score.obj@lig.slot
#     recep.slot <- score.obj@recep.slot
#     dge.lig <- object@imputation[[lig.slot]]@imp.data
#     dge.recep <- object@imputation[[recep.slot]]@imp.data

#     # sample null values for plotting
#     sample.index <- sample(ncol(score.obj@lig.null), ncol(dge.lig))
#     null.dge.gau <- score.obj@lig.null[, sample.index]
#     null.dge.dt <- score.obj@recep.null[, sample.index]
#     colnames(null.dge.gau) <- colnames(null.dge.dt) <- colnames(dge.lig)
#     rownames(null.dge.gau) <- rownames(null.dge.dt) <- rownames(dge.lig)

#     # null.dge.dt <- dge.recep

#     # sample null scores for plotting
#     score.mtx <- score.obj@score
#     null.lrscore.mtx <- score.obj@score.null
#     null.lrscore.mtx <- null.lrscore.mtx[sample(nrow(null.lrscore.mtx), nrow(score.mtx)), ]
#     rownames(null.lrscore.mtx) <- rownames(score.mtx)

#     res.list <- score.obj@res.list[[signif.use]]
#     res.intr.list = names(res.list)
#     cells.loc = as.data.frame(object@cells.loc)
#     cells.loc = cells.loc[colnames(dge.lig), ]

#     if (is.null(colors.list)){
#         levels.use <- levels(clusters)
#         colors.list <- uniqueColors(length(levels.use))
#         # as.character(paletteer::paletteer_d("ggsci::default_igv",
#                         #n = length(levels.use)))
#         names(colors.list) = levels.use
#     }

#     if (nrow(cells.loc) != ncol(dge.lig)){
#         stop("Incorrect beads positions!\n")
#     }

#     cat("Now plotting INTRs...\nRank No.")

#     plots.list <- lapply(num.plot, function(i){ # for each intr in the res.list
#         cat(i, ", ", sep = "")
#         sub.df = cells.loc
#         null.sub.df = cells.loc

#         pair.index = which(object@intr.valid$intr.index$id_cp_interaction == res.intr.list[i])

#         ligand.name = object@intr.valid$intr.index[pair.index, 4]
#         if (ligand.name == ""){ # if complex
#             ligand.name = object@intr.valid$intr.index[pair.index, 2]
#         }
#         ligand.name = gsub("_HUMAN", "", ligand.name)

#         receptor.name = object@intr.valid$intr.index[pair.index, 5]
#         if (receptor.name == ""){
#             receptor.name = object@intr.valid$intr.index[pair.index, 3]
#         }
#         receptor.name = gsub("_HUMAN", "", receptor.name)

#         # generate interaction names
#         # intr.name = paste0(gsub("_HUMAN", "", ligand.name), "-", gsub("_HUMAN", "", receptor.name))
#         intr.name = paste0(ligand.name, "-", receptor.name)

#         colnames(sub.df) = c("x", "y")
#         colnames(null.sub.df) = c("x", "y")

#         ligands = dge.lig[names( intr.db[[2]][intr.db[[2]] == res.intr.list[i]] ), ]
#         ligands.ori = dge.raw[names(intr.db[[2]][intr.db[[2]] == res.intr.list[i]]), ]
#         ligands.null = null.dge.gau[names(intr.db[[2]][intr.db[[2]] == res.intr.list[i]]), ]

#         receps = dge.recep[names(intr.db[[3]][intr.db[[3]] == res.intr.list[i]]), ]
#         receps.ori = dge.raw[names(intr.db[[3]][intr.db[[3]] == res.intr.list[i]]), ]
#         receps.null = null.dge.dt[names(intr.db[[3]][intr.db[[3]] == res.intr.list[i]]), ]

#         if (!is.null(nrow(ligands))){ligands = Matrix::colSums(ligands)}
#         if (!is.null(nrow(receps))){receps = Matrix::colSums(receps)}
#         if (!is.null(nrow(ligands.ori))){ligands.ori = Matrix::colSums(ligands.ori)}
#         if (!is.null(nrow(ligands.null))){ligands.null = Matrix::colSums(ligands.null)}
#         if (!is.null(nrow(receps.ori))){receps.ori = Matrix::colSums(receps.ori)}
#         if (!is.null(nrow(receps.null))){receps.null = Matrix::colSums(receps.null)}

#         sub.df$ligands = ligands
#         sub.df$ligands_ori = ligands.ori
#         sub.df$receptors = receps
#         sub.df$receptors_ori = receps.ori
#         null.sub.df$ligands_null = ligands.null
#         null.sub.df$receptors_null = receps.null

#         p.lig = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = ligands))+
#             ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#             ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
#             cowplot::theme_cowplot(12)+
#             labs(color = paste0("Ligand \n", ligand.name), x = NULL, y = NULL)

#         p.lig.ori = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = ligands_ori))+
#             ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#             ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
#             cowplot::theme_cowplot(12)+
#             labs(color = paste0("Ligand-ori \n", ligand.name), x = NULL, y = NULL)

#         p.lig.null = ggplot2::ggplot(null.sub.df, ggplot2::aes(x, y, color = ligands_null))+
#             ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#             ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
#             cowplot::theme_cowplot(12)+
#             labs(color = paste0("Ligand-null \n", ligand.name), x = NULL, y = NULL)

#         p.recep = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = receptors))+
#             ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#             ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
#             cowplot::theme_cowplot(12)+
#             labs(color = paste0("Receptor \n", receptor.name), x = NULL, y = NULL)

#         p.recep.ori = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = receptors_ori))+
#             ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#             ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
#             cowplot::theme_cowplot(12)+
#             labs(color = paste0("Receptor-ori \n", receptor.name), x = NULL, y = NULL)

#         p.recep.null = ggplot2::ggplot(null.sub.df, ggplot2::aes(x, y, color = receptors_null))+
#             ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#             ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
#             cowplot::theme_cowplot(12)+
#             labs(color = paste0("Receptor-null \n", receptor.name), x = NULL, y = NULL)

#         p.list = list(p.lig, p.lig.ori, p.lig.null, p.recep, p.recep.ori, p.recep.null)
#         names(p.list) = c("lig", "lig.ori", "lig.null", "recep", "recep.ori", "recep.null")

#         ##### plot scores

#         if (ncol(dge.lig) != nrow(score.mtx)){
#             cat("Setting LR-scores of cells with no neighbors to empty.\n")
#             sub.df$scores = 0
#             sub.df[Matrix::rownames(score.mtx), "scores"] = score.mtx[, res.intr.list[i]]

#             null.sub.df$null_scores = 0
#             null.sub.df[Matrix::rownames(null.lrscore.mtx), "null_scores"] = null.lrscore.mtx[, res.intr.list[i]]
#         } else{
#             sub.df$scores = score.mtx[, res.intr.list[i]]
#             null.sub.df$null_scores = null.lrscore.mtx[, res.intr.list[i]]

#         }

#         p.score = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = scores))+
#             ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#             ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
#             cowplot::theme_cowplot(12)+
#             labs(x = NULL, y = NULL)

#         p.list[["score"]] = p.score

#         p.score.null = ggplot2::ggplot(null.sub.df, ggplot2::aes(x, y, color = null_scores))+
#             ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#             ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
#             cowplot::theme_cowplot(12)+
#             labs(x = NULL, y = NULL)

#         p.list[["score.null"]] = p.score.null

#         ##### plot significant beads

#         sub.df$group = "beads"
#         # sub.df$group[res.list[[i]]] = "significant"
#         sub.df[res.list[[i]], "group"] = "significant"

#         p.sig = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = group))+
#             ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#             # geom_segment(ggplot2::aes(xend = xend, yend = yend), data = edges.df, color = "red", size = 0.2)+
#             scale_colour_manual(values = c("beads" = "grey", "significant" = "red"))+
#             cowplot::theme_cowplot(12)+
#             labs(x = NULL, y = NULL)

#         p.list[["sig"]] = p.sig

#         if (plot.clusters){
#             clusters = clusters[rownames(sub.df)]
#             sub.df$clusters = as.character(clusters)

#             sig.sub.df = sub.df[sub.df$group == "significant", ]

#             p.clust = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = clusters))+
#                 ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#                 # scale_colour_viridis_d(direction = -1)+
#                 # ggsci::scale_color_nejm()+
#                 # ggsci::scale_color_igv(palette = "default")+
#                 scale_colour_discrete(type = colors.list, breaks = levels.use) +
#                 ggplot2::geom_point(data = sig.sub.df, ggplot2::aes(x, y), color = "red", size = pt.size*3, alpha = 0.15)+
#                 cowplot::theme_cowplot(12)+
#                 guides(color = guide_legend(override.aes = list(size = 4)))+
#                 labs(x = NULL, y = NULL)

#             p.list[["clust"]] = p.clust
#         }

#         if (plot.velo){
#             sub.df$velo = velo.mtx[, res.intr.list[i]]

#             p.velo = ggplot2::ggplot(sub.df, ggplot2::aes(x, y, color = velo))+
#                 ggplot2::geom_point(size = pt.size, stroke = pt.stroke)+
#                 ggplot2::scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
#                 cowplot::theme_cowplot(12)+
#                 labs(x = NULL, y = NULL)

#             p.list[["velo"]] = p.velo
#         }

#         if (return.plot){ return(p.list) }

#         plot.len = length(p.list)
#         # cat(paste0("Length of plot list: ", plot.len, "\n"))
#         intr.name = paste0(intr.name, "-n_", length(res.list[[i]]))

#         if (plot.details){
#             p.list <- p.list[c(1, 4, 7, 9, 3, 2, 5, 8, 10, 6)]
#             png(paste0(res_dir, "/Rank_", i, "_", intr.name, ".png"),
#                         res = set.res, width = u_width*4, height = u_hgt*2, unit="in")
#             print(cowplot::plot_grid(plotlist = p.list, nrow = 2))
#             dev.off()

#         } else{
#             p.list <- p.list[c(1, 4, 7, 9, 2, 5, 8, 10)]
#             png(paste0(res_dir, "/Rank_", i, "_", intr.name, ".png"),
#                         res = set.res, width = u_width*4, height = u_hgt*2, unit="in")
#             print(cowplot::plot_grid(plotlist = p.list, nrow = 2))
#             dev.off()
#         }

#     })

#     cat("End.\nFinished!\n")

#     if (return.plot){return(plots.list)}
# }



# #' Permute Imputation Results of specific imputation method
# #' 
# #' This function is a follow-up function of imputeNiche, which permutes the default or user-defined
# #' imputation method and stored the results in the ImpData object.
# #' 
# #' @param object A Cytosignal object
# #' @param nn.type The imputation method to use
# #' @param perm.size Size of the permutation test, by defualt NULL
# #' @param times Number of times to permute the whole dataset, by default 10
# #' 
# #' @return A Cytosignal object

# permuteImpLR <- function(
# 	object,
# 	nn.type = NULL,
# 	perm.size = NULL,
# 	times = 10
# ){
# 	if (is.null(nn.type)){
# 		nn.type <- object@imputation[["default"]]
# 	}

# 	if (!nn.type %in% names(object@imputation)){
# 		stop("Cannot find corresponding imputation method.")
# 	}

# 	use.gau <- grepl("GauEps", nn.type)
# 	use.dt <- grepl("DT", nn.type)
# 	use.raw <- grepl("Raw", nn.type)

# 	message("Permuting imputation data on method ", nn.type, "...")

# 	dimnames.list <- dimnames(object@imputation[[nn.type]]@imp.data)
# 	cells.loc <- object@cells.loc
# 	dge.raw <- object@counts

# 	if (!is.null(times)) {
# 		if (!is.numeric(times)){
# 			stop("times must be a numeric value.")
# 		}

# 		times <- ceiling(times)

# 		if (times < 1) {
# 			stop("times must be larger than 1.")
# 		}

# 	} else {
# 		# decide how many times to permute according to the size of the data
# 		if (is.null(perm.size)){
# 			message("Permutation size not specified, using 100000 instead.")
# 			perm.size <- 100000
# 		}

# 		if (!is.numeric(perm.size)){
# 			stop("perm.size must be a numeric value.")
# 		}

# 		if (perm.size < ncol(dge.raw)){
# 			message("Permutation size too small, using ", ncol(dge.raw), " instead.")
# 			perm.size <- ncol(dge.raw)
# 		}
		
# 		times <- ceiling(perm.size / ncol(dge.raw))
# 	}

# 	message("Permuting whole dataset ", times, " times...")

# 	null.cells.loc = cells.loc[sample(nrow(cells.loc)), ]
# 	rownames(null.cells.loc) = rownames(cells.loc)

# 	# get a random graph for permuting
# 	if (use.gau) {
# 		param.eps <- object@parameters$r.diffuse.scale
# 		param.sigma <- object@parameters$sigma.scale

# 		null.eps.nb <- findNNGauEB.matrix(null.cells.loc, eps = param.eps, 
# 						sigma = param.sigma, weight.sum = 2); gc()

# 		null.graph <- imputeNiche.dgCMatrix(dge.raw, null.eps.nb$id, null.eps.nb$dist, 
# 							weights = "none"); gc()
# 	} else if (use.dt) {
# 		null.nb.fac <- findNNDT.matrix(null.cells.loc)
# 		null.graph <- imputeNiche.dgCMatrix(dge.raw, null.nb.fac$id, null.nb.fac$dist, 
# 							weights = "none"); gc()
# 	} else if (use.raw) {
# 		cat("No need to generate a random graph for raw imputation.\n")
# 	} else {
# 		stop("Cannot find corresponding imputation method.")
# 	}

# 	cat("Permuting NULL expressions...\nTimes No.")

# 	null.dge.list <- lapply(1:times, function(i){
# 		## MUST shuffule raw dge!!
# 		cat(i, ", ", sep = "")
# 		perm.index <- sample(ncol(dge.raw))
# 		null.dge.raw <- dge.raw[, perm.index]
# 		scale.fac <- object@parameters[["lib.size"]][perm.index]
		
# 		colnames(null.dge.raw) <- colnames(dge.raw)
# 		names(scale.fac) <- colnames(dge.raw)

# 		scale.fac <- Matrix::Matrix(scale.fac, nrow = 1, byrow = T, sparse = T) # computing scale factor
# 		# null.dge.raw <- changeUniprot.matrix_like(null.dge.raw.all, object@intr.valid[["gene_to_uniprot"]])[[1]]

# 		if (use.raw) {
# 			null.dge.eps <- null.dge.raw
# 			scale.fac.imp <- scale.fac
# 		} else {
# 			null.dge.eps <- null.dge.raw %*% null.graph; gc()
# 			scale.fac.imp <- scale.fac %*% null.graph; gc()
# 		}
		
# 		null.dge.eps <- normCounts.list(list(mat=null.dge.eps, scale.fac=as.numeric(scale.fac.imp)), "default"); gc()
# 		rm(null.dge.raw); gc()

# 		return(null.dge.eps)
# 	})

# 	cat("End.\nFinished!\n")


# 	null.dge.eps.unpt = meanMat_cpp(null.dge.list, nrow(null.dge.list[[1]]), ncol(null.dge.list[[1]]))
# 	dimnames(null.dge.eps.unpt) <- dimnames.list
# 	rm(null.dge.list); gc()

#     object@imputation[[nn.type]]@imp.data.null <- null.dge.eps.unpt

# 	return(object)

# }


# #' Permute LR score for specific ligand-receptor imputation obj pairs
# #' 
# #' This function is a follow-up function of inferScoreLR. It computes the NUL LRscores 
# #' using the NULL imputation results and stores the results in the LR score object.
# #' The null distribution of the LR scores can be used to test the significance of the LR scores.
# #' 
# #' @param object A Cytosignal object
# #' @param slot.use The imputation method to use
# #' 
# #' @return A Cytosignal object

# permuteScoreLR <- function(
# 	object,
# 	slot.use = NULL
# ){
# 	if (is.null(slot.use)){
# 		slot.use <- object@lrscore[["default"]]
# 	}

# 	message("Permuting scores on Score slot: ", slot.use, "...")

# 	# score.obj <- object@lrscore[[slot.use]]
# 	lig.slot <- object@lrscore[[slot.use]]@lig.slot
# 	recep.slot <- object@lrscore[[slot.use]]@recep.slot
# 	cells.loc <- object@cells.loc
# 	intr.valid <- object@lrscore[[slot.use]]@intr.list

# 	null.dge.eps.unpt <- object@imputation[[lig.slot]]@imp.data.null
# 	null.dge.dt.unpt <- object@imputation[[recep.slot]]@imp.data.null

#     # perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
#     # null.cells.loc = cells.loc[perm.index,]
#     # rownames(null.cells.loc) = rownames(cells.loc)
#     # null.nb.fac = findNNDT.matrix(null.cells.loc); gc() # Mean num of neighbors: 44, median: 36

# 	### test permuting NULL mtx again
# 	# null.dge.eps.unpt = null.dge.eps.unpt[, sample(ncol(null.dge.eps.unpt))]
# 	# null.dge.dt.unpt = null.dge.dt.unpt[, sample(ncol(null.dge.dt.unpt))]

#     null.lrscore.mtx = inferScoreLR.dgCMatrix(null.dge.eps.unpt, null.dge.dt.unpt,
# 						intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()

# 	# null.lrscore.mtx <- null.dge.eps.unpt * null.dge.dt.unpt

#     if (sum(Matrix::colSums(null.lrscore.mtx) == 0) != 0){
#         message("A total of ", sum(Matrix::colSums(null.lrscore.mtx) == 0), " intr are empty in NULL scores.")
#         null.lrscore.mtx = null.lrscore.mtx[, !Matrix::colSums(null.lrscore.mtx) == 0]
#     }

# 	intr.both <- intersect(colnames(null.lrscore.mtx), colnames(object@lrscore[[slot.use]]@score))

# 	if (length(intr.both) != ncol(null.lrscore.mtx)) {
# 		message("Removing ", ncol(null.lrscore.mtx) - length(intr.both), " more intr from NULL scores.")
# 		null.lrscore.mtx = null.lrscore.mtx[, intr.both]
# 	}

# 	if (length(intr.both) != ncol(object@lrscore[[slot.use]]@score)) {
# 		message("Removing ", ncol(object@lrscore[[slot.use]]@score) - length(intr.both), " corresponding intr from REAL scores.")
# 		object@lrscore[[slot.use]]@score = object@lrscore[[slot.use]]@score[, intr.both]
# 	}

# 	object@lrscore[[slot.use]]@score.null <- null.lrscore.mtx

#     return(object)
# }





# #' Permute LR score for specific ligand-receptor imputation obj pairs
# #' 
# #' This function is a follow-up function of inferScoreLR. It computes the NUL LRscores 
# #' using the NULL imputation results and stores the results in the LR score object.
# #' The null distribution of the LR scores can be used to test the significance of the LR scores.
# #' 
# #' @param object A Cytosignal object
# #' @param slot.use The imputation method to use
# #' 
# #' @return A Cytosignal object
# #' @export
# permuteNicheScoreLR <- function(
# 	object,
# 	slot.use = NULL
# ){
# 	if (is.null(slot.use)){
# 		slot.use <- object@lrscore[["default"]]
# 	}

# 	message("Permuting scores on Score slot: ", slot.use, "...")

# 	# score.obj <- object@lrscore[[slot.use]]
# 	lig.slot <- object@lrscore[[slot.use]]@lig.slot
# 	recep.slot <- object@lrscore[[slot.use]]@recep.slot
# 	cells.loc <- object@cells.loc
# 	intr.valid <- object@lrscore[[slot.use]]@intr.list

# 	null.dge.eps.unpt <- object@imputation[[lig.slot]]@imp.data.null
# 	null.dge.dt.unpt <- object@imputation[[recep.slot]]@imp.data.null

#     perm.index = sample(seq_len(nrow(cells.loc)), nrow(cells.loc), replace = F)
#     null.cells.loc = cells.loc[perm.index,]
#     rownames(null.cells.loc) = rownames(cells.loc)
#     null.nb.fac = findNNDT.matrix(null.cells.loc); gc() # Mean num of neighbors: 44, median: 36

#     null.lrscore.mtx = graphNicheLR.dgCMatrix(null.dge.eps.unpt, null.dge.dt.unpt, null.nb.fac[["id"]],
# 						intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()

# 	# null.lrscore.mtx <- null.dge.eps.unpt * null.dge.dt.unpt
 
#     if (sum(Matrix::colSums(null.lrscore.mtx) == 0) != 0){
#         message("A total of ", sum(Matrix::colSums(null.lrscore.mtx) == 0), " intr are empty in NULL scores.")
#         null.lrscore.mtx = null.lrscore.mtx[, !Matrix::colSums(null.lrscore.mtx) == 0]
#     }

# 	intr.both <- intersect(colnames(null.lrscore.mtx), colnames(object@lrscore[[slot.use]]@score))

# 	if (length(intr.both) != ncol(null.lrscore.mtx)) {
# 		message("Removing ", ncol(null.lrscore.mtx) - length(intr.both), " more intr from NULL scores.")
# 		null.lrscore.mtx = null.lrscore.mtx[, intr.both]
# 	}

# 	if (length(intr.both) != ncol(object@lrscore[[slot.use]]@score)) {
# 		message("Removing ", ncol(object@lrscore[[slot.use]]@score) - length(intr.both), " corresponding intr from REAL scores.")
# 		object@lrscore[[slot.use]]@score = object@lrscore[[slot.use]]@score[, intr.both]
# 	}

# 	object@lrscore[[slot.use]]@score.null <- null.lrscore.mtx

#     return(object)
# }



# #' Permute Imputation Results of specific imputation method
# #' 
# #' This function is a follow-up function of imputeNiche, which permutes the default or user-defined
# #' imputation method and stored the results in the ImpData object.
# #' 
# #' @param object A Cytosignal object
# #' @param nn.type The imputation method to use
# #' @param perm.size Size of the permutation test
# #' 
# #' @return A Cytosignal object
# #' @export
# permuteRandomNB <- function(
# 	object,
# 	nn.type = NULL,
# 	perm.size = 100000
# ){
# 	if (is.null(nn.type)){
# 		nn.type <- object@imputation[["default"]]
# 	}

# 	if (!nn.type %in% names(object@imputation)){
# 		stop("Cannot find corresponding imputation method.")
# 	}

# 	use.gau <- grepl("GauEps", nn.type)
# 	use.dt <- grepl("DT", nn.type)
# 	use.raw <- grepl("Raw", nn.type)

# 	# if (!use.gau & !use.dt){
# 	# 	stop("Cannot find corresponding imputation method.")
# 	# }

# 	message("Permuting imputation data on method ", nn.type, "...")

# 	dimnames.list <- dimnames(object@imputation[[nn.type]]@imp.data)
# 	cells.loc <- object@cells.loc
# 	dge.raw <- object@counts

# 	# decide how many times to permute according to the size of the data
# 	if (!is.numeric(perm.size)){
# 		stop("perm.size must be a numeric value.")
# 	}

# 	if (perm.size < ncol(dge.raw)){
# 		message("Permutation size too small, using ", ncol(dge.raw), " instead.")
# 		perm.size <- ncol(dge.raw)
# 	}
	
# 	# reason for computing by batch is to save costs
# 	times <- ceiling(perm.size / ncol(dge.raw))
# 	each.size <- ceiling(perm.size / times)

# 	message("Permuting whole dataset ", times, " times...")

# 	# permute cell locations and remember to change the rownames!!
# 	null.cells.loc = cells.loc[sample(nrow(cells.loc)), ]
# 	rownames(null.cells.loc) = rownames(cells.loc)

# 	# use the same imputation graph for permutation
# 	graph.nn <- object@imputation[[nn.type]]@nn.graph

# 	cat("Permuting NULL expressions...\nTimes No.")

# 	null.dge.list <- lapply(1:times, function(i){
# 		## shuffule dge, remember to change the colnames!!
# 		## could also shuffle the null graph, mathematically the same
# 		cat(i, ", ", sep = "")

# 		# # if shuffle dge.raw
# 		# perm.index <- sample(ncol(dge.raw))
# 		# null.dge.raw <- dge.raw[, perm.index]
# 		# scale.fac <- object@parameters[["lib.size"]][perm.index]
# 		# colnames(null.dge.raw) <- colnames(dge.raw)
# 		# names(scale.fac) <- colnames(dge.raw)

# 		scale.fac <- object@parameters[["lib.size"]]
# 		scale.fac <- Matrix::Matrix(scale.fac, nrow = 1, byrow = T, sparse = T) # computing scale factor
# 		# null.dge.raw <- changeUniprot.matrix_like(null.dge.raw.all, object@intr.valid[["gene_to_uniprot"]])[[1]]

# 		if (use.raw) {
# 			# sample the dge.raw to control the size of the permutation
# 			# not really necessary, but put here for consistency
# 			sample.index <- sample(ncol(dge.raw), each.size)
# 			null.dge.eps <- dge.raw[, sample.index]
# 			scale.fac.imp <- scale.fac[sample.index]
# 		} else {
# 			##### core permutation step #####
# 			# null.graph <- shuffle_sp_mat_col(graph.nn)
# 			null.graph <- shuffleEdgeRandomNB(graph.nn)
# 			##### core permutation step #####

# 			# sample the null graph to control the size of the permutation
# 			null.graph <- null.graph[, sample(ncol(null.graph), each.size)]

# 			null.dge.eps <- dge.raw %*% null.graph
# 			scale.fac.imp <- scale.fac %*% null.graph
# 		}
		
# 		null.dge.eps <- normCounts.list(list(mat=null.dge.eps, scale.fac=as.numeric(scale.fac.imp)), "default")
# 		# rm(null.dge.raw)

# 		return(null.dge.eps)
# 	})

# 	null.dge.imp <- cbind_list(null.dge.list)
# 	rm(null.dge.list); gc()

# 	cat("End.\nFinished!\n")

# 	object@imputation[[nn.type]]@imp.data.null <- null.dge.imp

# 	return(object)

# }




# #' Permute LR score for specific ligand-receptor imputation obj pairs
# #' 
# #' This function is a follow-up function of inferScoreLR. It computes the NUL LRscores 
# #' using the NULL imputation results and stores the results in the LR score object.
# #' The null distribution of the LR scores can be used to test the significance of the LR scores.
# #' 
# #' @param object A Cytosignal object
# #' @param slot.use The imputation method to use
# #' 
# #' @return A Cytosignal object
# #' @export
# permuteRandomScore <- function(
# 	object,
# 	slot.use = NULL
# ){
# 	if (is.null(slot.use)){
# 		slot.use <- object@lrscore[["default"]]
# 	}

# 	message("Permuting scores on Score slot: ", slot.use, "...")

# 	# score.obj <- object@lrscore[[slot.use]]
# 	lig.slot <- object@lrscore[[slot.use]]@lig.slot
# 	recep.slot <- object@lrscore[[slot.use]]@recep.slot
# 	cells.loc <- object@cells.loc
# 	intr.valid <- object@lrscore[[slot.use]]@intr.list

# 	null.dge.lig <- object@imputation[[lig.slot]]@imp.data.null
# 	null.dge.recep <- object@imputation[[recep.slot]]@imp.data

# 	if (ncol(null.dge.lig) != ncol(null.dge.recep)){
# 		null.dge.recep <- null.dge.recep[, sample(ncol(null.dge.recep), 
# 							ncol(null.dge.lig), replace = T)]
# 	}

# 	null.lrscore.mtx <- inferScoreLR.dgCMatrix(null.dge.lig, null.dge.recep,
# 					intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()

# 	colnames(null.lrscore.mtx) <- colnames(object@lrscore[[slot.use]]@score)
# 	rownames(null.lrscore.mtx) <- paste0("perm_", 1:nrow(null.lrscore.mtx))

#     if (sum(Matrix::colSums(null.lrscore.mtx) == 0) != 0){
#         message("A total of ", sum(Matrix::colSums(null.lrscore.mtx) == 0), " intr are empty in NULL scores.")
#         null.lrscore.mtx = null.lrscore.mtx[, !Matrix::colSums(null.lrscore.mtx) == 0]
#     }

# 	intr.both <- intersect(colnames(null.lrscore.mtx), colnames(object@lrscore[[slot.use]]@score))

# 	if (length(intr.both) != ncol(null.lrscore.mtx)) {
# 		message("Removing ", ncol(null.lrscore.mtx) - length(intr.both), " more intr from NULL scores.")
# 		null.lrscore.mtx = null.lrscore.mtx[, intr.both]
# 	}

# 	if (length(intr.both) != ncol(object@lrscore[[slot.use]]@score)) {
# 		message("Removing ", ncol(object@lrscore[[slot.use]]@score) - length(intr.both), " corresponding intr from REAL scores.")
# 		object@lrscore[[slot.use]]@score = object@lrscore[[slot.use]]@score[, intr.both]
# 	}

# 	object@lrscore[[slot.use]]@score.null <- null.lrscore.mtx

#     return(object)
# }







# #' Permute Imputation Results of specific imputation method
# #' 
# #' This function is a follow-up function of imputeNiche, which permutes the default or user-defined
# #' imputation method and stored the results in the ImpData object.
# #' 
# #' @param object A Cytosignal object
# #' @param nn.type The imputation method to use
# #' @param perm.size Size of the permutation test
# #' 
# #' @return A Cytosignal object
# #' @export
# permuteListImpLR <- function(
# 	object,
# 	nn.type = NULL,
# 	perm.size = 100000
# ){
# 	if (is.null(nn.type)){
# 		nn.type <- object@imputation[["default"]]
# 	}

# 	if (!nn.type %in% names(object@imputation)){
# 		stop("Cannot find corresponding imputation method.")
# 	}

# 	use.gau <- grepl("GauEps", nn.type)
# 	use.dt <- grepl("DT", nn.type)
# 	use.raw <- grepl("Raw", nn.type)

# 	# if (!use.gau & !use.dt){
# 	# 	stop("Cannot find corresponding imputation method.")
# 	# }

# 	message("Permuting imputation data on method ", nn.type, "...")

# 	dimnames.list <- dimnames(object@imputation[[nn.type]]@imp.data)
# 	cells.loc <- object@cells.loc
# 	dge.raw <- object@counts

# 	# decide how many times to permute according to the size of the data
# 	if (!is.numeric(perm.size)){
# 		stop("perm.size must be a numeric value.")
# 	}

# 	if (perm.size < ncol(dge.raw)){
# 		message("Permutation size too small, using ", ncol(dge.raw), " instead.")
# 		perm.size <- ncol(dge.raw)
# 	}
	
# 	times <- ceiling(perm.size / ncol(dge.raw))
# 	each.size <- ceiling(perm.size / times)

# 	message("Permuting whole dataset ", times, " times...")

# 	# permute cell locations and remember to change the rownames!!
# 	null.cells.loc = cells.loc[sample(nrow(cells.loc)), ]
# 	rownames(null.cells.loc) = rownames(cells.loc)

# 	# use the same imputation graph for permutation
# 	null.graph <- object@imputation[[nn.type]]@nn.graph

# 	cat("Permuting NULL expressions...\nTimes No.")

# 	null.dge.list <- lapply(1:times, function(i){
# 		## shuffule dge, remember to change the colnames!!
# 		## could also shuffle the null graph, mathematically the same
# 		cat(i, ", ", sep = "")
# 		perm.index <- sample(ncol(dge.raw))
# 		null.dge.raw <- dge.raw[, perm.index]
# 		scale.fac <- object@parameters[["lib.size"]][perm.index]
# 		# scale.fac <- object@parameters[["lib.size"]]
		
# 		colnames(null.dge.raw) <- colnames(dge.raw)
# 		names(scale.fac) <- colnames(dge.raw)

# 		scale.fac <- Matrix::Matrix(scale.fac, nrow = 1, byrow = T, sparse = T) # computing scale factor
# 		# null.dge.raw <- changeUniprot.matrix_like(null.dge.raw.all, object@intr.valid[["gene_to_uniprot"]])[[1]]

# 		if (use.raw) {
# 			# sample the dge.raw to control the size of the permutation
# 			sample.index <- sample(ncol(null.dge.raw), each.size)

# 			null.dge.eps <- null.dge.raw[, sample.index]
# 			scale.fac.imp <- scale.fac[sample.index]
# 		} else {
# 			# try to permute the values of each column
# 			null.graph <- shuffle_sp_mat_col(null.graph)
			
# 			# sample the null graph to control the size of the permutation
# 			null.graph <- null.graph[, sample(ncol(null.graph), each.size)]

# 			null.dge.eps <- null.dge.raw %*% null.graph
# 			scale.fac.imp <- scale.fac %*% null.graph
# 		}
		
# 		null.dge.eps <- normCounts.list(list(mat=null.dge.eps, scale.fac=as.numeric(scale.fac.imp)), "default")
# 		rm(null.dge.raw)

# 		return(null.dge.eps)
# 	})

# 	null.dge.imp <- cbind_list(null.dge.list)
# 	rm(null.dge.list); gc()

# 	cat("End.\nFinished!\n")

# 	object@imputation[[nn.type]]@imp.data.null <- null.dge.imp

# 	return(object)

# }



# #' Permute LR score for specific ligand-receptor imputation obj pairs
# #' 
# #' This function is a follow-up function of inferScoreLR. It computes the NUL LRscores 
# #' using the NULL imputation results and stores the results in the LR score object.
# #' The null distribution of the LR scores can be used to test the significance of the LR scores.
# #' 
# #' @param object A Cytosignal object
# #' @param slot.use The imputation method to use
# #' 
# #' @return A Cytosignal object
# #' @export
# permuteListScoreLR <- function(
# 	object,
# 	slot.use = NULL
# ){
# 	if (is.null(slot.use)){
# 		slot.use <- object@lrscore[["default"]]
# 	}

# 	message("Permuting scores on Score slot: ", slot.use, "...")

# 	# score.obj <- object@lrscore[[slot.use]]
# 	lig.slot <- object@lrscore[[slot.use]]@lig.slot
# 	recep.slot <- object@lrscore[[slot.use]]@recep.slot
# 	cells.loc <- object@cells.loc
# 	intr.valid <- object@lrscore[[slot.use]]@intr.list

# 	null.dge.eps <- object@imputation[[lig.slot]]@imp.data.null
# 	null.dge.dt <- object@imputation[[recep.slot]]@imp.data.null

# 	null.lrscore.mtx <- inferScoreLR.dgCMatrix(null.dge.eps, null.dge.dt,
# 					intr.valid[["ligands"]], intr.valid[["receptors"]]); gc()

# 	colnames(null.lrscore.mtx) <- colnames(object@lrscore[[slot.use]]@score)
# 	rownames(null.lrscore.mtx) <- paste0("perm_", 1:nrow(null.lrscore.mtx))

#     if (sum(Matrix::colSums(null.lrscore.mtx) == 0) != 0){
#         message("A total of ", sum(Matrix::colSums(null.lrscore.mtx) == 0), " intr are empty in NULL scores.")
#         null.lrscore.mtx = null.lrscore.mtx[, !Matrix::colSums(null.lrscore.mtx) == 0]
#     }

# 	intr.both <- intersect(colnames(null.lrscore.mtx), colnames(object@lrscore[[slot.use]]@score))

# 	if (length(intr.both) != ncol(null.lrscore.mtx)) {
# 		message("Removing ", ncol(null.lrscore.mtx) - length(intr.both), " more intr from NULL scores.")
# 		null.lrscore.mtx = null.lrscore.mtx[, intr.both]
# 	}

# 	if (length(intr.both) != ncol(object@lrscore[[slot.use]]@score)) {
# 		message("Removing ", ncol(object@lrscore[[slot.use]]@score) - length(intr.both), " corresponding intr from REAL scores.")
# 		object@lrscore[[slot.use]]@score = object@lrscore[[slot.use]]@score[, intr.both]
# 	}

# 	object@lrscore[[slot.use]]@score.null <- null.lrscore.mtx

#     return(object)
# }




#' Sub function for graphNicheLR, input a CytoSignal object
#' 
#' @param object A Cytosignal object
#' @param lig.slot The ligand slot to use
#' @param recep.slot The receptor slot to use
#' @param intr.db.name The intr database name to use
#' @param nn.use slot that the neighbor index should be taken from, by default is the same as
#' 			the recep.slot. For example, if score.obj = GauEps-DT, then nn.use = "DT".
#' 			nn.use could also be a user-defind factor.
#' 
#' @return A Cytosignal object
#' @export
graphNicheLR.CytoSignal <- function(
	object,
	lig.slot,
	recep.slot,
	intr.db.name,
	nn.use = NULL
){
	if (!lig.slot %in% names(object@imputation)){
		stop("Ligand slot not found.")
	}

	if (!recep.slot %in% names(object@imputation)){
		stop("Receptor slot not found.")
	}

	if (!intr.db.name %in% c("diff_dep", "cont_dep")) {
		stop("intr.db.name must be either 'diff_dep' or 'cont_dep'.")
	}

	if (is.null(nn.use)) {
		nn.use <- recep.slot
	}
	
	if (is.character(nn.use)) {
		if (!nn.use %in% names(object@imputation)){
			stop("Imputation slot not found.")
		}
		nb.id.fac <- object@imputation[[nn.use]]@nn.id
	} else if (is.factor(nn.use)) {
		if (length(nn.use) != ncol(object@imputation[[lig.slot]]@imp.data))
			stop("nn.use must have the same length as the number of cells.")
		nb.id.fac <- nn.use
	} else {
		stop("nn.use must be either a factor or a character.")
	}

	dge.lig <- object@imputation[[lig.slot]]@imp.data
	dge.recep <- object@imputation[[recep.slot]]@imp.data

	if (!all.equal(dim(dge.lig), dim(dge.recep))){
		stop("dge.lig and dge.recep must have the same dimension.")
	}

	if (ncol(dge.lig) < length(levels(nb.id.fac))){
		stop("Number of index beads larger than the number of beads in DGE.")
	}

	# if (!all.equal(object@intr.valid[["symbols"]][[lig.slot]], 
	# 				object@intr.valid[["symbols"]][[recep.slot]])){
	# 	stop("Unpt symbols generated from imputations not equal!")
	# }

	message("Computing scores using ", intr.db.name, " database.")
	message("Ligand: ", lig.slot, ", Receptor: ", recep.slot, ".")

	intr.db.list <- checkIntr(unname(object@intr.valid[["symbols"]][["intr"]]), 
							object@intr.valid[[intr.db.name]])

	res.mtx <- graphNicheLR.dgCMatrix(dge.lig, dge.recep, nb.id.fac, 
				intr.db.list[["ligands"]], intr.db.list[["receptors"]])

	score.obj <- methods::new(
		"lrScores",
		lig.slot = lig.slot,
		recep.slot = recep.slot,
		intr.slot = intr.db.name,
		intr.list = intr.db.list,
		score = res.mtx,
		score.null = methods::new("matrix"),
		res.list = list(),
		log = list(
			"Parmaeters:" = c(lig.slot, recep.slot),
			"Date:" = Sys.time()
		)
	)

	tag <- paste0(lig.slot, "-", recep.slot)
	object@lrscore[["default"]] <- tag
	object@lrscore[[tag]] <- score.obj

	return(object)
}


#' Compute the LR score for specific ligand-receptor imputation obj pairs
#' 
#' @param object A Cytosignal object
#' @param lig.slot The ligand slot to use
#' @param recep.slot The receptor slot to use
#' @param intr.db.name The intr database name to use
#' @param nn.use The neighbor index as niche
#' 
#' @return A Cytosignal object
#' @export
#' 
graphNicheLR <- function(
  object,
  ...
) {
  UseMethod(generic = 'graphNicheLR', object = object)
}

#' Sub function for graphNicheLR, input a sparse matrix
#' 
#' @param dge.lig A sparse matrix for ligand
#' @param dge.recep A sparse matrix for receptor
#' @param nb.id.fac A factor of neighbor indices
#' @param lig.fac A factor of ligand indices
#' @param recep.fac A factor of receptor indices
#' 
#' @return A sparse matrix
#' @export
#' 
graphNicheLR.dgCMatrix <- function(
	dge.lig,
    dge.recep,
    nb.id.fac,
	lig.fac,
	recep.fac
){

	### cavaet: remember to convert the Uniprot ids to indices!
	# convert nb fac
	nb.id.fac = sort(nb.id.fac)
	nb.index = facToIndex(nb.id.fac)

	# lig.fac = intr.valid[["ligands"]]
	# recep.fac = intr.valid[["receptors"]]

	if (max(as.integer(names(lig.fac))) > nrow(dge.lig)){
		stop("Intr index out of dge bounds.")
	}

	lig.index = facToIndex(lig.fac)
	recep.index = facToIndex(recep.fac)

	# compute scores
	res.mtx = graphNicheLR_cpp(
		unname(as.matrix(dge.lig)),
		unname(as.matrix(dge.recep)),
		lig.index[[1]], lig.index[[2]], 
		recep.index[[1]], recep.index[[2]],
		nb.index[[1]], nb.index[[2]]
	)

	dimnames(res.mtx) = list(colnames(dge.lig)[as.integer(levels(nb.id.fac))], levels(lig.fac))
	# dimnames(res.mtx) = list(colnames(dge.lig), levels(lig.fac))
	# res.mtx = Matrix(res.mtx, sparse = T)

	return(res.mtx)
}


#' Sub function for graphNicheLR, input a CytoSignal object
#' 
#' @param object A Cytosignal object
#' @param lig.slot The ligand slot to use
#' @param recep.slot The receptor slot to use
#' @param intr.db.name The intr database name to use
#' @param nn.use slot that the neighbor index should be taken from, by default is the same as
#' 			the recep.slot. For example, if score.obj = GauEps-DT, then nn.use = "DT".
#' 			nn.use could also be a user-defind factor.
#' 
#' @return A Cytosignal object
#' @export
graphNicheLR.CytoSignal <- function(
	object,
	lig.slot,
	recep.slot,
	intr.db.name,
	nn.use = NULL
){
	if (!lig.slot %in% names(object@imputation)){
		stop("Ligand slot not found.")
	}

	if (!recep.slot %in% names(object@imputation)){
		stop("Receptor slot not found.")
	}

	if (!intr.db.name %in% c("diff_dep", "cont_dep")) {
		stop("intr.db.name must be either 'diff_dep' or 'cont_dep'.")
	}

	dge.lig <- object@imputation[[lig.slot]]@imp.data
	dge.recep <- object@imputation[[recep.slot]]@imp.data

	if (!all.equal(dim(dge.lig), dim(dge.recep))){
		stop("dge.lig and dge.recep must have the same dimension.")
	}

	message("Computing scores using ", intr.db.name, " database.")
	message("Ligand: ", lig.slot, ", Receptor: ", recep.slot, ".")

	intr.db.list <- checkIntr(unname(object@intr.valid[["symbols"]][["intr"]]), 
							object@intr.valid[[intr.db.name]])

	# res.mtx <- graphNicheLR.dgCMatrix(dge.lig, dge.recep, nb.id.fac, 
	# 			intr.db.list[["ligands"]], intr.db.list[["receptors"]])

	res.mtx <- dge.lig * dge.recep

	score.obj <- methods::new(
		"lrScores",
		lig.slot = lig.slot,
		recep.slot = recep.slot,
		intr.slot = intr.db.name,
		intr.list = intr.db.list,
		score = t(res.mtx),
		score.null = methods::new("matrix"),
		res.list = list(),
		log = list(
			"Parmaeters:" = c(lig.slot, recep.slot),
			"Date:" = Sys.time()
		)
	)

	tag <- paste0(lig.slot, "-", recep.slot)
	object@lrscore[["default"]] <- tag
	object@lrscore[[tag]] <- score.obj

	return(object)
}



# #' Sub function for changeUniprot, input a cytosignal object
# #' 
# #' @param object A Cytosignal object
# #' @param slot.use The slot to use for subsetting
# #' @param verbose Whether to print out the progress
# #' 
# #' @return A Cytosignal object
# #' @export
# changeUniprot.CytoSignal <- function(
# 	object,
# 	slot.use = NULL, 
# 	# mode = "unpt",
# 	verbose = T
# ){
# 	if (is.null(slot.use)){
# 		slot.use <- object@imputation[["default"]]
# 	}

# 	if (!slot.use %in% names(object@imputation)){
# 		stop("Data not found.")
# 	}

# 	gene_to_uniprot <- object@intr.valid[["gene_to_uniprot"]]
# 	dge.raw <- object@counts
# 	# check duplicate items exists in database

# 	unpt.list <- changeUniprot.matrix_like(dge.raw, gene_to_uniprot, verbose = verbose)
	
# 	object@imputation[[slot.use]]@intr.data <- unpt.list[[1]]
# 	object@intr.valid[["symbols"]][[slot.use]] <- unpt.list[[2]]

# 	scale.fac <- object@imputation[[slot.use]]@scale.fac$raw
# 	scale.fac <- scale.fac[colnames(unpt.list[[1]])]
# 	object@imputation[[slot.use]]@scale.fac$raw <- scale.fac

# 	# if (!all.equal(names(slot(object, slot.use)), names(object@intr.data))) {
# 	# 	warning("Names of imp.norm.data and intr.data not equal!")
# 	# }

# 	return(object)
# }
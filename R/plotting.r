plotSignif <- function(object, num.plot, res_dir, plot.details = T, slot.use = NULL, signif.use = NULL, plot.clusters = T,
                    plot.velo = F, colors.list = NULL, pt.size=0.5, pt.stroke = 0.2, u_width = 6, u_hgt = 5, set.res = 200,
                    return.plot = F
){
    if (length(object@lrscore) == 0){
        stop("No score matrix found!\n")
    }

    if (is.null(slot.use)){
        slot.use = object@lrscore[["default"]]
    }

    score.obj = object@lrscore[[slot.use]]
    intr.db <- object@intr.valid[[score.obj@intr.slot]]

    if (is.null(signif.use)){
        signif.use = "result.hq.pear"
    }

    if (!signif.use %in% names(score.obj@res.list)){
        stop("No such significance level found!\n")
    }

    lig.slot <- score.obj@lig.slot
    recep.slot <- score.obj@recep.slot
    dge.lig <- object@imputation[[lig.slot]]@intr.data
    dge.recep <- object@imputation[[recep.slot]]@intr.data
    null.dge.gau <- object@imputation[[lig.slot]]@intr.data.null
    null.dge.dt <- object@imputation[[recep.slot]]@intr.data.null
    dge.raw <- changeUniprot.dgCMatrix(object@counts, object@intr.valid[["gene_to_uniprot"]])[[1]]
    clusters <- object@clusters

    score.mtx <- score.obj@score
    null.lrscore.mtx <- score.obj@score.null

    res.list <- score.obj@res.list[[signif.use]]
    res.intr.list = names(res.list)
    cells.loc = as.data.frame(object@cells.loc)
    cells.loc = cells.loc[colnames(dge.lig), ]

    if (is.null(colors.list)){
        levels.use <- levels(clusters)
        colors.list = as.character(paletteer::paletteer_d("ggsci::default_igv",
                        n = length(levels.use)))
        names(colors.list) = levels.use
    }

    if (nrow(cells.loc) != ncol(dge.lig)){
        stop("Incorrect beads positions!\n")
    }

    cat("Now processing: ")

    plots.list <- lapply(num.plot, function(i){ # for each intr in the res.list
        cat("No.", i, ", ", sep = "")
        sub.df = cells.loc
        null.sub.df = cells.loc

        pair.index = which(object@intr.valid$intr.index$id_cp_interaction == res.intr.list[i])

        ligand.name = object@intr.valid$intr.index[pair.index, 4]
        if (ligand.name == ""){ # if complex
            ligand.name = object@intr.valid$intr.index[pair.index, 2]
        }
        ligand.name = gsub("_HUMAN", "", ligand.name)

        receptor.name = object@intr.valid$intr.index[pair.index, 5]
        if (receptor.name == ""){
            receptor.name = object@intr.valid$intr.index[pair.index, 3]
        }
        receptor.name = gsub("_HUMAN", "", receptor.name)

        # generate interaction names
        # intr.name = paste0(gsub("_HUMAN", "", ligand.name), "-", gsub("_HUMAN", "", receptor.name))
        intr.name = paste0(ligand.name, "-", receptor.name)

        colnames(sub.df) = c("x", "y")
        colnames(null.sub.df) = c("x", "y")

        ligands = dge.lig[names( intr.db[[2]][intr.db[[2]] == res.intr.list[i]] ), ]
        ligands.ori = dge.raw[names(intr.db[[2]][intr.db[[2]] == res.intr.list[i]]), ]
        ligands.null = null.dge.gau[names(intr.db[[2]][intr.db[[2]] == res.intr.list[i]]), ]

        receps = dge.recep[names(intr.db[[3]][intr.db[[3]] == res.intr.list[i]]), ]
        receps.ori = dge.raw[names(intr.db[[3]][intr.db[[3]] == res.intr.list[i]]), ]
        receps.null = null.dge.dt[names(intr.db[[3]][intr.db[[3]] == res.intr.list[i]]), ]

        if (!is.null(nrow(ligands))){ligands = colSums(ligands)}
        if (!is.null(nrow(receps))){receps = colSums(receps)}
        if (!is.null(nrow(ligands.ori))){ligands.ori = colSums(ligands.ori)}
        if (!is.null(nrow(ligands.null))){ligands.null = colSums(ligands.null)}
        if (!is.null(nrow(receps.ori))){receps.ori = colSums(receps.ori)}
        if (!is.null(nrow(receps.null))){receps.null = colSums(receps.null)}

        sub.df$ligands = ligands
        sub.df$ligands_ori = ligands.ori
        sub.df$receptors = receps
        sub.df$receptors_ori = receps.ori
        null.sub.df$ligands_null = ligands.null
        null.sub.df$receptors_null = receps.null

        p.lig = ggplot(sub.df, aes(x, y, color = ligands))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Ligand \n", ligand.name), x = NULL, y = NULL)

        p.lig.ori = ggplot(sub.df, aes(x, y, color = ligands_ori))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Ligand-ori \n", ligand.name), x = NULL, y = NULL)

        p.lig.null = ggplot(null.sub.df, aes(x, y, color = ligands_null))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Ligand-null \n", ligand.name), x = NULL, y = NULL)

        p.recep = ggplot(sub.df, aes(x, y, color = receptors))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Receptor \n", receptor.name), x = NULL, y = NULL)

        p.recep.ori = ggplot(sub.df, aes(x, y, color = receptors_ori))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Receptor-ori \n", receptor.name), x = NULL, y = NULL)

        p.recep.null = ggplot(null.sub.df, aes(x, y, color = receptors_null))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(color = paste0("Receptor-null \n", receptor.name), x = NULL, y = NULL)

        p.list = list(p.lig, p.lig.ori, p.lig.null, p.recep, p.recep.ori, p.recep.null)
        names(p.list) = c("lig", "lig.ori", "lig.null", "recep", "recep.ori", "recep.null")

        ##### plot scores

        if (ncol(dge.lig) != nrow(score.mtx)){
            cat("Setting LR-scores of cells with no neighbors to empty.\n")
            sub.df$scores = 0
            sub.df[rownames(score.mtx), "scores"] = score.mtx[, res.intr.list[i]]

            null.sub.df$null_scores = 0
            null.sub.df[rownames(null.lrscore.mtx), "null_scores"] = null.lrscore.mtx[, res.intr.list[i]]
        } else{
            sub.df$scores = score.mtx[, res.intr.list[i]]
            null.sub.df$null_scores = null.lrscore.mtx[, res.intr.list[i]]

        }

        p.score = ggplot(sub.df, aes(x, y, color = scores))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(x = NULL, y = NULL)

        p.list[["score"]] = p.score

        p.score.null = ggplot(null.sub.df, aes(x, y, color = null_scores))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
            theme_cowplot(12)+
            labs(x = NULL, y = NULL)

        p.list[["score.null"]] = p.score.null

        ##### plot significant beads

        sub.df$group = "beads"
        # sub.df$group[res.list[[i]]] = "significant"
        sub.df[res.list[[i]], "group"] = "significant"

        p.sig = ggplot(sub.df, aes(x, y, color = group))+
            geom_point(size = pt.size, stroke = pt.stroke)+
            # geom_segment(aes(xend = xend, yend = yend), data = edges.df, color = "red", size = 0.2)+
            scale_colour_manual(values = c("beads" = "grey", "significant" = "red"))+
            theme_cowplot(12)+
            labs(x = NULL, y = NULL)

        p.list[["sig"]] = p.sig

        if (plot.clusters){
            clusters = clusters[rownames(sub.df)]
            sub.df$clusters = as.character(clusters)

            sig.sub.df = sub.df[sub.df$group == "significant", ]

            p.clust = ggplot(sub.df, aes(x, y, color = clusters))+
                geom_point(size = pt.size, stroke = pt.stroke)+
                # scale_colour_viridis_d(direction = -1)+
                # ggsci::scale_color_nejm()+
                # ggsci::scale_color_igv(palette = "default")+
                scale_colour_discrete(type = colors.list, breaks = levels.use) +
                geom_point(data = sig.sub.df, aes(x, y), color = "red", size = pt.size*3, alpha = 0.15)+
                theme_cowplot(12)+
                guides(color = guide_legend(override.aes = list(size = 4)))+
                labs(x = NULL, y = NULL)

            p.list[["clust"]] = p.clust
        }

        if (plot.velo){
            sub.df$velo = velo.mtx[, res.intr.list[i]]

            p.velo = ggplot(sub.df, aes(x, y, color = velo))+
                geom_point(size = pt.size, stroke = pt.stroke)+
                scale_color_viridis_c(option = "plasma", direction = -1, na.value = '#F5F5F5')+
                theme_cowplot(12)+
                labs(x = NULL, y = NULL)

            p.list[["velo"]] = p.velo
        }

        if (return.plot){ return(p.list) }

        plot.len = length(p.list)
        # cat(paste0("Length of plot list: ", plot.len, "\n"))
        intr.name = paste0(intr.name, "-n_", length(res.list[[i]]))

        if (plot.details){
            p.list <- p.list[c(1, 4, 7, 9, 3, 2, 5, 8, 10, 6)]
            png(paste0(res_dir, "/Rank_", i, "_", intr.name, ".png"),
                        res = set.res, width = u_width*4, height = u_hgt*2, unit="in")
            print(cowplot::plot_grid(plotlist = p.list, nrow = 2))
            dev.off()

        } else{
            p.list <- p.list[c(1, 4, 7, 9, 2, 5, 8, 10)]
            png(paste0(res_dir, "/Rank_", i, "_", intr.name, ".png"),
                        res = set.res, width = u_width*4, height = u_hgt*2, unit="in")
            print(cowplot::plot_grid(plotlist = p.list, nrow = 2))
            dev.off()
        }

    })

    cat("End.\nFinished!\n")

    if (return.plot){return(plots.list)}
}




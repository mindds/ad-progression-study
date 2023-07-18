#' subcluster filter: identify donor specific sub-clusters
#' donor specific cluster:
#' * ncell/mean(ncell in that subcluster) > 5 in more than 2 
#' * ncell/mean(ncell in that subcluster) > 10 in one 
#' @param seurat.obj seurat object
#' @param res numeric; resolution to use
find_donor_subc = function(seurat.obj, res){
  cat("resolution",res, "...\n")
  #nsubjects = length(unique(seurat.obj$Donor.ID))
  subc_ncell = table(seurat.obj[[paste0("integrated_snn_res.",res)]][,1], seurat.obj$Donor.ID)
  subc_pctcell = subc_ncell/rowMeans(subc_ncell)
  # check:
  
  is_donor_subc = rowSums(subc_pctcell > 5) == 2 | rowSums(subc_pctcell > 10)
  return(list(
    sel.clust = as.numeric(names(is_donor_subc[!is_donor_subc])),
    remove.clust = as.numeric(names(is_donor_subc[is_donor_subc]))
  ))
}


#' remove donor specific clusters
#' @param seurat.obj seurat object;
#' @param res numeric; resolution
removd_donor_clust = function(seurat.obj, res){
  Idents(seurat.obj) = paste0("integrated_snn_res.", res)
  ## donor specific clusters
  remove_c = find_donor_subc(seurat.obj, res)
  subc_seled = !(Idents(seurat.obj) %in% remove_c$remove.clust)
  obj_removed = seurat.obj[, subc_seled]
  
  DefaultAssay(obj_removed) = "RNA"
  Idents(obj_removed) = paste0("integrated_snn_res.", res)
  cat("remove donor specific subcluster:\n")
  print(remove_c$remove.clust)
  
  return(obj_removed)
}


#' find all markers
#' @param seurat.obj seurat object;
#' @param res numeric; resolution
#' @param ncore numeric; number of core to use
#' @param known_marker_genes vector; list of known marker genes
find_all_markers = function(seurat.obj, res, ncore=1, known_marker_genes){
  DefaultAssay(seurat.obj) = "RNA"
  Idents(seurat.obj) = paste0("integrated_snn_res.", res)
  
  cat("identify subcluster markers\n")
  if(ncore!=1){
    suppressMessages({
      plan("multiprocess", workers = ncore)
      options(future.globals.maxSize= 891289600)
      options(future.seed = TRUE)
      plan()
    })
  }
  ### Warning message:
  ### glm.fit: fitted probabilities numerically 0 or 1 occurred
  ### https://www.biostars.org/p/404577/
  markers = FindAllMarkers(seurat.obj, test.use = "LR",
                           latent.vars = c("Donor.ID"))
  plan("sequential")
  
  return(markers)
}

#' Create heatmap for given cluster resolution and set of cell type markers
#' @param seurat.obj seurat object;
#' @param brain_region string; name of brain region (eg BA20)
#' @param res numeric; cluster resolution
#' @param plot_dir string; directory to output heatmap image
#' @param marker_dir string;
#' @param known_marker_df dataframe; known marker genes and cell types
known_marker_heatmap = function(seurat.obj, data_name, brain_region, 
                                res, output_image, known_marker_df){
  DefaultAssay(seurat.obj) = "RNA"
  Idents(seurat.obj) = paste0("integrated_snn_res.", res)
  
  ## remove donor-specific subclusters
  ## defined as those with 10x greater number of cells than donor average
  seurat.obj = removd_donor_clust(seurat.obj, res)
  
  ## heatmaps
  DefaultAssay(seurat.obj) = "RNA"
  cat("scale data...\n")
  known_marker_df = known_marker_df[known_marker_df$genes %in% intersect(known_marker_df$genes, rownames(seurat.obj)), ]
  known_marker_genes = known_marker_df$genes
  seurat.obj = ScaleData(seurat.obj, features = known_marker_genes)
  
  cat("heatmap based on 500 downsampled data...\n")
  plot_seurat.obj = seurat.obj[, WhichCells(seurat.obj,downsample = 500)]
  rm(seurat.obj)
  
  ### heatmap of known marker genes
  do_heatmap_chunk(plot_seurat.obj, known_marker_df$genes, known_marker_df$group)
  n_genes = length(known_marker_genes)
  n_subclust = length(unique(Idents(plot_seurat.obj)))
  ggsave(output_image,
         width = n_subclust*2, height = (n_genes * 0.1 + length(unique(known_marker_df$group))*1), 
         limitsize = FALSE, units="in", dpi=300)
}


#' anaysis resolution by marker
#' @param ctype_obj seurat object;
#' @param res numeric; resolution
#' @param plot_dir string;
#' @param marker_dir string;
#' @param known_marker
evaluate_markers = function(ctype_obj, data_name, brain_region, 
                            res, plot_dir, marker_dir, known_marker){
  DefaultAssay(ctype_obj) = "RNA"
  Idents(ctype_obj) = paste0("integrated_snn_res.", res)
  
  ## remove donor-specific subclusters
  ## defined as those with 10x greater number of cells than donor average
  ctype_obj_removed = removd_donor_clust(ctype_obj, res)
  rm(ctype_obj)
  
  ## identify subcluster markers
  marker_file = paste0(marker_dir, sprintf("/%s_%s_subcluster-res%s-LR.regress.DonorID-removed.donorsclust.csv",
                                           brain_region, data_name, res))
  if(file.exists(marker_file)){
    cat("load markers from preprocessed file: ", marker_file, "\n")
    markers = fread(marker_file)
  }else{
    cat("identify subcluster markers\n")
    markers = find_all_markers(ctype_obj_removed, res, ncore)
    
    ## save markers
    setDT(markers)
    markers[, resolution := res]
    fwrite(markers, marker_file)
    cat("markers saved...\n")
  }
  
  ## top20 markers per cluster
  markers[, cluster := as.numeric(cluster)]
  top20_marker = markers[p_val_adj < 0.05, ] %>%
    .[order(cluster, -avg_log2FC, p_val_adj), head(.SD, 20), by = .(cluster)]
  cat("number of top20 markers", nrow(top20_marker), "\n")
  
  ## heatmaps
  DefaultAssay(ctype_obj_removed) = "RNA"
  cat("scale data...\n")
  known_marker = known_marker[genes %in% intersect(known_marker$genes, rownames(ctype_obj_removed)), ]
  known_marker_genes = known_marker$genes
  plot_genes = c(top20_marker$gene, known_marker_genes)
  ctype_obj_removed = ScaleData(ctype_obj_removed, features = plot_genes)
  
  cat("heatmap based on 500 downsampled data...\n")
  plot_ctype_obj_removed = ctype_obj_removed[, WhichCells(ctype_obj_removed,downsample = 500)]
  rm(ctype_obj_removed)
  
  ### heatmaps based on top 20 genes
  cat("heatmaps on top20 genes...\n")
  #### warnings: because of heatmap color
  #### Scale for 'fill' is already present. Adding another scale for 'fill',
  #### which will replace the existing scale.
  do_heatmap_chunk(plot_ctype_obj_removed, top20_marker$gene, top20_marker$cluster)
  n_subclust = length(unique(Idents(plot_ctype_obj_removed)))
  n_genes = nrow(top20_marker)
  ggsave(paste0(plot_dir, sprintf("/%s_heatmap-subcluster-%s-LR-top20.genes-removed.donorclust.png",
                                  data_name, res)),
         width = n_subclust*2, height = (n_genes * 0.1 + n_subclust*0.5), limitsize = FALSE, units="in", dpi=300)
  
  ### heatmap of known marker genes
  do_heatmap_chunk(plot_ctype_obj_removed, known_marker$genes, known_marker$group)
  n_genes = length(known_marker_genes)
  ggsave(paste0(plot_dir, sprintf("/%s_heatmap-subcluster-%s-LR-known.marker.genes-removed.donorclust.png",
                                  data_name, res)),
         width = n_subclust*2, height = (n_genes * 0.1 + length(unique(known_marker$group))*1), 
         limitsize = FALSE, units="in", dpi=300)
}


#' identify B group subclusters
#' @param seurat.obj
#' @param res
#' @param fc
B_subc = function(seurat.obj, res, fc){
  ind_meta = unique(seurat.obj[, c("Donor.ID", "B")])
  nBgroup_donor = table(ind_meta$B)
  nBgroup_pct = nBgroup_donor/sum(nBgroup_donor)
  n_table = table(seurat.obj[[paste0("integrated_snn_res.", res)]], seurat.obj$B)
  return(list(B1 = which(n_table[, 1]) > nBgroup_pct[1] * fc,
              B2 = which(n_table[, 2]) > nBgroup_pct[2] * fc,
              B3 = which(n_table[, 3]) > nBgroup_pct[3] * fc))
}


#' return number of expressed genes in seurat obj
#' @param obj seurat object;
#' @param min_pct express genes = count > 0 in at least min_pct cells
expr_in_obj = function(obj, min_pct){
  n_cell = ncol(obj)
  # genes express in how many cells
  expr_in_donor = rowSums(obj@assays$RNA@counts > 0)  > 0.01 * n_cell
  return(list(expr_genes = rownames(expr_in_donor),
              n_expr_genes = sum(expr_in_donor),
              n_cell = n_cell))
}

#' return HVFinfo in seurat obj
#' @param obj
get_HVFinfo = function(obj, donor){
  dt = HVFInfo(obj)
  dt["Donor.ID"] = donor
  dt["B"] = as.character(unique(obj$B))
  dt["is_HVG_in_donor"] = rownames(dt) %in% obj@assays$RNA@var.features
  dt = df_to_dt(dt)
  return(dt)
}


#' chi-square test
#' @param n_table matrix; proportion of Braak stage
#' [B, subc)]
chisq_test = function(n_table){
  chisq = chisq.test(n_table)
  obs_gt_exp = ifelse(chisq$observed > chisq$expected, "greater", "less")
  return(obs_gt_exp)
}


#' binominal test
#' @param n_table matrix; proportion of Braak stage 
#' [B, subc)]
binomial_test = function(n_table){
  p_B = rowSums(n_table)
  p_B = p_B/sum(p_B)
  vals = sapply(names(p_B), function(b){
    ## count table
    n_B = rbind(n_table[b, ], colSums(n_table))
    rownames(n_B) = c("n_B", "n_total")
    binom = sapply(colnames(n_B), function(subc){
      binom = binom.test(n_B["n_B", subc], n_B["n_total", subc], p = p_B[b], alternative = "greater")
      p_binom = binom$p.value
      p_fc = round(binom$estimate / p_B[b], digits = 2)
      return(c(p_binom, p_fc))
    })
    return(list(p = binom[1,], fc = binom[2,]))
  })
  pvals = do.call(rbind, vals["p", ])
  fcs = do.call(rbind, vals["fc", ])
  return(list(p = pvals, fc = fcs))
}


#' correlation analysis
#' @param 
#' 
cor_test = function(x, y){
  corr = cor.test(x=x, y=y, method = 'spearman')
  return(corr$pval)
}

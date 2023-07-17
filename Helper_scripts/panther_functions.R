
#' Adapt from Annie's script
#' @param jason_f string; path to the panther json file
#' @param result_f string; path to the result file
#' @param query_size numeric; length of PANTHER tested genes
#' @param effective_domain_size numeric; length of reference gene list; default = 20589 
load_panther_jason = function(
  json_f, result_f, query_size, genes, 
  effective_domain_size = 20589
){
  if (!file.exists(result_f)) {
    require(clusterProfiler)
    message("read json from:", json_f)
    panther_data = fromJSON(file=json_f)
    
    res_list = list()
    results_parsed = unlist(panther_data$overrepresentation$group,recursive=FALSE)
    
    # load GO data
    go_f = "/autofs/space/mindds_001/projects/AbbvieSnRNASeq/data/Rdata/annie/All_GO_Data.Rdata"
    if (!file.exists(go_f)) {
      require(org.Hs.eg.db)
      go_data <- get_GO_data(org.Hs.eg.db, "BP", "SYMBOL")
      #save(go_data, file=go_f)
    } else {
      load(go_f)
    }
    
    for (i in 1:length(results_parsed)) {
      result_group = results_parsed[[i]]
      list_depth = purrr::vec_depth(result_group)
      
      if (list_depth==3) {
        result <- result_group
        pathway_name <- result$term$label
        pathway_id <- result$term$id
        FDR <- result$input_list$fdr
        pValue <- result$input_list$pValue
        plus_minus <- result$input_list$plus_minus
        overlapping_genes <- vector(mode="character")
        
        if (pathway_name != "UNCLASSIFIED") {
          
          all_genes_in_pathway <- go_data$PATHID2EXTID[[pathway_id]]
          geneID <- paste0(overlapping_genes, collapse="/")
          
          term_size <- length(all_genes_in_pathway)
          intersection_size <- length(overlapping_genes)
          
          this_res_df <- data.frame(Description = pathway_name,
                                    ID = pathway_id,
                                    pValue = pValue,
                                    qvalue = FDR,
                                    geneID = geneID,
                                    GeneRatio = paste0(intersection_size, "/", query_size),
                                    BgRatio = paste0(term_size, "/", effective_domain_size),
                                    GeneSetSize = term_size,
                                    Count = intersection_size)
          
          if (plus_minus=="+") {
            res_list <- rlist::list.append(res_list, this_res_df)
          }
        }
        
      }
      else if (list_depth==5) {
        for (j in 1:length(result_group)) {
          result <- result_group[[j]]
          pathway_name <- result$term$label
          pathway_id <- result$term$id
          FDR <- result$input_list$fdr
          pValue <- result$input_list$pValue
          plus_minus <- result$input_list$plus_minus
          
          overlapping_genes <- result$input_list$mapped_id_list[[1]]
          
          if (pathway_name != "UNCLASSIFIED") {
            all_genes_in_pathway <- go_data$PATHID2EXTID[[pathway_id]]
            geneID <- paste0(overlapping_genes, collapse="/")
            
            term_size <- length(all_genes_in_pathway)
            intersection_size <- length(overlapping_genes)
            
            this_res_df <- data.frame(Description = pathway_name,
                                      ID = pathway_id,
                                      pValue = pValue,
                                      qvalue = FDR,
                                      geneID = geneID,
                                      GeneRatio = paste0(intersection_size, "/", query_size),
                                      BgRatio = paste0(term_size, "/", effective_domain_size),
                                      GeneSetSize = term_size,
                                      Count = intersection_size)
            
            if (plus_minus=="+") {
              res_list <- rlist::list.append(res_list, this_res_df)
            }
          }
        } 
      }
      else if (any(str_detect(result_group, "mapped_id_list"))) {
        result <- result_group
        pathway_name <- result$term$label
        pathway_id <- result$term$id
        FDR <- result$input_list$fdr
        pValue <- result$input_list$pValue
        plus_minus <- result$input_list$plus_minus
        overlapping_genes <- result$input_list$mapped_id_list[[1]]
        
        if (pathway_name != "UNCLASSIFIED") {
          all_genes_in_pathway <- go_data$PATHID2EXTID[[pathway_id]]
          geneID <- paste0(overlapping_genes, collapse="/")
          
          term_size <- length(all_genes_in_pathway)
          intersection_size <- length(overlapping_genes)
          
          this_res_df <- data.frame(Description = pathway_name,
                                    ID = pathway_id,
                                    pValue = pValue,
                                    qvalue = FDR,
                                    geneID = geneID,
                                    GeneRatio = paste0(intersection_size, "/", query_size),
                                    BgRatio = paste0(term_size, "/", effective_domain_size),
                                    GeneSetSize = term_size,
                                    Count = intersection_size)
          
          if (plus_minus=="+") {
            res_list <- rlist::list.append(res_list, this_res_df)
          }
        } 
      }
      else {
        for (j in 1:length(result_group)) {
          result <- result_group[[j]]
          pathway_name <- result$term$label
          pathway_id <- result$term$id
          FDR <- result$input_list$fdr
          pValue <- result$input_list$pValue
          plus_minus <- result$input_list$plus_minus
          overlapping_genes <- vector(mode="character")
          
          if (pathway_name != "UNCLASSIFIED") {
            all_genes_in_pathway <- go_data$PATHID2EXTID[[pathway_id]]
            geneID <- paste0(overlapping_genes, collapse="/")
            
            term_size <- length(all_genes_in_pathway)
            intersection_size <- length(overlapping_genes)
            
            this_res_df <- data.frame(Description = pathway_name,
                                      ID = pathway_id,
                                      pValue = pValue,
                                      qvalue = FDR,
                                      geneID = geneID,
                                      GeneRatio = paste0(intersection_size, "/", query_size),
                                      BgRatio = paste0(term_size, "/", effective_domain_size),
                                      GeneSetSize = term_size,
                                      Count = intersection_size)
            
            if (plus_minus=="+") {
              region_res_list <- rlist::list.append(region_res_list, this_res_df)
            }
          } 
        }
      }
    }
    
    DEG_enriched <- do.call(plyr::rbind.fill, res_list)
    DEG_enriched_sig <- DEG_enriched %>%
      # filter(GeneSetSize <= 500 & GeneSetSize >= 15) %>%
      filter(GeneSetSize <= 10000 & GeneSetSize >= 0) %>%
      mutate(p.adjust = p.adjust(pValue, method="BH")) %>%
      filter(pValue < 0.05, p.adjust < 0.1, qvalue < 0.1) %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pValue, p.adjust, qvalue, geneID, Count)
    
    custom_df <- DEG_enriched_sig
    
    DEG_enrichResult <- new("enrichResult", result=custom_df,
                                        pvalueCutoff=0.1, pAdjustMethod="BH",
                                        qvalueCutoff=0.1, organism="Homo sapiens",
                                        ontology="BP", gene=genes,
                                        keytype="SYMBOL", geneSets=go_data$PATHID2EXTID,
                                        readable=T)
    
    DEG_enrichResult_res <- DEG_enrichResult@result
    
    r = list(
      DEG_enriched = DEG_enriched, 
      DEG_enriched_sig = DEG_enriched_sig,
      DEG_enrichResult = DEG_enrichResult, 
      DEG_enrichResult_res = DEG_enrichResult_res
    )
    
    saveRDS(r, file=result_f)
    message("Rds file saved at:", result_f)
  } else {
    message("load data from:", result_f)
    r = readRDS(result_f)
  }
  return(r)
}


hclust = function(
  enriched_simple, 
  showCategory = 30,
  color = "p.adjust",
  nCluster = 12,
  cex_category = 1,
  label_format = 50,
  fontsize = 3.8,
  offset = 0.5,
  offset_tiplab = 0.5,
  hclust_method = "ward.D2",
  group_color = NULL,
  extend = 0.3,
  hilight = TRUE, 
  hexpand = .1,
  align = "both"
){
  n <- showCategory
  keep <- seq_len(n)
  
  ## Fill the upper triangular matrix completely
  termsim2 <- fill_termsim(enriched_simple, keep)
  
  # Write this similarity matrix to a csv (Table 13)
  # write.csv(termsim2, paste0(annie_path, "results/manuscript_tables/Table13_JC_Matrix.csv"), row.names = T)
  
  ## Use the ward.D method to avoid overlapping ancestor nodes of each group
  hc <- stats::hclust(stats::as.dist(1- termsim2),
                      method = "ward.D2")
  # clus <- stats::cutree(hc, nCluster)
  clus <- stats::cutree(hc, h=1)
  d <- data.frame(label = names(clus),
                  color = enriched_simple@result[keep, as.character(color)],
                  count = enriched_simple@result$Count[keep]) %>%
    mutate(count_fac = case_when(count <= 10 ~ "<= 10",
                                 count > 11 & count <= 25 ~ "11-25",
                                 count > 26 & count <= 50 ~ "26-50",
                                 count > 50 & count <= 100 ~ "51-100",
                                 count > 100 & count <= 200 ~ "101-200",
                                 T ~ "> 200"
    )) %>%
    mutate(count_fac = factor(count_fac, levels=c("<= 10",
                                                  "11-25", 
                                                  "26-50",
                                                  "51-100",
                                                  "101-200",
                                                  "> 200")))
  
  my_plot <- group_tree(hc = hc, clus = clus, d = d, offset_tiplab = offset_tiplab, 
                        nWords = nWords, label_format = label_format, offset = offset, 
                        fontsize = fontsize, group_color = group_color, extend = extend, 
                        hilight = hilight, cex_category = cex_category, align = align,
                        wrap_length=60) +
    theme_cowplot() +
    ggnewscale::new_scale_colour() +
    geom_tippoint(aes(color = -1*log10(color), size = count_fac)) +
    scale_colour_continuous(low="blue", high="red", name = "-log10 padj",
                            limits=c(1, 4),
                            breaks=c(1:4)) +
    scale_size_discrete(name = "# Genes", drop=F)
  
  # Using the cowplot package
  my_plot_legend <- my_plot +
    theme(legend.position="bottom",
          legend.direction = "horizontal") + 
    guides(color = guide_colorbar(order = 2, barwidth=15, barheight=1),
           size = guide_legend(nrow=1, byrow=T))
  legend <- cowplot::get_legend(my_plot_legend)
  
  # grid.newpage()
  # png(paste0(plot_dir, "Figure3E_legend.png"), width=10.5, height=1, units="in", res=300)
  # grid.draw(legend)
  # dev.off()
  
  
  my_plot = my_plot + 
    ggtree::hexpand(ratio = 2) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position="none")
  return(my_plot)
}


source("/autofs/space/mindds_001/projects/AbbvieSnRNASeq/scripts/annie/scripts/manuscript_figures/Figure4_Helper_Functions/enrichplot_custom_funcs.R")

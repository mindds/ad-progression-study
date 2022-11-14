get_ave_exp = function(obj, data_dir, sel_genes){
  cat(systime(), "averaging expression\n")
  DefaultAssay(obj) = "RNA"
  if(is.null(sel_genes)) obj = obj[sel_genes, ]
  
  obj$sample_ids = paste(obj$Donor.ID, obj$Region, sep = "_")
  Idents(obj) = "sample_ids"
  
  ave.exp = AverageExpression(
    obj, assays = "RNA", slot = "data")$RNA
  
  ave.exp_p = ave.exp %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as.data.table() %>%
    melt(., id.vars = "gene", variable.name = "sample_id", value.name = "ave.exp")
  
  
  ind_meta = obj@meta.data[, c("Donor.ID", "Region", "Ptau.Total.Tau..A.U..", "Path..Group.", "A", "B", "C", "Sex", "Age")] %>%
    unique() %>%
    setDT()
  ind_meta[, sample_id := paste(Donor.ID, Region, sep = "_")]
  
  ## load abeta
  abeta = fread(file.path(data_dir, "Hyman_AD_CTR_samples_Halo_data_3D6quant_plot.csv"))
  abeta = abeta %>%
    .[!is.na(ADRC), ] %>%
    .[, abeta_3D6 := mean(percent_DAB_positive_tissue, na.rm = T), by = .(ADRC, Region)]
  abeta[, Region := gsub("BA46/9", "BA46", Region)]
  abeta$Region %>% unique()
  abeta[, sample_id := paste(ADRC, Region, sep = "_")]
  
  ind_meta[abeta[, .(sample_id, abeta_3D6)] %>% unique(),
           on = .(sample_id),
           abeta_3D6 := i.abeta_3D6]
  
  ave.exp_p[ind_meta, on = .(sample_id), ptau := i.Ptau.Total.Tau..A.U..]
  ave.exp_p[ind_meta, on = .(sample_id), donor_id := i.Donor.ID]
  ave.exp_p[ind_meta, on = .(sample_id), path_group := i.Path..Group.]
  ave.exp_p[ind_meta, on = .(sample_id), A := i.A]
  ave.exp_p[ind_meta, on = .(sample_id), B := i.B]
  ave.exp_p[ind_meta, on = .(sample_id), C := i.C]
  ave.exp_p[ind_meta, on = .(sample_id), age := i.Age]
  ave.exp_p[ind_meta, on = .(sample_id), sex := i.Sex]
  ave.exp_p[ind_meta, on = .(sample_id), abeta_3D6 := i.abeta_3D6]
  ave.exp_p[ind_meta, on = .(sample_id), region := i.Region]
  
  ave.exp_p[, region := factor(region, levels = c("EC", "BA20", "BA46", "V2", "V1"))]
  ave.exp_p[, log2.ave.exp := log2(ave.exp + 0.01)]
  
  # cell expression levels
  cat(systime(), "\tgene expressed in cell if count > 0\n")
  cat(systime(), "\tcalculating pct of cells express given gene\n")
  cell_pct = data.table(
    gene = rownames(obj@assays$RNA@counts),
    cell_exp_pct = rowMeans(obj@assays$RNA@counts > 0)
  )
  ave.exp_p[cell_pct, on = .(gene), cell_exp_pct := i.cell_exp_pct]
  return(ave.exp_p)
}

# lmer model based on given 
library(nlme)
lmer_gene = function(g, ave.exp_p, f, coef_vars = "ptau", use_lmer = T, return_summary = F){
  if(!g %in% ave.exp_p$gene){
    return(cat(systime(), "\t", g, ":skiped because gene is low expressed\n"))
  }
  avexp = ave.exp_p[gene == g,]
  if(use_lmer){
    m1 = lmer(as.formula(f), avexp)
  }else{
    m1 = lm(as.formula(f), avexp)
  }
  
  v = summary(m1)
  if(return_summary){
    return(v)
  }else{
    if(coef_vars == "all"){
      coef_vars = rownames(v$coefficients)
    }
    r = data.table(
      primerid = g,
      test = "sample level validation",
      f = f,
      levels = coef_vars,
      coef = v$coefficients[coef_vars, "Estimate"],
      p = v$coefficients[coef_vars, "Pr(>|t|)"]
    )
    return(r)
  }
  
}

lmer_genes = function(gs, ave.exp_p, pct_cut, f, coef_vars = "ptau", use_lmer = T){
  cat(systime(), "sample level validation\n")
  cat(systime(), "\tformular:", f,"\n")
  r = pblapply(gs, function(g){
    lmer_gene(g, ave.exp_p[cell_exp_pct > pct_cut, ], f, coef_vars, use_lmer)
  })
  dt = Reduce(rbind, r)
  dt[, adjp := p.adjust(p, method = "BH"), by = .(levels)]
  
  # direction adj.pvalue
  dt[, direction_adjp := "ns"]
  dt[adjp < 0.05 & coef > 0, direction_adjp := "up"]
  dt[adjp < 0.05 & coef < 0, direction_adjp := "down"]
  
  # direction pvalue
  dt[, direction_p := "ns"]
  dt[p < 0.05 & coef > 0, direction_p := "up"]
  dt[p < 0.05 & coef < 0, direction_p := "down"]
  cat(systime(), "done\n")
  return(dt)
}


venn_smith = function(lmer_dt, smith_dt){
  smith_up_gs = smith_dt[padj < 0.05 & logFC > 0.25, gene]
  smith_down_gs = smith_dt[padj < 0.05 & logFC < -0.25, gene]
  up_gs = lmer_dt[direction == "up", primerid]
  down_gs = lmer_dt[direction == "down", primerid]
  up_venn = ggVennDiagram(
    list(model = up_gs, smith = smith_up_gs),
    label_alpha=0,
    label = "count"
  ) + 
    scale_fill_gradient(low="white",high = "#FF1818") +
    labs(title = "pos genes")
  down_venn = ggVennDiagram(
    list(model = down_gs, smith = smith_down_gs),
    label_alpha=0,
    label = "count"
  ) + 
    scale_fill_gradient(low="white",high = "#4D96FF") +
    labs(title = "neg genes")
  venn = ggarrange(
    up_venn,
    down_venn,
    nrow = 1,
    ncol = 2
  )
  return(venn)
}


p_mast_g = function(g, ave.exp_p, cov = "pTau/Tau", show_text = F){
  cov = ifelse(cov == "pTau/Tau", "ptau", "ptau_scaled")
  avexp = ave.exp_p[gene == g,]
  
  p = ggplot(
    data = avexp,
    aes_string(x = cov, y = "log2.ave.exp", color = "region",
               donor = "donor_id", sex = "sex", age = "age", 
               A = "A", B = "B", C = "C", 
               path_group = "path_group", 
               pTau = "ptau",
               abeta_3D6 = "abeta_3D6")
  ) +
    labs(x = "pTau/Tau", y = "Average expression",
         title = g) +
    geom_line(aes(group = donor_id), color = "gray") +
    geom_point() +
    ggthemes::theme_few(base_size = 14)
  
  if(show_text){
    p = p + geom_text(aes(label = region))
  }
  return(p)
}

p_density_filter_genes = function(dt, pct_cut){
  p = ggplot(
    dt,
    aes(x = log2.ave.exp, group = "sample_id")
  ) + 
    geom_density() +
    ggthemes::theme_few() +
    labs(title = "all genes") +
    ggplot(
      dt[cell_exp_pct > pct_cut, ],
      aes(x = log2.ave.exp, group = "sample_id")
    ) + 
    geom_density() +
    ggthemes::theme_few() +
    labs(title = "filtered genes")
  return(p)
}


p_cov = function(g, dt, cov = "path_group"){
  p = ggplot(
    dt[gene == g, ],
    aes_string(x = cov, 
               y = "log2.ave.exp",
               fill = cov)
  ) +
    geom_boxplot() +
    geom_point() +
    facet_grid(.~ region) +
    labs(x = "Path group", title = g, y = "log2 average expression") +
    ggthemes::theme_few()
  return(p)
}

p_lm = function(g, dt, cov = "ptau"){
  ggplot(
    dt[gene == g, ],
    aes_string(x = cov, 
               y = "log2.ave.exp")
  ) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) + 
    facet_grid(.~ region) +
    labs(x = cov, title = g) +
    theme_bw()
}


p_volcano = function(
  dt, choose_level, level_col = "levels",
  abslogFC_cut = 0, p_cut = 0.05, pcol = "p"
){
  dt2 = dt[get(level_col) == choose_level,]
  top10_up = dt2[direction == "up", ] %>%
    .[order(coef, decreasing = T), primerid] %>%
    head(10)
  top10_down = dt2[direction == "down", ] %>%
    .[order(coef, decreasing = T), primerid] %>%
    tail(10)
  dt2[primerid %in% c(top10_down, top10_up), label_genes := primerid]
  dt2[, y := get(pcol)]
  ggplot(
    dt2,
    aes(x = coef, y = -log10(y), color = direction)
  ) +
    geom_point() +
    geom_text_repel(aes(label = label_genes),
                    colour = "black") +
    geom_vline(aes(xintercept = abslogFC_cut), linetype="dashed") +
    geom_vline(aes(xintercept = -abslogFC_cut), linetype="dashed") +
    geom_hline(aes(yintercept = -log10(p_cut)), linetype="dashed") +
    scale_color_manual(
      values = c(
        "up" = "#FF1818",
        "down" = "#3DB2FF",
        "ns" = "gray"
      )
    ) +
    ggthemes::theme_few()
}


view_tT = function(dt, test_genes = NULL){
  if(!is.null(test_genes)){
    dt = dt[primerid %in% test_genes, ]
  }
  DT::datatable(dt)
}


lm_logfc_heatmap = function(sub_logfc_dt, row_split = NULL){
  sub_logfc_mtx = sub_logfc_dt[, .(primerid, path_group2, path_group3, path_group4)] %>%
    as.data.frame() %>%
    column_to_rownames("primerid") %>% 
    as.matrix()
  if(is.null(row_split)) row_split = sub_logfc_dt$cluster
  
  p = Heatmap(
    matrix = sub_logfc_mtx,
    row_split = row_split,
    cluster_columns = F,
    cluster_rows = F,
    show_row_names = F,
    name = "logFC",
    row_title_rot = 0
  )
  return(p)
}

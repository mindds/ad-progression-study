
#-------------------------------------------------------------------------------
# custom functions used to perform cell level analysis using MAST
# 
# History:
#   2022-04-29 Zhaozhi_Li initial

# reference:
#   1. https://github.com/RGLab/MAST/blob/master/vignettes/MAST-intro.pdf
#.  2. seurat: https://github.com/satijalab/seurat/blob/f1b2593ea72f2e6d6b16470dc7b9e9b645179923/R/differential_expression.R
#-------------------------------------------------------------------------------


#' perform cell level analysis given formula and sca object using MAST
#' @param sca sca object
#' @param f string; formula
#' @param use_genes vector; a vector of gene names to be selected
#' @return sca object
mast_zlm = function(sca, f, use_glmer = T){
  # reference from smith et al. paper:
  # zlm(~ histopath_marker + (1|sample)
  #     + cngeneson + pc_mito + sex + brain_region, sca,
  #     method = "glmer", ebayes = F).
  cat(systime(), "runing mast model\n")
  cat(systime(), "\tformula:", f, "\n")
  # reference: http://rglab.github.io/MAST/articles/MAST-Intro.html
  if(use_glmer){
    de = zlm(as.formula(f), sca, method = "glmer", ebayes = F)
  }else{
    de = zlm(as.formula(f), sca)
  }
  
  cat(systime(), "\tmodel done")
  return(de)
}


#' convert seurat object to sca object (for MAST)
#' @param obj seurat object
#' @param pct_cut numeric; if not null then it will be used as cut off select 
#' expressed genes. expressed genes: expressed in at leat pct_cut of the cells
#' @param use_genes vector; a vector of gene names to be selected
#' @return sca object
seurat_to_sca = function(obj, pct_cut = NULL, use_genes = NULL){
  DefaultAssay(obj) = "RNA"
  cdata = obj@meta.data
  cdata = cbind(
    data.frame(wellKey = rownames(cdata)),
    cdata
  )
  fdata = data.table(primerid = rownames(obj))
  rownames(fdata) = rownames(obj)
  
  cat(systime(), "normalization\n")
  obj = NormalizeData(obj)
  counts_n = as(obj@assays$RNA@data, "dgTMatrix")
  counts_n = as.matrix(counts_n)
  
  cat(systime(), "convert seurat object to sca object\n")
  sca = FromMatrix(
    exprsArray = counts_n,
    cData = cdata,
    fData = fdata
  )
  
  if(!is.null(pct_cut)){
    cat(systime(), "filter low expressed genes\n")
    cat(systime(), "\texpressed genes: expressed in at least", pct_cut * 100, "% of the cells\n")
    sca = MAST::filterLowExpressedGenes(sca, pct_cut)
  }
  
  if(!is.null(use_genes)){
    cat(systime(), "use selected genes\n")
    sca = sca[use_genes, ]
  }

  return(sca)
}


#' summarize MAST results
#' @param zlm_de mast object; returned by zlm()
#' @param cont string
#' @return a data.table with{
#' primerid: gene
#' Pr(>Chisq): hurdle P values
#' coef: logFC coefficients
#' adj.pvalue: adjusted pvalues by BH
#' }
sum_zlm = function(zlm_de, cont = 'ptau'){
  sum_dt = summary(zlm_de, doLRT=cont)$datatable
  fcHurdle = merge(
    # hurdle P values
    sum_dt[contrast %in% cont & component=='H',
          .(contrast, primerid, `Pr(>Chisq)`)],
    # logFC coefficients
    sum_dt[contrast %in% cont & component=='logFC',
          .(contrast, primerid, coef, ci.hi, ci.lo)], 
    by=c('primerid', 'contrast')
  )
    
  fcHurdle[,adj.pvalue:=p.adjust(`Pr(>Chisq)`, method = "BH"), by = .(contrast)]
    
  return(list(dt = fcHurdle, summary = sum_dt))
}


#' plot density of percentage of expression of genes given seurat object
#' @param obj seurat object
#' @param pct_cut numeric; if not null then it will be used as cut off select 
#' expressed genes. expressed genes: expressed in at leat pct_cut of the cells
#' @returns {
#' p: a density plot (ggplot)
#' pct_exp_dt: a data.table/
#' gene: gene name
#' pct_exp: percentage of express
#' /
#' }
p_distribution_of_gene = function(obj, pct_cut){
  DefaultAssay(obj) = "RNA"
  pct_exp_dt = data.table(
    gene = rownames(obj),
    pct_exp = rowMeans(obj@assays$RNA@counts > 0)
  )
  p = ggplot(
    pct_exp_dt,
    aes(x = pct_exp)
  ) +
    geom_density() +
    geom_vline(aes(xintercept=pct_cut),
               color="red", linetype="dashed", size=1) +
    labs(
      x = "percentage of expressed cells (%)",
      title = "dirstribution of percentage expression of genes across all genes"
    ) +
    ggthemes::theme_few() + 
    scale_x_continuous(labels=scales::percent) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  return(list(p = p, pct_exp_dt = pct_exp_dt))
}



match_cov = function(obj, meta, cov){
  meta = meta[, c("Donor.ID", "Region", cov), with = F]
  meta[, Donor.ID := as.character(Donor.ID)]
  test = obj@meta.data[, c("Donor.ID", "Region", cov)]
  setDT(test)
  test_cov = meta[test, on = .(Donor.ID, Region), cov, with = F][[1]]
  return(test_cov)
}

prepare_variable = function(obj, has_abeta = F, fix_apoe = F){
  calc.cdr = function(obj) colMeans(obj@assays$RNA@counts > 0)
  obj$cngeneson = calc.cdr(obj)
  obj$ptau = obj$Ptau.Total.Tau..A.U..
  obj$ptau_scaled = scale(obj$Ptau.Total.Tau..A.U..)
  obj$Braak_stage = ifelse(obj$Braak %in% c("VI", "V"), "B3",
                    ifelse(obj$Braak %in% c("III", "IV"), "B2",
                    ifelse(obj$Braak %in% c("0", "I", "II"), "B1", NA)))
  obj$Braak_stage = factor(obj$Braak_stage, levels = c("B1", "B2", "B3"))
  obj$Cerad = factor(obj$C, levels = c(0, 1, 2, 3))
  obj$percent.mito_scaled = scale(obj$percent.mito)
  obj$ratio.mito = obj$percent.mito/100
  if(has_abeta){
    obj$abeta_ratio = obj$abeta_3D6/100
    obj$abeta_scaled = scale(obj$abeta_3D6)
  }
  
  if(fix_apoe){
    apoe_meta = obj@meta.data[!is.na(obj$APOE), c("APOE", "Donor.ID")] %>% 
      unique() %>%
      as.data.table()
    obj$APOE = apoe_meta[
      data.table(Donor.ID = obj$Donor.ID),
      on = .(Donor.ID),
      APOE
    ]
    obj$has_e4 = ifelse(grepl("4", obj$APOE), T, F)
  }
  
  return(obj)
  
}


get_abeta = function(obj, meta){
  meta[, Donor.ID := as.character(Donor.ID)]
  obj$abeta_3D6 = meta[data.table(
    Donor.ID = obj$Donor.ID,
    Region = obj$Region
  ), on = .(Donor.ID, Region), abeta_3D6]
  return(obj)
}


add_direction = function(dt, p_cut = 0.05, abs_logFC_cut = 0, use_adjp = F){
  dt[, direction := "ns"]
  if(use_adjp){
    dt[adj.pvalue < p_cut & coef > abs_logFC_cut, direction := "up"]
    dt[adj.pvalue < p_cut & coef < -abs_logFC_cut, direction := "down"]
  }else{
    dt[`Pr(>Chisq)` < p_cut & coef > abs_logFC_cut, direction := "up"]
    dt[`Pr(>Chisq)` < p_cut & coef < -abs_logFC_cut, direction := "down"]
  }
  return(dt)
}

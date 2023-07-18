#' reformat snRNA-seq counts: sum given celltype for individuals 
#' @param counts   dgTMatrix; shape(n_genes, n_samples)
#' @param projids  vector; vector of string: individual ids of each cell
ind_expr = function(counts, projids, projids.u) {
  ind_counts = matrix(0, nrow = nrow(counts), ncol = length(projids.u))
  rownames(ind_counts) = rownames(counts)
  colnames(ind_counts) = as.character(projids.u)
  for(projid in unique(projids)){
    sel = projids == projid
    if(sum(sel) == 0){
      ind_counts[, as.character(projid)] = 0
    }else{
      sub_counts = counts[, sel]
      if(is.vector(sub_counts)){
        ind_counts[, as.character(projid)] = sub_counts
      }else{
        ind_counts[, as.character(projid)] = rowSums(sub_counts)
      }
    }
  }
  return(ind_counts)
}


#' RUV: remove unwanted variation
#' @param ind.expr matrix; individual level gene expression; 
#' shape(n_gene, n_samples)
#' @param ruvn integer; the number of factors of unwanted variation to be estimated from the data
ind_RUV_expr = function(indExpr, use.factor = NULL, ruvn = 10){
  cat("number of unwanted variance:", ruvn, "\n")
  d_e = DGEList(indExpr, genes=rownames(indExpr))
  keep = rowSums(edgeR::cpm(d_e)>1) >= 3
  cat("number of genes keep: ", sum(keep), "\n")
  d_e = d_e[keep, , keep.lib.sizes=FALSE]
  d_e = calcNormFactors(d_e, method="TMM")
  if(length(use.factor) == 0){
    design = model.matrix(~B)
  }else{
    design = model.matrix(~ use.factor)
  }
  cat("design matrix:\n")
  print(design)
  d_e = estimateGLMCommonDisp(d_e, design)
  d_e = estimateGLMTagwiseDisp(d_e, design)
  fit1 = glmFit(d_e, design)
  res1 = residuals(fit1, type="deviance")
  ruv_cov = RUVr(round(d_e$counts), as.character(rownames(d_e$counts)), k=ruvn, res1)
  return(list(RUV = ruv_cov, Counts = d_e$counts))
}


#' calculate RUV
#' @param ctype.obj seurat.object; seurat object of a given cell type
#' @param ctype.meata data.table; meta data of given cell type
#' @param ind_meta data.table; meta data of individuals
#' @param projids.u vector; vector of strings; unqiue projids for all samples
ctype_ruv = function(ctype.obj, ctype_meta, ind_meta, projids.u, use.factor = NULL){
  # individual level expression
  ind_ctype_expr = ind_expr(counts = ctype.obj@assays$RNA@counts, 
                            projids = ctype.obj@meta.data$Donor.ID,
                            projids.u = projids.u)
  # RUV
  if(!all(colnames(ind_ctype_expr) == ind_meta$Donor.ID)){
    ind_ctype_expr = ind_ctype_expr[, as.character(ind_meta$Donor.ID)]
  }
  cat("identify RUVr other than", use.factor, "\n")
  ctype_ruv = ind_RUV_expr(ind_ctype_expr, 
                             use.factor = unlist(ind_meta[, use.factor, with = FALSE]))
  
  # glm 
  ctype_ruv_C = data.table(Donor.ID = unique(ctype_meta$Donor.ID))
  ctype_ruv_C = cbind(ctype_ruv_C, ctype_ruv$RUV$W)
  ctype_meta = ctype_meta[ctype_ruv_C, on = .(Donor.ID)]
  
  return(list(meta = ctype_meta,
               ruv = ctype_ruv))
}


#' poisson mixed model based on given gene & RUV value
#' @param counts matrix; raw counts of given cell type
#' @param genes vector; vector of gene symbol of interested genes
#' @param ctype_meta data.table; meta data of given celltype: cell level
ruv_re = function(counts, genes, ctype_meta, use.factors){
  genes = gsub("-", "_", genes)
  rownames(counts) = gsub("-", "_", rownames(counts))
  ctype_meta.cp = copy(ctype_meta)
  if(length(genes)>1){
    sub_counts = t(counts[genes, ])
    ctype_meta.cp = cbind(ctype_meta.cp, as.matrix(sub_counts))
    chunk_results = sapply(genes, function(gene) ruv_re_genes(gene, ctype_meta.cp, use.factors))
    chunk_results = df_to_dt(as.data.frame(t(chunk_results)))
  }else if(length(genes) == 1){
    ctype_meta.cp[, (genes) := counts[genes, ]]
    chunk_results = as.data.frame(as.list(ruv_re_genes(genes, ctype_meta.cp, use.factors = use.factors)))
    rownames(chunk_results) = genes
    chunk_results = df_to_dt(chunk_results)
  }else{
    stop("no gene input\n")
  }
  return(chunk_results)
}


#' poisson mixed model based on given gene & RUV value
#' @param genes vector; gene symbol of interested genes
#' @param ctype_meta.cp data.table; meta data of given celltype: cell level
ruv_re_genes = function(gene, ctype_meta.cp, nAGQ = 10, use.factors = NULL){
  formular = paste(gene, "~ ", use.factors, "+ W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9",
                     "+ W_10 + offset(log(nCount_RNA)) + (1 | Donor.ID)")
  
  ctype_re_Ruv = glmer(formula = as.formula(formular), 
                       family = poisson, data = ctype_meta.cp, nAGQ=nAGQ)
  gene_sum_re_Ruv = summary(ctype_re_Ruv)
  gene_sum = c(gene_sum_re_Ruv$coefficients[use.factors,])
  names(gene_sum) = paste0(use.factors, c("estimate", "stderr", "z", "pval"))
  return(gene_sum)
}


#' pmm model for given seurat objects
#' @param ctype.obj seurat object
#' @param meta data.table; meta data 
#' @param ctype_meta_file string; path to the cell level meta file with RUVs
#' @param use.factors vector; name of variable in the model
#' @param model_genes vector; list of genes for pmm model
#' @param ncore numeric; number of cores to use when parallel
#' @param chunk_length numeric; length of each chunk
pmm_model = function(ctype.obj, use.factors, ctype_meta_file, model_genes, 
                     processed_file, ncore = 4, 
                     chunk_length = 1){
  
  # individual level meta data
  ind_meta = ctype.obj@meta.data
  setDT(ind_meta)
  ind_meta = unique(ind_meta[, c("Donor.ID", use.factors), with = FALSE])
  
  # unique projids
  projids.u = unique(ctype.obj$Donor.ID)
  
  # RUV
  counts = ctype.obj@assays$RNA@counts
  if(!file.exists(ctype_meta_file)){
    cat("remove unwanted variance by RUVr...\n")
    ctype_meta = df_to_dt(ctype.obj@meta.data, rowcol = "barcode")
    ctype_meta = ctype_meta[, c("barcode", "Donor.ID", use.factors, "nCount_RNA"), with = FALSE]
    ruv = ctype_ruv(ctype.obj, ctype_meta, ind_meta, projids.u, use.factors)
    ctype_meta = ruv$meta
    fwrite(ctype_meta, ctype_meta_file)
  }else{
    cat("load previous calculated RUV from", ctype_meta_file, "...\n")
    ctype_meta = fread(ctype_meta_file)
    if("B" %in% use.factors){
      ctype_meta[, B := as.factor(B)]
      cat("Braak levels (in order):", levels(ctype_meta$B), "\n")
    }
  }
  
  rm(ctype.obj)
  
  # poisson mixed model
  cat("mixed linear model..\n")
  cat(length(model_genes), "genes will do.\n")
  
  # parallel
  chunks_genes = split(model_genes,         
                       ceiling(seq_along(model_genes) / chunk_length))
  
  # parallel
  plan("multiprocess", workers = ncore)
  options(future.globals.maxSize= 8912896000)
  plan()
  
  full_results = future_lapply(seq_along(chunks_genes), function(idx){
    tryCatch({
      genes = chunks_genes[[idx]]
      cat("process gene: ", genes)
      chunk_re_Ruv = ruv_re(counts, genes, ctype_meta, use.factors)
      chunk_re_Ruv[, genes := gsub("_", "-", genes)]
      fwrite(chunk_re_Ruv, file = processed_file, append = T)
    },
    error = function(e){
      cat("error genes: ")
      cat(genes, "\n")
      error_chunk_files = paste0("pmm_error_genes_", Sys.Date(),".txt")
      cat("store error genes in", error_chunk_files, "\n")
      write(genes, file = error_chunk_files, append = T)
    })
  }, future.seed = TRUE)
}

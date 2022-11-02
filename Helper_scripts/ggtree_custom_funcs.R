fill_termsim <- function(x, keep) {
  termsim <- x@termsim[keep, keep]
  termsim[which(is.na(termsim))] <- 0
  termsim2 <- termsim + t(termsim)
  for ( i in seq_len(nrow(termsim2)))
    termsim2[i, i] <- 1
  return(termsim2)
}

add_cladelab <- function(p, nWords, label_format, offset, roots, 
                         fontsize, group_color, cluster_color, 
                         pdata, extend, hilight, align) {
  
  n_color <- length(levels(cluster_color)) - length(cluster_color)
  if (is.null(group_color)) {
    color2 <- scales::hue_pal()(length(roots) + n_color)
    if (n_color > 0) color2 <- color2[-seq_len(n_color)]
  } else {
    color2 <- group_color
  }
  df <- data.frame(node = as.numeric(roots),
                   cluster=cluster_color,
                   color = color2
  )
  
  if (hilight) {
    p <- p + ggtree::geom_hilight(
      data = df,
      mapping = aes_(node =~ node, fill =~ cluster),
      show.legend = FALSE, 
      align = align) + 
      scale_fill_manual(values = df$color, 
                        guide = 'none')
    
  }
  
  return(p)
  
}

group_tree <- function(hc, clus, d, offset_tiplab, nWords, 
                       label_format, offset, fontsize, group_color, 
                       extend, hilight, cex_category, 
                       ID_Cluster_mat = NULL, geneClusterPanel = NULL,
                       align, wrap_length=25) {
  group <- NULL
  # cluster data
  dat <- data.frame(name = names(clus), cls=paste0("cluster_", as.numeric(clus)))
  grp <- apply(table(dat), 2, function(x) names(x[x == 1]))  
  p <- ggtree(hc, branch.length = "none", show.legend=FALSE)
  # extract the most recent common ancestor
  noids <- lapply(grp, function(x) unlist(lapply(x, function(i) ggtree::nodeid(p, i))))
  roots <- unlist(lapply(noids, function(x) ggtree::MRCA(p, x)))
  # cluster data
  p <- ggtree::groupOTU(p, grp, "group") + aes_(color =~ group)
  
  
  if (geneClusterPanel == "pie" || is.null(geneClusterPanel)) {
    ## 1.5 * max(radius_of_pie)
    offset_tiplab <- offset_tiplab * 1.5 * max(sqrt(d$count / sum(d$count) * cex_category))
  }  else if (geneClusterPanel == "heatMap") {
    ## Close to the width of the tree
    offset_tiplab <- offset_tiplab * 0.16 * ncol(ID_Cluster_mat) * max(p$data$x)
  } else if (geneClusterPanel == "dotplot") {
    ## Close to the width of the tree
    offset_tiplab <- offset_tiplab * 0.09 * ncol(ID_Cluster_mat) * max(p$data$x)
  }
  # max_nchar <- max(nchar(p$data$label), na.rm = TRUE)
  offset <- offset * (max(p$data$x) * 1.2 + offset_tiplab)    
  pdata <- data.frame(name = p$data$label, color2 = p$data$group)
  pdata <- pdata[!is.na(pdata$name), ]
  cluster_color <- unique(pdata$color2)
  n_color <- length(levels(cluster_color)) - length(cluster_color)
  if (!is.null(group_color)) {
    color2 <- c(rep("black", n_color), group_color)
    p <- p + scale_color_manual(values = color2, guide = 'none')
  }
  
  p <- p %<+% d +
    geom_tiplab(offset = offset_tiplab, hjust = 0, 
                aes(label = stringr::str_wrap(label, wrap_length)),
                show.legend = FALSE, align = TRUE, linesize = 0,
                size=fontsize,
                lineheight = .75)
  
  p <- add_cladelab(p = p, nWords = nWords, label_format = label_format, 
                    offset = offset, roots = roots, fontsize = fontsize, 
                    group_color = group_color, cluster_color = cluster_color, 
                    pdata = pdata, extend = extend, hilight = hilight, align = align) 
  return(p)
}


#-------------------------------------------------------------------------------
# ggplot themes for figures
#-------------------------------------------------------------------------------

#' theme with border
#' used for volcano plots 
my_border_theme = function(){
  theme(
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background =element_rect(fill="white"),
    axis.title.x = element_text(
      margin = margin(0.1, 0, 0, 0, unit = "cm"),
      size = 15
    ),
    axis.title.y = element_text(
      margin = margin(0, 0.1, 0, 0, unit = "cm"),
      angle =90, 
      size = 15
    ),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    legend.key = element_rect(fill = "white"),
    legend.background = element_rect(fill="white"),
    panel.border = element_rect(
      colour = "black", 
      fill=NA, 
      size=1
    ),
    plot.title = element_text(hjust = 0.5, size = 15)
  )
}

#' theme for UMAP plot
umap_theme = function(){
  theme(
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
    NoLegend()
}
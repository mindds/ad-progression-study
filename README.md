# AD Progression Study
A single-nucleus transcriptomics study of 629,755 astrocytes from five brain regions of 32 donors representing the hierarchical spatiotemporal progression of Alzheimer’s disease (AD).

## Integration of the data
The Integration.Rmd file has the scripts for integrating all the single cell RNA sequencing data from the 32 donors for the 5 brain regions (EC, ITG, PFC V1 and V2). The input data are 5 Seurat objects saved as a .Rdata file for each brain region. The script also performs clustering. The output is an integrated Seurat objects. 

## Annotation of the integrated data
Cell_type_annotation.Rmd script script was used to identify the cell types. Cell_type_annotation-figure1b.Rmd has the code for Figure 1b.

## Regional heterogeneity of astrocyte transcriptome
The regional_astrocyte_heterogeneity_normal_ageing_brain.Rmd script performs differential expression on rPCA integrated astrocytes data (input is a Seurat object). The script also has the code for visualizations of Figure 2 in the manuscript. 

## Spatio-temporal progression of Alzheimer's disease. 
The Spatial trajectories.Rmd contains the scripts that perform clustering of the 504 DEGs between any two “adjacent” nodes of the AD network from EC to V1. It also has the code to generate Fig 3 in the manuscript


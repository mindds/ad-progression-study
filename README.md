# AD Progression Study
A single-nucleus transcriptomics study of 629,755 astrocytes from five brain regions of 32 donors representing the hierarchical spatiotemporal progression of Alzheimer’s disease (AD).
Link to manuscript: 

Note:  
Data required to run the scripts is located in the ./Example_data directory of this github repository. If the data is not available, we will give the required data upon request. 


## Annotation of cell types
The script cell_type_annotation.Rmd has the required code to plot figure 1b in the manuscrip, which is the UMAP clustering of the cell types and violin plots to visualize the expression of marker genes of the cell types.

##### Input data:
A subset of the Seurat object for each brain region (eg: EC_subset) in the form .Rdata is given in the ./Example_data directory


## Regional heterogeneity of astrocyte transcriptome
The regional_astrocyte_heterogeneity.Rmd script performs differential expression on rPCA integrated astrocytes data and the correlation analysis between pTau/Tau and gene expression levels in pathology stage 1 astrocytes. The script also has the code for visualizations of Figure 2 in the manuscript. 

##### Input data:
rPCA integrated astrocytes as a Seurat object (.Rdata) for differential expression (data will be available upon request).


## Spatio-temporal progression of Alzheimer's disease
The spatial_trajectory.Rmd perform clustering of the 504 DEGs between any two “adjacent” nodes of the AD network from EC to V1. It also has the code to generate Fig 3 in the manuscriptThe temporal_trajectory.Rmd performs clustering the n=798 DEGs between any two “adjacent” pathology stages from early to end-stage. 

##### Input data:
Merged Seurat object (.Rdata) of all brain regions - differential expression between adjacent brain regions (spatial) or pathology group(temporal) and average expression for each gene (step 1). This data will be available upon request.
Results from previous step - clustering (step2). Data is available in the ./Example_data directory

Pathways analysis was done on gsea web tool https://www.gsea-msigdb.org/gsea/msigdb/annotate.jsp  and the .csv files were generated. Selected pathways were plotted in Figure 3b from the above files


## Correlation between different astrocyte subclusters
The astrocyte_integration_subclustering_and_DE.Rmd has a example script for our EC brain region data. It was used to perform CCA integration, subclustering of the astrocytes and to find the marker genes of each astrocyte subcluster. 

##### Input data: 
Seurat objects of each brain region. Data will be available upon request.

The gene_marker_data.Rmd processes the marker genes results of each brain region, to be used as input for the correlation analysis of astrocyte subclusters.

##### Input data: 
Results from running astrocyte_integration_subclustering_and_DE.Rmd. The data is available in the ./Example_data directory. 

The correlation_analysis_ast_subclusters.Rmd performs spectral clustering mapping astrocytes subclusters across brain regions based on the strength of the correlation of their gene expression results in 10 main astrocyte clusters. The results are visualized as a heatmap in figure 5b

##### Input data: 
Gene expression results of all brain regions. (results from running gene_marker_data.Rmd)  The data is available in the ./Example_data directory.

## AB plaque and ptau correlation

Results of Ab plaque load (% immunoreactive area fraction) and pTau/Tau ratio (measured by ELISA) in adjacent samples to those used for snRNA-seq across brain regions and pathology stages

##### Input data: 
Metadata from the rPCA integrated Suerat objects. The metadata is available in the ./Example_data directory. 





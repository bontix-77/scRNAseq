# scRNAseq
This is my first single-cell RNA-seq (scRNA-seq) analysis using the Seurat package in R. It's a learning project where I explore how to process, analyze, and visualize scRNA-seq data step by step.
This script is a comprehensive collection of capabilities provided by the Seurat library.
The result are presented as a HTML knit of a RMarkdown file (also available in the "Code" folder).

## The three is composed of 3 folders: Data, Results and Code.

#### 📂 What’s Inside the script:

**Basic QC** (filtering by features, counts, mitochondrial content)
 **Normalization** with SCTransform
 **Dimensionality reduction** (PCA, UMAP)
 **Clustering** of cells
 **Marker gene identification** Markers are calculated using specific features or between all the clusters.
 **Visualizations** like UMAP, violin plots, feature plots among others
 **Anottation of cell type** This has been performed using the markers and litterature or by using references profiles, cell by cell using SingleR and celldex.




## 📦 Requirements

This project was done in R. Main packages used:

```r
library(Seurat) # Seurat toolkit
install.packages("tidyverse") # For some graphs in ggplot2 an dyplr
library(hdf5r) # For data import
library(patchwork) # For plotting
library(glmGamPoi) # For SCTtransform using poisson gamma distribution
library(presto)
library(SingleR) # For cell type annotation; Bioconductor
library(celldex) # For cell type annotation reference; Bioconductor
library(MAST) # For differential expression; Bioconductor

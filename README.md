# scRNAseq using Seurat
This is my first single-cell RNA-seq (scRNA-seq) analysis using the Seurat package in R. It's a learning project where I explore how to process, analyze, and visualize scRNA-seq data step by step.
This script is a comprehensive collection of capabilities provided by the Seurat library.<br>
The result are presented as a [HTML](https://github.com/bontix-77/scRNAseq/blob/f5051972b5a61c77fd46fed45d6fde9530551114/Results/final-Seurat-pipeline.zip) knit of a [RMarkdown](https://github.com/bontix-77/scRNAseq/blob/80790a4a9c137109b1e5bd1726b6ffef2d15d276/Code/Seurat%20RMarkdown.Rmd) file (also available in the "Code" folder). The HTML is a .zip and due to size need to be downloaded.

## The tree is composed of 2 folders: Data and Results.

Database used for this demonstration can be found [here](https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/GettingStarted_scRNASeq.zip).


#### 📂 What’s Inside the [script](https://github.com/bontix-77/scRNAseq/blob/73dbdc6cad8a683d4c84116f58426034b654f272/Code/final%20Seurat%20script.R):

 -**Basic QC** (filtering by features, counts, mitochondrial content).<br>
 -**Normalization** with SCTransform.<br>
 -**Dimensionality reduction** (PCA, UMAP).<br>
 -**Clustering** of cells using knn.<br>
 -**Marker gene identification** Markers are calculated using specific features or between all the clusters.<br>
 -**Visualizations** like UMAP, violin plots, feature plots among others.<br>
 -**Anottation of cell type** This has been performed using the markers and litterature or by using references profiles, cell by cell using SingleR and celldex.<br>




## 📦 Requirements

This project was done in R. Main packages used:

```r
library(Seurat) # Seurat toolkit
library(tidyverse) # For some graphs in ggplot2 an dyplr
library(hdf5r) # For data import
library(patchwork) # For plotting
library(glmGamPoi) # For SCTtransform using poisson gamma distribution
library(presto)
library(SingleR) # For cell type annotation; Bioconductor
library(celldex) # For cell type annotation reference; Bioconductor
library(MAST) # For differential expression; Bioconductor
```

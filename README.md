# scRNAseq based on Seurat package

This is my independent single-cell RNA-seq (scRNA-seq) analysis using the Seurat package in R and a number of algorithms to prepare the data (SoupX, scrublet) and post processing such as CellChat to analyze cell-cell communication pathways, and enrichGO for the gene ontology pathway enrichment analysis. 
It's a learning repository where I show apossible flow for the  processing, analysis, and visualization of scRNA-seq data step by step. This repository is a comprehensive collection of state-of-the art tools in the field on scRNASeq.

## 📂 What’s Inside the repository

\-**Code folder** - containing all the algorithm each as module in its own sub-folder. here there is also a sub-directory called "python" containing a couple of scripts to download the GSM files from GEO and for the gene ID conversion for example from Ensembl ID to gene symbol.

\-**Results folder** - containing a series of results and graphs. The final result document is a knit of a RMarkdown in .html, consequentely it needs to be dowloaded in order to visualize it.

### The analysis has been performed on public data sets found at GEO, details can be found in the relative README file in the specif project folders (in construction)

Database used for the adipocite project can be found [here](https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/GettingStarted_scRNASeq.zip).

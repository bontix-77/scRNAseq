# nf-scRNA: Nextflow Pipeline for 10x scRNA-seq Processing (In Progress)

This repository hosts a **Nextflow** pipeline under active development for processing **10x Genomics single-cell RNA-seq** data.\
At this stage, the workflow includes two functional processes:

1.  **Reading and merging** multiple 10x Chromium output folders into a unified Seurat object.\
2.  **Filtering** low-quality cells based on customizable thresholds.

------------------------------------------------------------------------

## 🧩 Current Workflow Overview

```         
Code/ 
   └── nf-scRNA/ 
        ├── main.nf # Nextflow main workflow 
        ├── nextflow.config # Config parameters (paths, results) 
        ├── r-scripts/  
        │      ├── read_h5.R          # Read & merge 10x data  
        │      ├── cell_filter.R      # Filter the cell based on countings and mitocondrial RNAs
        │      ├── RDS_to_scanpy.R    # convert the seurat object to a scanpy object .h5ad
        │      ├── SCTransform.R      # Seurat normalization using negative binomila distribution
        │      └── seurat_PCA_UMAP.R  # Implementation of PCA y UMAP. Opcional batch effect                    │                               harmonization.
        ├── README 
        └── results/ # visualizations and images
```

### Processes

| Process | Description | Output |
|------------------|------------------------------------|------------------|
| **RUN_readh5** | Reads multiple 10x directories (e.g. GSM folders) and merges them into one Seurat object | no images |
| **RUN_cell_filter** | Filters cells by quality metrics (features, counts, mitochondrial content) | `5 .png analytics immages` |
| **RUN_SCTransform** | Normalize the count | no images |
| **RUN_seuratDisk** | Convert Seurat object .rds in a scampy object .h5ad |  |

------------------------------------------------------------------------

## ⚙️ Requirements

-   **Nextflow ≥ 23.10**
-   **R ≥ 4.3**
-   R packages:
    -   `Seurat`, `tidyverse`, `ggplot2`, `Matrix`

------------------------------------------------------------------------

📜 Status

🚧 Work in progress Current modules:

✅ RUN_readh5 — merge multiple 10x datasets

✅ RUN_cell_filter — apply QC-based cell filtering

✅ RUN_SCTransform — normalization and scaling using negative binomial

✅ RUN_seuratDisk — save Seurat object in h5Seurat file than convert to .h5ad to be transfered to a scanpy script

🔜 RUN_plot_QC, RUN_cluster, RUN_annotation (coming soon)

## 🧪 Running the Pipeline

\`\`\`bash git clone https://github.com/bontix-77/scRNAseq.git cd scRNAseq/Code/nf-scRNA nextflow run main.nf
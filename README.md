------------------------------------------------------------------------

# scRNAseq: A Modular Single-Cell RNA-seq Analysis Framework

**Author:** Alexander Bontempo **Repository:** [bontix-77/scRNAseq](https://github.com/bontix-77/scRNAseq) **Keywords:** scRNA-seq, Seurat, Nextflow, R, Python, Bioinformatics

------------------------------------------------------------------------

## 🧠 Overview

This repository provides a **complete and modular single-cell RNA-seq (scRNA-seq) analysis framework**, built to handle every step from raw data to biological interpretation. It combines **R-based analysis** (Seurat, SCTransform, QC, visualization) with **Nextflow orchestration** for reproducibility and scalability.

The project’s goal is to serve both as a **research-ready workflow** and a **learning platform** for integrating R, Python, and workflow languages in real scRNA-seq pipelines. The pipeline permits also to pull the T cell receptors diverisy into 4 main groups (alfa, beta, delta and gama) to prevent miss identification of the T cells population as discussed here :https://doi.org/10.1016/j.immuno.2025.100063. The option can be accessed switching parameter in the nextflow.config file.

|  |  |  |
|------------------------|------------------------|------------------------|
| ![](Code/nf-scRNA/results/cell_filter/nCount2_RNA_per_sample.png) | ![](Code/nf-scRNA/results/seurat_cellType/CellType_final.png) | ![](Code/nf-scRNA/results/seurat_PCA_UMAP/DimHeatmap.png) |

------------------------------------------------------------------------

## 🧬 Structure

```         
scRNAseq/
│
├── Code/
│   ├── nf-scRNA/             # Nextflow workflow and process definitions
│   │   ├── main.nf           # Main workflow (Seurat/Scanpy modes)
│   │   ├── nexflow.config    Parameter files and environment setup
│   │   ├── Docker_files/
│   │   │    ├── Seurat_base/        # folder containing all documents for the sc_rna:v1.0 docker image
│   │   │    │   ├── building.log         # complet log of the image building step
│   │   │    │   ├── dockerfile           # source dockerfile
│   │   │    │   ├── packages_list.csv    # complete list of the packages included in the image
│   │   │    │   └── README.md            # documentation
│   │   │    └── Seurat+celldex      # sc_rna:v1.1 expandiding v1.0 with seurat fix/v.5.3.1 celldex and
│   │   │        └── dockerfile        SingleR
│   │   │         
│   │   │
│   │   ├── r-scripts/          # Reusable process blocks (optional)
│   │   │     ├── read_h5.R         # Import and merge 10x datasets
│   │   │     ├── cell_filter.R     # QC and cell-level filtering
│   │   │     ├── SCTransform.R     # Normalization using Seurat SCTransform
│   │   │     ├── RDS_to_h5ad.R     # Conversion Seurat → h5Seurat → h5ad
│   │   │     ├── seurat_celldex.R  # Autmatic cell annotation using celldex and SingleR packages
│   │   │     ├── seurat_markers.R
│   │   │     └── Harmony_PCA.R     # Perform PCA with or whitout harmony 
│   │   │                                 (in nextflow.config harmony parameter)
│   │   ├── HIV scanpy/
│   │   │       ├── seurat_to_scanpy.py # script to produce clusters and cell markers after seuratDisk
│   │   │       └── results/            # contains some sample images
│   │   │                                
│   │   └──  results/                # contains images and other deliverables
│   │   
│   ├── CellChat/                    # For cellular comunication analysis\
│   ├── GO analysis/                 # gene ontology script
│   ├── Scrublet/                    # dubplets removal
│   ├── SoupX/                       # removal of environmental RNA
│   ├── tools/
│   │   ├── biomaRt_script.R        # to change between gene IDs
│   │   └── download_GSM.py         # snipet to automatize GSM downloading from GEO repositories
│   │
│   ├── seurat/               # contains the HIV seurat pipe used for nf-scRNA and a analysis of adipocytes
│   │   ├── Adipocytes
│   │   └── HIV PBMC  
│   │
│   └── Docs/                 # Notebooks and explanations
│       └── workflow_diagram.html / .png
│
├── Results/                  # Output structure (auto-created)
└── README.md                 # This file
```

------------------------------------------------------------------------

## ⚙️ Pipeline Flow

The project can be executed in modular R mode or as a complete **Nextflow pipeline**.

**Main steps:**

1.  **Data import and merging** (`read_h5.R`) Reads 10x Genomics HDF5 or matrix folders, merges across samples.

2.  **Quality control and filtering** (`cell_filter.R`) Filters cells by counts, mitochondrial genes, or customizable metrics.

3.  **Normalization** (`SCTransform.R`) Performs variance stabilization using Seurat’s SCTransform.

4.  **Dimensional reduction and clustering** PCA/UMAP/Leiden analysis in Seurat using (`seurat_PACA_UMAP.R`)or Scanpy (`RDS_to_h5ad.R`).

5.  **Cell type annotation** automated using (`seurat_celldex.R`) or manual

6.  **Cross-platform export** (`RDS_to_h5ad.R`) Converts `.rds` objects into `.h5Seurat` and `.h5ad` formats for Python.

7.  **Parallel execution** (`main.nf`) Nextflow handles process dependencies, caching, and reproducibility.

------------------------------------------------------------------------

## 🧩 Key Features

-   **Hybrid R–Nextflow design:** combine scripting flexibility with pipeline reproducibility.
-   **Scanpy compatibility:** use SeuratDisk for export into the Python ecosystem.
-   **QC visualization:** automatic PNG reports of filtering thresholds and cell statistics.
-   **Modular structure:** every stage can be run independently or within the main pipeline.
-   **Local container support:** Docker/Singularity optional, no internet required.

------------------------------------------------------------------------

## 🧰 Requirements

| Component | Version | Purpose |
|-----------------------------|-------------------|------------------------|
| **Nextflow** | ≥ 23.04 | Pipeline orchestration |
| **R** | ≥ 4.2 | Core analysis environment |
| **Seurat** | ≥ 5.0 | Data handling and normalization |
| **SeuratDisk** | ≥ 1.1 | Conversion to `.h5Seurat` and `.h5ad` |
| **celldex** | \- | Markers to annotate cell type |
| **SingleR** | \- |  |
| **SoupX**, **dplyr**, **ggplot2**, **patchwork** | – | QC and visualization |
| *(Optional)* **Harmony** | \- | batch effect harmonization |
| *(Optional)* **Scanpy**, **Anndata** | – | Downstream analysis in Python |

------------------------------------------------------------------------

## 🚀 Quick Start

1.  **Clone the repository**

    ``` bash
    git clone https://github.com/bontix-77/scRNAseq.git
    cd scRNAseq/Code/nf-scRNA
    ```

2.  **Run the pipeline**

    ``` bash
    nextflow run main.nf \
      --inputDir ../../Data \
      --resultDir ../../Results \
      --scriptFile ../R \
      --analysis scanpy
    ```

3.  **Inspect results**

    -   Outputs are organized under `Results/`
    -   Check the `.rds`, `.h5Seurat`, `.h5ad`, and `.png` files for processed data and QC plots.

------------------------------------------------------------------------

## 📊 Outputs

| Folder         | Description      |
|----------------|------------------|
| `merged_h5/`   | no visuals       |
| `cell_filter/` | 5 resuming .png  |
| `SCTransform/` | no visuals       |
| `seuratDisk/`  | no visuals       |
| seurat_PCA/    | 2 resulting .png |
| seura_cellType | 2 resulting .png |

------------------------------------------------------------------------

## 🧭 Purpose & Vision

The repository is meant to bridge the gap between **bench experience** and **bioinformatic automation**. It shows how a researcher can transform individual R analysis scripts into a full, portable workflow.

Future directions include:

-   Integrating cell-cell communication and trajectory analysis modules (e.g., CellChat, Monocle3)
-   Adding automated report generation via **Quarto** or **RMarkdown**
-   Providing Docker/Singularity containers for reproducible environments

------------------------------------------------------------------------

## 📄 License

MIT License © 2025 Alexander Bontempo

------------------------------------------------------------------------

Would you like me to make this version **optimized for public recruiters/labs** (a bit more narrative about your transition and design thinking), or keep it **strictly technical/documentation-style** for GitHub users and contributors?

Database used for the adipocite project can be found [here](https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/GettingStarted_scRNASeq.zip).
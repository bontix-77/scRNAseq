library(Seurat)
library(SeuratDisk)

args <- commandArgs(trailingOnly = TRUE)

HIV_SCT <- readRDS(args[1]) # Load the Seurat object
# during the SeuratData i have a problem as the gene names are lost so i save the vector in a file
# Get feature names from SCT assay (or RNA)
gene_names <- rownames(HIV_SCT@assays$SCT@counts) # or SCT@counts if using SCTransform

# Save to a text file
writeLines(gene_names, "gene_names.txt")

# Convert Seurat object (.rds) to h5Seurat
SaveH5Seurat(HIV_SCT, filename = "data.h5Seurat")

# Convert .h5Seurat to .h5ad
Convert("data.h5Seurat", dest = "scanpy.h5ad")

library(harmony)
library(Seurat)
library(uwot)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
path <- args[1]

HIV_SCT <- readRDS(path)


################################################################
# HIV_SCT <- readRDS('/home/alexander-bontempo/Desktop/GitHub/scRNAseq/Code/nf-scRNA/work/e1/73c3b5d7a7582cfd2a9fbac16c7fe1/HIV_SCTransform.rds')

DefaultAssay(object = HIV_SCT) <- "SCT"
# run PCA

HIV_SCT <- RunPCA(HIV_SCT, verbose = FALSE, assay = "SCT")

if (args[2] == "yes") {
  library(harmony)
  HIV_SCT <- RunHarmony(HIV_SCT, c("orig.ident"))
  print("harmony selected")
}


# visualizethte first 9 PC
png("DimHeatmap.png", width = 1200, height = 900, res = 150)
DimHeatmap(
  HIV_SCT,
  dims = 1:9,
  cells = 500,
  balanced = TRUE,
  ncol = 3
)
dev.off()
# perform the elbow analisys
elbow <- ElbowPlot(HIV_SCT, ndims = 40)

png("elbowplot.png", width = 1200, height = 900, res = 150)
print(elbow)
dev.off()
print("control print")

reduction_type <- if ("harmony" %in% names(HIV_SCT@reductions)) {
  "harmony"
} else {
  "pca"
}

# performing the neighbors and clusters
resolution  <-  as.numeric(args[3])
HIV_SCT <- HIV_SCT %>%
  FindNeighbors(reduction = reduction_type) %>%
  FindClusters(resolution = resolution)

png(paste0(reduction_type, ".png"), width = 1200, height = 900, res = 150)
DimPlot(
  HIV_SCT,
  reduction = reduction_type,
  group.by = c("orig.ident", "seurat_clusters"),
  alpha = 0.2,
  ncol = 2
)
dev.off()
# save the final RDS file
print(paste0("argomento3=",args[3]))
HIV_SCT <- HIV_SCT %>% RunUMAP(reduction = reduction_type, dims = 1:20)
png("UMAP_DimPlot.png", width = 1200, height = 900, res = 150)
DimPlot(
  HIV_SCT,
  reduction = "umap",
  group.by = c("orig.ident", "seurat_clusters"),
  alpha = 0.2,
  ncol = 2
)
dev.off()
# save the final RDS file

saveRDS(HIV_SCT, "HIV_filt_PCA.rds")

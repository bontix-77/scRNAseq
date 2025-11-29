# # 1. Define a cache folder in your Scratch space
# my_cache_dir <- "/home/abontemp/.cache/gypsum/bucket/blueprint_encode/2024-02-26"

# # 2. Create it if it doesn't exist
# if (!dir.exists(my_cache_dir)) {
#   dir.create(my_cache_dir, recursive = TRUE)
# }

# # 3. Tell R (and Bioconductor/celldex) to use this folder
# Sys.setenv("XDG_CACHE_HOME" = my_cache_dir)
# Sys.setenv("R_USER_CACHE_DIR" = my_cache_dir)




library(celldex)
library(Seurat)
library(SingleR)
# library(MAST)

# HIV_pp_mks <- readRDS("/home/alexander-bontempo/Desktop/GitHub/scRNAseq/Code/nf-scRNA/work/78/7b6dab1202a195fd238b385c79d9c5/markers_cluster.rds")

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
HIV_pp_mks <- readRDS(path)

HIV.sce <- as.SingleCellExperiment(HIV_pp_mks, assay = "SCT") # This selects *only* the SCT assay
#mouseRNASeq <- celldex::MouseRNAseqData()

cellde_path <- paste0(args[2],"/celldex.rds")
humanRNA <- readRDS(cellde_path)
head(humanRNA)
table(humanRNA$label.main) # principal cell types
table(humanRNA$label.fine) # subtypes

annot <- SingleR::SingleR(test = HIV.sce, ref = humanRNA, labels = humanRNA$label.fine)
head(annot)

# if a label is to weak during SingleR the cell is taged as NA. Now we add the labels determine to the metadata of the Seurat object

table(annot$pruned.labels, useNA = "always") # useNA can be turned on in the `table` function

HIV_pp_mks$humanRNASeq.fine <- annot$pruned.labels
print('complete here before ------------------------------')
## let's visualize the final result

# annotFig1 <- DimPlot(HIV_pp_mks, group.by = "cell_type", label = T) + NoLegend()
annotFig2 <- DimPlot(HIV_pp_mks, group.by = "humanRNASeq.fine", label = T, repel = TRUE, label.size = 2)
print('complete here after ------------------------------')

png("CellType_celldex.png",width = 1200,height = 900,res=150)
print(annotFig2)
dev.off()


## annotation by cluster

Idents(HIV_pp_mks) <- "seurat_clusters" # Assign clusters as the identities
avgExp <- AggregateExpression(HIV_pp_mks, assays = "SCT")$SCT # Run AverageExpression on the SCT assay and return only SCT
clustAnnot <- SingleR::SingleR(test = avgExp, ref = humanRNA, labels = humanRNA$label.fine) # Run SingleR on the averaged expression matrix
clustAnnot


clustLabels <- as.vector(clustAnnot$pruned.labels) # retrieve only the cluster-derived annotations
number_clusters  <- length(levels(HIV_pp_mks@meta.data$seurat_clusters))-1
names(clustLabels) <- c(0:number_clusters) # assign the cluster numbers as the annotations
clustLabels.vect <- clustLabels[match(HIV_pp_mks$seurat_clusters, names(clustLabels))] # match the cluster identities per cell in the Seurat data to the cluster labels
names(clustLabels.vect) <- colnames(HIV_pp_mks) # ensure that the cluster identities are assigned the cell names
HIV_pp_mks$humanRNA.fine.clust <- clustLabels.vect # add the cluster annotations to the vector

clustAnnotFig1 <- DimPlot(HIV_pp_mks, group.by = "seurat_clusters", label = T,repel= T, label.size=2) + NoLegend()
# clustAnnotFig2 <- DimPlot(HIV_pp_mks, group.by = "cell_type", label = T) + NoLegend()
clustAnnotFig3 <- DimPlot(HIV_pp_mks, group.by = "humanRNASeq.fine",label= T,repel= T, label.size=2) + NoLegend()
clustAnnotFig4 <- DimPlot(HIV_pp_mks, group.by = "humanRNA.fine.clust", label = T,repel= T, label.size=2)
png("CellType_final.png",width = 1200,height = 900,res=150)
print(clustAnnotFig1/( clustAnnotFig3| clustAnnotFig4))
dev.off()
# (clustAnnotFig1 | clustAnnotFig2) / (clustAnnotFig3 | clustAnnotFig4)
saveRDS(HIV_pp_mks,"celldex.rds")
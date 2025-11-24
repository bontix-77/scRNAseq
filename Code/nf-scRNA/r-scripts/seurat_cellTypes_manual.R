library(Seurat)
HIV_pp_mks <- readRDS("/home/alexander-bontempo/Desktop/GitHub/scRNAseq/Code/nf-scRNA/work/78/7b6dab1202a195fd238b385c79d9c5/markers_cluster.rds")

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
HIV_pp_mks <- readRDS(path)

cells <- vector(length = ncol(HIV_pp_mks))

#assign at each cell in the cluster the cell type tag
cells[which(HIV_pp_mks$seurat_clusters%in% c(0))] <- "T cells naive"
cells[which(HIV_pp_mks$seurat_clusters %in% c(1))] <- "T cells activated"
cells[which(HIV_pp_mks$seurat_clusters%in% c(2))] <- "citotoxic T CD8 and NK cells"
cells[which(HIV_pp_mks$seurat_clusters%in% c(3))] <- "monocytes/mieloid"
cells[which(HIV_pp_mks$seurat_clusters%in% c(4))] <- "B cells"
cells[which(HIV_pp_mks$seurat_clusters%in% c(5))] <- "plasma cells"
cells[which(HIV_pp_mks$seurat_clusters%in% c(6))] <- "megakaryocytes/platelets"


HIV_pp_mks$cell_type <- cells


f1 <- DimPlot(HIV_pp_mks, group.by = "seurat_clusters", label = T)
#f2 <- DimPlot(HIV_pp, group.by = "time_point") + NoLegend()
f3 <- DimPlot(HIV_pp_mks, group.by = "cell_type", label = T) + NoLegend()
f3
png("CellType.png",width = 1200,height = 900,res=150)
print( f1 | f3)
dev.off()

saveRDS(HIV_pp_mks,"Cell_type.rds")
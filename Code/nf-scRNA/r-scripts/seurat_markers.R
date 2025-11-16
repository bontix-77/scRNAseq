library(Seurat)
# find markers (many method available, default a Wilcox rank sum. Look up test.use parameters for options)
# requires presto installation to speed up the next step
# devtools::install_github('immunogenomics/presto')
library(dplyr)
library(presto)
args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
HIV_PCA_cluster<- readRDS(path)

# HIV_PCA_cluster <- readRDS("/home/alexander-bontempo/Desktop/GitHub/scRNAseq/Code/nf-scRNA/work/3c/b2db263f999da2c26c2d99201d2011/HIV_filt_PCA.rds")

HIV_PCA_cluster <- PrepSCTFindMarkers(HIV_PCA_cluster, verbose = T)



# find all markers
HIV_PCA_cluster_markers <- FindAllMarkers(HIV_PCA_cluster, only.pos = TRUE)
# ordering the results
HIV_PCA_cluster_markers <- HIV_PCA_cluster_markers %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj)

# examine a small subset
HIV_PCA_cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
# markers visualization
top20 <- HIV_PCA_cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup()
top20_df <- data.frame(
     cluster=top20$cluster,
     gene=top20$gene
)

write.csv(top20_df,"Top20_markers.csv")

DoHeatmap(HIV_PCA_cluster, features = top20$gene) + NoLegend()

# associate pubblished markers with cell type
# contaminants platelets
platelets <- c("Ptprc", "Plp1", "Pecam1")

Naive_B <- c("CD19","IGHD")

Memory_B<- c("CD19","CD27","IGHG1","MS4A1")

CD14_monocytes <- c("CD14","S100A8","CD4")

CD16_monocytes  <- c("FCGR3A","CD4")

NK <- c("GZMA","GZMB","GNLY","NKG7","FCGR3A","CD8A")
#negative markers TRDC, GZMK and GATA3
CD56brightNK <- c("NCAM1", "TRDC","GZMK","GATA3") 
CD4_naiveT <- c("CD3E","CD4","CCR7","IL7R","TCF7","PTPRC","CD27")
#negative markers PTPRC
CD4_Tcm <- c("CD3E","CD4","CD27","SELL","CCR7","PTPRC")
#negative markers CD27, CCR7, PTPRC, SELL, TCF7
CD4_Tem <- c("CD3","CD4","CD27", "CCR7", "PTPRC", "SELL", "TCF7")
CD4_Trm <- c("CD3E","CD4","CD69","ITGAE")
CD4_Tscm <- c("CD3E", "CD4", "CD27", "SELL", "PTRPC", "IL2RB", "FAS")

CD8_niveT  <- c("CD3E","CD8A","CCR7","TCF7","PTPRC","CD27")
# negative marker PTPRC
CD8_Tcm <- c("CD3E","CD8A","CD27","SELL","CCR7","PTPRC")
#negative markers CD27, CCR7, PTPRC, SELL, TCF7
CD8_Tem <- c("CD3","CD8A","CD27", "CCR7", "PTPRC", "SELL", "TCF7")
CD8_Trm <- c("CD3E","CD8A","CD69","ITGAE")
CD8_Tscm <- c("CD3E", "CD8A", "CD27", "SELL", "PTRPC", "IL2RB", "FAS")
#mucosal-associated invariant T
MAIT <- c("CD3E", "TRAV1-2", "IL7R", "GZMK", "CCR6")

Vγ9Vδ2T	 <- c("CD3E", "TRGV9", "TRDV2", "GZMA", "CCL5", "TRDC")
#negative markers TRGV9, TRDV2
gdT <- c("CD3E", "TRGC1","TRGC2", "TRDC",	"TRGV9", "TRDV2")

Treg <- c("CD3E", "FOXP3", "IL2RA")
T_prolif <- c("CD3E","MKI67")

mDC <- c("CST3","CD1C", "HLA-DRA", "CD4")

pDC <- c("CST3", "CLEC4C", "CXCR3", "IL3RA", "GZMB", "CD4")
Plasmoblasts <- c("JCHAIN", "CD27", "MKI67", "IGHD", "IGHG1")
#hematopoietic stem cells
HSC <- c("PPBP", "ITGA2B")
#Red cells
RBC <- c("HBB", "HBA1", "HBA2")

VlnPlot(HIV_PCA_cluster, features = Plasmoblasts)

# visualize in the cluster the cel type corresponding to the  feature= parameter above (this case contam)
FeaturePlot(HIV_PCA_cluster, features = CD4_Tcm)

png("Clusters.png",width = 1200,height = 900,res=150)
DimPlot(HIV_PCA_cluster,reduction="umap",group.by="seurat_clusters")
dev.off()

saveRDS(HIV_PCA_cluster,"markers_cluster.rds")
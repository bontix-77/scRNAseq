library(Seurat)
library (DESeq2)


args <- commandArgs(trailingOnly=TRUE)
path <- args[1]

HIV <- readRDS(path)

HIV_pp_mks <- readRDS(path)
pseudo_HIV<- AggregateExpression(HIV_pp_mks, assays = "SCT", return.seurat = T, group.by = c("orig.ident","group"))#, "time_point", "condition", "condition_tp"))


# Assuming 'pb' is your aggregated Seurat object
# Run standard PCA on the pseudo-bulk samples
pb <- NormalizeData(pseudo_HIV)
pb <- FindVariableFeatures(pb)
pb <- ScaleData(pb)
pcs <- length(pb@meta.data$orig.ident) - 1

pb <- RunPCA(pb, npcs = pcs)

DimPlot(pb, group.by = "group", pt.size = 3)


pseudo_HIV@assays$SCT@layers$counts@Dimnames[1] <- HIV_pp_mks@assays$SCT@data@Dimnames[1]
head(pseudo_HIV@assays$SCT$counts)
pseudo_HIV@meta.data
glimpse(pseudo_HIV)

# just to clean up the look a little bit
pseudo_HIV <- RenameCells(pseudo_HIV, new.names = gsub("_.*", "", pseudo_HIV$orig.ident))
pseudo_HIV$orig.ident <- gsub("_.*", "", pseudo_HIV$orig.ident)
head(pseudo_HIV@assays$SCT$counts)
pseudo_HIV@meta.data

## performin bulk DE

Idents(HIV_pp_mks) <- "group"
Idents(pseudo_HIV) <- "group"
table(Idents(HIV_pp_mks))
scDE <- FindMarkers(HIV_pp_mks, ident.1 = "undetectable",ident.2 = "detectable", test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1)


bulk_HIV_de <- FindMarkers(pseudo_HIV, ident.1 = "undetectable", ident.2 = "detectable", test.use = "DESeq2")
head(bulk_HIV_de)

# comparing how many differentially expressed genes between SC and bulk analisys comparing the conditions

scDE.genes <- rownames(scDE)[which(scDE$p_val_adj < 0.05)]
bulkDE.genes <- rownames(bulk_HIV_de)[which(bulk_HIV_de$p_val_adj < 0.2)]
length(scDE.genes)
length(bulkDE.genes)

## heck the common features between sc and bulk
length(intersect(scDE.genes, bulkDE.genes))
head(intersect(scDE.genes, bulkDE.genes), 30)

## to chech spefic features

bulk_HIV_de[c("TCR_A", "TRC_B"), ]

## visualize the DE genes

Idents(HIV_pp) <- "seurat_clusters"
DotPlot(HIV_pp, features = unique(top5PerCluster$gene), dot.scale = 3) + coord_flip()

# violine as alternative visualization

Idents(HIV_pp) <- "time_point"
VlnPlot(HIV_pp, features = c("Acta2", "Cd36"), alpha = 0.1)


## anotate the differential genes

markers <- c("CCR7", "NKG7", "LYZ", "MS4A1", "MZB1", "PPBP", "TCF7")
Idents(HIV_pp_mks) <- "SCT_snn_res.0.13"
avgExp <- AverageExpression(HIV_pp_mks, markers, assay = "SCT")$SCT
avgExp
DimPlot(HIV_pp_mks, label = T)
FeaturePlot(HIV_pp_mks, features = markers, ncol = 3, order = T)
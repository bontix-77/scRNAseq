library(tidyverse)
library(Seurat)

##############
#explanation of the S4 structure
# i  is the index associated to the features listed in Dimnames[[1]]
# p is the pointer to show where the cells expression start and finish. It start to 0.  
#    i 34 234 534 454
#    p 0 1754
#    x 1 2 4 3 
#    this meand that the first 1754 indexes found in the i vector belong to the expression of the first cell 
# Dimnames
#      [[1]] is the name of the genes indicized in the vector i
#      [[2]] is the name/tag of each cell   
# x is the gene expression 
#     in the exemple above the gen 34 in the cell 1 has expression 1. looking into Dimnames[[1]] i can find the name of the gene 34
setwd("/home/alexander-bontempo/Desktop/HIV GSM/GSM1")

# List of the samples files. In this case we have 3  files for each sample: matrix.mtx, barcodes.tsv.gz and features.tsv.gz used to map raw reads in X10 Genomics Chromium systems.
dirs <- list.dirs()
dirs_name <- basename(dirs[dirs !="./"])
dirs_name <- dirs_name[-1]

# Create a list of count matrices
paste( "/home/alexander-bontempo/Desktop/HIV GSM/GSM1/",dirs[2])

reads <- lapply(paste0("/home/alexander-bontempo/Desktop/HIV GSM/GSM1/",dirs_name,"/")
  , Read10X)
# Assign names manually GEO serie GSE220790
#names(reads) <- c("6817423", "6817424", "6817425", "6817426", "6817427", "6817428")
names(reads) <- c("6817423", "6817431", "6817432", "6817433", "6817434", "6817435","6817436")
#########################################################################
# Single sample
##   W10 <- CreateSeuratObject(counts=W10, project="W10", min.cells = 3, min.features = 200)
# View W10
## W10
#########################################################################

# here i filter the cells with at least 200 features (genes) detected and genes detected in at least 3 cells

  
# HIV <-  mapply(CreateSeuratObject,
 # counts = reads,
  #project = names(reads),
  #MoreArgs = list(min.cells = 20, min.features = 200)
#)
#but since HIV expression in PLWH is very low lets remoove min cell filtering
 
HIV <-  mapply(CreateSeuratObject,
  counts = reads1,
  project = names(reads),
  MoreArgs = list( min.cells = 2,min.features = 200)
)

# remove the original sparse matrices
rm(reads)

# assign metadata group
HIV$`6817423`$group <- "detectable"
HIV$`6817431`$group <- "undetectable"
HIV$`6817432`$group <- "undetectable"
HIV$`6817433`$group <- "detectable"
HIV$`6817434`$group <- "undetectable"
HIV$`6817435`$group <- "undetectable"
HIV$`6817436`$group <- "detectable"


HIV <- merge(HIV[[1]],
  y = HIV[2:length(HIV)],
  add.cell.ids = names(HIV), 
  project = "HIV PBMC on patient under ART"
)

#############################################################################################################
# If you want to apply QC metrics independently for each sample, you can use this for
# loop to create an individual object for each sample.
## for (file in files){
##  seurat_data <- Read10X_h5(paste0("C:/Users/Owner/Documents/github/Seurat test script/seurat 2 basic script/R project dowloaded from manual/GettingStarted_scRNASeq/data/",
##                                   file))
##  seurat_obj <- CreateSeuratObject(counts = seurat_data,
##                                   min.features = 200,
##                                   min.cells = 3,
##                                   project = file)
##  assign(file, seurat_obj)
## }
## The metadata in the Seurat object is located in HIV@metadata and contains the
## information associated with each cell.
################################################################################



head(HIV@meta.data) # using head to return only the first 6 rows


# Access a single column.
head(HIV$orig.ident)
# or
head(HIV[["orig.ident"]])
# acess multiple columns
# Access multiple columns, rows.
head(HIV[[c("orig.ident", "nCount_RNA")]])[1:3, ]

# to save the Seurat object


#saveRDS(HIV, "/home/alexander-bontempo/Desktop/HIV GSM/RDS/HIV_merged.rds")


library(tidyverse) # dplyr and ggplot2
library(Seurat) # Seurat toolkit
library(hdf5r) # for data import
library(patchwork) # for plotting
library(presto) # for differential expression
library(glmGamPoi) # for sctransform
library(ggplot2)

# read the object
## HIV <- readRDS("C:/Users/Owner/Documents/github/Seurat test script/seurat 2 basic script/outputs/HIV_merged.rds")


glimpse(HIV)


# calculate mitocondrial genes porcentage for quality control and normalization. Mitocondrial RNA increase need to be assesd 
#to determine if related to poor cell viability or in response of biological relevant processes.

HIV[["percent.mt"]] <- PercentageFeatureSet(HIV, pattern = "^MT-")
# set colors
cnames <- setNames(c("green","blue","orange","yellow","black0","red",'purple'), levels(factor(HIV@meta.data$orig.ident)))
cnames
#HIV <- NormalizeData(HIV)

# plot total counts per sample
VlnPlot(HIV, features = "nCount_RNA",  group.by = "orig.ident", raster = FALSE, alpha = 0.2) 
VlnPlot(HIV, features = "nCount_RNA",  group.by = "group", raster = FALSE, alpha = 0.2) 
# or using ggplot2
HIV@meta.data %>%
  ggplot(aes( color = orig.ident,x = nCount_RNA, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 650, color = "red", linetype = "dotted")

HIV@meta.data %>%
  ggplot(aes( color = group,x = nCount_RNA, fill = group)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 650, color = "red", linetype = "dotted")

# list all assays present
Assays(HIV)
# check what slots are in my RNA assay
slotNames(HIV[["RNA"]])



# determine number of features (genes)

plot_featur_i<- VlnPlot(HIV, features = "nFeature_RNA", group.by = "orig.ident") +
  geom_hline(yintercept = 600, color = "red")
show(plot_featur_i)
plot_featur_g <- VlnPlot(HIV, features = "nFeature_RNA", group.by = "group") +
  geom_hline(yintercept = 600, color = "red")
show(plot_featur_g)
# visualize % mitocondrial genes

plot_mt_i <- VlnPlot(HIV, features = "percent.mt", group.by = "orig.ident") +
  geom_hline(yintercept = 10, color = "red")
show(plot_mt_i)
plot_mt_g <- VlnPlot(HIV, features = "percent.mt", group.by = "group") +
  geom_hline(yintercept = 10, color = "red")



plot_count <- VlnPlot(HIV, features = "nCount_RNA", group.by = "orig.ident") +
  geom_hline(yintercept = 10, color = "red")
(plot_featur_i/plot_mt_i)|(plot_featur_g/plot_mt_g)

# scatter the point by number of features and mitocondrial genes

FeatureScatter(HIV, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by = "orig.ident")

# and by number of fetures and total counts
FeatureScatter(HIV, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", log = TRUE)


# filter
# Set one set of parameters for Day 0 samples;
# keep the rownames (Cell barcodes)
# F_30 <- HIV@meta.data |>
#   filter(orig.ident=="6817430",
#      nFeature_RNA > 500,
#      nFeature_RNA < 4000,
#     percent.mt < 10
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)
# F_31 <- HIV@meta.data |>
#   filter(orig.ident=="6817431",
#          nFeature_RNA > 500,
#          percent.mt < 10,
#          nFeature_RNA < 4000
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)
# F_32 <- HIV@meta.data |>
#   filter(orig.ident=="6817432",
#          percent.mt < 10,
#          nFeature_RNA > 500,
#         nFeature_RNA < 4000
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)

# F_33 <- HIV@meta.data |>
#   filter(orig.ident=="6817433",
#          nFeature_RNA > 500,
#          percent.mt < 10
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)
# F_34 <- HIV@meta.data |>
#   filter(orig.ident=="6817434",
#          nFeature_RNA > 700,
#          percent.mt < 10
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)
# F_35 <- HIV@meta.data |>
#   filter(orig.ident=="6817435",
#          nFeature_RNA > 500,
#          nFeature_RNA < 4000,
#          percent.mt < 10
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)
# F_36 <- HIV@meta.data |>
#   filter(orig.ident=="6817436",
#          nFeature_RNA > 500,
#          nFeature_RNA < 4000,
#          percent.mt < 10
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)

# keep <- c(F_30,F_31,F_32,F_33,F_34,F_35,F_36)

to_keep <- HIV@meta.data |>
  filter(nFeature_RNA > 600,
     nFeature_RNA < 2500,
    percent.mt < 10
  ) |>
  tibble::rownames_to_column("Cell") |>
  pull(Cell)


# use different parameters; established above
HIV_filt <- subset(HIV, cells = to_keep)
FeatureScatter(HIV_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", log = TRUE)
FeatureScatter(HIV_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "group", log = TRUE)


# save file after filtering

saveRDS(HIV_filt, "/home/alexander-bontempo/Desktop/HIV GSM/RDS/HIV_merge_filt.rds")


##############################################################
## Normalization, variable feature  and scale               ##
##############################################################


# run sctransform
#SCTransform usaually keep just 3000 features to increase variable.features.rv.th = 1.0 means to keep all the features. 
#0.9 means to keep 90% of the features based on variance. 
#alternively you can set variable.features.n to a specific number of features to keep.
HIV_filt <- SCTransform(HIV_filt, vars.to.regress = "percent.mt", verbose = FALSE)

# Check default assay
############ DefaultAssay(object = HIV_filt)

# Set default assay
########### DefaultAssay(object = HIV_filt) <- "RNA"
#########################################################################################
######################## save seurat object to h5ad to transfer to scanpy###############
library(SeuratDisk)                                                                     #
# during the SeuratData i have a problem as the gene names are lost so i save the vector in a file 
# Get feature names from SCT assay (or RNA)
gene_names <- rownames(HIV_filt@assays$SCT@scale.data)  # or SCT@counts if using SCTransform

# Save to a text file
writeLines(gene_names, "/home/alexander-bontempo/Desktop/HIV GSM/h5ad/gene_names.txt")

# Convert Seurat object (.rds) to h5Seurat
SaveH5Seurat(HIV_filt, filename = "/home/alexander-bontempo/Desktop/HIV GSM/h5ad/data.h5Seurat")

# Convert .h5Seurat to .h5ad
Convert("/home/alexander-bontempo/Desktop/HIV GSM/h5ad/data.h5Seurat", dest = "/home/alexander-bontempo/Desktop/HIV GSM/h5ad/data.h5ad")
#########################################################################################
###########################################################################################


# run PCA

HIV_filt <- RunPCA(HIV_filt, verbose = FALSE, assay = "SCT")


library(harmony)

HIV_filt <- RunHarmony(HIV_filt,"orig.ident")
# visualizethte first 9 PC
DimHeatmap(HIV_filt, dims = 1:9, cells = 500, balanced = TRUE, ncol = 3)
# perform the elbow analisys
ElbowPlot(HIV_filt, reduction = "harmony",ndims = 40)

# performing the neighbors analisys to prepara for the clusterin
HIV_filt <- FindNeighbors(HIV_filt,reduction="harmony", dims = 1:30)

# calculate the clusterin (suggessted redolution 0.4-1.2)

HIV_filt <- FindClusters(HIV_filt, reduction = "harmony",resolution = 0.8)

# UMAP for the viasualization of the clusters
HIV_filt <- RunUMAP(HIV_filt, dims = 1:20,reduction = "harmony")
DimPlot(HIV_filt,
  reduction = "harmony", group.by = c("orig.ident", "SCT_snn_res.0.8","group"),
  alpha = 0.2, ncol = 2,
  label=F
)
figumap_clusters <- DimPlot(HIV_filt,
  reduction = "umap", group.by =  "SCT_snn_res.0.8",
  alpha = 0.2, ncol = 1,
  label=F
  
)#+ NoLegend()
figumap_orig <- DimPlot(HIV_filt,
  reduction = "umap", group.by = "orig.ident" ,
  alpha = 0.2, ncol = 1,
  label=F

)
figumap_orig|figumap_clusters
# save the HIV filtered file

#saveRDS(HIV_filt, "/home/alexander-bontempo/Desktop/HIV GSM/RDS/HIV_merge_filt_sctran_clust0.1.rds")

# prepare the SCT data for the search of markers among clusters
HIV_filt <- PrepSCTFindMarkers(HIV_filt, verbose = T)


# find markers (many method available, default a Wilcox rank sum. Look up test.use parameters for options)
# requires presto installation to speed up the next step
# devtools::install_github('immunogenomics/presto')
library(presto)
# find all markers
HIV_filt_markers <- FindAllMarkers(HIV_filt, only.pos = F)
# ordering the results
HIV_filt_markers <- HIV_filt_markers %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj)

# examine a small subset
HIV_filt_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
# markers visualization
top20 <- HIV_filt_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup()
DoHeatmap(HIV_filt, features = top20$gene) + NoLegend()

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
CD4_Tem <- c("CD3E","CD4","CD27", "CCR7", "PTPRC", "SELL", "TCF7")
CD4_Trm <- c("CD3E","CD4","CD69","ITGAE")
CD4_Tscm <- c("CD3E", "CD4", "CD27", "SELL", "PTPRC", "IL2RB", "FAS")

CD8_niveT  <- c("CD3E","CD8A","CCR7","TCF7","PTPRC","CD27")
# negative marker PTPRC
CD8_Tcm <- c("CD3E","CD8A","CD27","SELL","CCR7","PTPRC")
#negative markers CD27, CCR7, PTPRC, SELL, TCF7
CD8_Tem <- c("CD3E","CD8A","CD27", "CCR7", "PTPRC", "SELL", "TCF7")
CD8_Trm <- c("CD3E","CD8A","CD69","ITGAE")
CD8_Tscm <- c("CD3E", "CD8A", "CD27", "SELL", "PTPRC", "IL2RB", "FAS")
#mucosal-associated invariant T
MAIT <- c("CD3E", "TRAV1-2", "IL7R", "GZMK", "CCR6")

Vγ9Vδ2T	 <- c("CD3E", "TCR","TRGV9", "TRDV2", "GZMA", "CCL5", "TRDC")
#negative markers TRGV9, TRDV2
gdT <- c("CD3E", "TCR","TRGC1","TRGC2", "TRDC",	"TRGV9", "TRDV2")

Treg <- c("CD3E", "FOXP3", "IL2RA")
T_prolif <- c("CD3E","MKI67")

mDC <- c("CST3","CD1C", "HLA-DRA", "CD4")

pDC <- c("CST3", "CLEC4C", "CXCR3", "IL3RA", "GZMB", "CD4")
Plasmoblasts <- c("JCHAIN", "CD27", "MKI67", "IGHD", "IGHG1")
#hematopoietic stem cells
HSC <- c("PPBP", "ITGA2B")
#Red cells
RBC <- c("HBB", "HBA1", "HBA2")


cd8 <- VlnPlot(HIV_filt, features = "CD8A")
cd4 <- VlnPlot(HIV_filt, features = "CD4")
cd4/cd8

# visualize in the cluster the cel type corresponding to the  feature= parameter above (this case contam)
FeaturePlot(HIV_filt, features = "CD4")

#saveRDS(HIV_filt_markers, "/home/alexander-bontempo/Desktop/HIV GSM/RDS/HIV_merge_filt_markers.rds")


## to find different expression between a label and the others in a metadata idents (if between two specific labels "idents.2=" can be used to specify.)
Idents(HIV_filt) <- "group"
# here if you want to look differential expression between specific clusters
to_keep <- HIV_filt@meta.data|>
  filter(seurat_clusters%in%c(2,3)
  ) |>
  tibble::rownames_to_column("Cell") |>
  pull(Cell)
HIV_subset <- subset(HIV_filt, cells = to_keep)

det_vs_unde <- FindMarkers(HIV_subset, ident.1 = "undetectable",ident.2="detectable", test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1)

##############################################
##    Differential Expression analisys      ##
##############################################
DefaultAssay(HIV_filt)


library(SingleR) # for cell type annotation; Bioconductor
library(celldex) # for cell type annotation reference; Bioconductor
library(MAST) # for differential expression; Bioconductor

mem.maxVSize(vsize = 15000)
glimpse(HIV_filt)
#table(HIV_filt$condition_tp)
Idents(HIV_filt) <- "SCT_snn_res.0.13"
table(Idents(HIV_filt))
######################################################################################
#         Createa new metadata with the expression of gene (example CD3D)           ##
#                                                                                   ##
## gene_name <- "CD3D"                                                               ##
#                                                                                   ##
#    Extract gene expression values from the SCT assay                              ##
## seurat_obj[[gene_name]] <- FetchData(seurat_obj, vars = gene_name, assay = "SCT") ##
######################################################################################


# Viewing the code, JoinLayers was called to link all the data sets together. This is
# a new feature of Seurat5, and is required for analyzing data after integration and batch correction, in this case fue all the samples (D10,D16 , etc.)
plot(density(sample(JoinLayers(HIV_filt@assays$RNA)$counts["GAPDH", ], 2500)), cex = 0, lty = 1, main = "Density of Gapdh in 2500 random cells")

hist(sample(JoinLayers(HIV@assays$RNA)$counts["GAPDH", ], 2500), breaks = 52, main = "Histogram of Gapdh in 2500 random cells", ylab = "Frequency", xlab = "Gene counts")

###############################################################################
##        how estract expression matrixes                                    ##
##  expr<-GetAssayData(HIV_filt,assay = "SCT",slot = "data")                 ##
##  expr_df <- as.data.frame(as.matrix(expr)) %>%                            ##
##   rownames_to_column(var = "gene") %>%                                    ##
##   pivot_longer(-gene, names_to = "cell", values_to = "expression")        ##
###############################################################################


HIV_pp <- HIV_filt

Idents(HIV_pp) <- "SCT"
HIV_pp_mks <- PrepSCTFindMarkers(HIV_pp)
#Idents(HIV_pp_mks) <- "time_point"
#table(Idents(HIV_pp_mks))

## to find different expression betwe a label and the others in a metadata idents (if between two specific labels "idents.2=" can be used to specify.)
#Day0_Day6_DE <- FindMarkers(HIV_pp_mks, ident.1 = "Day 0", test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1)


Idents(HIV_pp_mks) <- "SCT_snn_res.0.13"
#DefaultAssay(HIV_pp_mks) <- "SCT"

## find all markers between clusters

de_allClusters <- FindAllMarkers(HIV_pp_mks, test.use = "wilcox", min.pct = 0.1, only.pos = TRUE)

head(de_allClusters)

## to check the most positivelly different markers per cluster
top5PerCluster <- matrix(ncol = 7)
colnames(top5PerCluster) <- colnames(de_allClusters)
for (i in 0:7) {
  top5PerCluster <- rbind(top5PerCluster, head(de_allClusters[which(de_allClusters$cluster == i), ], 5))
}
top5PerCluster <- top5PerCluster[-1, ]
top5PerCluster

DoHeatmap(HIV_pp_mks, features = top5PerCluster$gene, slot = "scale.data")
#| Cluster | Top markers                              | Likely cell type                 |
#| ------- | ---------------------------------------- | -------------------------------- |
#| 0       | NKG7, CST7, GZMH, CCL5, GNLY             | **NK cells / cytotoxic T cells** |
#| 1       | LEF1, TCF7, CCR7, IL7R, PIK3IP1          | **Naive CD4+ T cells**           |
#| 2       | IL7R, LTB, IL32, AQP3, CD3E              | **Memory/activated T cells**     |
#| 3       | CST3, LYZ, IFI30, FCN1, S100A9           | **Monocytes**                    |
#| 4       | CD79A, MS4A1, HLA-DQA1, BANK1, LINC00926 | **B cells**                      |
#| 5       | PPBP, TUBB1, CAVIN2, SPARC, NRGN         | **Platelets / megakaryocytes**   |
#| 6       | MZB1, TXNDC5, JCHAIN, ITM2C, TNFRSF17    | **Plasma cells**                 |


#fig1 <- DimPlot(HIV_pp_mks, group.by = "time_point")

# to find a feature with a specific pattern
grep("HIV", rownames(HIV_pp_mks), value = TRUE,ignore.case = TRUE)

#changinf the assay and looking in the RNA reading before SCTransform:

DefaultAssay(HIV_pp_mks) <- "SCT"
grep("HIV", rownames(HIV_pp_mks), value = TRUE,ignore.case = TRUE)

sum(FetchData(HIV_pp_mks, vars = "HIV", assay = "RNA")>0)



fig1 <- FeaturePlot(HIV_pp_mks, features = "CD3D", order = T)
fig2 <- FeaturePlot(HIV_pp_mks, features = "CD4", order = T)
fig3 <- FeaturePlot(HIV_pp_mks, features = "CD8A", order = T)
fig2|fig3
fig1 / (fig2 | fig3)
# CD14 monocites marker
fig4 <- FeaturePlot(HIV_pp_mks, features = "CD14", order = T)
fig4
fig6 <- FeaturePlot(HIV_pp_mks, features = "CXCR4", order = T)
fig6
fig5 <- FeaturePlot(HIV_pp_mks, features = "CCR5", order = T)
fig5

fig5|fig6

FeaturePlot(
  HIV_filt,
  features = c("CD3E", "CD4","CD8A", "NKG7"),
  order = TRUE,
  #split.by = "seurat_clusters",
  cols = c("lightgrey","yellow", "red"),
  reduction = "umap"
)


#differentiate CD8 from NK cells

fig7 <- FeaturePlot(HIV_pp_mks, features = "NCAM1", order = T)
fig8 <- FeaturePlot(HIV_pp_mks, features = "NKG7", order = T)
fig7|fig8/fig3
DoHeatmap(HIV_filt, features = c("CD3E", "CD8A", "NCAM1", "NKG7")) 

## is possible to join all cells as a seample and do a pseudobulk analisys

pseudo_HIV <- AggregateExpression(HIV_pp_mks, assays = "SCT", return.seurat = T, group.by = c("orig.ident"))#, "time_point", "condition", "condition_tp"))

head(pseudo_HIV@assays$SCT$counts)
pseudo_HIV@meta.data
glimpse(pseudo_HIV)

# just to clean up the look a little bit
pseudo_HIV <- RenameCells(pseudo_HIV, new.names = gsub("_.*", "", pseudo_HIV$orig.ident))
pseudo_HIV$orig.ident <- gsub("_.*", "", pseudo_HIV$orig.ident)
head(pseudo_HIV@assays$RNA$counts)
pseudo_HIV@meta.data

## performin bulk DE

Idents(pseudo_HIV) <- "time_point"

## DESeq2 if not istalled:  BiocManager::install("DESeq2")
bulk_HIV_de <- FindMarkers(pseudo_HIV, ident.1 = "Day 0", ident.2 = "Day 6", test.use = "DESeq2")
head(bulk_HIV_de)

# comparing how many differentially expressed genes between SC and bulk analisys comparing the conditions

scDE.genes <- rownames(Day0_Day6_DE)[which(Day0_Day6_DE$p_val_adj < 0.05)]
bulkDE.genes <- rownames(bulk_HIV_de)[which(bulk_HIV_de$p_val_adj < 0.05)]
length(scDE.genes)
length(bulkDE.genes)

## heck the common features between sc and bulk
length(intersect(scDE.genes, bulkDE.genes))
head(intersect(scDE.genes, bulkDE.genes), 30)

## to chech spefic features

bulk_HIV_de[c("Acta2", "Cd36"), ]

## visualize the DE genes

Idents(HIV_pp) <- "SCT_snn_res.0.13"
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
## after find markers characteristic for a particular cluster we can annotate
# cluster 0 T cells naive
#cluster 1 T cells activated
#cluster 2 citotoxic  T CD8 and NK cells
#cluster 3  monocytes/mieloid
#cluster 4 B cells
#cluster 5 plasma cells
#cluster 6 megakaryocytes/platelets
## visualization of features

cells <- vector(length = ncol(HIV_pp_mks))
cells[which(HIV_pp_mks$SCT_snn_res.0.13 %in% c(0))] <- "T cells naive"
cells[which(HIV_pp_mks$SCT_snn_res.0.13 %in% c(1))] <- "T cells activated"
cells[which(HIV_pp_mks$SCT_snn_res.0.13 %in% c(2))] <- "citotoxic T CD8 and NK cells"
cells[which(HIV_pp_mks$SCT_snn_res.0.13 %in% c(3))] <- "monocytes/mieloid"
cells[which(HIV_pp_mks$SCT_snn_res.0.13 %in% c(4))] <- "B cells"
cells[which(HIV_pp_mks$SCT_snn_res.0.13 %in% c(5))] <- "plasma cells"
cells[which(HIV_pp_mks$SCT_snn_res.0.13 %in% c(6))] <- "megakaryocytes/platelets"

HIV_pp_mks$cell_type <- cells

sum(JoinLayers(HIV_filt@assays$RNA)$counts["CD4", ]!=0)
  

# Subset only cells expressing CD4 (expression > 0)
CD4_positive <- subset(HIV_pp_mks, subset = CD4 > 0)

# Visualize CD4 expression in UMAP only for those cells
FeaturePlot(CD4_positive, features = "CD4")
CD4 <- DimPlot(CD4_positive, group.by = "SCT_snn_res.0.13", label = T) + NoLegend()+ ggtitle("CD4+")



# Subset only cells expressing CD4 (expression > 0)
CD8_positive <- subset(HIV_pp_mks, subset = CD8A > 0 & CD8B > 0)

# Visualize CD4 expression in UMAP only for those cells
FeaturePlot(CD8_positive, features = "CD4")
CD8 <- DimPlot(CD8_positive, group.by = "SCT_snn_res.0.13", label = T) + NoLegend()+ggtitle("CD8+")

CD4|CD8

f1 <- DimPlot(HIV_pp_mks, group.by = "SCT_snn_res.0.13", label = T) + NoLegend()
#f2 <- DimPlot(HIV_pp, group.by = "time_point") + NoLegend()
f3 <- DimPlot(HIV_pp_mks, group.by = "cell_type", label = T) + NoLegend()
f3
( f1 | f3)
## SingleR annotation. this method annotate each cell in the dataset against a reference dataset
library(celldex)

HIV.sce <- as.SingleCellExperiment(HIV_pp_mks, assay = "SCT") # This selects *only* the SCT assay
#mouseRNASeq <- celldex::MouseRNAseqData()
humanRNA <- BlueprintEncodeData()
head(humanRNA)
table(humanRNA$label.main)
table(humanRNA$label.fine)

annot <- SingleR::SingleR(test = HIV.sce, ref = humanRNA, labels = humanRNA$label.main)
head(annot)

# if a label is to weak during SingleR the cell is taged as NA. Now we add the labels determine to the metadata of the Seurat object

table(annot$pruned.labels, useNA = "ifany") # useNA can be used turned on in the `table` function

HIV_pp_mks$humanRNASeq.main <- annot$pruned.labels

## let's visualize the final result

annotFig1 <- DimPlot(HIV_pp_mks, group.by = "cell_type", label = T) + NoLegend()
annotFig2 <- DimPlot(HIV_pp_mks, group.by = "humanRNASeq.main", label = T)

annotFig1 | annotFig2

## annotation by cluster

Idents(HIV_pp_mks) <- "SCT_snn_res.0.13" # Assign clusters as the identities
avgExp <- AverageExpression(HIV_pp_mks, assays = "SCT")$SCT # Run AverageExpression on the SCT assay and return only SCT
clustAnnot <- SingleR::SingleR(test = avgExp, ref = humanRNA, labels = humanRNA$label.main) # Run SingleR on the averaged expression matrix
clustAnnot


clustLabels <- as.vector(clustAnnot$pruned.labels) # retrieve only the cluster-derived annotations
names(clustLabels) <- c(0:6) # assign the cluster numbers as the annotations
clustLabels.vect <- clustLabels[match(HIV_pp_mks$SCT_snn_res.0.13, names(clustLabels))] # match the cluster identities per cell in the Seurat data to the cluster labels
names(clustLabels.vect) <- colnames(HIV_pp_mks) # ensure that the cluster identities are assigned the cell names
HIV_pp_mks$humanRNA.main.clust <- clustLabels.vect # add the cluster annotations to the vector

clustAnnotFig1 <- DimPlot(HIV_pp_mks, group.by = "SCT_snn_res.0.13", label = T) + NoLegend()
clustAnnotFig2 <- DimPlot(HIV_pp_mks, group.by = "cell_type", label = T) + NoLegend()
clustAnnotFig3 <- DimPlot(HIV_pp_mks, group.by = "humanRNASeq.main",label= T) + NoLegend()
clustAnnotFig4 <- DimPlot(HIV_pp_mks, group.by = "humanRNA.main.clust", label = T)

(clustAnnotFig1 | clustAnnotFig2) / (clustAnnotFig3 | clustAnnotFig4)

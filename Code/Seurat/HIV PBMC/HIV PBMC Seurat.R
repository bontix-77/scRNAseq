library(Seurat)
library(tidyverse)


setwd("C:/Users/Owner/Documents/scRNAseq/Data/GSM")

# List of the .h5 files. h5 files are outpur of CellRanger, used to map raw reads in X10 Genomics Chromium systems.
dirs <- list.dirs()
dirs_name <- basename(dirs[dirs !="./"])
dirs_name <- dirs_name[-1]

# Create a list of count matrices
paste( "C:/Users/Owner/Documents/scRNAseq/",dirs[2])

reads <- lapply(paste0("C:/Users/Owner/Documents/scRNAseq/Data/GSM/",dirs_name,"/")
  
, Read10X)

# Automated way to assign names; modify for your purposes
# names(h5_read)<- sapply(files,
#                        function(x){str_split_1(x,"_")[1]},
#                        USE.NAMES = FALSE)

# Assign names manually
             names(reads) <- c("6817423", "6817424", "6817425", "6817426", "6817427", "6817428")
#########################################################################
# Single sample
##   W10 <- CreateSeuratObject(counts=W10, project="W10", min.cells = 3, min.features = 200)
# View W10
## W10
#########################################################################

# All samples

  
 HIV <-  mapply(CreateSeuratObject,
  counts = reads,
  project = names(reads),
  MoreArgs = list(min.cells = 3, min.features = 200)
)



# remove the original sparse matrices
rm(reads)
HIV <- merge(HIV[[1]],
  y = HIV[2:length(HIV)],
  HIV.cell.ids = names(HIV), 
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


# You can see the primary slots using:
## glimpse(W10_filtered_feature_bc_matrix.h5)
# to go back to the names layer
## Layers(W10_filtered_feature_bc_matrix.h5)

## W10_filtered_feature_bc_matrix.h5[["RNA"]]$counts |> head()
# or
## W10_filtered_feature_bc_matrix.h5@assays$RNA$counts  |> head()
# or
## LayerData(W10_filtered_feature_bc_matrix.h5, assay="RNA", layer='counts') |> head()


# The metadata in the Seurat object is located in adp@metadata and contains the
# information associated with each cell.
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

saveRDS(HIV, "C:/Users/Owner/Documents/scRNAseq/Code/Seurat/HIV PBMC/intermediate files/HIV_merged.rds")


library(tidyverse) # dplyr and ggplot2
library(Seurat) # Seurat toolkit
library(hdf5r) # for data import
library(patchwork) # for plotting
library(presto) # for differential expression
library(glmGamPoi) # for sctransform
library(ggplot2)

# read the object
## adp <- readRDS("C:/Users/Owner/Documents/github/Seurat test script/seurat 2 basic script/outputs/adp_merged.rds")


glimpse(HIV)


# calculate mitocondrial genes porcentage for quality control and normalization. Mitocondrial RNA increase need to be assesd 
#to determine if related to poor cell viability or in response of biological relevant processes.

HIV[["percent.mt"]] <- PercentageFeatureSet(HIV, pattern = "^MT-")
# set colors
cnames <- setNames(c("green","blue","orange","yellow","black0"), levels(factor(HIV@meta.data$orig.ident)))
cnames

# plot total counts per sample
VlnPlot(HIV, features = "nCount_RNA", layer = "counts", group.by = "orig.ident", raster = FALSE, alpha = 0.2) 
# or using ggplot2
HIV@meta.data %>%
  ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 650, color = "red", linetype = "dotted")

# determine number of features (genes)

plot_featur <- VlnPlot(HIV, features = "nFeature_RNA", group.by = "orig.ident") 
 

# visualize % mitocondrial genes

VlnPlot(HIV, features = "percent.mt", group.by = "orig.ident") +
  geom_hline(yintercept = 10, color = "red")

plot_count <- VlnPlot(HIV, features = "nCount_RNA", group.by = "orig.ident") +
  geom_hline(yintercept = 10, color = "red")
plot_featur|plot_count

# scatter the point by number of features and mitocondrial genes

FeatureScatter(HIV, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by = "orig.ident")

# and by number of fetures and total counts
FeatureScatter(HIV, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", log = TRUE)


# filter
# Set one set of parameters for Day 0 samples;
# keep the rownames (Cell barcodes)
F_23 <- HIV@meta.data |>
  filter(orig.ident=="6817423",
     nFeature_RNA > 600,
    percent.mt < 10
  ) |>
  tibble::rownames_to_column("Cell") |>
  pull(Cell)
F_24 <- HIV@meta.data |>
  filter(orig.ident=="6817424",
         nFeature_RNA > 300,
         percent.mt < 10
  ) |>
  tibble::rownames_to_column("Cell") |>
  pull(Cell)
F_25 <- HIV@meta.data |>
  filter(orig.ident=="6817425",
         percent.mt < 10
  ) |>
  tibble::rownames_to_column("Cell") |>
  pull(Cell)

F_26 <- HIV@meta.data |>
  filter(orig.ident=="6817426",
         nFeature_RNA > 300,
         percent.mt < 10
  ) |>
  tibble::rownames_to_column("Cell") |>
  pull(Cell)
F_27 <- HIV@meta.data |>
  filter(orig.ident=="6817427",
         nFeature_RNA > 700,
         percent.mt < 10
  ) |>
  tibble::rownames_to_column("Cell") |>
  pull(Cell)
F_28 <- HIV@meta.data |>
  filter(orig.ident=="6817428",
         nFeature_RNA > 900,
         percent.mt < 10
  ) |>
  tibble::rownames_to_column("Cell") |>
  pull(Cell)

keep <- c(F_23,F_24,F_25,F_26,F_27,F_28)
# use different parameters; established above
HIV_filt <- subset(HIV, cells = keep)
FeatureScatter(HIV_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", log = TRUE)
# save file after filtering

saveRDS(HIV_filt, "C:/Users/Owner/Documents/scRNAseq/Code/Seurat/HIV PBMC/intermediate files/intermediate files/HIV_merge_filt.rds")


##############################################################
## Normalization, variable feature  and scale               ##
##############################################################


# run sctransform
HIV_filt <- SCTransform(HIV_filt, vars.to.regress = "percent.mt", verbose = FALSE)

# Check default assay
############ DefaultAssay(object = adp_filt)

# Set default assay
########### DefaultAssay(object = adp_filt) <- "RNA"


# run PCA

HIV_filt <- RunPCA(HIV_filt, verbose = FALSE, assay = "SCT")

# visualizethte first 9 PC


DimHeatmap(HIV_filt, dims = 1:9, cells = 500, balanced = TRUE, ncol = 3)
# perform the elbow analisys
ElbowPlot(HIV_filt, ndims = 40)

# performing the neighbors analisys to prepara for the clusterin
HIV_filt <- FindNeighbors(HIV_filt, dims = 1:10)

# calculate the clusterin (suggessted redolution 0.4-1.2)

HIV_filt <- FindClusters(HIV_filt, resolution = 0.1)

# UMAP for the viasualization of the clusters

HIV_filt <- RunUMAP(HIV_filt, dims = 1:15)
DimPlot(HIV_filt,
  reduction = "pca", group.by = c("orig.ident", "seurat_clusters"),
  alpha = 0.2, ncol = 2
)
DimPlot(HIV_filt,
  reduction = "umap", group.by = c("orig.ident", "seurat_clusters"),
  alpha = 0.2, ncol = 2
)
# save the adp filtered file

saveRDS(HIV_filt, "C:/Users/Owner/Documents/scRNAseq/Code/Seurat/HIV PBMC/intermediate files//HIV_merge_filt_sctran_clust0.1.rds")

# prepare the SCT data for the search of markers among clusters
HIV_filt <- PrepSCTFindMarkers(HIV_filt, verbose = T)


# find markers (many method available, default a Wilcox rank sum. Look up test.use parameters for options)
# requires presto installation to speed up the next step
# devtools::install_github('immunogenomics/presto')
library(presto)
# find all markers
adp_filt_markers <- FindAllMarkers(adp_filt, only.pos = TRUE)
# ordering the results
adp_filt_markers <- adp_filt_markers %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj)

# examine a small subset
adp_filt_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
# markers visualization
top20 <- adp_filt_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup()
DoHeatmap(adp_filt, features = top20$gene) + NoLegend()

# associate pubblished markers with cell type
# contaminants
contam <- c("Ptprc", "Plp1", "Pecam1")
# beige adipocytes
beige <- c("Ucp1", "Ppargc1a", "Elovl3", "Cidea")
# preadipocytes
preadip <- c("Mmp3", "Cd142", "Itgb1")
# proliferating
prolif <- "Mki67"
# differentiating adipocytes
diffadip <- c("Col5a3", "Serpina3n")

VlnPlot(adp_filt, features = contam)

# visualize in the cluster the cel type corresponding to the  feature= parameter above (this case contam)
FeaturePlot(adp_filt, features = contam)

saveRDS(adp_filt_markers, "C:/Users/Owner/Documents/scRNAseq/Code/Seurat/HIV PBMC/intermediate files/adp_merge_filt_markers.rds")

##############################################
##    Differential Expression analisys      ##
##############################################
DefaultAssay(adp_filt)


library(SingleR) # for cell type annotation; Bioconductor
library(celldex) # for cell type annotation reference; Bioconductor
library(MAST) # for differential expression; Bioconductor

mem.maxVSize(vsize = 15000)
glimpse(adp_filt)
table(adp_filt$condition_tp)
Idents(adp_filt) <- "SCT_snn_res.0.1"
table(Idents(adp_filt))
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
plot(density(sample(JoinLayers(adp_filt@assays$RNA)$count["Gapdh", ], 2500)), cex = 0, lty = 1, main = "Density of Gapdh in 2500 random cells")

hist(sample(JoinLayers(adp@assays$RNA)$count["Gapdh", ], 2500), breaks = 52, main = "Histogram of Gapdh in 2500 random cells", ylab = "Frequency", xlab = "Gene counts")

###############################################################################
##        how estract expression matrixes                                    ##
##  expr<-GetAssayData(adp_filt,assay = "SCT",slot = "data")                 ##
##  expr_df <- as.data.frame(as.matrix(expr)) %>%                            ##
##   rownames_to_column(var = "gene") %>%                                    ##
##   pivot_longer(-gene, names_to = "cell", values_to = "expression")        ##
###############################################################################


adp_pp <- adp_filt

Idents(adp_pp) <- "SCT"
adp_pp_mks <- PrepSCTFindMarkers(adp_pp)
Idents(adp_pp_mks) <- "time_point"
table(Idents(adp_pp_mks))

## to find different expression betwe a label and the others in a metadata idents (if between two specific labels "idents.2=" can be used to specify.)
Day0_Day6_DE <- FindMarkers(adp_pp_mks, ident.1 = "Day 0", test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1)


Idents(adp_pp_mks) <- "SCT_snn_res.0.1"
DefaultAssay(adp_pp_mks) <- "SCT"

## find all markers between clusters

de_allClusters <- FindAllMarkers(adp_pp_mks, test.use = "wilcox", min.pct = 0.1, only.pos = TRUE)

head(de_allClusters)

## to check the most positivelly different markers per cluster
top5PerCluster <- matrix(ncol = 7)
colnames(top5PerCluster) <- colnames(de_allClusters)
for (i in 0:7) {
  top5PerCluster <- rbind(top5PerCluster, head(de_allClusters[which(de_allClusters$cluster == i), ], 5))
}
top5PerCluster <- top5PerCluster[-1, ]
top5PerCluster

DoHeatmap(adp_pp_mks, features = top5PerCluster$gene, slot = "scale.data")
## visualization of features

fig1 <- DimPlot(adp_pp_mks, group.by = "time_point")
fig2 <- FeaturePlot(adp_pp_mks, features = "Acta2", order = T)
fig3 <- FeaturePlot(adp_pp_mks, features = "Cd36", order = T)

fig1 / (fig2 | fig3)


## is possible to join all cells as a seample and do a pseudobulk analisys

pseudo_adp <- AggregateExpression(adp_pp_mks, assays = "SCT", return.seurat = T, group.by = c("orig.ident", "time_point", "condition", "condition_tp"))

head(pseudo_adp@assays$SCT$counts)
pseudo_adp@meta.data
glimpse(pseudo_adp)

# just to clean up the look a little bit
pseudo_adp <- RenameCells(pseudo_adp, new.names = gsub("_.*", "", pseudo_adp$orig.ident))
pseudo_adp$orig.ident <- gsub("_.*", "", pseudo_adp$orig.ident)
head(pseudo_adp@assays$RNA$counts)
pseudo_adp@meta.data

## performin bulk DE

Idents(pseudo_adp) <- "time_point"

## DESeq2 if not istalled:  BiocManager::install("DESeq2")
bulk_adp_de <- FindMarkers(pseudo_adp, ident.1 = "Day 0", ident.2 = "Day 6", test.use = "DESeq2")
head(bulk_adp_de)

# comparing how many differentially expressed genes between SC and bulk analisys comparing the conditions

scDE.genes <- rownames(Day0_Day6_DE)[which(Day0_Day6_DE$p_val_adj < 0.05)]
bulkDE.genes <- rownames(bulk_adp_de)[which(bulk_adp_de$p_val_adj < 0.05)]
length(scDE.genes)
length(bulkDE.genes)

## heck the common features between sc and bulk
length(intersect(scDE.genes, bulkDE.genes))
head(intersect(scDE.genes, bulkDE.genes), 30)

## to chech spefic features

bulk_adp_de[c("Acta2", "Cd36"), ]

## visualize the DE genes

Idents(adp_pp) <- "SCT_snn_res.0.1"
DotPlot(adp_pp, features = unique(top5PerCluster$gene), dot.scale = 3) + coord_flip()

# violine as alternative visualization

Idents(adp_pp) <- "time_point"
VlnPlot(adp_pp, features = c("Acta2", "Cd36"), alpha = 0.1)


## anotate the differential genes

markers <- c("Mmp3", "Mki67", "Fabp4", "Scd1", "Ucp1", "Ppargc1a", "Elovl3", "Cidea")
Idents(adp_pp) <- "SCT_snn_res.0.1"
avgExp <- AverageExpression(adp_pp, markers, assay = "SCT")$SCT
avgExp
DimPlot(adp_pp, label = T)
FeaturePlot(adp_pp, features = markers, ncol = 3, order = T)
## after find markers characteristic for a particular cluster we can annotate

adipocyte <- vector(length = ncol(adp_pp))
adipocyte[which(adp_pp$SCT_snn_res.0.1 %in% c(0, 5))] <- "Preadipocytes"
adipocyte[which(adp_pp$SCT_snn_res.0.1 %in% c(2, 6))] <- "Proliferating cells"
adipocyte[which(adp_pp$SCT_snn_res.0.1 %in% c(1, 3))] <- "Differentiating beige adipocytes"
adipocyte[which(adp_pp$SCT_snn_res.0.1 %in% c(4))] <- "Differentiated beige adipocytes"
adipocyte[which(adp_pp$SCT_snn_res.0.1 %in% c(7))] <- "Unclassified"
adp_pp$adipocyte <- adipocyte

f1 <- DimPlot(adp_pp, group.by = "SCT_snn_res.0.1", label = T) + NoLegend()
f2 <- DimPlot(adp_pp, group.by = "time_point") + NoLegend()
f3 <- DimPlot(adp_pp, group.by = "adipocyte", label = T) + NoLegend()
(f1 / f2 | f3)
## SingleR annotation. this method annotate each cell in the dataset against a reference dataset
library(celldex)

adp.sce <- as.SingleCellExperiment(adp_pp, assay = "SCT") # This selects *only* the SCT assay
mouseRNASeq <- celldex::MouseRNAseqData()
head(mouseRNASeq)
table(mouseRNASeq$label.main)
table(mouseRNASeq$label.fine)

annot <- SingleR::SingleR(test = adp.sce, ref = mouseRNASeq, labels = mouseRNASeq$label.main)
head(annot)

# if a label is to weak during SingleR the cell is taged as NA. Now we add the labels determine to the metadata of the Seurat object

table(annot$pruned.labels, useNA = "ifany") # useNA can be used turned on in the `table` function

adp_pp$mouseRNASeq.main <- annot$pruned.labels

## let's visualize the final result

annotFig1 <- DimPlot(adp_pp, group.by = "adipocyte", label = T) + NoLegend()
annotFig2 <- DimPlot(adp_pp, group.by = "mouseRNASeq.main", label = T)

annotFig1 | annotFig2

## annotation by cluster

Idents(adp_pp) <- "SCT_snn_res.0.1" # Assign clusters as the identities
avgExp <- AverageExpression(adp_pp, assays = "SCT")$SCT # Run AverageExpression on the SCT assay and return only SCT
clustAnnot <- SingleR::SingleR(test = avgExp, ref = mouseRNASeq, labels = mouseRNASeq$label.main) # Run SingleR on the averaged expression matrix
clustAnnot


clustLabels <- as.vector(clustAnnot$pruned.labels) # retrieve only the cluster-derived annotations
names(clustLabels) <- c(0:7) # assign the cluster numbers as the annotations
clustLabels.vect <- clustLabels[match(adp_pp$SCT_snn_res.0.1, names(clustLabels))] # match the cluster identities per cell in the Seurat data to the cluster labels
names(clustLabels.vect) <- colnames(adp_pp) # ensure that the cluster identities are assigned the cell names
adp_pp$mouseRNASeq.main.clust <- clustLabels.vect # add the cluster annotations to the vector

clustAnnotFig1 <- DimPlot(adp_pp, group.by = "SCT_snn_res.0.1", label = T) + NoLegend()
clustAnnotFig2 <- DimPlot(adp_pp, group.by = "adipocyte", label = T) + NoLegend()
clustAnnotFig3 <- DimPlot(adp_pp, group.by = "mouseRNASeq.main")
clustAnnotFig4 <- DimPlot(adp_pp, group.by = "mouseRNASeq.main.clust")

(clustAnnotFig1 | clustAnnotFig2) / (clustAnnotFig3 | clustAnnotFig4)

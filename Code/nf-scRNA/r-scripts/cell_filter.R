
library(dplyr) # dplyr and ggplot2
library(Seurat) # Seurat toolkit
library(hdf5r) # for data import
library(patchwork) # for plotting
library(presto) # for differential expression
library(glmGamPoi) # for sctransform
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
path <- args[1]

HIV <- readRDS(path)





# to save the Seurat object


#saveRDS(HIV, "/home/alexander-bontempo/Desktop/HIV GSM/RDS/HIV_merged.rds")

# read the object
## HIV <- readRDS("C:/Users/Owner/Documents/github/Seurat test script/seurat 2 basic script/outputs/HIV_merged.rds")



# calculate mitocondrial genes porcentage for quality control and normalization. Mitocondrial RNA increase need to be assesd 
#to determine if related to poor cell viability or in response of biological relevant processes.

HIV[["percent.mt"]] <- PercentageFeatureSet(HIV, pattern = "^MT-")
# set colors
samples <- length(unique(HIV@active.ident))

cnames <- setNames(sample(grDevices::colors(),samples), levels(factor(HIV@meta.data$orig.ident)))
cnames
#HIV <- NormalizeData(HIV)

# plot total counts per sample
counts_plot <- VlnPlot(HIV, features = "nCount_RNA",  group.by = "orig.ident", raster = FALSE, alpha = 0.2) 
png("nCount_RNA_per_sample.png",width = 1200, height = 900, res = 150)
print(counts_plot)
dev.off()
# or using ggplot2
counts_plot2 <- HIV@meta.data %>%
  ggplot(aes( color = orig.ident,x = nCount_RNA, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 650, color = "red", linetype = "dotted")
 

png("nCount2_RNA_per_sample.png", width = 1200, height = 900, res = 150)
print(counts_plot2)
dev.off()
# list all assays present
Assays(HIV)
# check what slots are in my RNA assay
slotNames(HIV[["RNA"]])



# determine number of features (genes)

plot_featur <- VlnPlot(HIV, features = "nFeature_RNA", group.by = "orig.ident") +
  geom_hline(yintercept = 2500, color = "red")
show(plot_featur)
png("nFeature_RNA_per_sample.png", width = 1200, height = 900, res = 150)
print(plot_featur)
dev.off()

# visualize % mitocondrial genes

plot_mt <- VlnPlot(HIV, features = "percent.mt", group.by = "orig.ident") +
  geom_hline(yintercept = 10, color = "red")
show(plot_mt)
png("percent_mt_per_sample.png", width = 1200, height = 900, res = 150)
print(plot_mt)
dev.off()

plot_count <- VlnPlot(HIV, features = "nCount_RNA", group.by = "orig.ident") +
  geom_hline(yintercept = 10, color = "red")
plot_featur|plot_mt

# scatter the point by number of features and mitocondrial genes

FeatureScatter(HIV, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by = "orig.ident")

# and by number of fetures and total counts
FeatureScatter(HIV, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", log = TRUE)


# filter
# Set one set of parameters for Day 0 samples;
# keep the rownames (Cell barcodes)
# F_30 <- HIV@meta.data |>
#   filter(orig.ident=="GSM6817423",
#      nFeature_RNA > 500,
#      nFeature_RNA < 4000,
#     percent.mt < 10
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)
# F_31 <- HIV@meta.data |>
#   filter(orig.ident=="GSM6817431",
#          nFeature_RNA > 500,
#          percent.mt < 10,
#          nFeature_RNA < 4000
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)
# F_32 <- HIV@meta.data |>
#   filter(orig.ident=="GSM6817432",
#          percent.mt < 10,
#          nFeature_RNA > 500,
#         nFeature_RNA < 4000
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)

# F_33 <- HIV@meta.data |>
#   filter(orig.ident=="GSM6817433",
#          nFeature_RNA > 500,
#          percent.mt < 10
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)
# F_34 <- HIV@meta.data |>
#   filter(orig.ident=="GSM6817434",
#          nFeature_RNA > 700,
#          percent.mt < 10
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)
# F_35 <- HIV@meta.data |>
#   filter(orig.ident=="GSM6817435",
#          nFeature_RNA > 500,
#          nFeature_RNA < 4000,
#          percent.mt < 10
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)
# F_36 <- HIV@meta.data |>
#   filter(orig.ident=="GSM6817436",
#          nFeature_RNA > 500,
#          nFeature_RNA < 4000,
#          percent.mt < 10
#   ) |>
#   tibble::rownames_to_column("Cell") |>
#   pull(Cell)

cells_keep <- HIV@meta.data |>
  filter(nFeature_RNA > 500,
    #  nFeature_RNA < 2500,     # since scrublet is applied, upper limit for nFeature_RNA can be removed
    percent.mt < 10
  ) |>
  tibble::rownames_to_column("Cell") |>
  pull(Cell)

# keep <- c(F_30,F_31,F_32,F_33,F_34,F_35,F_36)
# use different parameters; established above
HIV_filt <- subset(HIV, cells = cells_keep)
features <- FeatureScatter(HIV_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", log = TRUE)
png("filtered_nCount_vs_nFeature.png", width = 1200, height = 900, res = 150)
print(features)
dev.off()
saveRDS(HIV_filt,"HIV_filtered.rds")

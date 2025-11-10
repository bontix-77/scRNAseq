library(Seurat)

args <- commandArgs(trailingOnly=TRUE)
path <- args[1]

HIV <- readRDS(path)
HIV_SCT <- SCTransform(HIV, vars.to.regress = "percent.mt", verbose = FALSE)


saveRDS(HIV_SCT, "HIV_SCTransform.rds")
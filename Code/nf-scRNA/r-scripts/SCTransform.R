library(Seurat)
#if memory problem:
library(future)
library(glmGamPoi)
plan("sequential") 
##################
set.seed(1997)
args <- commandArgs(trailingOnly=TRUE)
path <- args[1]

HIV <- readRDS(path)
HIV_SCT <- SCTransform(HIV, 
  vars.to.regress = "percent.mt",
  verbose = TRUE,
  # conserve.memory = TRUE,   #to save memory activate this 
  # variable.features.n = 1500,   #to save memory activate this 
  return.only.var.genes = TRUE,
  method = "glmGamPoi",
  do.correct.umi = FALSE           
      
)
saveRDS(HIV_SCT, "HIV_SCTransform.rds")
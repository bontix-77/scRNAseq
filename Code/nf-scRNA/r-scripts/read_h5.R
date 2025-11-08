library(tidyverse)
library(Seurat)
args <- commandArgs(trailingOnly=TRUE)


folders <- strsplit(args[1], " ")[[1]]

if (!grepl("/$", args[2])) args[2] <- paste0(args[2], "/")
setwd(args[2])

reads <- lapply(folders, Read10X)
names(reads) <- folders
#writeLines(names(reads), con="results.txt")

# sink("results.txt")
# cat("Working dir:", getwd(), "\n")
# cat("Folders:\n")
# print(folders)
# cat("\nFiles in first folder:\n")
# print(list.files(folders[1]))
#sink()

HIV <-  mapply(CreateSeuratObject,
  counts = reads,
  project = names(reads),
  MoreArgs = list( min.cells = 2,min.features = 200)
)


rm(reads)
HIV <- merge(HIV[[1]],
  y = HIV[2:length(HIV)],
  HIV.cell.ids = names(HIV), 
  project = "HIV PBMC on patient under ART"
)

saveRDS(HIV, "HIV_merged.rds")

library(Seurat)
args <- commandArgs(trailingOnly = TRUE)


folders <- strsplit(args[1], " ")[[1]]
base_dir <- args[2]
if (!grepl("/$", base_dir)) {
  base_dir <- paste0(base_dir, "/")
}


reads <- lapply(folders, function(x) Read10X(file.path(base_dir, x)))
names(reads) <- folders
#writeLines(names(reads), con="results.txt")

# sink("results.txt")
# cat("Working dir:", getwd(), "\n")
# cat("Folders:\n")
# print(folders)
# cat("\nFiles in first folder:\n")
# print(list.files(folders[1]))
#sink()

HIV <- mapply(
  CreateSeuratObject,
  counts = reads,
  project = names(reads),
  MoreArgs = list(min.cells = 2, min.features = 200)
)


rm(reads)
# assign metadata group
HIV$`6817423`$group <- "detectable"
HIV$`6817431`$group <- "undetectable"
HIV$`6817432`$group <- "undetectable"
HIV$`6817433`$group <- "detectable"
HIV$`6817434`$group <- "undetectable"
HIV$`6817435`$group <- "undetectable"
HIV$`6817436`$group <- "detectable"

HIV <- merge(
  HIV[[1]],
  y = HIV[2:length(HIV)],
  add.cell.ids = names(HIV),
  project = "HIV PBMC on patient under ART"
)

saveRDS(HIV, "HIV_merged.rds")

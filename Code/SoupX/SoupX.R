library(SoupX)
library(Seurat)

library(SoupX)

# Use the path to your .h5 file (likely filtered)
sc <- load10X_h5("C:/Users/Owner/Documents/github/Seurat test script/seurat 2 basic script/R project dowloaded from manual/GettingStarted_scRNASeq/data")

# Proceed to estimate contamination
sc <- autoEstCont(sc)

# Get decontaminated counts
out <- adjustCounts(sc)










# Create a list of count matrices
files <- list.files(path = "C:/Users/Owner/Documents/github/Seurat test script/seurat 2 basic script/R project dowloaded from manual/GettingStarted_scRNASeq/data", recursive = T, pattern = "*.h5")


h5_read <- lapply(paste0(
  "C:/Users/Owner/Documents/github/Seurat test script/seurat 2 basic script/R project dowloaded from manual/GettingStarted_scRNASeq/data/",
  files
), Read10X_h5)
names(h5_read) <- c("D10", "D16", "D20", "D26", "W10", "W16", "W20", "W26")

adp <- mapply(CreateSeuratObject,
              counts = h5_read,
              project = names(h5_read),
              MoreArgs = list(min.cells = 3, min.features = 200)
)

adp <- merge(adp[[1]],
             y = adp[2:length(adp)],
             add.cell.ids = names(adp), project = "Adipose"
)

sc <- SoupChannel(tod = adp, toc = adp)

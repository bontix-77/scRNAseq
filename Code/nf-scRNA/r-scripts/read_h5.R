library(Seurat)
args <- commandArgs(trailingOnly = TRUE)


folders <- strsplit(args[1], " ")[[1]]
base_dir <- args[2]
if (!grepl("/$", base_dir)) {
  base_dir <- paste0(base_dir, "/")
}


reads <- lapply(folders, function(x) Read10X(file.path(base_dir, x)))
names(reads) <- folders

# function to colapse TCR alfa bete delta and gama in 4 groups
if (args[3] =="yes"){
                library(Matrix)
                marker_T <- c()
                Alfa <- c(
                  "TRAV1-1",
                  #"TRAV1-2",
                  "TRAV2",
                  "TRAV3",
                  "TRAV4",
                  "TRAV5",
                  "TRAV6",
                  "TRAV8-1",
                  "TRAV10",
                  "TRAV12-1",
                  "TRAV8-2",
                  "TRAV8-3",
                  "TRAV13-1",
                  "TRAV12-2",
                  "TRAV8-4",
                  "TRAV8-5",
                  "TRAV13-2",
                  "TRAV14DV4",
                  "TRAV9-2",
                  "TRAV12-3",
                  "TRAV8-6",
                  "TRAV16",
                  "TRAV17",
                  "TRAV18",
                  "TRAV19",
                  "TRAV20",
                  "TRAV21",
                  "TRAV22",
                  "TRAV23DV6",
                  "TRAV24",
                  "TRAV25",
                  "TRAV26-1",
                  "TRAV27",
                  "TRAV29DV5",
                  "TRAV30",
                  "TRAV26-2",
                  "TRAV34",
                  "TRAV35",
                  "TRAJ53",
                  "TRAJ52",
                  "TRAJ51",
                  "TRAJ50",
                  "TRAJ49",
                  "TRAJ48",
                  "TRAJ47",
                  "TRAJ46",
                  "TRAJ45",
                  "TRAJ44",
                  "TRAJ43",
                  "TRAJ42",
                  "TRAJ41",
                  "TRAJ40",
                  "TRAJ39",
                  "TRAJ38",
                  "TRAJ37",
                  "TRAJ36",
                  "TRAJ35",
                  "TRAJ34",
                  "TRAJ33",
                  "TRAJ32",
                  "TRAJ31",
                  "TRAJ30",
                  "TRAJ29",
                  "TRAJ28",
                  "TRAJ27",
                  "TRAJ26",
                  "TRAJ25",
                  "TRAJ24",
                  "TRAJ23",
                  "TRAJ22",
                  "TRAJ21",
                  "TRAJ20",
                  "TRAJ19",
                  "TRAJ18",
                  "TRAJ17",
                  "TRAJ16",
                  "TRAJ14",
                  "TRAJ13",
                  "TRAJ12",
                  "TRAJ11",
                  "TRAJ10",
                  "TRAJ9",
                  "TRAJ8",
                  "TRAJ1",
                  "TRAC",
                  "TRAV7",
                  "TRAV9-1",
                  "TRAV11",
                  "TRAV15",
                  "TRAV8-7",
                  "TRAV28",
                  "TRAV31",
                  "TRAV32",
                  "TRAV33",
                  "TRAV37")
                Beta <- c(
                  "TRBV1",
                  "TRBV2",
                  "TRBV3-1",
                  "TRBV4-1",
                  "TRBV5-1",
                  "TRBV6-1",
                  "TRBV4-2",
                  "TRBV6-2",
                  "TRBV7-2",
                  "TRBV6-4",
                  "TRBV7-3",
                  "TRBV5-3",
                  "TRBV9",
                  "TRBV10-1",
                  "TRBV11-1",
                  "TRBV10-2",
                  "TRBV11-2",
                  "TRBV6-5",
                  "TRBV7-4",
                  "TRBV5-4",
                  "TRBV6-6",
                  "TRBV5-5",
                  "TRBV7-6",
                  "TRBV5-6",
                  "TRBV7-7",
                  "TRBV5-7",
                  "TRBV7-9",
                  "TRBV13",
                  "TRBV10-3",
                  "TRBV11-3",
                  "TRBV12-3",
                  "TRBV12-4",
                  "TRBV12-5",
                  "TRBV14",
                  "TRBV15",
                  "TRBV18",
                  "TRBV19",
                  "TRBV20-1",
                  "TRBV21-1",
                  "TRBV23-1",
                  "TRBV24-1",
                  "TRBV25-1",
                  "TRBV27",
                  "TRBV28",
                  "TRBV29-1",
                  "TRBC1",
                  "TRBC2",
                  "TRBV30",
                  "TRBV16",
                  "TRBV7-5",
                  "TRBV6-8",
                  "TRBV7-1",
                  "TRBV8-1",
                  "TRBV5-2",
                  "TRBV8-2",
                  "TRBV12-1",
                  "TRBV12-2",
                  "TRBV6-7",
                  "TRBV17",
                  "TRBV22-1",
                  "TRBVA",
                  "TRBV26",
                  "TRBVB",
                  "TRBV20OR9-2",
                  "TRBV21OR9-2",
                  "TRBV22OR9-2",
                  "TRBV23OR9-2",
                  "TRBV24OR9-2",
                  "TRBV25OR9-2",
                  "TRBV26OR9-2",
                  "TRBV29OR9-2",
                  "TRBJ1-1",
                  "TRBJ1-2",
                  "TRBJ1-3",
                  "TRBJ1-4",
                  "TRBJ1-5",
                  "TRBJ1-6",
                  "TRBJ2-1",
                  "TRBJ2-2",
                  "TRBJ2-2P",
                  "TRBJ2-3",
                  "TRBJ2-4",
                  "TRBJ2-5",
                  "TRBJ2-6",
                  "TRBJ2-7")
                  Gamma <- c(
                  "TRGC2",
                  "TRGJ2",
                  "TRGJP2",
                  "TRGC1",
                  "TRGJ1",
                  "TRGJP",
                  "TRGJP1",
                  "TRGV11",
                  "TRGVB",
                  "TRGV10",
                  #"TRGV9",
                  "TRGVA",
                  "TRGV8",
                  "TRGV7",
                  "TRGV6",
                  "TRG-AS1",
                  "TRGV5P",
                  "TRGV5",
                  "TRGV4",
                  "TRGV3",
                  "TRGV2",
                  "TRGV1")
                Delta <- c(
                  "TRDV1",
                  "TRDV2",
                  "TRDV3",
                  "TRDD1",
                  "TRDD2",
                  "TRDD3",
                  "TRDJ1",
                  "TRDJ2",
                  "TRDJ3",
                  #  "TRDC",
                  "TRDJ4"
                
                )

                grep("TRD", reads[[1]]@Dimnames[[1]], value = T)


                group_genes_and_collapse <- function(
                  mtx,
                  remove,
                  new_feature_name
                ) {
                  # extracting the matrix 
                  # mtx <- GetAssayData(seu)#, assay = assay, slot = slot_use)

                  # tengo solo i geni che esistono davvero nella matrice
                  genes_present <- intersect(remove, rownames(mtx))
                  if (length(genes_present) == 0) {
                    stop("Nessuno dei geni indicati esiste nell assay")
                  }

                  # somma per colonna dei geni da accorpare
                  grouped_counts <- Matrix::colSums(mtx[genes_present, , drop = FALSE])

                  # aggiungo la nuova "feature" come ultima riga
                  # costruisco una matrice una riga con i counts raggruppati
                  new_row <- Matrix(
                    grouped_counts,
                    nrow = 1,
                    sparse = T,
                    dimnames = list(new_feature_name, colnames(mtx))
                  )

                  mtx_new <- rbind(mtx, new_row)

                  # elimino i geni originali
                  keep_features <- setdiff(rownames(mtx_new), genes_present)
                  mtx_new <- mtx_new[keep_features, ]

                  
                  return(mtx_new)
                }

                TCR <- list(`TCR-A`=Alfa,`TCR-B`=Beta,`TCR-D`=Delta,`TCR-G`=Gamma)

                reads1 <- list()

                for (i in names(reads)) {
                  mtx_i <- reads[[i]]   # parto dalla matrice originale
                  
                  for (g in names(TCR)) {
                    mtx_i <- group_genes_and_collapse(
                      mtx = mtx_i,
                      remove = TCR[[g]],
                      new_feature_name = g
                    )
                  }
                  
                  reads[[i]] <- mtx_i  # salvo la versione finale collassata
                }
                  print ("TCR colapsed")
                }


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
# Include the group metadata to samples as determined in the nextflow.config parameter samples.
samples <- args[4]
samples <- paste0("list(",samples,")")
samples <- eval(parse(text=samples))
for (item in names(samples)){
HIV[[item]]$group <- samples[[item]]
print (paste0("sample ",item, " is assigned to the group ",samples[[item]]))

}

HIV <- merge(
  HIV[[1]],
  y = HIV[2:length(HIV)],
  add.cell.ids = names(HIV),
  project = "HIV PBMC on patient under ART"
)

saveRDS(HIV,"HIV_merged.rds")

library(CellChat)# Load the CellChat library
library(Seurat)
Data_set <- readRDS( "/home/alexander-bontempo/Desktop/GitHub/scRNAseq/adp_merge_filt_sctran_clust_harmony.rds") # Load the dataset containing adipocytes")
#Data.cc <- GetAssayData(Data_set, assay = "SCT", layer = "data") #Data_set[["SCT"]]$data
#labels <- Idents(Data_set)
#meta <- data.frame(samples = labels, row.names = names(labels)) # create a dataframe of the cell labels
# Subset the whole Seurat object in 4 subset defined by the treatment status an 2 time point using the condition_tp metadata 
Idents(Data_set) <- "condition_tp"
for (i in unique(Data_set$condition_tp)) {
  sub_set <-   subset(Data_set, subset = condition_tp== i)
  assign(i,sub_set)

}


#performing the CellChat analysis on the DKO group at day 6.CellChat permits to compare more groups if a differential analysis is required.
`DKO Day 6`$samples <- as.factor(`DKO Day 6`$orig.ident)
`DKO Day 6`$seurat_clusters <- factor( as.numeric(as.character(`DKO Day 6`$seurat_clusters)) + 1)
cellchat <- createCellChat(object = `DKO Day 6`, group.by = "seurat_clusters",assay="SCT") # create a CellChat object

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
 CellChatDB.use <- subsetDB(CellChatDB)

 # use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat@DB <- CellChatDB.use



# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#> The number of highly variable ligand-receptor pairs used for signaling inference is 692


cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways


#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# there are evidenced for cell-cell communication between cluste 1 vs clusters 2 and 3.

#At this point we can anlyze specific interacctions or pathways.
pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network","in DKO gropu at 6 days"))

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair

netAnalysis_contribution(cellchat, signaling = pathways.show)


#We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. We provide a function `extractEnrichedLR` to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.  

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord",save=T)

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways

# check the order of cell identity to set suitable vertex.receiver
# the following code doesn't work. it return an error i coudn't figure out why.
#levels(cellchat@idents)
#vertex.rece = rep(TRUE, 5)      # all TRUE initially
#vertex.rece[1] = FALSE 
#for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
#  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.rece, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
#}

## Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
CellChat can also show all the significant interactions mediated by L-R pairs and signaling pathways, and interactions provided by users from some cell groups to other cell groups using the function netVisual_bubble (option A) and netVisual_chord_gene (option B).

### (A)  Bubble plot
#We can also show all the significant interactions (L-R pairs) from some cell groups to other cell groups using `netVisual_bubble`.
# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
# (2) show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
# (3) show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(1,2), pairLR.use = pairLR.use, remove.isolate = TRUE)


# set the order of interacting cell pairs on x-axis
# (4) Default: first sort cell pairs based on the appearance of sources in levels(object@idents), and then based on the appearance of targets in levels(object@idents)
# (5) sort cell pairs based on the targets.use defined by users
#netVisual_bubble(cellchat, targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.target = T)
# (6) sort cell pairs based on the sources.use defined by users
#netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T)
# (7) sort cell pairs based on the sources.use and then targets.use defined by users
#netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T)
# (8) sort cell pairs based on the targets.use and then sources.use defined by users
#netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE)

#Similar to Bubble plot, CellChat provides a function `netVisual_chord_gene` for drawing Chord diagram to

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3), targets.use = c(4,5), legend.pos.x = 15)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3), targets.use = c(4,5), signaling = c("CCL","CXCL"),legend.pos.x = 8)
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5), slot.name = "netP", legend.pos.x = 10)CellChat can plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathways using a Seurat wrapper function `plotGeneExpression` if the Seurat R package has been installed. This function provides three types of visualizztion, including "violin", "dot", "bar". Alternatively, users can extract the signaling genes related to the inferred L-R pairs or signaling pathway using `extractEnrichedLR`, and then plot gene expression using Seurat or other packages.

#CellChat can plot the gene expression distribution of signaling genes related to L-R pairs or signaling 
# pathways using a Seurat wrapper function `plotGeneExpression` if the Seurat R package has been installed.
#  This function provides three types of visualizztion, including "violin", "dot", "bar". 
# Alternatively, users can extract the signaling genes related to the inferred L-R pairs or signaling pathway
#  using `extractEnrichedLR`, and then plot gene expression using Seurat or other packages.


plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = TRUE, type = "violin")

plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)

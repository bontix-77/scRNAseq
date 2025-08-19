
#Preparation: if neccesary run the following:
#devtools::install_github("ctlab/fgsea")
#BiocManager::install("clusterProfiler")
#BiocManager::install("corg.Mm.eg.er")  for mouse use: org.Hs.eg.db
#BiocManager::install("enrichplot")
#1. Load necessary packages and data:

library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db) # Or appropriate organism package human: org.Hs.eg.db
library(enrichplot) # For visualization

#2. Load your Seurat object:

# Assuming your Seurat object is named 'adp'
adp <- readRDS("/home/alexander-bontempo/Desktop/GitHub/scRNAseq/adp_merge_filt_sctran_clust_harmony.rds")

#3. Identify marker genes:

 #   Use FindAllMarkers or FindMarkers to identify differentially expressed genes for each cluster.
  #  For example: 
     adp@active.assay <- "SCT"
     Idents(adp) <- "condition_tp"
     # Find markers for all clusters
     all.markers <- FindAllMarkers(adp, #only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
     # Find markers for cluster 1 compared to all other cells
     adp$condition_tp <- as.factor(adp$condition_tp)
     cluster1.markers <- FindMarkers(adp, ident.1 = "DKO Day 0",ident.2 = "DKO Day 6", min.pct = 0.25)

#4. Prepare gene lists for enrichment analysis: 

  #  Extract the gene names of interest from the marker gene results. For example, if you used FindAllMarkers, you can filter for the top markers:

    

      top500<- all.markers %>% group_by(cluster) %>% dplyr::top_n(n = 100, wt = avg_log2FC)
    #Convert gene symbols to Entrez IDs using bitr from clusterProfiler. This requires specifying the organism database (e.g., org.Hs.eg.db for human). 



     # Example using top10 from above
     gene.df <- bitr(top500$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = TRUE)
     all.genes <- rownames(adp@assays$SCT)
     universe  <- bitr(all.genes,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = TRUE)
    
#5. Perform GO enrichment analysis: 

    #se enrichGO with the converted Entrez IDs, background genes (optional, but recommended), and the desired ontology (BP, CC, or MF). 



     enrichGO_results <- enrichGO(gene          = gene.df$ENTREZID,
                                  OrgDb         = org.Mm.eg.db,
                                  ont           = "BP", # Or CC or MF or BP
                                  pAdjustMethod = "BH",
                                  #keyType       = "ENTREZID",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.05,
                                  universe      = universe$ENTREZID, # Optional: background genes
                                  readable      = TRUE)

#6. Visualize the results:

  #  Use various plotting functions from enrichplot to visualize the enriched GO terms. 



     # Bar plot
     barplot(enrichGO_results, showCategory=20)
     # Dot plot
     dotplot(enrichGO_results, showCategory=20)
     # Enrichment map
     emapplot(enrichGO_results, showCategory = 20)
     # Cnet plot
     rnd_genes <- universe$ENTREZID[sample(1:1000, 30,replace = FALSE)]
     cnetplot(enrichGO_results, categorySize="pvalue", max.ovelaps=200,color.params = list(foldChange = rnd_genes))

#Example using a pre-defined gene list:
#If you have a list of genes you want to analyze, you can skip steps 3 and 4 (gene identification and conversion) and directly use the gene list in enrichGO.

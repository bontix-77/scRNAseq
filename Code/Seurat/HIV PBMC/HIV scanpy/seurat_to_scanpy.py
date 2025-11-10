
'''
tutorial:
https://scanpy.readthedocs.io/en/latest/tutorials/index.htmlhttps://scanpy.readthedocs.io/en/latest/tutorials/index.html
'''

import scanpy as sc
#import louvain    #library can't install: missing dependencices i can't solve

adata = sc.read_h5ad("/home/alexander-bontempo/Desktop/HIV GSM/h5ad/data.h5ad")

######## change gene names
# Load gene namessc.pl.dotplot(
#    adata,
#    marker_genes,
#    groupby="leiden_res_0.02",
#    standard_scale="var",
#    use_raw=False
#)
with open("/home/alexander-bontempo/Desktop/HIV GSM/h5ad/gene_names.txt") as f:
    gene_names = [line.strip() for line in f]

# Replace raw.var_names
adata.raw.var_names = gene_names

# Optional: replace var_names too
adata.var_names = gene_names


########################################################################
# Find highly variable genes (recommended before PCA)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# PCA
sc.tl.pca(adata, svd_solver='arpack')



# Compute neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

#louvain clustering
sc.tl.louvain(adata, resolution=0.5)


# Leiden clustering
sc.tl.leiden(adata, resolution=0.5)
# 7️⃣ Compute UMAP embedding
sc.tl.umap(adata)
# Plot UMAP
sc.pl.umap(adata, color='leiden')



adata.obs['leiden'].value_counts()


#cluster markers
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')
# Get raw array
names_array = adata.uns['rank_genes_groups']['names']






#check top markers
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

# Top 5 genes for each cluster
adata.uns['rank_genes_groups']['names']

# Scores
adata.uns['rank_genes_groups']['scores']

# P-values
adata.uns['rank_genes_groups']['pvals']
adata.uns['rank_genes_groups']['logfoldchanges']

adata.var_names[:10]



########## grid UMAP#####################
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph")


sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)


#cluster markers
sc.tl.rank_genes_groups(adata, groupby='leiden_res_1', method='wilcoxon')

sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res_1", standard_scale="var", n_genes=5,use_raw=False)
############ cell type#######################3
# to look if a gene is present

# to look if a gene is present

[x for x in adata.var_names if "GYPA" in x]

adata.rank_genes_groups['names']

marker_genes = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    # Note: DMXL2 should be negative
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Erythroblast": ["MKI67",  "HBB"],
    # Note HBM and GYPA are negative markers
    "Proerythroblast": ["CDK6", "SYNGR1"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
    "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
    "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
    # Note IGHD and IGHM are negative markers
    "B cells": [
        "MS4A1",
        "ITGB1",
        "PRDM1",
        "IRF4",
        "PAX5",
        "BCL11A",
        "BLK",
        "IGHD",
        "IGHM",
    ],
    "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
    # Note PAX5 is a negative marker
    "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
}



sc.pl.dotplot(
    adata,
    marker_genes,
    groupby="leiden_res_1",
    standard_scale="var",
    use_raw=False
)


adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_1"].map(
    {
        "0": "T cells double positive",
        "1": "T naive",
        "2": "T cells",
        "3": "T cells",
        "4": "T CD8+",
        "5": "no determined",
        "6": "naive CD20+ B",
        "7" : "CD14+ Mono",
        "8": "CD16+ Mono",
        "9":"plasmacells",
        "10": "no determined",
        "11":"no detemined"

    }
)
# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(adata, groupby="leiden_res_1", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res_1", standard_scale="var", n_genes=5, use_raw=False)

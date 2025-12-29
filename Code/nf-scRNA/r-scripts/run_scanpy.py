"""
tutorial:
https://scanpy.readthedocs.io/en/latest/tutorials/index.htmlhttps://scanpy.readthedocs.io/en/latest/tutorials/index.html
"""

import scanpy as sc
import sys
# import louvain    #library can't install: missing dependencices i can't solve

path_data= sys.argv[1]
adata = sc.read_h5ad(path_data)


# for some reason the indexes in the var is not the gene name but a numeric. Thefollowin
adata.raw.var.set_index(adata.raw.var["_index"], inplace=True)

########################################################################
# Find highly variable genes (recommended before PCA)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# PCA
sc.tl.pca(adata, svd_solver="arpack")


# Compute neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# louvain clustering
# sc.tl.louvain(adata, resolution=0.5)


# Leiden clustering
sc.tl.leiden(adata, resolution=1.4, flavor="igraph")
# Compute UMAP embedding
sc.tl.umap(adata)
# Plot UMAP
sc.pl.umap(adata, color="leiden", show=False, save="_leiden_res_1.4.png")



# how many cells per cluster?
adata.obs["leiden"].value_counts()


# cluster markers
sc.tl.rank_genes_groups(adata, groupby="leiden", method="t-test")
"""
UPDATE: indexing solution below was partial. complete solution has been introduced by line 15
For some reason adata.uns['rank_genes_groups']['names'] derived from runing seuratData
doesn't contain gene name but the adata.raw.var['_index'].
To convert the indexes to gene name i had to recreate a rec assigning dtype (one for cluster).

"""
# # names_raw=adata.uns['rank_genes_groups']['names']
# names_raw=adata.raw.var_names


# names_array=[]

# for i,v in enumerate(names_raw):
#     partial=[]
#     for j in range(len(v)):
#         if names_raw[i][j].isdigit():
#             partial.append(adata.raw.var['_index'].iloc[int(names_raw[i][j]),])
#         else:
#             partial.append(names_raw[i][j])
#     names_array.append(partial)


# # 1. i took the group names from anothe dtype in the rank_genes_groups array
# group_names = adata.uns['rank_genes_groups']['scores'].dtype.names
# print("group_names:", group_names)

# # 2. list to numpy array
# names_np = np.array(names_array)
# print("names_array shape:", names_np.shape)
# print("n groups in dtype:", len(group_names))

# # 3. build column based on shape
# if names_np.ndim == 1:

#     raise ValueError("names_array ha una sola dimensione, mi aspetto 2D (genes x groups o groups x genes)")

# if names_np.shape[1] == len(group_names):

#     cols = [names_np[:, i] for i in range(names_np.shape[1])]
# elif names_np.shape[0] == len(group_names):

#     cols = [names_np[i, :] for i in range(names_np.shape[0])]
# else:
#     raise ValueError(
#         f"Shape di names_array {names_np.shape} non compatibile con numero gruppi {len(group_names)}"
#     )

# # 4. create the rec usin cols as arrays and group names
# names_rec = np.rec.fromarrays(cols, names=group_names)

# print(type(names_rec), names_rec.shape, names_rec.dtype)

# # 5. change to adata
# adata.uns['rank_genes_groups']['names'] = names_rec
# ############################################


# check top markers
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

# Top 5 genes for each cluster
adata.uns["rank_genes_groups"]["names"]

# Scores
adata.uns["rank_genes_groups"]["scores"]

# P-values
adata.uns["rank_genes_groups"]["pvals"]
adata.uns["rank_genes_groups"]["logfoldchanges"]

adata.var_names[:10]


########## grid UMAP#####################
for res in [0.02, 0.5, 1.4]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )


sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_1.40"],
    legend_loc="on data",
)


# cluster markers
sc.tl.rank_genes_groups(adata, groupby="leiden_res_1.40", method="wilcoxon")

sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_res_1.40", standard_scale="var", n_genes=5, use_raw=True
)
############ cell type#######################3
# to look if a gene is present

[x for x in adata.var_names if "GYPA" in x]

adata.uns["rank_genes_groups"]["names"]

marker_genes = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    # Note: DMXL2 should be negative
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Erythroblast": [ "HBB"],#"MKI67"], gene not found in the set
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
    "CD4+ T": ["CD4", "IL7R"],#, "TRBC2"], # TRBC2 not found in raw
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
}

marker_genes_extend = {
    "platelets" : ["PTPRC", "PECAM1"], #, "Plp1"   no present in DF

    "Naive_B" : ["CD19", "IGHD"],

    "Memory_B" : ["CD19", "CD27", "IGHG1", "MS4A1"],

    "CD14_monocytes" : ["CD14", "S100A8", "CD4"],

    "CD16_monocytes" : ["FCGR3A", "CD4"],

    "NK" : ["GZMA", "GZMB", "GNLY", "NKG7", "FCGR3A", "CD8A"],

    "CD56bright_NK" : ["NCAM1", "TCR-D", "GZMK", "GATA3"],

    "CD4_naive_T" : ["CD3E", "CD4", "CCR7", "IL7R", "TCF7", "PTPRC", "CD27"],

    "CD4_Tcm" : ["CD3E", "CD4", "CD27", "SELL", "CCR7", "PTPRC"],

    "CD4_Tem" : ["CD3E", "CD4", "CD27", "CCR7", "PTPRC", "SELL", "TCF7"],

    "CD4_Trm" : ["CD3E", "CD4", "CD69", "ITGAE"],

    "CD4_Tscm" : ["CD3E", "CD4", "CD27", "SELL", "PTPRC", "IL2RB", "FAS"],

    "CD8_naive_T" : ["CD3E", "CD8A", "CCR7", "TCF7", "PTPRC", "CD27"],

    "CD8_Tcm" : ["CD3E", "CD8A", "CD27", "SELL", "CCR7", "PTPRC"],

    "CD8_Tem" : ["CD3E", "CD8A", "CD27", "CCR7", "PTPRC", "SELL", "TCF7"],

    "CD8_Trm" : ["CD3E", "CD8A", "CD69", "ITGAE"],

    "CD8_Tscm" : ["CD3E", "CD8A", "CD27", "SELL", "PTPRC", "IL2RB", "FAS"],

    "MAIT" : ["CD3E", "TRAV1-2", "IL7R", "GZMK", "CCR6"],

    "Vγ9Vδ2_T" : ["CD3E",  "TCR-G", "TCR-D", "GZMA", "CCL5", "TRDC"],

    "gd_T" : ["CD3E", "TCR-G",  "TRDC", "TCR-D"],

    "Treg" : ["CD3E", "FOXP3", "IL2RA"],

    "T_proliferating" : ["CD3E", "MKI67"],

    "mDC" : ["CST3", "CD1C", "HLA-DRA", "CD4"],

    "pDC" : ["CST3", "CLEC4C", "CXCR3", "IL3RA", "GZMB", "CD4"],

    "Plasmablasts" : ["JCHAIN", "CD27", "MKI67", "IGHD", "IGHG1"],

    "HSC" : ["PPBP", "ITGA2B"],

    "RBC" : ["HBB", "HBA1", "HBA2"]
}

sc.pl.dotplot(
    adata, marker_genes_extend, groupby="leiden_res_1.40", standard_scale="var", use_raw=True)

adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_1.40"].map(
    {
        "0": "T cells double positive",
        "1": "T naive",
        "2": "T cells",
        "3": "T cells",
        "4": "T CD8+",
        "5": "no determined",
        "6": "naive CD20+ B",
        "7": "CD14+ Mono",
        "8": "CD16+ Mono",
        "9": "plasmacells",
        "10": "no determined",
        "11": "no detemined",
    }
)
# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(adata, groupby="leiden_res_1.40", method="wilcoxon")
sc.pl.rank_genes_groups(
    adata,
    groupby="leiden_res_1.40",
    standard_scale="var",
    n_genes=5,
    use_raw=True,
    show=False,
    save=".png"
)
sc.pl.umap(adata, color="cell_type_lvl1", show=False, save="_cell_type_lvl1.png")

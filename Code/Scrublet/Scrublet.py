#to install: pip install scrublet in terminal

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import gzip
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42



#read the matix file


root_dir = "/home/alexander-bontempo/Desktop/GitHub/scRNAseq/Code/nf-scRNA/data"

directories =  [d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))]
matrixes = {}
genes={}
for d in directories:
    dir_path = os.path.join(root_dir, d)
    matrix_path = os.path.join(dir_path, 'matrix.mtx.gz')
    genes_path = os.path.join(dir_path, 'features.tsv.gz')
    with gzip.open(matrix_path, 'rt') as tar:
        matrixes[d] = scipy.io.mmread(tar).T.tocsc()
    with gzip.open(genes_path, 'rt') as tar1:
         genes_temp = [line.strip().split('\t')[1] for line in tar1]
         genes[d]= np.array(genes_temp)
   
for i in range(len(directories)):
    print('DireSample: {}'.format(directories[i]))
    print('Number of cells: {}'.format(matrixes[directories[i]].shape[0]))
    print('Number of genes: {}'.format(matrixes[directories[i]].shape[1]))
    print('Gene names: {}'.format(genes[directories[i]][:10]))  # Print first 10 gene names

scrub = scr.Scrublet(matrixes[directories[0]], expected_doublet_rate=0.06)


print(root_dir +directories[0]+ 'genes.tsv')


doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

scrub.call_doublets(threshold=0.25)                                                          
scrub.plot_histogram()


print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

# # Uncomment to run tSNE - slow
# print('Running tSNE...')
# scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))

# # Uncomment to run force layout - slow
# print('Running ForceAtlas2...')
# scrub.set_embedding('FA', scr.get_force_layout(scrub.manifold_obs_, n_neighbors=5. n_iter=1000))
    
print('Done.')

scrub.plot_embedding('UMAP', order_points=True);

# scrub.plot_embedding('tSNE', order_points=True);
# scrub.plot_embedding('FA', order_points=True);

#the predicted doublets are :
scrub.predicted_doublets_

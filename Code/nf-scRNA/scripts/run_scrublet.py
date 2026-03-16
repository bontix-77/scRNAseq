import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import gzip
import sys
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'DejaVu Sans' # if in windows 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42
# Define your seed for random operations
RANDOM_SEED = 1997


#read the argument pointing to the data folder
root_dir= sys.argv[1]
print("#################################################root_dir: "+root_dir)
# root_dir = "/home/alexander-bontempo/Desktop/HIV GSM/GSM1"

directories =  [root_dir+p for p in sys.argv[2] ]
matrixes = {}
genes={}

print(sys.argv[2]+" :argv[2]################################")
directories= str(sys.argv[2]).strip().split()

for d in directories:
    print (f"###########################################d : {d}")
    dir_path = os.path.join(root_dir, d)
    print (f"###########################################d : {dir_path}")
    matrix_path = os.path.join(dir_path,"matrix.mtx.gz")
    print(matrix_path)
    genes_path = os.path.join(dir_path, 'features.tsv.gz')
    print(genes_path)
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



def scrublet_ (matrix_):
    scrub = scr.Scrublet(matrix_, expected_doublet_rate=0.06)
  
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
    scrub.plot_embedding('UMAP', order_points=True)

# scrub.plot_embedding('tSNE', order_points=True);
# scrub.plot_embedding('FA', order_points=True);
#the predicted doublets are in scrub.predicted_doublets_:
    return( np.where(scrub.predicted_doublets_==True)[0])

dublets_vector={}
for m in range(len(matrixes)):
    matrix=matrixes[directories[m]]
    matrix=scrublet_(matrix)
    #adding 1 because the idexes will be used in R where the index start from 1
    dublets_vector[directories[m]]=(matrix+1).tolist()

dublets_vector=pd.Series(dublets_vector)
dublets_vector=dublets_vector.to_frame().reset_index()
dublets_vector.to_csv(('dublets_vectors.csv'), index=False)





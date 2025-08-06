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


root_dir = "C:/Users/Owner/Documents/scRNAseq/Data/GSM"

directories =  [d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))]
matrices = {}
genes={}
for d in directories:
    dir_path = os.path.join(root_dir, d)
    matrix_path = os.path.join(dir_path, 'matrix.mtx.gz')
    genes_path = os.path.join(dir_path, 'features.tsv.gz')
    with gzip.open(matrix_path, 'rt') as tar:
        matrices[d] = scipy.io.mmread(tar).T.tocsc()
    with gzip.open(genes_path, 'rt') as tar1:
         genes_temp = [line.strip().split('\t')[1] for line in tar1]
         genes[d]= np.array(genes_temp)
   
for i in range(len(directories)):
    print('DireSample: {}'.format(directories[i]))
    print('Number of cells: {}'.format(matrices[directories[i]].shape[0]))
    print('Number of genes: {}'.format(matrices[directories[i]].shape[1]))
    print('Gene names: {}'.format(genes[directories[i]][:10]))  # Print first 10 gene names

scrub = scr.Scrublet(matrices[directories[1]], expected_doublet_rate=0.06)


print(root_dir +directories[0]+ 'genes.tsv')

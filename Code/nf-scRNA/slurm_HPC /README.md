The sbatch file permits to execute the seurat nextflow on HPC using singularity.

Please pull the docker image in the nf-scRNA folder on the HPC using :

``` bash
singularity pull docker://bontix77/sc_rna:v1.1
```

ATTENTION: to run the HPC version of the nextflow please use the `nextflow.config` and the `seurat_celldex.R` versions contained in this folder. Include to the `r-scripts` folder the file `celldex.rds`.

run the slum file by:

``` bash
sbatch seurat_slurm.sh
```

## Expected folder structure



├── nf_scRNA
    ├── r-scripts/ 
    │     ├── read_h5.R 
    │     ├── cell_filter.R 
    │     ├── SCTransform.R 
    │     ├── RDS_to_h5ad.R 
    │     ├── seurat_celldex.R     # Make sure to use the version found in slurm_HPC/
    │     ├── seurat_markers.R 
    │     ├── seurat_PCA_UMAP.R 
    │     ├── seurat_pseudobulk.R
    │     └── celldex.rds          # This file is found in slurm_HPC/
    │     
    │
    ├── data/
    │     ├── sample1/   #Change `sample1` for the name of your sample. Use the same names in nextflow.config
    │     │    ├── matrix.mtx.gz
    │     │    ├── features.tsv.gz
    │     │    └── barcodes.tsv.gz
    │     ├── sapmple2/
    │     │    ├── matrix.mtx.gz
    │     │    ├── features.tsv.gz
    │     │    └── barcodes.tsv.gz│     
    │     .
    │     . 
    │
    ├── main_HPC.nf
    ├── nextflow.config              # Make sure to use the version found in slurm_HPC/
    └── seurat_slurm.sh
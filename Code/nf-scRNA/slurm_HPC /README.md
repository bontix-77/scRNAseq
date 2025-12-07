The sbatch file permits to execute the seurat nextflow on HPC using singularity.

Please pull the docker image in the nf-scRNA folder on the HPC using :

``` bash
singularity pull docker://bontix77/sc_rna:v1.1
```
For runnung the scanpy branch you also need to compile the python_modules docker file into a docker image and pull it to the HPC as singularity image. All the files are found in the Docker folder.


ATTENTION: to run the HPC version of the nextflow please use the `nextflow.config` and the `seurat_celldex.R` versions contained in this folder. Include to the `r-scripts` folder the file `celldex.rds`.

modules required in the HPC: slurm, nextflow and singularity.

run the slum file by:

``` bash
sbatch seurat_slurm.sh
```

## Expected folder structure
(All the R script for HPC are same of those found in nf-scRNA/r-script/ except seurat_celldex.R)

```         
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
    ├── sc_rna_v1.1.sif   # obtained by : singularity pull docker://bontix77/sc_rna:v1.1
    ├── pytnon_modules.sif  # obtained by compiling the docker file into a docker image then
    │                         and then by: singularity pull docker://<your_dockerhub_username>/python_modules:latest
    │                         Docker files available in the `docker_files/pyton_moldules` folder
    │
    ├── main_HPC.nf
    ├── nextflow.config            # Make sure to use the version found in slurm_HPC/
    └── seurat_slurm.sh
```

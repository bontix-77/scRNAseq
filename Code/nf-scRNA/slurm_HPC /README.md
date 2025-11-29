The sbatch file permits to execute the seurat nextflow on HPC using singularity.

Please pull the docker image in the nf-scRNA folder on the HPC using :

``` bash
singularity pull docker://bontix77/sc_rna:v1.1
```

ATTENTION: to run the HPC version of the nextflow please use the `nextflow.config` and the `seurat_celldex.R` versions contained in this folder. Include to the `r-scripts` folder the file `celldex.rds`.

run the slum file by:
```bash 
sbatch seurat_slurm.sh
```

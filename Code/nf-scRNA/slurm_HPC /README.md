The sbatch file permits to execute the seurat nextflow on HPC using singularity.

Please pull the docker image in the nf-scRNA folder on the HPC using :

``` bash
singularity pull docker://bontix77/sc_rna:v1.1
```

ATTENTION: At the moment the seurat_cellType.R is not available on HPC since celldex is higly dependant to internet access, and compute nodes usually don't have access to the web. For practicity the script will be adapted using Azimuth soon....STAY TUNED.
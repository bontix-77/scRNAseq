## docker image sc-rna:v1.1


include the image v1.0 adding the Seurat fix necesary to run properly SingleR and solve othe warnings:
devtools::install_github("satijalab/Seurat", ref ="fix/v.5.3.1")

It include the packages celldex and SingleR
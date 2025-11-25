# Seurat v.5.3.1 Base Docker Image

This image provides a complete environment for single cell RNA analysis in R using Seurat v.5.3.1. It includes all system dependencies and the most common R packages used for processing, annotation and visualization of scRNAseq data. The goal is to offer a ready to use container that works out of the box on any system, including local machines, HPC clusters and workflow managers such as Nextflow.

Base image

This container is built on:

rocker/r-ver:4.5.1

The image contains a full R installation together with the required system libraries for building and running Seurat and its extensions.

Included system libraries

The container installs a wide set of libraries needed by Seurat, Azimuth, BPCells, presto, Signac, hdf5r and other tools. These include:

SSL, XML, CURL, LZMA, BZ2 and general compression libraries

PNG, JPEG, TIFF, font and graphic libraries

Cairo, fontconfig, freetype, harfbuzz and related components

Mathematical and mapping libraries such as GSL, PROJ, GDAL, GEOS and UDUNITS

HDF5 and hdf5 tools

FFTW for numerical transforms

Build tools such as pkg-config, cmake, patch, gfortran and git2

This ensures that all R packages compile correctly without manual intervention.

Installed R packages

The container installs a large collection of R packages that cover:

Core components

Matrix

SeuratObject

Seurat v.5.3.1

BPCells

presto

glmGamPoi

hdf5r

Bioinformatics

Signac

Rsamtools

Rhtslib

BSgenome.Hsapiens.UCSC.hg38

Statistical and enrichment tools

multtest

mutoss

metap

qqconf

DirichletMultinomial

Transcription factor and motif analysis

TFBSTools

Spatial and geometry support

sf

Satija Lab ecosystem

seurat-data

seurat-wrappers

azimuth

GitHub installations

Installed using the remotes package where needed.

The CRAN mirror is configured to cloud.r-project.org for reproducibility.

Purpose

This image provides a stable and reproducible environment for:

scRNAseq preprocessing

clustering and dimensionality reduction

annotation with Azimuth

multimodal analysis with Signac

large matrix handling with BPCells

enrichment and statistical workflows

next generation container based pipelines

It is ideal for integration into Nextflow pipelines that require consistency and portability.

Usage

Pull the image directly:

``` bash
docker pull bontix77/sc_rna:v1.0
```

Run it interactively:

``` bash
docker run -it bontix77/sc_rna:v1.0 R
```

Use it in Nextflow:

``` groovy
process {
    container = 'bontix77/sc_rna:v1.0'
}
```

Repository

Maintainer

Alexander Bontempo Email: bontix77\@yahoo.it

GitHub: bontix77

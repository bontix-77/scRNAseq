#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// def mode = params.analysis.tostring().toLowerCase()
// assert  ['seurat', 'scanpy'].contains(mode) : "Invalid analysis mode: ${mode}. Choose 'seurat' or 'scanpy'" 
workflow{
folder_ch= channel.fromPath("${params.inputDir}/*", type:"dir")
def mode = params.analysis.toString().toLowerCase()
assert  ['seurat', 'scanpy'].contains(mode) : "Invalid analysis mode: ${mode}. Choose 'seurat' or 'scanpy'" 

//rscript= file("${params.scriptFile}/test.R")


RUN_PRINT_R(folder_ch.collect())

RUN_readh5(folder_ch.collect(), params.inputDir)
RUN_cell_filter( RUN_readh5.out.merged_rds ) 
RUN_SCTransform( RUN_cell_filter.out.filtered_rds ) 

if (mode == 'seurat'){
    RUN_Seurat_PCA( RUN_SCTransform.out.SCTransformed_rds )
    // RUN_Seurat_UMAP( RUN_Seurat_PCA.out.S_PCA_rds )
} else {
    RUN_seuratDisk( RUN_SCTransform.out.SCTransformed_rds )
}
}
process RUN_PRINT_R{
    publishDir params.resultDir, mode: 'copy', overwrite: true
input:
    path folders

    

    script:
    """

echo "Processing file: ${folders}" > print_log.txt
    """

}
// Define process to run R script
process RUN_readh5 {
    publishDir "${params.resultDir}/merged_h5", mode: 'copy', overwrite: true
    input:
    path folders
    val  base_dir
    
    
    output:
    path "HIV_merged.rds", emit: merged_rds

    script:
    """
  
    Rscript ${params.scriptFile}/read_h5.R "${folders}" ${base_dir}
    """
    
 }


process RUN_cell_filter {
    publishDir "${params.resultDir}/cell_filter", mode: 'copy', overwrite: true

    input:
    path rds_file 

    output:
    path "HIV_filtered.rds", emit: filtered_rds
    path "*.png"
    script:
    """
    Rscript ${params.scriptFile}/cell_filter.R "${rds_file}"
    """
} 


process RUN_SCTransform {
    publishDir "${params.resultDir}/SCTrasnform", mode: 'copy', overwrite: true

    input:
    path SCTransform

    output:
    path "HIV_SCTransform.rds", emit: SCTransformed_rds
    script:
    """
    Rscript ${params.scriptFile}/SCTransform.R "${SCTransform}"
    """
} 

process RUN_seuratDisk {
    publishDir "${params.resultDir}/seuratDisk", mode: 'copy', overwrite: true

    input:
    path SCTransform

    output:
    path "*.txt", emit: gene_names
    path "data.h5Seurat"
    path "scanpy.h5ad", emit: h5Seurat_file
    script:
    """
    Rscript ${params.scriptFile}/RDS_to_h5ad.R "${SCTransform}"
    """
}
#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Create a channel with all files in ./data/\

workflow{
folder_ch= channel.fromPath("${params.inputDir}/*", type:"dir")

//rscript= file("${params.scriptFile}/test.R")


RUN_PRINT_R(folder_ch.collect())

RUN_readh5(folder_ch.collect(), params.inputDir)
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
    publishDir params.resultDir, mode: 'copy', overwrite: true
    input:
    path folders
    val  base_dir
    
    
    output:
    path "HIV_merged.rds"

    script:
    """
  
    Rscript ${params.scriptFile}/read_h5.R "${folders}" ${base_dir}
    """
    
 }

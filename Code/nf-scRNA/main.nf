#!/usr/bin/env nextflow
//note: redundant comments and explanation have been added since the porpose of this is to be instructional.

// ───────────────────────────────────────────────────────────────────────────────
// Shebang: lets this file be executed directly if it has execute permissions.
// Nextflow will interpret everything below according to its DSL (domain-specific language).
// ───────────────────────────────────────────────────────────────────────────────

nextflow.enable.dsl = 2
// Enabling DSL2: activates modular Nextflow syntax with processes and workflows,
// allowing input/output channels to be composed more cleanly.

// def mode = params.analysis.tostring().toLowerCase()
// assert  ['seurat', 'scanpy'].contains(mode) : "Invalid analysis mode: ${mode}. Choose 'seurat' or 'scanpy'" 
// NOTE: The two lines above are intentionally commented out. They show an earlier
// validation approach that only allowed 'seurat' or 'scanpy'. In the current script,
// a third option 'parallel' is also supported and validated inside `workflow { }`.

// ───────────────────────────────────────────────────────────────────────────────
// The main workflow block orchestrates channel creation, parameter checking, and
// the order in which processes are executed. Think of it as the pipeline's control flow.
// ───────────────────────────────────────────────────────────────────────────────
workflow {
    // Create a channel that emits directory paths from the input directory.
    // The glob `*` selects all immediate entries under params.inputDir.
    // `type: "dir"` ensures only directories are matched (not files).
    folder_ch = channel.fromPath("${params.inputDir}/*", type: "dir")
    
    // Read an analysis mode from params, normalize to lowercase for robust matching.
    def mode = params.analysis.toString().toLowerCase()
     
    // Validate the mode against allowed values. If invalid, abort with a helpful message.
    assert ['seurat', 'scanpy', 'parallel'].contains(mode) : "Invalid analysis mode: ${mode}. Choose 'seurat', 'scanpy' or 'parallel params: --analysis"
    

    def samples =params.samples
    // Run the HDF5 reading/merging step, passing both the collected list of folders
    // and the base input directory (params.inputDir) as an additional value.
    RUN_scrublet(folder_ch.collect(),params.inputDir)
    def dublets=RUN_scrublet.out.dublets
    RUN_readh5(folder_ch.collect(), params.inputDir,samples,dublets)
    
    // Filter cells based on quality metrics or other criteria. This consumes the
    // merged RDS produced by the prior step.
    RUN_cell_filter(RUN_readh5.out.merged_rds)

    // Apply Seurat's SCTransform normalization. Consumes the filtered RDS and produces
    // a transformed object for downstream analysis or conversion.
    RUN_SCTransform(RUN_cell_filter.out.filtered_rds)

    // Branch the workflow based on the selected mode. If 'seurat', continue
    // with PCA in Seurat. Otherwise (for 'scanpy' or 'parallel'), convert the
    // Seurat object to h5Seurat/h5ad using SeuratDisk and proceed toward Scanpy.
    def pval = params.pseudoBulk_p
    if (mode == 'seurat') {
        RUN_Seurat_PCA_UMAP(RUN_SCTransform.out.SCTransformed_rds)
        RUN_Seurat_markers(RUN_Seurat_PCA_UMAP.out.PCA_UMAP_seurat)
        RUN_Seurat_pseudoBulk(RUN_Seurat_markers.out.markers_cluster,pval)
        // Read from params cellType_manual to include manual cell type annotation (stillin development)
        def man_annotation = params.cellType_manual.toString().toLowerCase()
        if (man_annotation == "yes"){
              RUN_Seurat_cellTypes_manual(RUN_Seurat_markers.out.markers_cluster)
        }
        RUN_Seurat_cellType_automatic(RUN_Seurat_markers.out.markers_cluster)
    }
    else if (mode=='scanpy'){
        RUN_seuratDisk(RUN_SCTransform.out.SCTransformed_rds)
        RUN_scanpy(RUN_seuratDisk.out.h5Seurat_file,RUN_seuratDisk.out.gene_names)
    }
    else {
        RUN_Seurat_PCA_UMAP(RUN_SCTransform.out.SCTransformed_rds)
        RUN_Seurat_markers(RUN_Seurat_PCA_UMAP.out.PCA_UMAP_seurat)
        RUN_Seurat_pseudoBulk(RUN_Seurat_markers.out.markers_cluster,pval)
        // Read from params cellType_manual to include manual cell type annotation (stillin development)
        def man_annotation = params.cellType_manual.toString().toLowerCase()
        if (man_annotation == "yes"){
              RUN_Seurat_cellTypes_manual(RUN_Seurat_markers.out.markers_cluster)
        }
        RUN_Seurat_cellType_automatic(RUN_Seurat_markers.out.markers_cluster)
        RUN_seuratDisk(RUN_SCTransform.out.SCTransformed_rds)
        RUN_scanpy(RUN_seuratDisk.out.h5Seurat_file,RUN_seuratDisk.out.gene_names)
    }
}
// ───────────────────────────────────────────────────────────────────────────────
// PROCESS: RUN_readh5
// Purpose: Run an R script that reads raw 10x/other HDF5s contained in the given
// directories, merges them, and outputs a single Seurat object (RDS).
// Output artifact: HIV_merged.rds
// Notes:
// - `folders` contains the list of directory paths to scan.
// - `base_dir` passes the original input directory string for the R script's internal logic.
// ───────────────────────────────────────────────────────────────────────────────
process RUN_scrublet {
    publishDir "${params.resultDir}/scrublet", mode: 'copy', overwrite: true
    container 'python:v1.0'

    input:

    
    path folders
    val base_dir

    output:
    path "dublets_vectors.csv", emit: dublets

    script:
    """
    python "${params.scriptFile}/run_scrublet.py" "${base_dir}" "${folders}"
    """
}

process RUN_readh5 {
    // Publish all outputs to results/merged_h5 for consistent artifact organization.
    publishDir "${params.resultDir}/read_h5", mode: 'copy', overwrite: true
    // container 'bontix77/sc_rna:v1.0'

    input:
    path folders
    val base_dir
    val samples
    val dublets

    output:
    path "HIV_merged.rds", emit: merged_rds

    script:
    
    """
    Rscript "${params.scriptFile}/read_h5.R" "${folders}" ${base_dir} "${params.colapse_T_eceptors}" "${samples}" "${dublets}"
    """
}


// ───────────────────────────────────────────────────────────────────────────────
// PROCESS: RUN_cell_filter
// Purpose: Apply cell-level QC filters (e.g., min genes, max mitochondrial %, etc.).
// Input artifact: HIV_merged.rds
// Output artifacts: HIV_filtered.rds plus any PNG plots (QC visualizations).
// ───────────────────────────────────────────────────────────────────────────────
process RUN_cell_filter {
    // Organize outputs into a dedicated subdirectory for clarity.
    publishDir "${params.resultDir}/cell_filter", mode: 'copy', overwrite: true
    // container 'bontix77/sc_rna:v1.0'

    input:

    path rds_file

    output:
    path "HIV_filtered.rds", emit: filtered_rds
    path "*.png"

    script:
    """
    Rscript "${params.scriptFile}/cell_filter.R" "${rds_file}"
    """
}


// ───────────────────────────────────────────────────────────────────────────────
// PROCESS: RUN_SCTransform
// Purpose: Normalize and variance-stabilize the filtered object using Seurat's
// SCTransform workflow. Produces a new RDS for downstream steps.
// Note: The publishDir path intentionally says "SCTrasnform" (typo mirrored as-is).
// ───────────────────────────────────────────────────────────────────────────────
process RUN_SCTransform {
    // Keep outputs grouped under SCTrasnform for traceability with prior runs.
    publishDir "${params.resultDir}/SCTrasnform", mode: 'copy', overwrite: true
    // container 'bontix77/sc_rna:v1.0'

    input:
    path SCTransform

    output:
    path "HIV_SCTransform.rds", emit: SCTransformed_rds

    script:
    """
    Rscript "${params.scriptFile}/SCTransform.R" "${SCTransform}"
    """
}

// ───────────────────────────────────────────────────────────────────────────────
// PROCESS: RUN_seuratDisk
// Purpose: Convert a Seurat RDS to h5Seurat and then to AnnData (.h5ad) using
// SeuratDisk (and usually SeuratDisk::Convert). This enables downstream analysis
// in Python/Scanpy while preserving metadata and layers as much as possible.
// Outputs:
// - One or more "*.txt" (e.g., gene lists or logs) [emitted as `gene_names`]
// - `data.h5Seurat`: intermediate SeuratDisk format
// - `scanpy.h5ad`: final AnnData for Scanpy [emitted as `h5Seurat_file`]
// ───────────────────────────────────────────────────────────────────────────────
process RUN_seuratDisk {
    // Publish all conversion outputs into a dedicated seuratDisk subfolder.
    publishDir "${params.resultDir}/seuratDisk", mode: 'copy', overwrite: true
    // container 'bontix77/sc_rna:v1.0'

    input:
    path SCTransform

    output:
    path "*.txt", emit: gene_names
    path "data.h5Seurat"
    path "scanpy.h5ad", emit: h5Seurat_file

    script:
    """
    Rscript "${params.scriptFile}/RDS_to_h5ad.R" "${SCTransform}"
    """
}
process RUN_scanpy {
    publishDir "${params.resultDir}/scanpy_analysis", mode: 'copy', overwrite: true
     container 'python:v1.0'
    input:
    path h5Seurat_file
    path genes_names

    output:
    path "figures/umap_leiden_res_1.4.png"
    path "figures/rank_genes_groups_leiden_res_1.40.png"
    // path "*.csv"
    // path "*.h5ad"

    script:
    """
    python "${params.scriptFile}/run_scanpy.py" "${h5Seurat_file}" "${genes_names}"
    """
}
process RUN_Seurat_PCA_UMAP {
    publishDir "${params.resultDir}/seurat_PCA_UMAP", mode: 'copy', overwrite: true
    // container 'bontix77/sc_rna:v1.0'

    input:
    path SCTransform

    output:
    path "*.rds", emit: PCA_UMAP_seurat
    path "elbowplot.png"
    path "DimHeatmap.png"
    path "*.png"
    path "UMAP_DimPlot.png"

    script:
    """
    Rscript "${params.scriptFile}/seurat_PCA_UMAP.R" "${SCTransform}" "${params.harmony}" "${params.resolution}"
    """
}
process RUN_Seurat_markers {

    publishDir "${params.resultDir}/seurat_Markers", mode: 'copy', overwrite: true
    // container 'bontix77/sc_rna:v1.0'

    input:
    path PCA_UMAP

    output:
    path "Top20_markers.csv"
    path "Clusters.png"
    path "*.rds", emit: markers_cluster

    script:
    """
    Rscript "${params.scriptFile}/seurat_markers.R" "${PCA_UMAP}"
    """
}
process RUN_Seurat_pseudoBulk {
    publishDir "${params.resultDir}/Pseudo_bulk", mode: 'copy', overwrite: true

    input:
    path Markers
    val pval

    output:

    path "differentialExpression_singlecell.csv"
    path "differentialExpression_bulk.csv"
    script:
    """
    Rscript "${params.scriptFile}/seurat_pseudoBulk.R" "${Markers}" ${pval}
    """


}
process RUN_Seurat_cellTypes_manual {

    publishDir "${params.resultDir}/seurat_cellType", mode: 'copy', overwrite: true
    // container 'bontix77/sc_rna:v1.0'

    input:
    path CellType

    output:
    path "CellType.png"
    path "*.rds"

    script:

    """
    Rscript "${params.scriptFile}/seurat_cellTypes.R" "${CellType}"
    """
}
process RUN_Seurat_cellType_automatic {

    publishDir "${params.resultDir}/seurat_cellType", mode: 'copy', overwrite: true
    // container 'sc_rna:v1.1'

    input:
    path Markers

    output:
    path "CellType_celldex.png"
    path "CellType_final.png"
    path "*.rds"

    script:

    """
    Rscript "${params.scriptFile}/seurat_celldex.R" "${Markers}"
    """
}

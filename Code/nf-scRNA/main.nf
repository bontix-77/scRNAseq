//note: redundant explanation and comments have been added as the porposue is to make it clear for new users.



#!/usr/bin/env nextflow
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
    assert ['seurat', 'scanpy','parallel'].contains(mode) : "Invalid analysis mode: ${mode}. Choose 'seurat'. 'scanpy' or 'parallel params: --analysis"

    //rscript= file("${params.scriptFile}/test.R")
    // (Commented out) Example of how you might have bound a script path to a file handle.

    // Kick off a quick diagnostic print to verify what folders were captured.
    // `.collect()` turns the multi-value channel into a single list value for the process.
    //in the final version this step is not necessary, so it is commented out.
    //RUN_PRINT_R(folder_ch.collect())

    // Run the HDF5 reading/merging step, passing both the collected list of folders
    // and the base input directory (params.inputDir) as an additional value.
    RUN_readh5(folder_ch.collect(), params.inputDir)

    // Filter cells based on quality metrics or other criteria. This consumes the
    // merged RDS produced by the prior step.
    RUN_cell_filter(RUN_readh5.out.merged_rds)

    // Apply Seurat's SCTransform normalization. Consumes the filtered RDS and produces
    // a transformed object for downstream analysis or conversion.
    RUN_SCTransform(RUN_cell_filter.out.filtered_rds)

    // Branch the workflow based on the selected mode. If 'seurat', continue
    // with PCA in Seurat. Otherwise (for 'scanpy' or 'parallel'), convert the
    // Seurat object to h5Seurat/h5ad using SeuratDisk and proceed toward Scanpy.
    if (mode == 'seurat') {
        RUN_Seurat_PCA(RUN_SCTransform.out.SCTransformed_rds)
    }
    else {
        RUN_seuratDisk(RUN_SCTransform.out.SCTransformed_rds)
    }
}

// ───────────────────────────────────────────────────────────────────────────────
// PROCESS: RUN_PRINT_R
// Purpose: Simple debugging/visibility step to print which folders are being processed.
// Publishes a small log file to the results directory for traceability.
// ───────────────────────────────────────────────────────────────────────────────
process RUN_PRINT_R {
    // publishDir: copies outputs to a user-defined results folder (idempotent with overwrite).
    publishDir params.resultDir, mode: 'copy', overwrite: true

    // INPUTS:
    // - `path folders`: a value representing a path or list of paths. Here it's the
    //   `collect()`'ed list of directories emitted by folder_ch.
    input:
    path folders

    // SCRIPT:
    // Writes a line with the resolved input to a text file for later inspection.
    script:
    """

echo "Processing file: ${folders}" > print_log.txt
    """
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
process RUN_readh5 {
    // Publish all outputs to results/merged_h5 for consistent artifact organization.
    publishDir "${params.resultDir}/merged_h5", mode: 'copy', overwrite: true

    // INPUTS:
    // - `path folders`: list of directories to process
    // - `val base_dir`: plain value (string) with the input directory root
    input:
    path folders
    val base_dir

    // OUTPUTS:
    // - A single RDS file containing the merged Seurat object.
    output:
    path "HIV_merged.rds", emit: merged_rds

    // SCRIPT:
    // Invoke the R script with positional arguments:
    //   1) the list (or Nextflow-resolved) of folders
    //   2) the base directory
    // The script must create HIV_merged.rds in the working directory for Nextflow to collect it.
    script:
    """
  
    Rscript ${params.scriptFile}/read_h5.R "${folders}" ${base_dir}
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

    // INPUTS:
    // - `path rds_file`: the merged Seurat object from RUN_readh5
    input:
    path rds_file

    // OUTPUTS:
    // - `HIV_filtered.rds`: the filtered Seurat object
    // - `*.png`: any QC plots the R script generates (e.g., violin plots, scatter plots)
    output:
    path "HIV_filtered.rds", emit: filtered_rds
    path "*.png"

    // SCRIPT:

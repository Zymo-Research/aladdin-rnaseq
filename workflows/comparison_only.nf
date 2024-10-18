#!/usr/bin/env nextflow
/*
Workflow of only sample/group comparison
*/

params.summary = [:]

// Load functions
include { setup_channel } from ('../libs/setup_channel')

/*
 * SET UP CONFIGURATION VARIABLES
 */
// Genome options
if (params.genome) {
    if (!params.genomes.containsKey(params.genome)) {
        exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
    }
} else {
    exit 1, "--genome is a required input!"
}
// Validate read quantification method
if (!(params.read_quant_method in ['STAR_featureCounts', 'STAR_Salmon'])) {
    exit 1, "Read quantification method only accepts STAR_featureCounts or STAR_Salmon, for now!"
}
if (params.dtu_analysis && params.read_quant_method!='STAR_Salmon') {
    exit 1, "When DTU analysis is requested, read quantification method must be STAR_Salmon!"
}
// DTU cannot be done for bacteria
if (params.dtu_analysis && (params.genomes[ params.genome ].bacteria ?: false)) {
    exit 1, "DTU analysis cannot be carried out for a bacteria genome (only one transcript per gene)!"
}

/*
 * SET & VALIDATE INPUT CHANNELS
 */
if (params.read_quant_method=='STAR_featureCounts') {
    merged_counts = setup_channel(params.merged_counts, "When read_quant_method==STAR_featureCounts, merged read counts", true, "")
    salmon_results = Channel.of([])
} else {
    if (params.salmon_results) {
        salmon_results = Channel.fromPath(params.salmon_results, type:'dir').collect()
    } else {
        exit 1, "When read_quant_method==STAR_Salmon, --salmon_results is a required input!"
    }
    merged_counts = Channel.of([])
}
design = setup_channel(params.design, "design CSV file", true, "")
comparisons = setup_channel(params.comparisons, "comparison CSV file", false, "all pairwise comparisons will be carried out.")
multiqc_config = setup_channel(params.multiqc_config, "MultiQC config", true, "")
summary_header = Channel.fromPath("$baseDir/assets/workflow_summary_header.txt")

/*
 * COLLECT SUMMARY & LOG
 */
// Create a channel for settings to show in MultiQC report
info_to_report = ['Genome', 'Read Quant. Method', 'DESeq2 FDR cutoff', 'DESeq2 Log2FC cutoff', 'gProfiler FDR cutoff', 'DTU analysis FDR cutoff']
summary_to_report = params.summary.findAll { it.key in info_to_report }
Channel.from( summary_to_report.collect{ [it.key, it.value] } )
       .map { k,v -> "            <dt>$k</dt><dd><samp>$v</samp></dd>" }
       .collectFile(name: 'summary.txt', newLine: true, sort: 'index')
       .set { workflow_summary_to_report }
// Save workflow summary plain text
Channel.from( params.summary.collect{ [it.key, it.value] } )
       .map { k,v -> "${k.padRight(21)} : $v" }
       .collectFile(name: "${params.outdir}/pipeline_info/workflow_summary.txt", newLine: true, sort: 'index')

// IMPORT SUBWORKFLOWS
include { COMPARISON } from '../subworkflows/comparison'

/*
 * PROCESS DEFINITION
 */
include { check_design } from '../processes/check_design'
include { software_versions } from "../processes/software_versions" addParams(
    publish_dir: "${params.outdir}/pipeline_info",
    pipeline_version: workflow.manifest.version,
    nextflow_version: workflow.nextflow.version
)
include { multiqc_comparison_only } from "../processes/multiqc" addParams(
    publish_dir: "${params.outdir}/MultiQC",
    skip_multiqc: params.skip_multiqc,
    run_name: params.summary["Run Name"],
    ensembl_web: params.genomes[ params.genome ].ensembl_web ?: false,
    deseq2_fdr: params.deseq2_fdr,
    gprofiler_fdr: params.gprofiler_fdr,
    dexseq_fdr: params.dexseq_fdr
)
include { summarize_downloads } from "../processes/summarize_downloads" addParams(
    publish_dir: "${params.outdir}/download_data"
)

/*
 * WORKFLOW DEFINITION
 */
workflow COMPARISON_ONLY {
    // Check design file and set up channels
    check_design(design, comparisons.ifEmpty([]))

    // Comparisons
    COMPARISON(
        merged_counts, 
        salmon_results,
        check_design.out.deseq2_design,
        comparisons.ifEmpty([])
    )
    
    // Software versions
    software_versions(COMPARISON.out.versions.flatten().collect())

    // Report
    multiqc_comparison_only(multiqc_config, \
                            COMPARISON.out.report.collect().ifEmpty([]), \
                            software_versions.out.report.collect(), \
                            summary_header, \
                            workflow_summary_to_report)
    report_locations = multiqc.out.report.map { "${params.outdir}/MultiQC/" + it.getName() }

    // Collect file locations
    locations = COMPARISON.out.file_locations
                    .mix(report_locations)
                    .collectFile(name: "${params.outdir}/download_data/file_locations.txt", newLine: true )

    // Parse the list of files for downloading into a JSON file
    summarize_downloads( locations, [], check_design.out.checked_design )
}

/*
 * LOG ON COMPLETION
 */
workflow.onComplete {
    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.warn "Warning, pipeline completed, but with errored process(es)"
      log.info "Number of ignored errored process(es) : ${workflow.stats.ignoredCount}"
      log.info "Number of successfully ran process(es) : ${workflow.stats.succeedCount}"
    }
    if (workflow.success) {
        log.info "[rnaseq] Pipeline completed successfully"
    } else {
        log.warn "[rnaseq] Pipeline completed with errors"
    }
}

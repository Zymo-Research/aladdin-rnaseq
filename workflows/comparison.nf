#!/usr/bin/env nextflow
/*
Workflow of only sample/group comparison
*/

params.summary = [:]

// Load functions
include { setup_channel } from ('../libs/setup_channel')
include { parse_design } from ('../libs/parse_design')

/*
 * SET UP CONFIGURATION VARIABLES
 */
// Genome options
// Check if genome exists in the config file
if (params.genome) {
    if (!params.genomes.containsKey(params.genome)) {
        exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
    }
} else {
    log.warn "Genome not provided, gProfiler analysis will be skipped."
}
// Reference index path configuration
params.gprofiler_organism = params.genome ? params.genomes[ params.genome ].gprofiler ?: false : false
params.ensembl_web = params.genome ? params.genomes[ params.genome ].ensembl_web ?: false : false

/*
 * SET & VALIDATE INPUT CHANNELS
 */
merged_counts = setup_channel(params.merged_counts, "merged read counts", true, "")
design = setup_channel(params.design, "design CSV file", true, "")
comparisons = setup_channel(params.comparisons, "comparison CSV file", false, "all pairwise comparisons will be carried out.")
multiqc_config = setup_channel(params.multiqc_config, "MultiQC config", true, "")
summary_header = Channel.fromPath("$baseDir/assets/workflow_summary_header.txt")

/*
 * COLLECT SUMMARY & LOG
 */
// Create a channel for settings to show in MultiQC report
info_to_report = ['Genome', 'DESeq2 FDR cutoff', 'DESeq2 Log2FC cutoff', 'gProfiler FDR cutoff']
summary_to_report = params.summary.findAll { it.key in info_to_report }
Channel.from( summary_to_report.collect{ [it.key, it.value] } )
       .map { k,v -> "            <dt>$k</dt><dd><samp>$v</samp></dd>" }
       .collectFile(name: 'summary.txt', newLine: true, sort: 'index')
       .set { workflow_summary_to_report }
// Save workflow summary plain text
Channel.from( params.summary.collect{ [it.key, it.value] } )
       .map { k,v -> "${k.padRight(21)} : $v" }
       .collectFile(name: "${params.outdir}/pipeline_info/workflow_summary.txt", newLine: true, sort: 'index')

/*
 * PROCESS DEFINITION
 */
include { check_design } from '../processes/check_design'
include { deseq2 } from '../processes/deseq2.nf' addParams(
    publish_dir: "${params.outdir}/DESeq2",
    skip_deseq2: params.skip_deseq2,
    deseq2_fdr: params.deseq2_fdr,
    deseq2_lfc: params.deseq2_lfc
)
include { gprofiler } from '../processes/gprofiler.nf' addParams(
    publish_dir: "${params.outdir}/gProfiler",
    gprofiler_organism: params.gprofiler_organism,
    deseq2_fdr: params.deseq2_fdr,
    gprofiler_fdr: params.gprofiler_fdr
)
include { software_versions } from "../processes/software_versions" addParams(
    publish_dir: "${params.outdir}/pipeline_info",
    pipeline_version: workflow.manifest.version,
    nextflow_version: workflow.nextflow.version
)
include { multiqc_comparison_only } from "../processes/multiqc" addParams(
    publish_dir: "${params.outdir}/MultiQC",
    skip_multiqc: params.skip_multiqc,
    run_name: params.summary["Run Name"],
    ensembl_web: params.ensembl_web,
    deseq2_fdr: params.deseq2_fdr,
    gprofiler_fdr: params.gprofiler_fdr
)

/*
 * WORKFLOW DEFINITION
 */
workflow COMPARISON {
    // Check design file and set up channels
    check_design(design, comparisons.ifEmpty([]))

    // Study level analysis
    deseq2(merged_counts, check_design.out.checked_design, comparisons.ifEmpty([]))
    gprofiler(deseq2.out.results)

    // Generate download links
    deseq2_locations = deseq2.out.download.flatten().map { "${params.outdir}/DESeq2/" + it.getName() }
    gprofiler_locations = gprofiler.out.download.flatten().map { "${params.outdir}/gProfiler/" + it.getName() }
    locations = deseq2_locations
                    .mix(gprofiler_locations)
                    .collectFile(name: "${params.outdir}/download_data/file_locations.txt", newLine: true )

    // Software versions
    versions = deseq2.out.version.mix(gprofiler.out.version)
    software_versions(versions.collect())

    // Report
    multiqc_comparison_only(multiqc_config, \
                            deseq2.out.report.mix(deseq2.out.results).collect().ifEmpty([]), \
                            gprofiler.out.report.collect().ifEmpty([]), \
                            software_versions.out.report.collect(), \
                            summary_header, \
                            workflow_summary_to_report)
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

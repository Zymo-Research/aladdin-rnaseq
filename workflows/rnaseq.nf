#!/usr/bin/env nextflow
/*
Workflow of standard RNAseq
*/

params.summary = [:]

// Load functions
include { setup_channel } from ('../libs/setup_channel')
include { parse_design } from ('../libs/parse_design')

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
// Parse protocol
if (params.protocol) {
    if (!params.protocols.containsKey(params.protocol)) {
        exit 1, "The provided protocol '${params.protocol}' is not available in the protocols file. Currently the available protocols are ${params.protocols.keySet().join(", ")}"
    } else {
        params.protocol_settings = params.protocols[params.protocol]
        params.summary['Trimming'] = params.protocol_settings['trimming_text']
        params.summary['Strandedness'] = params.protocol_settings['strandedness_text']
        params.summary['Library Prep'] = params.protocol_settings['description']
    }
} else {
    exit 1, "--protocol is a required input!"
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
design = setup_channel(params.design, "design CSV file", true, "")
comparisons = setup_channel(params.comparisons, "comparison CSV file", false, "all pairwise comparisons will be carried out.")
ercc_info = setup_channel(params.ercc_spikein ? "$baseDir/assets/ERCC92_info.tsv" : false, "ERCC info", false, "ERCC count will not be conducted.")
multiqc_config = setup_channel(params.multiqc_config, "MultiQC config", true, "")
summary_header = Channel.fromPath("$baseDir/assets/workflow_summary_header.txt", checkIfExists: true)
if (params.protocol_settings.umi) {
    extra_multiqc_config = Channel.fromPath("$baseDir/assets/multiqc_config_3mrna.yaml", checkIfExists: true)
} else {
    extra_multiqc_config = Channel.empty()
}

/*
 * COLLECT SUMMARY & LOG
 */
// Create a channel for settings to show in MultiQC report
info_to_report = ['Genome', 'ERCC spike-in', 'Library Prep', 'Strandedness', 'Trimming', 'Read Quant. Method',
                  'DESeq2 FDR cutoff', 'DESeq2 Log2FC cutoff', 'gProfiler FDR cutoff', 'DTU analysis FDR cutoff']
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
 * IMPORT SUBWORKFLOWS
 */
include { READS_QC_TRIMMING } from '../subworkflows/reads_qc_trimming'
include { ALIGN_DEDUP_QUANT } from '../subworkflows/align_dedup_quant'
include { ALIGN_DEDUP_QUANT as ALIGN_DEDUP_QUANT_ERCC } from '../subworkflows/align_dedup_quant' addParams(
    outdir: "${params.outdir}/ERCC",
    genome: 'ERCC92',
    percent_mapped_cutoff: 0,
    fc_count_type: 'exon',
    fc_group_features: 'gene_id',
    fc_extra_attributes: false
)
include { ALIGNMENT_QC } from '../subworkflows/alignment_qc'
include { COMPARISON } from '../subworkflows/comparison'

/*
 * PROCESS DEFINITION
 */
include { check_design } from '../processes/check_design'
include { cat_fastq } from '../processes/cat_fastq'
include { add_ercc_info } from '../processes/add_ercc_info' addParams(
    publish_dir: "${params.outdir}/ERCC"
)
include { software_versions } from "../processes/software_versions" addParams(
    publish_dir: "${params.outdir}/pipeline_info",
    pipeline_version: workflow.manifest.version,
    nextflow_version: workflow.nextflow.version
)
include { multiqc } from "../processes/multiqc" addParams(
    publish_dir: "${params.outdir}/MultiQC",
    skip_multiqc: params.skip_multiqc,
    run_name: params.summary["Run Name"],
    ensembl_web: params.genomes[ params.genome ].ensembl_web ?: false,
    deseq2_fdr: params.deseq2_fdr,
    gprofiler_fdr: params.gprofiler_fdr,
    dexseq_fdr: params.dexseq_fdr,
    trimming_module: params.protocol_settings.trimming_2step ? 'trimming_2step' : 'Trim_Galore'
)
include { summarize_downloads } from "../processes/summarize_downloads" addParams(
    publish_dir: "${params.outdir}/download_data"
)

/*
 * WORKFLOW DEFINITION
 */
workflow RNASEQ {

    // Check design file and set up channels
    check_design(design, comparisons.ifEmpty([]))
    check_design.out.checked_design
        .splitCsv( header: true )
        .map { parse_design(it) }
        .groupTuple()
        .map {
            meta, reads ->
                [ meta, reads.flatten() ]
        }
        .branch {
            meta, reads ->
            // Separate samples with multiple runs from those without
            multiple: ( meta.single_end && reads.size() > 1 ) || ( !meta.single_end && reads.size() > 2 )
            single: true
        }
        .set { reads }

    // Merge FASTQs for samples with mulitple runs
    cat_fastq(reads.multiple)
    ch_input_reads = cat_fastq.out.reads.mix(reads.single)

    // Read QC and trimming
    READS_QC_TRIMMING(ch_input_reads)

    // Align, dedup, quantify
    ALIGN_DEDUP_QUANT(READS_QC_TRIMMING.out.reads)
    ALIGN_DEDUP_QUANT.out.md5sum
        .collectFile( name: "${params.outdir}/download_data/md5sum.txt", sort: { it.getName() } )
        .set { md5sums }

    // Align and count ERCC reads
    if (params.ercc_spikein) {
        ALIGN_DEDUP_QUANT_ERCC(ALIGN_DEDUP_QUANT.out.unmapped)
        add_ercc_info(ALIGN_DEDUP_QUANT_ERCC.out.counts, ercc_info.collect())
        ch_ercc_counts = add_ercc_info.out.report  
    } else {
        ch_ercc_counts = Channel.empty()
    }
    
    // alignment QC
    ALIGNMENT_QC(ALIGN_DEDUP_QUANT.out.bam)

    // Comparisons
    COMPARISON(
        ALIGN_DEDUP_QUANT.out.counts.ifEmpty([]), 
        ALIGN_DEDUP_QUANT.out.salmon_results.collect().ifEmpty([]),
        check_design.out.deseq2_design,
        comparisons.ifEmpty([])
    )

    // Software versions
    versions = ALIGN_DEDUP_QUANT.out.versions
                   .mix(READS_QC_TRIMMING.out.versions, ALIGNMENT_QC.out.versions, COMPARISON.out.versions)
                   .flatten()
                   .unique{ it.getName() }
    software_versions(versions.collect())

    // Report
    multiqc(multiqc_config, \
            extra_multiqc_config.ifEmpty([]), \
            READS_QC_TRIMMING.out.report.collect().ifEmpty([]), \
            ALIGN_DEDUP_QUANT.out.report.collect().ifEmpty([]), \
            ch_ercc_counts.collect().ifEmpty([]), \
            ALIGNMENT_QC.out.report.collect().ifEmpty([]), \
            COMPARISON.out.report.collect().ifEmpty([]), \
            software_versions.out.report.collect(), \
            summary_header, \
            workflow_summary_to_report, \
            ALIGN_DEDUP_QUANT.out.warning.ifEmpty('')
    )
    report_locations = multiqc.out.report.map { "${params.outdir}/MultiQC/" + it.getName() }

    // Gather locations of files to download
    locations = ALIGN_DEDUP_QUANT.out.file_locations
                    .mix(report_locations, ALIGNMENT_QC.out.file_locations, COMPARISON.out.file_locations)
                    .collectFile( name: "${params.outdir}/download_data/file_locations.txt", newLine: true )
    
    // Parse the list of files for downloading into a JSON file
    summarize_downloads( locations, md5sums, check_design.out.checked_design )
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

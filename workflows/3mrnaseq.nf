#!/usr/bin/env nextflow
/*
Workflow of Zymo-Seq 3' mRNA-Seq
*/

params.summary = [:]

// Load functions
include { parse_presets } from ('../libs/parse_presets')
include { setup_channel } from ('../libs/setup_channel')
include { check_star_log } from ('../libs/check_star_log')
include { parse_design; parse_design as parse_design_original_fastq } from ('../libs/parse_design')

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
    exit 1, "--genome is a required input!"
}

// Reference index path configuration
// Define these here - after the profiles are loaded with the iGenomes paths
params.star_index = params.genomes[ params.genome ].star ?: false
params.gtf = params.genomes[ params.genome ].gtf ?: false
params.bed12 = params.genomes[ params.genome ].bed12 ?: false
params.gprofiler_organism = params.genomes[ params.genome ].gprofiler ?: false
params.ensembl_web = params.genomes[ params.genome ].ensembl_web ?: false
params.rRNA_gtf = params.genomes[ params.genome ].rRNA_gtf ?: false
params.bacteria = params.genomes[ params.genome ].bacteria ?: false
params.csi_index = params.genomes[ params.genome ].csi_index ?: false

// Parse presets.
params.protocol_settings = parse_presets(params.protocol)
// Ensure a value is set up for params.ignore_R1
params.ignore_R1 = (params.protocol == "zymo_3mrna_nodedup")
// Ensure original bams are saved
params.save_bam = true

/*
 * SET & VALIDATE INPUT CHANNELS
 */
star_index = setup_channel(params.star_index, "STAR index", true, "")
gtf = setup_channel(params.gtf, "GTF annotation file", true, "")
bed12 = setup_channel(params.bed12, "BED annotation file", true, "")
rRNA_gtf = setup_channel(params.rRNA_gtf, "rRNA GTF file", false, "will depend on genome rRNA annotation.")
ercc_star_index = setup_channel(params.ercc_spikein ? params.genomes['ERCC92'].star : false, "ERCC STAR index", false, "will not look for ERCC reads.")
ercc_gtf = setup_channel(params.ercc_spikein ? params.genomes['ERCC92'].gtf : false, "ERCC GTF file", false, "will not look for ERCC reads.")
design = setup_channel(params.design, "design CSV file", true, "")
comparisons = setup_channel(params.comparisons, "comparison CSV file", false, "all pairwise comparisons will be carried out.")
multiqc_config = setup_channel(params.multiqc_config_3mrna, "MultiQC config", true, "")
ercc_info = setup_channel(params.ercc_spikein ? "$baseDir/assets/ERCC92_info.tsv" : false, "ERCC info", false, "ERCC count will not be conducted.")
summary_header = Channel.fromPath("$baseDir/assets/workflow_summary_header.txt", checkIfExists: true)

/*
 * COLLECT SUMMARY & LOG
 */
// Create a channel for settings to show in MultiQC report
info_to_report = ['Genome', 'ERCC spike-in', 'Library Prep', 'Strandedness', 'Trimming', 'DESeq2 FDR cutoff', 'DESeq2 Log2FC cutoff', 'gProfiler FDR cutoff']
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
include { check_design } from '../processes/check_design' addParams(
    ignore_R1: params.ignore_R1
)
include { fastqc } from '../processes/fastqc' addParams(
    publish_dir: "${params.outdir}/FastQC",
    skip_qc: params.skip_qc,
    skip_fastqc: params.skip_fastqc
)
include { umi_extract } from '../processes/umitools_extract' addParams(
    publish_dir: "${params.outdir}/umiextract"
)
include { trim_galore } from '../processes/trim_galore' addParams(
    publish_dir: "${params.outdir}/Trim_Galore",
    save_trimmed: params.save_trimmed,
    protocol_settings: params.protocol_settings,
    trim_nextseq: params.trim_nextseq,
    min_read_length: params.min_read_length,
    adapter_overlap: params.adapter_overlap,
    skip_trimming: params.skip_trimming
)
include { star } from '../processes/star' addParams(
    publish_dir: "${params.outdir}/STAR",
    save_unaligned: (params.save_unaligned || params.ercc_spikein), // Also save unaligned if ERCC alignment is requested
    save_secondary_alignments: params.save_secondary_alignments,
    save_bam: params.save_bam,
    csi_index: params.csi_index,
    bacteria: params.bacteria,
    star_twopass: params.star_twopass,
    star_min_overlap: params.star_min_overlap,
    star_max_overlap_mismatch: params.star_max_overlap_mismatch
)
include { alignERCC } from '../processes/star' addParams(
    publish_dir: "${params.outdir}/count_ERCC",
    star_min_overlap: params.star_min_overlap,
    star_max_overlap_mismatch: params.star_max_overlap_mismatch
)
include { umi_dedup } from '../processes/umitools_dedup' addParams(
    publish_dir: "${params.outdir}/umidedup",
    csi_index: params.csi_index
)
include { umi_dedup as umi_dedup_ercc } from '../processes/umitools_dedup' addParams(
    publish_dir: "${params.outdir}/count_ERCC",
    csi_index: false
)
include { count_ercc } from '../processes/count_ercc' addParams(
    publish_dir: "${params.outdir}/count_ERCC",
    strandedness: params.protocol_settings.strandedness,
    ercc_spikein: params.ercc_spikein
)
include { rseqc } from '../processes/rseqc' addParams(
    publish_dir: "${params.outdir}/RSeQC",
    skip_qc: params.skip_qc,
    skip_rseqc: params.skip_rseqc,
    csi_index: params.csi_index,
    bacteria: params.bacteria
)
include { preseq } from '../processes/preseq' addParams(
    publish_dir: "${params.outdir}/preseq",
    skip_qc: params.skip_qc,
    skip_preseq: params.skip_preseq,
    preseq_defect_filesize: params.preseq_defect_filesize
)
include { qualimap } from '../processes/qualimap' addParams(
    publish_dir: "${params.outdir}/Qualimap",
    skip_qc: params.skip_qc,
    skip_qualimap: params.skip_qualimap,
    strandedness: params.protocol_settings.strandedness
)
include { featurecounts } from '../processes/featurecounts' addParams(
    publish_dir: "${params.outdir}/featureCounts",
    strandedness: params.protocol_settings.strandedness,
    fc_extra_attributes: params.fc_extra_attributes,
    fc_group_features: params.fc_group_features,
    fc_count_type: params.fc_count_type,
    fc_biotype_count_type: params.fc_biotype_count_type,
    fc_biotype_group_features: params.fc_biotype_group_features,
    filter_primary_alignments: params.save_secondary_alignments // When secondary alignments were saved, it needs to be filtered before featureCounts
)
include { merge_featurecounts } from '../processes/merge_featurecounts' addParams(
    publish_dir: "${params.outdir}/featureCounts",
    bam_suffix: ".dedupped.bam"
)
include { parse_biotype_qc } from '../processes/parse_biotype_qc' addParams(
    skip_biotype_qc: params.skip_biotype_qc
)
include { count_genes_detected } from '../processes/count_genes_detected' addParams(
    publish_dir: "${params.outdir}/featureCounts",
    gene_detection_method: params.gene_detection_method
)
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
include {parse_dedup_stat} from '../processes/parse_umitools_dedup'
include { software_versions } from "../processes/software_versions" addParams(
    publish_dir: "${params.outdir}/pipeline_info",
    pipeline_version: workflow.manifest.version,
    nextflow_version: workflow.nextflow.version
)
include { multiqc_3mrna } from "../processes/multiqc" addParams(
    publish_dir: "${params.outdir}/MultiQC",
    skip_multiqc: params.skip_multiqc,
    run_name: params.summary["Run Name"],
    ensembl_web: params.ensembl_web,
    deseq2_fdr: params.deseq2_fdr,
    gprofiler_fdr: params.gprofiler_fdr,
    ignore_R1: params.ignore_R1
)
include { summarize_downloads } from "../processes/summarize_downloads" addParams(
    publish_dir: "${params.outdir}/download_data"
)

/*
 * WORKFLOW DEFINITION
 */
def fail_percent_mapped = [:]
workflow RNASEQ3M {
    // Check design file and set up channels
    check_design(design, comparisons.ifEmpty([]))
    check_design.out.checked_design
        .splitCsv( header: true )
        .map { parse_design(it, params.ignore_R1) }
        .set { reads } // Channel of [meta, reads]

    // Read QC
    fastqc(reads)

    // Extract UMI from Read 1 and add to Read 2. The value for meta.single_end is also updated to be true
    umi_extract(reads)

    // Trim only UMI-added Read 2
    trim_galore(umi_extract.out.umiextracted_reads)

    // Alignment and filtering based on alignment rate
    if (params.skip_trimming) {
        star(umi_extract.out.umiextracted_reads, star_index.collect())
    } else {
        star(trim_galore.out.reads, star_index.collect())
    }
    star.out.bam
        .branch { meta, star_log, bam, bai ->
            pass: check_star_log(star_log, params.percent_mapped_cutoff, log)
                return [meta, bam, bai]
            fail: !check_star_log(star_log, params.percent_mapped_cutoff, log)
                fail_percent_mapped[meta.name] = true
                return meta.name
        }
        .set { star_out_filtered }

    // Dudup bams based on UMI information and collect dedup stats for multiqc report
    umi_dedup(star_out_filtered.pass)
    parse_dedup_stat(umi_dedup.out.umidedup_log.collect())

    // Align, dedup, and count ERCC
    alignERCC(star.out.unmapped, ercc_star_index.collect())
    umi_dedup_ercc(alignERCC.out.bam)
    count_ercc(umi_dedup_ercc.out.bam_dedupped, ercc_gtf.collect(), ercc_info.collect())

    // BAM QC
    rseqc(umi_dedup.out.bam_dedupped, bed12.collect())
    preseq(umi_dedup.out.bam_dedupped)
    qualimap(umi_dedup.out.bam_dedupped, gtf.collect())

    // Read counting
    featurecounts(umi_dedup.out.bam_dedupped, gtf.collect(), rRNA_gtf.collect().ifEmpty([]))
    merge_featurecounts(featurecounts.out.counts.collect())
    parse_biotype_qc(featurecounts.out.biotype_counts.collect(), featurecounts.out.biotype_log.collect())
    count_genes_detected(merge_featurecounts.out.merged_counts, merge_featurecounts.out.gene_lengths)

    // Comparisons
    deseq2(merge_featurecounts.out.merged_counts, check_design.out.checked_design, comparisons.ifEmpty([]))
    gprofiler(deseq2.out.results)

    // Software versions
    versions = fastqc.out.version.first()
                   .mix(umi_extract.out.version.first(), trim_galore.out.version.first(), star.out.version.first(),
                        preseq.out.version.first(), rseqc.out.version.first(), qualimap.out.version.first(), featurecounts.out.version.first(),
                        deseq2.out.version, gprofiler.out.version)
    software_versions(versions.flatten().collect())

    // Report
    multiqc_3mrna(multiqc_config, \
            fastqc.out.report.collect().ifEmpty([]), \
            trim_galore.out.report.collect().ifEmpty([]), \
            star.out.report.collect(), \
            parse_dedup_stat.out.dedup_stats.collect().ifEmpty([]), \
            count_ercc.out.report.collect().ifEmpty([]), \
            rseqc.out.report.collect().ifEmpty([]), \
            preseq.out.report.collect().ifEmpty([]), \
            qualimap.out.report.collect().ifEmpty([]), \
            featurecounts.out.report.collect(), \
            parse_biotype_qc.out.report.collect().ifEmpty([]), \
            count_genes_detected.out.report.collect(), \
            deseq2.out.report.mix(deseq2.out.results).collect().ifEmpty([]), \
            gprofiler.out.report.collect().ifEmpty([]), \
            software_versions.out.report.collect(), \
            summary_header, \
            workflow_summary_to_report)
    
    // Gather locations of files to download
    bam_locations = star.out.bam
                            .map { meta, star_log, bam, bai -> "${params.outdir}/STAR/" + bam.getName() }
    counts_locations = merge_featurecounts.out.merged_counts.map { "${params.outdir}/featureCounts/" + it.getName() }
    deseq2_locations = deseq2.out.download.flatten().map { "${params.outdir}/DESeq2/" + it.getName() }
    gprofiler_locations = gprofiler.out.download.flatten().map { "${params.outdir}/gProfiler/" + it.getName() }
    report_locations = multiqc.out.report.map { "${params.outdir}/MultiQC/" + it.getName() }
    bam_locations
        .mix(counts_locations, deseq2_locations, gprofiler_locations, report_locations)
        .collectFile(name: "${params.outdir}/download_data/file_locations.txt", newLine: true )
        .set { locations }
    
    // Save md5sum into one file
    star.out.bam_md5sum
        .collectFile(name: "${params.outdir}/download_data/md5sum.txt", sort: { it.getName() } )
        .set { md5sums }
    
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
    if (fail_percent_mapped.size() > 0) {
        log.warn "[rnaseq] WARNING - ${fail_percent_mapped.size()} samples skipped due to poor alignment scores!"
    }
    if (workflow.success) {
        log.info "[rnaseq] Pipeline completed successfully"
    } else {
        log.warn "[rnaseq] Pipeline completed with errors"
    }
}
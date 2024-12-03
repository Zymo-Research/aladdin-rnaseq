// Run alignment, deduplication, and quantification
// This can be used for the reference genome, ERCC spike-in, or others such as a 2nd reference genome

include { check_star_log } from ('../libs/check_star_log')
include { setup_channel } from ('../libs/setup_channel')
include { make_transcripts } from '../processes/make_transcripts'
include { star } from '../processes/star' addParams(
    publish_dir: "${params.outdir}/STAR",
    save_unaligned: params.save_unaligned || params.ercc_spikein,
    save_secondary_alignments: params.save_secondary_alignments,
    csi_index: params.genomes[ params.genome ].csi_index ?: false,
    bacteria: params.genomes[ params.genome ].bacteria ?: false,
    star_twopass: params.star_twopass,
    star_min_overlap: params.star_min_overlap,
    star_max_overlap_mismatch: params.star_max_overlap_mismatch,
    read_quant_method: params.read_quant_method
)
include { umi_dedup } from '../processes/umitools_dedup' addParams(
    publish_dir: "${params.outdir}/umidedup",
    csi_index: params.genomes[ params.genome ].csi_index ?: false
)
include {parse_dedup_stat} from '../processes/parse_umitools_dedup'
include { filter_transcriptome_bam } from '../processes/filter_transcriptome_bam' addParams(
    publish_dir: "${params.outdir}/umidedup"
)
include { featurecounts } from '../processes/featurecounts' addParams(
    publish_dir: "${params.outdir}/featureCounts",
    strandedness: params.protocol_settings.strandedness,
    fc_extra_attributes: params.fc_extra_attributes,
    fc_group_features: params.fc_group_features,
    fc_count_type: params.fc_count_type,
    read_quant_method: params.read_quant_method
)
include { salmon_quant } from '../processes/salmon_quant' addParams(
    publish_dir: "${params.outdir}/Salmon",
    read_quant_method: params.read_quant_method,
    strandedness: params.protocol_settings.strandedness
)
include { summarize_samples } from '../processes/summarize_samples' addParams(
    publish_dir: "${params.outdir}/summarize_samples",
    read_quant_method: params.read_quant_method,
    gene_detection_method: params.gene_detection_method,
    fc_extra_attributes: params.fc_extra_attributes
)

workflow ALIGN_DEDUP_QUANT {
    take:
    reads // [meta, reads]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_file_locations = Channel.empty()

    // Setup channels related to the genome
    ch_index = setup_channel(params.genomes[ params.genome ].star ?: false, "STAR index", true, "")
    ch_gtf = setup_channel(params.genomes[ params.genome ].gtf ?: false, "GTF annotation file", true, "")

    // Prepare transcripts file if using Salmon
    if (params.read_quant_method=='STAR_Salmon') {
        if (params.genomes[ params.genome ].transcripts) {
            ch_transcripts = setup_channel(params.genomes[ params.genome ].transcripts, "transcripts FASTA file", true, "")
        } else {
            ch_genome = setup_channel(params.genomes[ params.genome ].fasta, "genome FASTA file", true, "")
            make_transcripts(ch_genome, ch_gtf)
            ch_transcripts = make_transcripts.out.transcripts
        }
    } else {
        ch_transcripts = Channel.empty()
    }
    
    // Align to reference
    star(reads, ch_index.collect())
    ch_versions = ch_versions.mix(star.out.version.first())
    ch_multiqc_files = ch_multiqc_files.mix(star.out.report)

    // Filter STAR results
    star.out.bam
        .branch { meta, star_log, bam, bai ->
            pass: check_star_log(star_log, params.percent_mapped_cutoff, log)
                return [meta, bam, bai]
            fail: !check_star_log(star_log, params.percent_mapped_cutoff, log)
                return meta.name
        }
        .set { star_out_filtered }
    // Also filter transcriptome BAM
    star.out.transcriptome_bam
        .join(star_out_filtered.pass)
        .map { [it[0], it[1], it[2]] }
        .set { star_out_filtered_transcriptome }
    // Get passed BAM locations
    bam_locations = star_out_filtered.pass
                        .map { meta, bam, bai -> "${params.outdir}/STAR/" + bam.getName() }
    ch_file_locations = ch_file_locations.mix(bam_locations)
    // Get a warning for failed samples
    star_out_filtered.fail
        .collect()
        .map {
            "The following samples were excluded from analysis due to alignment rate < ${params.percent_mapped_cutoff}%:\n${it.join("; ")}"
        }
        .set { ch_warning_message }
    ch_warning_message.subscribe { log.warn "$it" }

    if (params.protocol_settings.umi) {
        // Dedup bams based on UMI information and collect dedup stats for multiqc report
        umi_dedup(star_out_filtered.pass)
        parse_dedup_stat(umi_dedup.out.umidedup_log.map{ it[1] }.collect())
        ch_versions = ch_versions.mix(umi_dedup.out.version.first())
        ch_multiqc_files = ch_multiqc_files.mix(parse_dedup_stat.out.dedup_stats)
        ch_bam = umi_dedup.out.bam_dedupped
        // Dedup transcriptome bam if Salmon is the quant method
        if (params.read_quant_method == 'STAR_Salmon') {
            star_out_filtered_transcriptome
                .join(umi_dedup.out.bam_dedupped)
                .map { [it[0], it[2], it[3]] }
                .set { filter_transcriptome_bam_input }
            filter_transcriptome_bam(filter_transcriptome_bam_input)
            ch_transcriptome_bam = 
                umi_dedup.out.umidedup_log
                    .join(filter_transcriptome_bam.out.filtered_bam)
            ch_versions = ch_versions.mix(filter_transcriptome_bam.out.version.first())
        } else {
            ch_transcriptome_bam = star_out_filtered_transcriptome
        }
    } else {
        ch_bam = star_out_filtered.pass
        ch_transcriptome_bam = star_out_filtered_transcriptome
    }

    // Use featureCounts to count reads, these steps will run when params.read_quant_method == 'STAR_featureCounts'
    featurecounts(ch_bam, ch_gtf.collect())
    // Use Salmon to count reads, these steps will run when params.read_quant_method == 'STAR_Salmon'
    salmon_quant(ch_transcriptome_bam, ch_gtf.collect(), ch_transcripts.collect())
    // Summarize all samples
    summarize_samples(featurecounts.out.counts.mix(salmon_quant.out.results).collect())
    counts_locations = summarize_samples.out.counts
                        .mix(summarize_samples.out.transcript_counts)
                        .map { "${params.outdir}/summarize_samples/" + it.getName() }

    ch_versions = ch_versions.mix(featurecounts.out.version.first(), salmon_quant.out.version.first())
    ch_multiqc_files = ch_multiqc_files.mix(featurecounts.out.report, salmon_quant.out.report, summarize_samples.out.report)
    ch_file_locations = ch_file_locations.mix(counts_locations)

    emit:
    versions       = ch_versions
    report         = ch_multiqc_files
    file_locations = ch_file_locations
    bam            = ch_bam
    unmapped       = star.out.unmapped
    md5sum         = star.out.md5sum
    counts         = summarize_samples.out.counts
    tpm            = summarize_samples.out.tpm
    salmon_results = salmon_quant.out.results
    warning        = ch_warning_message
}
// Conduct various library QC after the alignment step

// Load functions
include { setup_channel } from ('../libs/setup_channel')
include { rseqc } from '../processes/rseqc' addParams(
    publish_dir: "${params.outdir}/RSeQC",
    skip_rseqc: params.skip_rseqc,
    csi_index: params.genomes[ params.genome ].csi_index ?: false,
    bacteria: params.genomes[ params.genome ].bacteria ?: false
)
include { preseq } from '../processes/preseq' addParams(
    publish_dir: "${params.outdir}/preseq",
    skip_preseq: params.skip_preseq
)
include { mark_duplicates } from '../processes/mark_duplicates' addParams(
    publish_dir: "${params.outdir}/MarkDuplicates",
    skip_markduplicates: params.skip_markduplicates || params.protocol_settings.umi,
    csi_index: params.genomes[ params.genome ].csi_index ?: false
)
include { dupradar } from '../processes/dupradar' addParams(
    publish_dir: "${params.outdir}/dupRadar",
    skip_dupradar: params.skip_dupradar,
    strandedness: params.protocol_settings.strandedness
)
include { qualimap } from '../processes/qualimap' addParams(
    publish_dir: "${params.outdir}/Qualimap",
    skip_qualimap: params.skip_qualimap,
    strandedness: params.protocol_settings.strandedness
)
include { bamscale } from '../processes/bamscale' addParams(
    publish_dir: "${params.outdir}/BAMscale",
    strandedness: params.protocol_settings.strandedness,
    generate_bigwig: params.generate_bigwig
)
include { featurecounts_biotype } from '../processes/featurecounts' addParams(
    publish_dir: "${params.outdir}/featureCounts",
    strandedness: params.protocol_settings.strandedness,
    fc_biotype_count_type: params.fc_biotype_count_type,
    fc_biotype_group_features: params.fc_biotype_group_features,
    skip_biotype_qc: params.skip_biotype_qc
)
include { parse_biotype_qc } from '../processes/parse_biotype_qc' addParams(
    skip_biotype_qc: params.skip_biotype_qc
)

workflow ALIGNMENT_QC {
    take:
    bam

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_file_locations = Channel.empty()

    if (!params.skip_qc) {
        ch_gtf = setup_channel(params.genomes[ params.genome ].gtf ?: false, "GTF annotation file", true, "")
        ch_bed12 = setup_channel(params.genomes[ params.genome ].bed12 ?: false, "BED annotation file", true, "")
        ch_rRNA_gtf = setup_channel(params.genomes[ params.genome ].rRNA_gtf ?: false, "rRNA GTF file", false, "will depend on genome rRNA annotation.")
        
        // BAM QC
        rseqc(bam, ch_bed12.collect())
        preseq(bam)
        mark_duplicates(bam)
        dupradar(mark_duplicates.out.bam, ch_gtf.collect())
        qualimap(bam, ch_gtf.collect())
        bamscale(bam)

        // Use featureCounts for biotype QC regardless of read_quant_method
        featurecounts_biotype(bam, ch_gtf.collect(), ch_rRNA_gtf.collect().ifEmpty([]))
        parse_biotype_qc(featurecounts_biotype.out.biotype_counts.collect(), featurecounts_biotype.out.biotype_log.collect())

        ch_versions =  ch_versions.mix(
            rseqc.out.version.first(), preseq.out.version.first(), mark_duplicates.out.version.first(),
            dupradar.out.version.first(), qualimap.out.version.first(), bamscale.out.version.first(),
            featurecounts_biotype.out.version.first()
        )

        ch_multiqc_files = ch_multiqc_files.mix(
            rseqc.out.report, preseq.out.report, mark_duplicates.out.report,
            dupradar.out.report, qualimap.out.report, parse_biotype_qc.out.report
        )

        ch_file_locations = ch_file_locations.mix(
            bamscale.out.download.flatten().map { "${params.outdir}/BAMscale/" + it.getName() }
        )
    }

    emit:
    versions       = ch_versions
    report         = ch_multiqc_files
    file_locations = ch_file_locations
}
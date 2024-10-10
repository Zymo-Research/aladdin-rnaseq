// Compare samples and groups
// Run DEG, pathway enrichment, DTU analysis

include { setup_channel } from ('../libs/setup_channel')
include { tx2gene } from '../processes/tx2gene' addParams(
    fc_extra_attributes: params.fc_extra_attributes
)
include { deseq2 } from '../processes/deseq2' addParams(
    publish_dir: "${params.outdir}/DESeq2",
    skip_deseq2: params.skip_deseq2,
    read_quant_method: params.read_quant_method,
    deseq2_fdr: params.deseq2_fdr,
    deseq2_lfc: params.deseq2_lfc,
    heatmap_group_order: params.heatmap_group_order
)
include { gprofiler } from '../processes/gprofiler' addParams(
    publish_dir: "${params.outdir}/gProfiler",
    gprofiler_organism: params.genomes[ params.genome ].gprofiler ?: false,
    deseq2_fdr: params.deseq2_fdr,
    gprofiler_fdr: params.gprofiler_fdr
)
include { dtu_analysis } from '../processes/dtu_analysis' addParams(
    publish_dir: "${params.outdir}/DTU_analysis",
    dtu_analysis: params.dtu_analysis,
    dexseq_fdr: params.dexseq_fdr,
    prop_filter_transcript_counts: params.prop_filter_transcript_counts,
    prop_filter_transcript_props: params.prop_filter_transcript_props,
    prop_filter_gene_counts: params.prop_filter_gene_counts
)

workflow COMPARISON {
    take:
    counts
    salmon_results
    design
    comparisons

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_file_locations = Channel.empty()

    if (params.read_quant_method=='STAR_Salmon') {
        ch_gtf = setup_channel(params.genomes[ params.genome ].gtf ?: false, "GTF annotation file", true, "")
        tx2gene(ch_gtf)
        ch_tx2gene = tx2gene.out.mapping
    } else {
        ch_tx2gene = Channel.empty()
    }

    deseq2(counts.ifEmpty([]), salmon_results.collect().ifEmpty([]), design, comparisons.ifEmpty([]), ch_tx2gene.ifEmpty([]))
    gprofiler(deseq2.out.results)
    dtu_analysis(salmon_results.collect().ifEmpty([]), design, ch_tx2gene.ifEmpty([]), comparisons.ifEmpty([]))

    ch_versions = ch_versions.mix(
        deseq2.out.version, gprofiler.out.version, dtu_analysis.out.version
        )
    ch_multiqc_files = ch_multiqc_files.mix(
        deseq2.out.report, deseq2.out.results, gprofiler.out.report, dtu_analysis.out.report
        )
    deseq2_locations = deseq2.out.download.flatten().map { "${params.outdir}/DESeq2/" + it.getName() }
    gprofiler_locations = gprofiler.out.download.flatten().map { "${params.outdir}/gProfiler/" + it.getName() }
    dtu_locations = dtu_analysis.out.download.flatten().map { "${params.outdir}/DTU_analysis/" + it.getName() }
    ch_file_locations = ch_file_locations.mix(
        deseq2_locations, gprofiler_locations, dtu_locations
        )

    emit:
    versions       = ch_versions
    report         = ch_multiqc_files
    file_locations = ch_file_locations
}
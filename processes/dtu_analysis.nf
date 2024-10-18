// Run differential transcript usage (DTU) analysis using DEXSeq along with DRIMSeq and stageR
params.publish_dir = "DTU_analysis"
params.dexseq_fdr = 0.05
params.prop_filter_transcript_counts = 0.5
params.prop_filter_transcript_props = 0.5
params.prop_filter_gene_counts = 1

process dtu_analysis {
    label 'low_memory'
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it.endsWith('.xlsx') ? it : null }

    when:
    params.dtu_analysis

    input:
    path ('salmon/*')
    path design_csv
    path ('tx2gene.tsv')
    path comparisons

    output:
    path "*_DTU_analysis_DEXSeq_results.xlsx", emit: download
    path "*_DTU_analysis_DEXSeq_results.tsv", emit: report
    path "v_*.txt", emit: version

    script:
    """
    DTU.r salmon $design_csv ${params.dexseq_fdr} ${params.prop_filter_transcript_counts} ${params.prop_filter_transcript_props} ${params.prop_filter_gene_counts} $comparisons
    Rscript -e "write(x=as.character(packageVersion('DEXSeq')), file='v_DEXSeq.txt')"
    Rscript -e "write(x=as.character(packageVersion('DRIMSeq')), file='v_DRIMSeq.txt')"
    Rscript -e "write(x=as.character(packageVersion('stageR')), file='v_stageR.txt')"
    """
}
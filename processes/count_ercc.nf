// Count ERCC reads
params.publish_dir = "count_ERCC"
params.strandedness = 0
params.ercc_spikein = false

process count_ercc {
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
    tuple val(meta), path (bam), path(bai) // bam index files are not necessary for counting; just to fit with the input channel
    path gtf
    path info

    output:
    path "*.featureCounts.txt*"
    path "*.featureCounts.txt.summary", emit: ercc
    path "*plot_ERCC.csv", emit: report
    
    script:
    // Concentration for mix 1 is in 4th column, mix 2 in 5th column of ERCC92 info file
    conc_col = params.ercc_spikein.toInteger() + 3
    concentration = "<(tail -n +2 $info | cut -f2,$conc_col | sort | cut -f2)"
    paired_end = meta.single_end ? '' : '-p'
    """
    featureCounts -a $gtf -g gene_id -t exon -o ${meta.name}_ERCC.featureCounts.txt ${paired_end} -s ${params.strandedness} $bam
    echo 'ID\tlength\tcount\tconcentration' > ${meta.name}_plot_ERCC.csv
    paste <(tail -n +3 ${meta.name}_ERCC.featureCounts.txt | cut -f1,6,7 | sort) $concentration >> ${meta.name}_plot_ERCC.csv
    """
}
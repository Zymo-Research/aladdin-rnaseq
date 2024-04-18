// Run genome coverage by bamscale
params.strandedness = 0
params.publish_dir = "bamscale"

process bamscale {
    tag "${meta.name}"
    label 'low_memory'
    errorStrategy 'retry'
    maxRetries 1
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it.endsWith('.bw') ? it : null }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "*.bw"          , emit: download
    path "v_BAMscale.txt", emit: version

    when:
    params.generate_bigwig

    script:
    def strand = params.strandedness == 0 ? "rna" : "strandrnaR"
    """
    SCALE_FACTOR=\$(grep 'primary mapped (' <(samtools flagstat $bam) | awk '{print 1000000/\$1}')
    BAMscale scale --bam $bam --operation $strand --scale custom --factor \$SCALE_FACTOR --binsize 1 --maxfrag 1000000 -t $task.cpus
    BAMscale scale &> v_BAMscale.txt || true
    """
}
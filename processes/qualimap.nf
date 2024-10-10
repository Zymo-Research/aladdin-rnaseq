// Run Qualimap
params.publish_dir = "Qualimap"
params.skip_qualimap = false
params.strandedness = 0

process qualimap {
    label 'low_memory'
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it=='v_qualimap.txt' ? null : it }

    when:
    !params.skip_qualimap

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf

    output:
    path "${meta.name}", emit: report
    path "v_qualimap.txt", emit: version

    script:
    qualimap_direction_text = [ 0:'non-strand-specific', 1:'strand-specific-forward', 2:'strand-specific-reverse' ]
    qualimap_direction = qualimap_direction_text[params.strandedness]
    paired = meta.single_end ? '' : '-pe'
    memory = task.memory.toGiga() + "G"
    """
    frac=\$(samtools idxstats $bam | awk -F '\t' '{s+=\$3} END {frac=10000000/s; if (frac>1) {print "1.0"} else {print frac}}')
    if [ "\$frac" == "1.0" ]; then
        unset DISPLAY
        qualimap --java-mem-size=${memory} rnaseq -p $qualimap_direction $paired -bam $bam -gtf $gtf -outdir ${meta.name}
    else
        samtools view -bs \$frac $bam > subsample.bam
        unset DISPLAY
        qualimap --java-mem-size=${memory} rnaseq -p $qualimap_direction $paired -bam subsample.bam -gtf $gtf -outdir ${meta.name}
    fi
    qualimap rnaseq > v_qualimap.txt 2>&1 || true
    """
}
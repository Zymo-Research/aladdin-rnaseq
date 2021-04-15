// Run Preseq
params.publish_dir = "preseq"
params.skip_qc = false
params.skip_preseq = false
params.preseq_defect_filesize = 10000000 // Run preseq in Defect mode if BAM file is less than this size

process preseq {
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it=='v_preseq.txt' ? null : it }

    when:
    !params.skip_qc && !params.skip_preseq

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "${meta.name}.ccurve.txt", emit: report
    path "v_preseq.txt", emit: version

    script:
    // Run preseq in Defect mode if the BAM file is too small, only relevant if using small test datasets
    if (bam.size() > params.preseq_defect_filesize) {
        """
        preseq lc_extrap -v -B $bam -o ${meta.name}.ccurve.txt
        preseq &> v_preseq.txt
        """
    } else {
        """
        preseq lc_extrap -v -D -B $bam -o ${meta.name}.ccurve.txt
        preseq &> v_preseq.txt
        """
    }
}
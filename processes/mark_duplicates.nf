// Run Picard MarkDuplicates
params.publish_dir = "MarkDuplicates"
params.skip_qc = false
params.skip_markduplicates = false
params.markdup_java_options = '"-Xms4000m -Xmx7g"'
params.csi_index = false

process mark_duplicates {
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("_metrics.txt") > 0) "metrics/$filename"
            else if (filename.indexOf("markDups.bam") > 0) filename
            else null
        }

    when:
    !params.skip_qc && !params.skip_markduplicates

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.name}.markDups.bam"), path("${meta.name}.markDups.bam.{bai,csi}"), emit: bam
    path "${meta.name}.markDups_metrics.txt", emit: report
    path "${meta.name}_bam_md5sum.txt", emit: md5sum optional true
    path "v_markduplicates.txt", emit: version

    script:
    markdup_java_options = (task.memory.toGiga() > 8) ? ${params.markdup_java_options} : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""
    dupsam_index = params.csi_index ? "samtools index -c ${meta.name}.markDups.bam" : "samtools index ${meta.name}.markDups.bam"
    """
    picard ${markdup_java_options} MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${meta.name}.markDups.bam \\
        METRICS_FILE=${meta.name}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT
    $dupsam_index
    md5sum ${meta.name}.markDups.bam > ${meta.name}_bam_md5sum.txt
    picard MarkDuplicates --version &> v_markduplicates.txt || true
    """
}
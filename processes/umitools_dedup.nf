// UMI-tools to deduplicate reads
params.publish_dir = "./umidedup"
params.csi_index = false

process umi_dedup {
    label 'mid_memory'
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("dedup_LOGFILE") > 0) "logs/$filename"
            else if (filename.indexOf("dedupped.bam") > 0) filename
            else null
        }

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("*dedupped.bam"), path("*.dedupped.bam.{bai, csi}"), emit: bam_dedupped
    path "*dedup_LOGFILE", emit: umidedup_log
    path "${meta.name}.md5sum.txt", emit: bam_md5sum optional true

    script:
    sam_dedupped_index = params.csi_index ? "samtools index -c ${meta.name}.dedupped.bam" : "samtools index ${meta.name}.dedupped.bam"
    """
    umi_tools dedup --random-seed=100 --spliced-is-unique --multimapping-detection-method=NH --output-stats=${meta.name}.dedup --stdin=$bams --log=${meta.name}.dedup_LOGFILE >${meta.name}.dedupped.bam
    $sam_dedupped_index
    """
}
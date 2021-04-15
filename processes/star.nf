// STAR alignment
params.publish_dir = "./STAR"
params.save_unaligned = false
params.save_secondary_alignments = false
params.save_bam = false
params.csi_index = false
params.bacteria = false
params.star_twopass = false
params.star_min_overlap = 10
params.star_max_overlap_mismatch = 0.01

process star {
    label 'high_memory'
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else if (params.save_unaligned && (filename.indexOf('Unmapped') > -1)) "unmapped/$filename"
            else if (params.save_bam) filename
            else null
        }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*Log.final.out"), path("*.bam"), path("*.{bai,csi}"), emit: bam
    path "*.out", emit: report
    path "*SJ.out.tab"
    tuple val(meta), path("*Unmapped*"), emit: unmapped optional true
    path "${meta.name}_bam_md5sum.txt", emit: md5sum optional true
    path "v_*.txt", emit: version

    script:
    avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    sam_index = params.csi_index ? "samtools index -c ${meta.name}Aligned.sortedByCoord.out.bam" : "samtools index ${meta.name}Aligned.sortedByCoord.out.bam"
    unaligned = params.save_unaligned ? "--outReadsUnmapped Fastx" : ''
    no_intron = params.bacteria ? "--alignIntronMax 1" : ''
    secondary = params.save_secondary_alignments ? '' : "--outSAMmultNmax 1"
    two_pass = params.star_twopass ? '--twopassMode Basic' : ''
    md5sum_cmd = params.save_bam ? "md5sum ${meta.name}Aligned.sortedByCoord.out.bam > ${meta.name}_bam_md5sum.txt": ''
    """
    STAR --genomeDir $index \\
        --readFilesIn $reads  \\
        --runThreadN ${task.cpus} \\
        --outSAMtype BAM SortedByCoordinate $avail_mem \\
        --readFilesCommand zcat $unaligned $no_intron $secondary $two_pass \\
        --outFileNamePrefix ${meta.name} \\
        --peOverlapNbasesMin ${params.star_min_overlap} \\
        --peOverlapMMp ${params.star_max_overlap_mismatch}

    $sam_index
    $md5sum_cmd
    STAR --version &> v_star.txt
    samtools --version &> v_samtools.txt
    """
}

process alignERCC {
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path('*.bam'), path('*.bam.bai'), emit: bam
    
    script:
    avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    STAR --genomeDir $index \\
        --readFilesIn $reads \\
        --runThreadN ${task.cpus} \\
        --outSAMtype BAM SortedByCoordinate $avail_mem \\
        --alignIntronMax 1 \\
        --outFileNamePrefix ${meta.name} \\
        --peOverlapNbasesMin ${params.star_min_overlap} \\
        --peOverlapMMp ${params.star_max_overlap_mismatch}
    samtools index ${meta.name}Aligned.sortedByCoord.out.bam
    """
}
// STAR alignment
params.publish_dir = "./STAR"
params.save_unaligned = false
params.save_secondary_alignments = false
params.csi_index = false
params.bacteria = false
params.star_twopass = false
params.star_min_overlap = 10
params.star_max_overlap_mismatch = 0.01
params.read_quant_method = 'STAR_featureCounts'

process star {
    label 'high_memory'
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename ->
            if (params.save_unaligned && filename.indexOf('unmapped') >= 0) filename
            else if (filename.indexOf('.genome.primary.bam') > 0) filename
            else if (filename.indexOf('toTranscriptome.out.bam') > 0) "transcriptome_bam/$filename"
            else if (params.save_secondary_alignments && filename.indexOf('Aligned.sortedByCoord.out.bam') > 0) "secondary_alignments/$filename"
            else if (filename.indexOf('.out') > 0 && filename.indexOf('.bam') < 0) "logs/$filename"
            else null
        }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*Log.final.out"), path("*.genome.primary.bam"), path("*.{bai,csi}"), emit: bam
    path "*.out", emit: report
    path "*SJ.out.tab"
    path "*Aligned.sortedByCoord.out.bam"
    tuple val(meta), path("unmapped/*unmapped*.fastq.gz"), emit: unmapped optional true
    tuple val(meta), path("*Log.final.out"), path("*.toTranscriptome.out.bam"), emit: transcriptome_bam optional true
    path "${meta.name}_bam_md5sum.txt", emit: md5sum optional true
    path "v_*.txt", emit: version

    script:
    avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    csi_index = params.csi_index ? '-c' : ''
    unaligned = params.save_unaligned ? '--outReadsUnmapped Fastx' : ''
    no_intron = params.bacteria ? '--alignIntronMax 1' : ''
    two_pass = params.star_twopass ? '--twopassMode Basic' : ''
    transcriptome_bam = params.read_quant_method == 'STAR_Salmon' ? '--quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend' : ''
    if (params.save_unaligned) {
        mkdir_unaligned = 'mkdir unmapped'
        gzip_unaligned_r1 = "cat ${meta.name}Unmapped.out.mate1 | gzip > unmapped/${meta.name}.unmapped_R1.fastq.gz"
        if (!meta.single_end) {
            gzip_unaligned_r2 = "cat ${meta.name}Unmapped.out.mate2 | gzip > unmapped/${meta.name}.unmapped_R2.fastq.gz"
        } else {
            gzip_unaligned_r2 = ''
        }
    } else {
        mkdir_unaligned = ''
        gzip_unaligned_r1 = ''
        gzip_unaligned_r2 = ''
    }
    """
    STAR --genomeDir $index \\
        --readFilesIn $reads  \\
        --runThreadN ${task.cpus} \\
        --outSAMtype BAM SortedByCoordinate $avail_mem \\
        --readFilesCommand zcat $unaligned $no_intron $two_pass $transcriptome_bam \\
        --outFileNamePrefix ${meta.name} \\
        --peOverlapNbasesMin ${params.star_min_overlap} \\
        --peOverlapMMp ${params.star_max_overlap_mismatch}

    samtools view -F 256 -b ${meta.name}Aligned.sortedByCoord.out.bam > ${meta.name}.genome.primary.bam
    md5sum ${meta.name}.genome.primary.bam > ${meta.name}_bam_md5sum.txt
    samtools index ${csi_index} ${meta.name}.genome.primary.bam

    $mkdir_unaligned
    $gzip_unaligned_r1
    $gzip_unaligned_r2

    STAR --version &> v_star.txt
    samtools --version &> v_samtools.txt
    """
}

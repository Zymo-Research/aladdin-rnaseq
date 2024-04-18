// Filter transcriptome BAM based on reads kept in the dedupped genome BAM
params.publish_dir = "./umidedup"

process filter_transcriptome_bam {
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it=='v_filtersamreads.txt' ? null : it }

    input:
    tuple val(meta), path(transcriptome_bam), path(genome_bam)

    output:
    tuple val(meta), path("${meta.name}.filtered.transcriptome.bam"), emit: filtered_bam
    path "v_picard.txt", emit: version

    script:
    """
    samtools view $genome_bam | cut -f1 > read_names.txt
    picard FilterSamReads \\
        INPUT=$transcriptome_bam \\
        OUTPUT=${meta.name}.filtered.transcriptome.bam \\
        READ_LIST_FILE=read_names.txt \\
        FILTER=includeReadList
    picard FilterSamReads --version &> v_picard.txt || true
    """
}
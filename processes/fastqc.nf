// FastQC before trimming
params.publish_dir = "FastQC"
params.skip_fastqc = false
params.skip_qc = false

process fastqc {
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it=='v_fastqc.txt' ? null : it }

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    tuple val(meta), path(reads)

    output:
    path '*_fastqc.{zip,html}', emit: report
    path 'v_fastqc.txt', emit: version

    script:
    if (meta.single_end) {
        """
        fastqc --version &> v_fastqc.txt
        [ ! -f  ${meta.name}_R1.fastq.gz ] && ln -s $reads ${meta.name}_R1.fastq.gz
        fastqc --quiet --threads $task.cpus ${meta.name}_R1.fastq.gz
        """
    } else {
        """
        fastqc --version &> v_fastqc.txt
        [ ! -f  ${meta.name}_R1.fastq.gz ] && ln -s ${reads[0]} ${meta.name}_R1.fastq.gz
        [ ! -f  ${meta.name}_R2.fastq.gz ] && ln -s ${reads[1]} ${meta.name}_R2.fastq.gz
        fastqc --quiet --threads $task.cpus ${meta.name}_R1.fastq.gz
        fastqc --quiet --threads $task.cpus ${meta.name}_R2.fastq.gz
        """
    }
}
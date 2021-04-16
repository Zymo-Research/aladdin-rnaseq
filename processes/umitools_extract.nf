// UMI-tools to extract UMI from Read 1
// meta.single_end is modifed to true to activate correct settings for subsequent analysis
params.publish_dir = "./umiextract"

process umi_extract {
    errorStrategy 'retry'
    maxRetries 1
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { filename -> if (filename.indexOf("umiextract_report.txt") > 0) "logs/$filename" else null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(new_meta), path("*R2_umi.fastq.gz"), emit: umiextracted_reads
    path "*umiextract_report.txt", emit: umitools_extract_reports
    path "v_umi_tools.txt", emit: version

    script:
    new_meta = meta.clone()
    new_meta.single_end = true
    """
    [ ! -f  ${meta.name}_R1.fastq.gz ] && ln -s ${reads[0]} ${meta.name}_R1.fastq.gz
    [ ! -f  ${meta.name}_R2.fastq.gz ] && ln -s ${reads[1]} ${meta.name}_R2.fastq.gz
    umi_tools extract -I ${meta.name}_R1.fastq.gz --bc-pattern=NNNNNNNN --read2-in=${meta.name}_R2.fastq.gz --stdout=${meta.name}_R1_umi.fastq.gz --read2-out=${meta.name}_R2_umi.fastq.gz > ${meta.name}_umiextract_report.txt
    umi_tools --version &> v_umi_tools.txt
    """
}
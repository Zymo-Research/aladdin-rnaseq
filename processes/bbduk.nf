// Run BBDuk to remove low-complexity reads
params.publish_dir = "BBDuk"

process bbduk {
    tag "${meta.name}"
    label 'process_medium'
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename ->
            if (filename.endsWith('.log')) "logs/$filename"
            else if (params.save_trimmed && filename.endsWith(".fastq.gz")) "reads/$filename"
            else null
        }

    container 'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    path('*.log')                      , emit: report
    path('v_bbduk.txt')                , emit: version

    script:
    def raw     = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed = meta.single_end ? "out=${meta.name}_R1.fastq.gz" : "out1=${meta.name}_R1.fastq.gz out2=${meta.name}_R2.fastq.gz"
    def entropy = params.bbduk_entropy ? "entropy=${params.bbduk_entropy} entropywindow=${params.bbduk_entropy_window}" : ''
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    bbduk.sh \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed \\
        $entropy \\
        threads=$task.cpus \\
        &> ${meta.name}.bbduk.log
    bbduk.sh --version &> v_bbduk.txt
    """
}

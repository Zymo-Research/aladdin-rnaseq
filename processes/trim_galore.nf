// Trimming
params.publish_dir = "Trim_Galore"
params.save_trimmed = false
params.protocol_settings = [ "clip_r1":0, "clip_r2":0, "three_prime_clip_r1":0, "three_prime_clip_r2":0,  
                             "adapter":false, "adapter2":false, "strandedness":0 ]
params.trim_nextseq = 0
params.min_read_length = 20
params.adapter_overlap = 1
params.skip_trimming = false

process trim_galore {
    label 'mid_memory'
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (params.save_trimmed && filename.endsWith("fq.gz")) filename
            else null
        }

    when:
    !params.skip_trimming

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*fq.gz"), emit: reads
    path "*trimming_report.txt", emit: report
    path "*_fastqc.{zip,html}"
    path "v_trim_galore.txt", emit: version

    script:
    c_r1 = params.protocol_settings.clip_r1 > 0 ? "--clip_r1 ${params.protocol_settings.clip_r1}" : ''
    c_r2 = params.protocol_settings.clip_r2 > 0 ? "--clip_r2 ${params.protocol_settings.clip_r2}" : ''
    tpc_r1 = params.protocol_settings.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.protocol_settings.three_prime_clip_r1}" : ''
    tpc_r2 = params.protocol_settings.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.protocol_settings.three_prime_clip_r2}" : ''
    nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
    a1 = params.protocol_settings.adapter ? "-a ${params.protocol_settings.adapter}" : ''
    a2 = params.protocol_settings.adapter2 ? "-a2 ${params.protocol_settings.adapter2}" : ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    if (meta.single_end) {
        """
        trim_galore --version &> v_trim_galore.txt
        [ ! -f  ${meta.name}_R1.fastq.gz ] && ln -s $reads ${meta.name}_R1.fastq.gz
        trim_galore -j $task.cpus --fastqc --gzip $c_r1 $tpc_r1 --stringency ${params.adapter_overlap} --length ${params.min_read_length} $a1 $nextseq ${meta.name}_R1.fastq.gz
        """
    } else {
        """
        trim_galore --version &> v_trim_galore.txt
        [ ! -f  ${meta.name}_R1.fastq.gz ] && ln -s ${reads[0]} ${meta.name}_R1.fastq.gz
        [ ! -f  ${meta.name}_R2.fastq.gz ] && ln -s ${reads[1]} ${meta.name}_R2.fastq.gz
        trim_galore -j $task.cpus --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 --stringency ${params.adapter_overlap} --length ${params.min_read_length} $a1 $a2 $nextseq ${meta.name}_R1.fastq.gz ${meta.name}_R2.fastq.gz
        """
    }
}
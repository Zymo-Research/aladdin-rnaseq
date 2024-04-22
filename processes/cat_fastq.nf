// Merge FASTQ files from multiple runs for the same sample

process cat_fastq {
    tag "${meta.name}"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")  // avoid filename collision

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads

    script:
    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
    if (meta.single_end) {
        if (readList.size > 1) {
            """
            cat ${readList.join(' ')} > ${meta.name}.merged.fastq.gz
            """
        }
    } else {
        if (readList.size > 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
            """
            cat ${read1.join(' ')} > ${meta.name}_1.merged.fastq.gz
            cat ${read2.join(' ')} > ${meta.name}_2.merged.fastq.gz
            """
        }
    }
}
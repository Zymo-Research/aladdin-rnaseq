// Run Salmon quantification on transcriptome BAM
params.publish_dir = "Salmon"
params.read_quant_method = "STAR_featureCounts"
params.strandedness = 0

process salmon_quant {
    label 'low_memory'
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { (it=='v_Salmon.txt'||it.endsWith('mqc.json')) ? null : it }
    
    when:
    params.read_quant_method == 'STAR_Salmon'

    input:
    tuple val(meta), path(logfile), path(bam)
    path gtf
    path transcripts

    output:
    path "${meta.name}", emit: results
    path "*mqc.json", emit: report
    path "v_Salmon.txt", emit: version

    script:
    libtypes = [0:'U', 1:'SF', 2:'SR']
    libtype = libtypes[params.strandedness]
    if (!meta.single_end) { libtype = 'I' + libtype }
    log_type = logfile.getName().endsWith('dedup_LOGFILE') ? '--dedup_log' : '--star_log'
    """
    salmon quant --geneMap $gtf --threads ${task.cpus} --libType=$libtype -t $transcripts -a $bam -o ${meta.name}
    salmon_quant_mqc.py $log_type $logfile --salmon_log ${meta.name}/logs/salmon_quant.log 
    salmon --version &> v_Salmon.txt
    """
}
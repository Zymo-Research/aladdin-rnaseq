// Run featureCounts to assign and count reads
params.publish_dir = "featureCounts"
params.strandedness = 0
params.fc_extra_attributes = 'gene_name'
params.fc_group_features = 'gene_id'
params.fc_count_type = 'exon'
params.fc_biotype_count_type = 'exon'
params.fc_biotype_group_features = 'gene_biotype'
params.filter_primary = false
params.read_quant_method = 'STAR_featureCounts'

process featurecounts {
    label 'low_memory'
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename ->
            if (filename.endsWith("_gene.featureCounts.txt.summary")) "gene_counts/$filename"
            else if (filename.endsWith("_gene.featureCounts.txt")) "gene_counts/$filename"
            else null
        }
    
    when:
    params.read_quant_method == 'STAR_featureCounts'

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    
    output:
    path "${meta.name}_gene.featureCounts.txt", emit: counts
    path "${meta.name}_gene.featureCounts.txt.summary", emit: report
    path "v_featurecounts.txt", emit: version

    script:
    extraAttributes = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
    paired_end = meta.single_end ? '' : "-p"
    """
    featureCounts -a $gtf -g ${params.fc_group_features} -t ${params.fc_count_type} -o ${meta.name}_gene.featureCounts.txt $extraAttributes ${paired_end} -s ${params.strandedness} $bam
    featureCounts -v &> v_featurecounts.txt
    """
}

process featurecounts_biotype {
    label 'low_memory'
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_biotype.featureCounts.txt") > 0) "biotype_counts/$filename"
            else null
        }
    
    when:
    !params.skip_biotype_qc

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path rRNA_gtf

    output:
    path "${meta.name}_biotype.featureCounts.txt", emit: biotype_counts
    path "${meta.name}_biotype.featureCounts.txt.summary", emit: biotype_log
    path "v_featurecounts.txt", emit: version

    script:
    append_gtf = rRNA_gtf ? "cat $gtf $rRNA_gtf > biotype.gtf" : "cp $gtf biotype.gtf"
    paired_end = meta.single_end ? '' : "-p"
    """
    $append_gtf
    featureCounts -a biotype.gtf -g ${params.fc_biotype_group_features} -t ${params.fc_biotype_count_type} -o ${meta.name}_biotype.featureCounts.txt ${paired_end} -s ${params.strandedness} -M $bam
    cut -f 1,7 ${meta.name}_biotype.featureCounts.txt | tail -n +3 > temp.txt && mv temp.txt ${meta.name}_biotype.featureCounts.txt
    featureCounts -v &> v_featurecounts.txt
    """
}
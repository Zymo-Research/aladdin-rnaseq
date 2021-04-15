// Run featureCounts to assign and count reads
params.publish_dir = "featureCounts"
params.strandedness = 0
params.fc_extra_attributes = 'gene_name'
params.fc_group_features = 'gene_id'
params.fc_count_type = 'exon'
params.fc_biotype_count_type = 'exon'
params.fc_biotype_group_features = 'gene_biotype'
params.filter_primary = false

process featurecounts {
    label 'low_memory'
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_biotype.featureCounts.txt") > 0) "biotype_counts/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else if (filename=='v_featurecounts.txt') null
            else filename
        }

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path rRNA_gtf

    output:
    path "${meta.name}_gene.featureCounts.txt", emit: counts
    path "${meta.name}_gene.featureCounts.txt.summary", emit: report
    path "${meta.name}_biotype.featureCounts.txt", emit: biotype_counts optional true
    path "${meta.name}_biotype.featureCounts.txt.summary", emit: biotype_log optional true
    path "v_featurecounts.txt", emit: version

    script:
    extraAttributes = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
    append_gtf = params.skip_biotype_qc ? '' : rRNA_gtf ? "cat $rRNA_gtf >> $gtf" : ''
    paired_end = meta.single_end ? '' : "-p"
    biotype_qc = params.skip_biotype_qc ? '' : "featureCounts -a $gtf -g ${params.fc_biotype_group_features} -t ${params.fc_biotype_count_type} -o ${meta.name}_biotype.featureCounts.txt ${paired_end} -s ${params.strandedness} -M $bam"
    mod_biotype = params.skip_biotype_qc ? '' : "cut -f 1,7 ${meta.name}_biotype.featureCounts.txt | tail -n +3 > temp.txt && mv temp.txt ${meta.name}_biotype.featureCounts.txt"
    filter_primary = params.filter_primary ? "samtools view -F 256 -b $bam > primary.bam && mv primary.bam $bam" : ''
    """
    $filter_primary
    featureCounts -a $gtf -g ${params.fc_group_features} -t ${params.fc_count_type} -o ${meta.name}_gene.featureCounts.txt $extraAttributes ${paired_end} -s ${params.strandedness} $bam
    $append_gtf
    $biotype_qc
    $mod_biotype
    featureCounts -v &> v_featurecounts.txt
    """
}
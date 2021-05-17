// Merge featureCounts results into one
params.publish_dir = "featureCounts"
params.bam_suffix = "Aligned.sortedByCoord.out.bam"

process merge_featurecounts {
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
    path input_files

    output:
    path 'merged_gene_counts.txt', emit: merged_counts
    path 'gene_lengths.txt', emit: gene_lengths

    script:
    // Geneid in 1st column; gene length in 6th; gene_name in 7th
    gene_ids = "<(tail -n +2 ${input_files[0]} | cut -f1,7 )"
    counts = input_files.collect{filename ->
    // Remove first line and take 8th column (counts)
    "<(tail -n +2 ${filename} | sed 's:${params.bam_suffix}::' | cut -f8)"}.join(" ")
    """
    paste $gene_ids $counts > merged_gene_counts.txt
    tail -n +2 ${input_files[0]} | cut -f1,6 > gene_lengths.txt
    """
}
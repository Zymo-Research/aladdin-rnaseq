// Add ERCC concentration to the counts matrix
params.publish_dir = "./ERCC"
params.ercc_spikein = false

process add_ercc_info {
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
    path counts
    path info

    output:
    path "ERCC_counts_info.tsv", emit: report
    
    script:
    """
    add_ercc_info.py -c $counts -i $info -m $params.ercc_spikein
    """
}
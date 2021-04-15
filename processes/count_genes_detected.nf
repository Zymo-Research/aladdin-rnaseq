// Calculate FPKM, TPM; count numbers of genes detected and format for MultiQC
params.publish_dir = "featureCounts"
params.gene_detection_method = "reads"

process count_genes_detected {

    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it.endsWith('.json') ? null : it }

    input:
    path merged_counts
    path gene_lengths
    
    output:
    path 'FPKM.tsv', emit: fpkm
    path 'TPM.tsv', emit: tpm
    path '*mqc.json', emit: report

    script:
    def gene_detection_str = "$merged_counts -m reads -c 1"
    if (params.gene_detection_method == 'fpkm') {
        gene_detection_str = 'FPKM.tsv -m fpkm -c 0.1'
    } else if (params.gene_detection_method == 'tpm') {
        gene_detection_str = 'TPM.tsv -m tpm -c 0.1'
    }
    """
    calculate_fpkm_tpm.py $merged_counts -l $gene_lengths
    count_genes_detected.py $gene_detection_str
    """
}
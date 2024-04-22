// Combine featureCounts or Salmon results into one file
// Calculate FPKM, TPM values where applicable
// Count No. genes detected
params.publish_dir = "summarize_samples"
params.read_quant_method = "STAR_featureCounts"
params.gene_detection_method = "reads"
params.fc_extra_attributes = "gene_name"

process summarize_samples {
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it.endsWith('.json') ? null : it }
    
    input:
    path input_files

    output:
    path "merged_gene_counts.txt", emit: counts
    path "gene_TPM.txt", emit: tpm
    path "*mqc.json", emit: report
    path "merged_transcript_counts.txt", emit: transcript_counts optional true
    path "transcript_TPM.txt" optional true
    path "*FPKM.txt"

    script:
    def gene_detection_str = "-m reads -c 1"
    if (params.gene_detection_method == 'fpkm') {
        gene_detection_str = "-m fpkm -c 0.1"
    } else if (params.gene_detection_method == 'tpm') {
        gene_detection_str = "-m tpm -c 0.1"
    }
    extra_attr = params.fc_extra_attributes ? "-e ${params.fc_extra_attributes}" : ""
    """
    summarize_samples.py $input_files -q ${params.read_quant_method} $extra_attr $gene_detection_str 
    """
}
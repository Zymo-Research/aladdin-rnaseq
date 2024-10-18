// Combine DEG and DTU results
process combine_deg_dtu {
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it.endsWith('.xlsx') ? it : null }

    input:
    path deg_results
    path dtu_results

    output:
    path "*_combined_DEG_DTU_padj.tsv" , emit: report
    path "*_combined_DEG_DTU_padj.xlsx", emit: download
    
    script:
    """
    combine_deg_dtu.py -g ${deg_results} -t ${dtu_results}
    """
}
// Merge and parse biotype QC results for MultiQC
params.skip_biotype_qc = false

process parse_biotype_qc {
    when:
    !skip_biotype_qc
    
    input:
    path counts
    path logs

    output:
    path 'biotype*mqc.json', emit: report

    script:
    """
    parse_biotype_for_multiqc.py $counts -f rRNA
    """
}
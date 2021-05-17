// Run pathway enrichment analysis with gProfiler
params.publish_dir = "gProfiler"
params.gprofiler_organism = false
params.deseq2_fdr = 0.05
params.gprofiler_fdr = 0.05

process gprofiler {
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it.endsWith('xlsx') ? it : null }

    input:
    path deseq_results

    when:
    params.gprofiler_organism

    output:
    path "*_gProfiler_results.tsv", emit: report
    path "*_gProfiler_results.xlsx", emit: download
    path "v_gProfiler.txt", emit: version

    script:
    """
    gProfiler.py $deseq_results -o $params.gprofiler_organism -q $params.deseq2_fdr -p $params.gprofiler_fdr
    pip freeze | grep gprofiler > v_gProfiler.txt
    """
}
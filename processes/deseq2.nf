// Run study level analysis/plots and DEG analysis with DESeq2
params.publish_dir = "DESeq2"
params.skip_deseq2 = false
params.deseq2_fdr = 0.05
params.deseq2_lfc = 0.585

process deseq2 {
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename -> 
            if (filename.endsWith('.xlsx') || filename.endsWith('.jpg')) filename
            else null
        }

    when:
    !params.skip_deseq2

    input:
    path merged_counts
    path design_csv
    path comparisons

    output:
    path "*.{xlsx,jpg}", emit: download
    path "*_DESeq_results.tsv", emit: results optional true
    path "*{heatmap,plot,matrix}.tsv", emit: report optional true
    path "v_DESeq2.txt", emit: version

    script:
    """
    DESeq2.r $merged_counts $design_csv $params.deseq2_fdr $params.deseq2_lfc $comparisons
    Rscript -e "write(x=as.character(packageVersion('DESeq2')), file='v_DESeq2.txt')"
    """
}
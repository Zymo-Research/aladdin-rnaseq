// Generate report via MultiQC
params.publish_dir = "MultiQC"
params.skip_multiqc = false
params.run_name = false
params.ensembl_web = false
params.deseq2_fdr = 0.05
params.gprofiler_fdr = 0.05
params.dexseq_fdr = 0.05
params.trimming_module = 'Trim_Galore'

process multiqc {
    label 'no_cache'
    publishDir "${params.publish_dir}", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    path multiqc_config
    path extra_multiqc_config
    path ('reads_qc_trimming/*')
    path ('align_dedup_quant/*')
    path ('ercc/*')
    path ('alignment_qc/*')
    path ('comparison/*')
    path ('software_versions/*')
    path summary_header
    path workflow_summary
    val warning

    output:
    path "*multiqc_report.html", emit: report
    path "*report_data"
    //path "multiqc_plots"

    script:
    rtitle = params.run_name ? "--title \"RNAseq report for $params.run_name\"" : ''
    rfilename = params.run_name ? "--filename " + params.run_name.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    ensembl_link = params.ensembl_web ? "gene_link_prefix: $params.ensembl_web" : ''
    comment = warning ? "--comment \"$warning\"" : ''
    extra_config = extra_multiqc_config ? "--config ${extra_multiqc_config}" : ''
    """
    cat $summary_header $workflow_summary > workflow_summary_mqc.yaml
    echo '    </dl>' >> workflow_summary_mqc.yaml
    rm $summary_header $workflow_summary
    echo 'DESeq2_alpha: $params.deseq2_fdr' >> $multiqc_config
    echo 'gProfiler_alpha: $params.gprofiler_fdr' >> $multiqc_config
    echo 'DEXSeq_alpha: $params.dexseq_fdr' >> $multiqc_config
    echo '$ensembl_link' >> $multiqc_config
    multiqc . -f $rtitle $rfilename --config $multiqc_config $extra_config $comment \\
        -m fastqc -m ${params.trimming_module} -m bbduk -m star -m rseqc -m preseq -m picard -m qualimap -m featureCounts \\
        -m custom_content -m DESeq2 -m gProfiler -m DTU -m plot_sample_distance -m plot_gene_heatmap -m plot_ERCC
    """
}

process multiqc_comparison_only {
    label 'no_cache'
    publishDir "${params.publish_dir}", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    path multiqc_config
    path ('comparison/*')
    path ('software_versions/*')
    path summary_header
    path workflow_summary

    output:
    path "*multiqc_report.html", emit: report
    path "*report_data"
    //path "multiqc_plots"

    script:
    rtitle = params.run_name ? "--title \"RNAseq report for $params.run_name\"" : ''
    rfilename = params.run_name ? "--filename " + params.run_name.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    ensembl_link = params.ensembl_web ? "gene_link_prefix: $params.ensembl_web" : ''
    """
    cat $summary_header $workflow_summary > workflow_summary_mqc.yaml
    echo '    </dl>' >> workflow_summary_mqc.yaml
    rm $summary_header $workflow_summary
    echo 'DESeq2_alpha: $params.deseq2_fdr' >> $multiqc_config
    echo 'gProfiler_alpha: $params.gprofiler_fdr' >> $multiqc_config
    echo 'DEXSeq_alpha: $params.dexseq_fdr' >> $multiqc_config
    echo '$ensembl_link' >> $multiqc_config
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m DESeq2 -m gProfiler -m DTU -m plot_sample_distance -m plot_gene_heatmap
    """
}
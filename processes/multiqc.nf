// Generate report via MultiQC
params.publish_dir = "MultiQC"
params.skip_multiqc = false
params.run_name = false
params.ensembl_web = false
params.ignore_R1 = false
params.deseq2_fdr = 0.05
params.gprofiler_fdr = 0.05

process multiqc {
    label 'no_cache'
    publishDir "${params.publish_dir}", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    path multiqc_config
    path ('fastqc/*')
    path ('trimgalore/*')
    path ('star/*')
    path ('ercc/*')
    path ('rseqc/*')
    path ('preseq/*')
    path ('markduplicates/*')
    path ('dupradar/*')
    path ('qualimap/*')
    path ('featurecounts/*')
    path ('biotype/*')
    path ('counted_genes/*')
    path ('deseq2/*')
    path ('gprofiler/*')
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
    ignore_r1_comment = params.ignore_R1 ? "--comment 'This protocol uses only Read 2 for analysis. All report sections present Read 2 as R1, Read 1 or read 1, except the Download data section. If you have any question regarding this report, please contact the Zymo representative who sent you this report.'" : ''
    ensembl_link = params.ensembl_web ? "gene_link_prefix: $params.ensembl_web" : ''
    """
    cat $summary_header $workflow_summary > workflow_summary_mqc.yaml
    echo '    </dl>' >> workflow_summary_mqc.yaml
    rm $summary_header $workflow_summary
    echo 'DESeq2_alpha: $params.deseq2_fdr' >> $multiqc_config
    echo 'gProfiler_alpha: $params.gprofiler_fdr' >> $multiqc_config
    echo '$ensembl_link' >> $multiqc_config
    multiqc . -f $rtitle $rfilename --config $multiqc_config $ignore_r1_comment \\
        -m fastqc -m Trim_Galore -m star -m rseqc -m preseq -m picard -m qualimap -m featureCounts \\
        -m custom_content -m DESeq2 -m gProfiler -m plot_sample_distance -m plot_gene_heatmap -m plot_ERCC
    """
}

process multiqc_3mrna {
    label 'no_cache'
    publishDir "${params.publish_dir}", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    path multiqc_config
    path ('fastqc/*')
    path ('trimgalore/*')
    path ('star/*')
    path ('umitools/*')
    path ('ercc/*')
    path ('rseqc/*')
    path ('preseq/*')
    path ('qualimap/*')
    path ('featurecounts/*')
    path ('biotype/*')
    path ('counted_genes/*')
    path ('deseq2/*')
    path ('gprofiler/*')
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
    echo '$ensembl_link' >> $multiqc_config
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m fastqc -m Trim_Galore -m star -m rseqc -m preseq -m qualimap -m featureCounts \\
        -m custom_content -m DESeq2 -m gProfiler -m plot_sample_distance -m plot_gene_heatmap -m plot_ERCC
    """
}

process multiqc_comparison_only {
    label 'no_cache'
    publishDir "${params.publish_dir}", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    path multiqc_config
    path ('deseq2/*')
    path ('gprofiler/*')
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
    echo '$ensembl_link' >> $multiqc_config
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m DESeq2 -m gProfiler -m plot_sample_distance -m plot_gene_heatmap -m download_data
    """
}

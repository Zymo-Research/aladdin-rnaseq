def collect_summary(params, workflow) {
    def summary = [:]
    if (workflow.revision) {
        summary['Pipeline Release'] = workflow.revision
    }
    // This has the bonus effect of catching both -name and --name
    def run_name = params.name
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        run_name = workflow.runName
    }
    summary['Run Name']                  = run_name ?: workflow.runName
    summary['Design']                    = params.design
    summary['Group Comparisons']         = params.comparisons
    summary['Genome']                    = params.genome
    summary['Protocol']                  = params.protocol
    summary['Save Unaligned']            = params.save_unaligned
    summary['Max Resources']             = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
    summary['Output dir']                = params.outdir
    summary['Launch dir']                = workflow.launchDir
    summary['Working dir']               = workflow.profile == 'awsbatch' ? "s3:/${workflow.workDir}" : workflow.workDir
    summary['Config Profile']            = workflow.profile
    if (workflow.profile == 'awsbatch') {
        summary['AWS Region']            = params.awsregion
        summary['AWS Queue']             = params.awsqueue
    }
    summary['DESeq2 FDR cutoff']         = params.deseq2_fdr
    summary['DESeq2 Log2FC cutoff']      = params.deseq2_lfc
    summary['gProfiler FDR cutoff']      = params.gprofiler_fdr
    if (params.dtu_analysis) {
        summary['DTU analysis FDR cutoff']       = params.dexseq_fdr
        summary['prop_filter_transcript_counts'] = params.prop_filter_transcript_counts
        summary['prop_filter_transcript_props']  = params.prop_filter_transcript_props
        summary['prop_filter_gene_counts']       = params.prop_filter_gene_counts
    }
    if (!params.merged_counts && !params.salmon_results) {
        summary['Min Adapter Overlap']   = params.adapter_overlap
        summary['Min Trimmed Length']    = params.min_read_length
        summary['Save Trimmed']          = params.save_trimmed
        summary['Percent Mapped Cutoff'] = params.percent_mapped_cutoff
        summary['Read Quant. Method']    = params.read_quant_method
        summary['Gene Detection Method'] = params.gene_detection_method
        if (params.ercc_spikein) {
            summary['ERCC spike-in']     = "ERCC92 Mix $params.ercc_spikein"
        }
    } else {
        summary['Salmon results']        = params.salmon_results
        summary['Merged read counts']    = params.merged_counts
    }
    return summary
}
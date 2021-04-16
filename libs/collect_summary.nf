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
    if (!params.merged_counts) {
        summary['Min Adapter Overlap']   = params.adapter_overlap
        summary['Min Trimmed Length']    = params.min_read_length
        summary['Save Trimmed']          = params.save_trimmed
        if (params.ercc_spikein) {
            summary['ERCC spike-in']     = "ERCC92 Mix $params.ercc_spikein"
        }
        presets = [ "zymo_ribofree" : [ "trimming": "5'R2: 10bp / adapter1: NNNNNNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC / adapter2: AGATCGGAAGAGCGTCGTGTAGGGAAAGA",
                                        "strandedness":'Reverse', "kit_name": 'Zymo-Seq RiboFree Total RNA Library Kit'], 
                    "zymo_3mrna"    : ["trimming" : "5'R1: 10bp / adapter:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "strandedness": 'Forward',
                                        "kit_name": "Zymo-Seq 3' mRNA Library Kit"],
                    "zymo_3mrna_nodedup" : [ "trimming": "5'R1: 10bp / adapter:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "strandedness": 'Forward',
                                        "kit_name": "Zymo-Seq 3' mRNA Library Kit"],
                    "illumina"      : [ "trimming": "only using Illumina adatpers", "strandedness":'None', "kit_name": "Illumina RNA kits"],
                    "pico"          : [ "trimming": "5'R1: 3bp / 3'R2: 3bp / Illumina adapters", "strandedness":'Forward', "kit_name": "SMARTer Stranded Total RNA-Seq Kit - Pico Input"]
        ]
        summary['Trimming']              = presets[params.protocol]['trimming']
        summary['Strandedness']          = presets[params.protocol]['strandedness']
        summary['Library Prep']          = presets[params.protocol]['kit_name']
    }
    return summary
}
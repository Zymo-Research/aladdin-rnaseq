#!/usr/bin/env nextflow
/*
Zymo rnaseq pipeline, adapted from nf-core/rnaseq pipeline (https://github.com/nf-core/rnaseq) release 1.4.2
Many changes were made including migration to DSL 2
*/

// DSL 2
nextflow.enable.dsl=2

// Load functions
include { help_message } from ('./libs/help_message')
include { collect_summary } from ('./libs/collect_summary')

// Show help message
if (params.help){
    help_message()
    exit 0
}

// AWSBatch sanity checking
if( workflow.profile == 'awsbatch') {
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// --outdir Remove potential trailing slash at the end
outdir = params.outdir - ~/\/*$/

// Genome options
// Check if genome exists in the config file
if (params.genome) {
    if (!params.genomes.containsKey(params.genome)) {
        exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
    } else {
        params.genome_settings = params.genomes[params.genome]
        keys = ['star','gtf','bed12','gprofiler','ensembl_web','rRNA_gtf','bacteria','csi_index']
        for (key in keys) {
            if (!(key in params.genome_settings.keySet())) {
                params.genome_settings[key] = false
            }
        }
    }
} else {
    exit 1, "--genome is a required input!"
}

// ERCC settings
if (params.ercc_spikein) {
    params.ercc_settings = params.genomes['ERCC92']
}

// Collect summary
log.info "Zymo rnaseq v${workflow.manifest.version}"
def summary = collect_summary(params, workflow)
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")

// Import appropriate workflow 
if (params.merged_counts) { // when merged_counts is supplied, it is assumed that only comparison is requested using merged_counts
    include { COMPARISON } from './workflows/comparison' addParams( summary: summary, outdir: outdir )
} else if (params.protocol == "zymo_3mrna") {
    include { RNASEQ3M } from './workflows/3mrnaseq' addParams( summary: summary, outdir: outdir )
} else {
    include { RNASEQ } from './workflows/rnaseq' addParams( summary: summary, outdir: outdir )
}

// Workflow 
workflow {
    if (params.merged_counts) { // when merged_counts is supplied, it is assumed that only comparison is requested using merged_counts
        COMPARISON()
    } else if (params.protocol == "zymo_3mrna") {
        RNASEQ3M()
    } else {
        RNASEQ()
    }
}

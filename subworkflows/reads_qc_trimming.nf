// Conduct Reads QC and trimming

include { fastqc } from '../processes/fastqc' addParams(
    publish_dir: "${params.outdir}/FastQC",
    skip_qc: params.skip_qc,
    skip_fastqc: params.skip_fastqc
)
include { umi_extract } from '../processes/umitools_extract' addParams(
    publish_dir: "${params.outdir}/umiextract",
)
include { trim_galore } from '../processes/trim_galore' addParams(
    publish_dir: "${params.outdir}/Trim_Galore",
    save_trimmed: params.save_trimmed,
    protocol_settings: params.protocol_settings,
    trim_nextseq: params.trim_nextseq,
    min_read_length: params.min_read_length,
    adapter_overlap: params.adapter_overlap
)
include { bbduk } from '../processes/bbduk' addParams(
    publish_dir: "${params.outdir}/BBDuk",
    save_trimmed: params.save_trimmed,
    bbduk_entropy: params.bbduk_entropy,
    bbduk_entropy_window: params.bbduk_entropy_window
)


workflow READS_QC_TRIMMING {
    take:
    reads

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Read QC
    fastqc(reads)
    ch_versions = ch_versions.mix(fastqc.out.version.first())
    ch_multiqc_files = ch_multiqc_files.mix(fastqc.out.report)
    
    if (!params.skip_trimming) {
        // Trim adapters
        if (params.protocol_settings.umi) {
            // Extract UMI from Read 1 and add to Read 2. The value for meta.single_end is also updated to be true
            umi_extract(reads)
            trim_galore(umi_extract.out.umiextracted_reads)
        } else {
            trim_galore(reads)
        }
        ch_versions = ch_versions.mix(trim_galore.out.version.first())
        ch_multiqc_files = ch_multiqc_files.mix(trim_galore.out.report)
        ch_reads = trim_galore.out.reads

        // Remove low complexity reads
        if (!params.skip_bbduk) {
            bbduk(trim_galore.out.reads)
            ch_versions = ch_versions.mix(bbduk.out.version.first())
            ch_multiqc_files = ch_multiqc_files.mix(bbduk.out.report)
            ch_reads = bbduk.out.reads
        }
    } else {
        ch_reads = reads
    }

    emit:
    versions = ch_versions
    report   = ch_multiqc_files
    reads    = ch_reads
}
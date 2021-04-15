def help_message() {
    log.info "Zymo rnaseq v${workflow.manifest.version}"
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run Zymo-Research/nxf-rnaseq --genome GRCh37 --design "<design CSV file on S3>" --protocol zymo_ribofree \
                                          -profile awsbatch -work-dir "<work dir on S3>" --awsqueue "<SQS ARN>" --outdir "<output dir on S3>"

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
    
    Mandatory arguments (full pipeline):
      --design                      Path to a CSV file with sample labels, input data locations, and sample grouping
      --genome                      Name of iGenomes reference

    Mandatory arguments (comparison only):
      --design                      Path to a CSV file with sample labels and sample grouping
      --merged_counts               Path to a TSV file with read counts of all samples
      --genome                      Name of iGenomes reference (Only required with pathway analysis)

    Mandatory arguments (generate URL only):
      --file_locations              Path to a file containing locations of files to generate download URL for

    Trimming options:
      --protocol                    Library prep protocol(kit) used. This will determine trimming parameters and strandedness.
                                    Available: illumina, zymo_ribofree, zymo_3mrna, zymo_3mrna_nodedup, pico. Default: illumina
      --adapter_overlap [int]       Instructs Trim Galore to require at least this many bp of overlap with adapters to trim a read. Default is 1
      --min_read_length [int]       Instructs Trim Galore to discard reads shorter than this. Default is 20
      --save_trimmed                Save trimmed FastQ file intermediates
      --skip_trimming               Skips trimming step and use the reads directly for alignment

    Alignment:
      --save_unaligned              Save unaligned reads from either STAR or HISAT2 to extra output files.
      --save_secondary_alignments   Save secondary alignments. By default only primary alignment were saved, resulting in only one alignment per read in BAM files.
      --star_min_overlap [int]      Minimum overlap bases for merging overlapping paired-end reads in STAR (default: 10)
      --star_max_overlap_mismatch   Maximum mismatch ratio in the overlap region for merging overlapping paired-end reads in STAR (default: 0.01)
      --star_twopass                Use two-pass mode in STAR alignment. This is useful for discovering novel junctions.
      --percent_mapped_cutoff       The least percentage of reads mapped required for a sample to proceed past alignment (default: 5)     

    Read Counting:
      --fc_extra_attributes         Define which extra parameters should also be included in featureCounts (default: 'gene_name')
      --fc_group_features           Define the attribute used to group features (default: 'gene_id')
      --fc_count_type               Define the type used to assign reads (default: 'exon')
      --fc_biotype_group_features   Define the attribute used to group features in the biotype analysis (default: 'gene_biotype')
      --fc_biotype_count_type       Define the type used to assign reads in the biotype analysis (default: 'exon')
      --skip_biotype_qc             Skip Biotype QC
      --gene_detection_method [str] Define which method is used to calculate numbers of genes detected.
                                    Available: reads, fpkm, tpm. Default: reads

    QC options:
      --skip_qc                     Skip all QC steps apart from MultiQC
      --skip_fastqc                 Skip FastQC
      --skip_rseqc                  Skip RSeQC
      --skip_preseq                 Skip Preseq
      --skip_markduplicates         Skip Picard MarkDuplicates (and also dupRadar)
      --skip_dupradar               Skip dupRadar
      --skip_qualimap               Skip Qualimap
      --skip_multiqc                Skip MultiQC

    Comparison options:
      --deseq2_fdr [float]          FDR cutoff for DESeq2, default is 0.05
      --deseq2_lfc [float]          Log2FoldChange cutoff for DESeq2, please use log2 of desired fold change cutoff, default is 0.585, which is log2(1.5)
      --comparisons                 Path to a CSV file stating the sample groups you want to compare. If not provided, all pairwise comparisons will be carried out.
      --gprofiler_fdr [float]       FDR cutoff for gProfiler, default is 0.05

    Other options:
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --ercc_spikein [int]          Whether ERCC spike-in was used. Use 1 for Mix 1, 2 for Mix 2 (default: false)
      --kit_QC                      Use this parameter to supply a path to a CSV file specifying dataset pairs for scatterplotting if running the pipeline on QC libraries of Zymo-Seq RiboFree Total RNA Libary Kit. Must be used with --protocol zymo_ribofree. Default is false.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on

    Cloudfront options:
      --cloudfront_origin_path      The origin path for Cloudfront
      --cloudfront_domain_name      The domain name for Cloudfront
      --cloudfront_private_key_ID   The private key ID for Cloudfront
      --cloudfront_privateKey       The location of the private key file for Cloudfront, if not provided, download links will not be included in the report
      --cloudfront_link_duration    The number of days Cloudfront URLs will stay active. Default: 60
      --deliver_fastqs              Whether to include links to download original FASTQ files. Default: false
    """.stripIndent()
}
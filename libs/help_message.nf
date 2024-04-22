def help_message() {
    log.info "Aladdin rnaseq v${workflow.manifest.version}"
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run Zymo-Research/aladdin-rnaseq --genome GRCh37 --design "<design CSV file on S3>" --protocol Zymo_RiboFree_PROTv1 \
                                              -profile awsbatch -work-dir "<work dir on S3>" --awsqueue "<SQS ARN>" --outdir "<output dir on S3>"

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: docker, awsbatch.
    
    Mandatory arguments (full pipeline):
      --design                      Path to a CSV file with sample labels, input data locations, and sample grouping
      --genome                      Name of iGenomes reference

    Mandatory arguments (comparison only):
      --design                      Path to a CSV file with sample labels and sample grouping
      --genome                      Name of iGenomes reference (Only required with pathway analysis)
      --merged_counts               Path to a TSV file with read counts of all samples
      OR
      --salmon_results              Path to a folder with Salmon quant results, usually one subfolder per sample

    Trimming options:
      --protocol                    Library prep protocol(kit) used. This will determine trimming parameters and strandedness
                                    Available: Zymo_RiboFree_PROTv2, Zymo_RiboFree_PROTv1, Zymo_SwitchFree_PROTv1
      --trim_nextseq                Enable trimming for erroneous Gs at 3' end that are common on NextSeq or NovaSeq platforms
      --adapter_overlap             Instructs Trim Galore to require at least this many bp of overlap with adapters to trim a read (default: 1)
      --min_read_length             Instructs Trim Galore to discard reads shorter than this (default: 20)
      --save_trimmed                Save trimmed FastQ file intermediates
      --skip_trimming               Skips trimming step and use the reads directly for alignment
      --skip_bbduk                  Skips BBDuk step which discards low complexity reads
      --bbduk_entropy               The entropy threshold for BBDuk to designate a read as low complexity (default: 0.5)
      --bbduk_entropy_window        The window size in bp in which BBDuk calculates entropy (default: 50)

    Alignment:
      --save_unaligned              Save unaligned reads from STAR to extra output files.
      --save_secondary_alignments   Save secondary alignments. By default only primary alignment were saved, resulting in only one alignment per read in BAM files.
      --star_min_overlap            Minimum overlap bases for merging overlapping paired-end reads in STAR (default: 10)
      --star_max_overlap_mismatch   Maximum mismatch ratio in the overlap region for merging overlapping paired-end reads in STAR (default: 0.01)
      --star_twopass                Use two-pass mode in STAR alignment. This is useful for discovering novel junctions.
      --percent_mapped_cutoff       The least percentage of reads mapped required for a sample to proceed past alignment (default: 5)     

    Read Counting:
      --read_quant_method           The method to use for read quantification 
                                    Available: STAR_featureCounts, STAR_Salmon (default: 'STAR_featureCounts')
      --fc_extra_attributes         Define which extra parameters should also be included in featureCounts (default: 'gene_name')
      --fc_group_features           Define the attribute used to group features (default: 'gene_id')
      --fc_count_type               Define the type used to assign reads (default: 'exon')
      --fc_biotype_group_features   Define the attribute used to group features in the biotype analysis (default: 'gene_biotype')
      --fc_biotype_count_type       Define the type used to assign reads in the biotype analysis (default: 'exon')
      --skip_biotype_qc             Skip Biotype QC
      --gene_detection_method       Define which method is used to calculate numbers of genes detected.
                                    Available: reads, fpkm, tpm (default: 'reads')

    QC options:
      --skip_qc                     Skip all QC steps apart from MultiQC
      --skip_fastqc                 Skip FastQC
      --skip_rseqc                  Skip RSeQC
      --skip_preseq                 Skip Preseq
      --skip_markduplicates         Skip Picard MarkDuplicates (and also dupRadar)
      --skip_dupradar               Skip dupRadar
      --skip_qualimap               Skip Qualimap
      --skip_multiqc                Skip MultiQC
      --generate_bigwig             Generate genome coverage bigWig files

    Comparison options:
      --deseq2_fdr                  FDR cutoff for DESeq2 (default: 0.05)
      --deseq2_lfc                  Log2FoldChange cutoff for DESeq2, please use log2 of desired fold change cutoff (default: 0)
      --comparisons                 Path to a CSV file stating the sample groups you want to compare. If not provided, all pairwise comparisons will be carried out.
      --gprofiler_fdr               FDR cutoff for gProfiler (default: 0.05)
      --dtu_analysis                Run Differential Transcript Usage analysis
      --dexseq_fdr                  FDR cutoff for DEXSeq (default: 0.05)
      --prop_filter_transcript_counts
                                    Filtering criteria for DTU analysis. Minimum proportion of samples that a transcript must have >=10 reads (default: 0.5)
      --prop_filter_transcript_props
                                    Filtering criteria for DTU analysis. Minimum proportion of samples that a transcript must have >=10% of the reads of its gene (default: 0.5)
      --prop_filter_gene_counts     Filtering criteria for DTU analysis. Minimum proportion of samples that the gene of a transcript must have >=10 reads (default: 1)
      --heatmap_group_order         Orders the Top Gene Expression Patterns Heatmap by group label instead of similarity between samples. (default: false)

    Other options:
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --ercc_spikein [int]          Whether ERCC spike-in was used. Use 1 for Mix 1, 2 for Mix 2 (default: false)

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on

    """.stripIndent()
}
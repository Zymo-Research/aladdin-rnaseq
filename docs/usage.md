# Zymo-Research/aladdin-rnaseq: Usage

## Table of contents

<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:true -->
* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Running full RNAseq pipeline](#running-full-rnaseq-pipeline)
  * [Running comparison only](#running-comparison-only)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
  * [`-profile`](#-profile)
  * [`--design`](#--design)
* [Reference genomes](#reference-genomes)
  * [`--genome`](#--genome)
  * [`--ercc_spikein`](#--ercc_spikein)
* [Adapter Trimming](#adapter-trimming)
  * [`--protocol`](#--protocol)
  * [`--trim_nextseq`](#--trim_nextseq)
  * [`--adapter_overlap [int]`](#--adapter_overlap-int)
  * [`--min_read_length [int]`](#--min_read_length-int)
  * [`--save_trimmed`](#--save_trimmed)
  * [`--skip_trimming`](#--skip_trimming)
  * [`--skip_bbduk`](#--skip_bbduk)
  * [`--bbduk_entropy [float]`](#--bbduk_entropy-float)
  * [`--bbduk_entropy_window [int]`](#--bbduk_entropy_window-int)
* [Alignment](#--alignment)
  * [`--save_unaligned`](#--save_unaligned)
  * [`--save_secondary_alignments`](#--save_secondary_alignments)
  * [`--star_min_overlap [int]`](#--star_min_overlap-int)
  * [`--star_max_overlap_mismatch [float]`](#--star_max_overlap_mismatch-float)
  * [`--star_twopass`](#--star_twopass)
  * [`--percent_mapped_cutoff [float]`](#--percent_mapped_cutoff-float)
* [Read Counting](#read-counting)
  * [`--read_quant_method [str]`](--read_quant_method-str)
  * [`--fc_extra_attributes [str]`](#--fc_extra_attributes-str)
  * [`--fc_group_features [str]`](#--fc_group_features-str)
  * [`--fc_count_type [str]`](#--fc_count_type-str)
  * [`--fc_biotype_group_features [str]`](#--fc_biotype_group_features-str)
  * [`--fc_biotype_count_type [str]`](#--fc_biotype_count_type-str)
  * [`--skip_biotype_qc`](#--skip_biotype_qc)
  * [`--gene_detection_method [str]`](#--gene_detection_method-str)
* [Alignment QC](#alignment-qc)
  * [Skipping QC steps](#skipping-qc-steps)
  * [`--generate_bigwig](#--generate_bigwig)
* [Comparisons](#comparisons)
  * [`--skip_deseq2`](#--skip_deseq2)
  * [`--deseq2_fdr [float]`](#--deseq2_fdr-float)
  * [`--deseq2_lfc [float]`](#--deseq2_lfc-float)
  * [`--comparisons`](#--comparisons)
  * [`--gprofiler_fdr [float]`](#--gprofiler_fdr-float)
  * [`--dtu_analysis`](#--dtu_analysis)
  * [`--dexseq_fdr [float]`](#--dexseq_fdr-float)
  * [`--prop_filter_transcript_counts [float]`](#--prop_filter_transcript_counts-float)
  * [`--prop_filter_transcript_props [float]`](#--prop_filter_transcript_props-float)
  * [`--prop_filter_gene_counts [float]`](#--prop_filter_gene_counts-float)
  * [`--heatmap_group_order`](#--heatmap_group_order)
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Resource allocation](#resource-allocation)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`-work-dir`](#-work-dir)
  * [`--name`](#--name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--multiqc_config`](#--multiqc_config)
<!-- TOC END -->

## Introduction
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

### Running full RNAseq pipeline
The typical command for running the full RNAseq pipeline is as follows:

```bash
nextflow run Zymo-Research/aladdin-rnaseq \
  --genome GRCh38 \
  --protocol zymo_ribofree \
  --design "<design CSV file on S3>" \
  -profile awsbatch \
	-work-dir "<work dir on S3>" \
	--awsqueue "<SQS ARN>" \
	--outdir "<output dir on S3>" \
	-r "1.0.0" \
	--name "<study name>"
```

This command will retrieve the pipeline code from GitHub. If you have downloaded the pipeline code, you can subsitute `Zymo-Research/aladdin-rnaseq` with `main.nf`.

### Running comparison only
Sometimes, one would need to rerun the comparison part of the pipeline without reruning sample processing, for example, with a different FDR cutoff or discarding some bad samples. These is an alternative entrypoint to the pipeline starting with combined read counts of genes, running only `DESeq2`, `gProfiler` and report. To run the pipeline via this entrypoint, the command is as follows:

```bash
nextflow run Zymo-Research/aladdin-rnaseq \
  --genome GRCh38 \
  --design "<design CSV file on S3>" \
  -profile awsbatch \
	-work-dir "<work dir on S3>" \
	--awsqueue "<SQS ARN>" \
	--outdir "<output dir on S3>" \
	-r "1.0.0" \
	--name "<study name>" \
  --merged_counts "<merged counts file on S3>"
```

If you have run the full pipeline on the same data before, you can use the merged counts file generated by the pipeline, usually in the `featureCounts` folder, as the input. If not, you must provide a TSV file similar to that. Note, if you don't need pathway enrichment analysis, then `--genome` can be omitted.

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull Zymo-Research/nxf-rnaseq
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [Zymo-Research/aladdin-rnaseq releases page](https://github.com/Zymo-Research/aladdin-rnaseq/releases) and find the latest version number - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Main arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. While there are several profiles listed below, Zymo strongly recommends using the `awsbatch` profile, as most of our tests are conducted using this profile. We do not guarantee other profiles would work.

* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/rnaseq`](http://hub.docker.com/r/nfcore/rnaseq/)

When using `awsbatch` profile, one must supply [other options related to AWS batch](#aws-batch-specific-parameters), and supply the locations of [work directory](#-work-dir) and [output directory](#--outdir) on AWS S3.

### `--design`
You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--design 'path/to/data/sample_sheet.csv'
```

#### Full analysis
The `group` identifier should be identical when you have multiple replicates from the same experimental group.

The `sample` identifier should be used to name your samples. That's how they will appear in the report.

The contents of these two columns must contain only alphanumerical characters or underscores, must start with a letter, and cannot start with "R1" or "R2". Sample and group labels must also exclude phrases that will be automatically removed by MultiQC. Sample and group label terms unavailable for use in this pipeline can be found [here](https://github.com/ewels/MultiQC/blob/b936a7a6d7050f3edc1ceefe8ae6ecd93865bf66/multiqc/utils/config_defaults.yaml#L150-L284), in the MultiQC source code.

The `read_1` and `read_2` identifiers should be used to specify the locations of FASTQ files for read 1 and 2.

The example below shows a simple experiment with two groups and two replicates each group. As you can see, mixing of single-end and paired-end data is allowed. The pipeline will infer that from your design file automatically. 

```
group,sample,read_1,read_2
Control,Sample1,s3://mybucket/this_is_s1_R1.fastq.gz,s3://mybucket/this_is_s1_R2.fastq.gz
Control,Sample2,s3://mybucket/this_is_s2_R1.fastq.gz,s3://mybucket/this_is_s2_R2.fastq.gz
Experiment,Sample3,s3://mybucket/that_is_s3_R1.fastq.gz,
Experiment,Sample4,s3://mybucket/that_be_s4_R1.fastq.gz,
```

#### Simple analysis
If your experiment does not have replicates, or you don't want comparison between groups, you can simply leave the group identifiers empty. In this case, comparison steps `DESeq2` and `gProfiler` will be skipped.

```
group,sample,read_1,read_2
,Sample1,s3://mybucket/this_is_s1_R1.fastq.gz,s3://mybucket/this_is_s1_R2.fastq.gz
,Sample2,s3://mybucket/this_is_s2_R1.fastq.gz,s3://mybucket/this_is_s2_R2.fastq.gz
,Sample3,s3://mybucket/that_is_s3_R1.fastq.gz,
,Sample4,s3://mybucket/that_be_s4_R1.fastq.gz,
```

| Column      | Description                                                                                                               |
|-------------|---------------------------------------------------------------------------------------------------------------------------|
| `group`     | Group/condition identifier for sample. This will be identical for replicate samples from the same experimental group.     |
| `sample`    | Sample label. How you want your samples to be identified in the report.                                                   |
| `read_1`    | Full path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".                 |
| `read_2`    | Full path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".                 |

## Reference genomes

The reference genomes information are stored in [genome config file](../conf/igenomes.config). A explanation of attributes of each genome follows:

| Setting        | Description                                                                    |
|----------------|--------------------------------------------------------------------------------|
| `star`         | Path to the STAR index folder                                                  |
| `gtf`          | Path to the annotation GTF file                                                |
| `bed12`        | Path to the annotation BED file                                                |
| `gprofiler`    | Organism ID used in gProfiler, if applicable                                   |
| `bacteria`     | Indicates this is a bacteria, which will disallow introns in STAR alignment    |
| `ensembl_web`  | The domain of Ensembl website the organism belongs to, if applicable           |
| `rRNA_gtf`     | Path to a GTF file containing additional rRNA regions                          |
| `csi_index`    | Indicates building BAM indices using csi, applicable to very long chromosomes  |

We supplemented [AWS iGenomes](https://ewels.github.io/AWS-iGenomes/) with our own genome indices. 

### `--genome`
To run the pipeline, you must specify which to use with the `--genome` flag. The genome you specify must be one of the options available in the genomes file, either the default one, or the one you specify.

### `--ercc_spikein`
Use this option to indicate that ERCC spike-in was added to the samples. Possible values are `1` and `2`, which indicates Mix 1 and Mix 2 of [ERCC spike-in](https://www.thermofisher.com/order/catalog/product/4456739#/4456739), respectively. When this option is activated, reads that are not aligned to the main genome are aligned again to ERCC92 reference genome and counted. The counts and concentrations of ERCC spike-in transcripts will be displayed in scatter plots in the report.

## Adapter Trimming

### `--protocol`
Use this to specify which library kit or protocol was used to make the RNAseq library. The possible options can be found in [protocol config file](../conf/protocols.config). We have chosen the appropriate trimming parameters such as whether to trim bps off 5' or 3' of R1 or R2 and adapter sequences for these protocols. This is a mandatory option for full RNAseq analysis, but not for comparison only analysis.

### `--trim_nextseq`
Instructs Trim Galore to apply the --nextseq=20 option, to trim based on quality after removing poly-G tails.

### `--adapter_overlap [int]`
Instructs Trim Galore to require at least this many bp of overlap with adapters to trim a read. Default is 1.

### `--min_read_length [int]`
Instructs Trim Galore to discard reads shorter than this. Default is 20.

### `--save_trimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this flag to copy these files when complete.

### `--skip_trimming`
You can use this option to skip the trimming step, for example, if your input sequences are already trimmed elsewhere. This also skips the BBDuk step mentioned below.

### `--skip_bbduk`
By default, the pipeline removes low complexity reads such as those containing mostly just polyAs using BBDuk. You can use this option to skip this step.

### `--bbduk_entropy [float]`
The entropy threshold used by BBDuk to remove low complexity reads. Default is 0.5. For details on how BBDuk works, see [here](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)

### `--bbduk_entropy_window [int]`
The window size in bp that BBDuk uses to calculate entropy. Default is 50bp.

## Alignment

Unlike the nfcore/rnaseq pipeline, this pipeline only supports STAR as the aligner. This decision was based on internal research at Zymo.

### `--save_unaligned`
By default reads that don't align to the reference genome are not saved. Set to true to output those reads in FASTQ format.

### `--save_secondary_alignments`
By default we ask STAR to only output one alignment per read (pair) using `--outSAMmultNmax 1` option in STAR. Secondary alignments are not counted in read counting step and therefore have no effect on the final results. However, they may have unintended effect on QC results, because not all tools filter secondary alignments. Therefore, we decided to just eliminate them at alignment. If you want those secondary alignment for some reason, you can use this option to include them in the BAM file.

### `--star_min_overlap [int]`
We use STAR to merge paired-end reads into one read if they overlap with each other, because reads that dovetail may cause a problem in STAR alignment (See [STAR Github issue](https://github.com/alexdobin/STAR/issues/484)). Those merged reads will be converted to two paired reads in the output BAM. By default, this merging will only be performed if two paired reads overlap each other by 10 bp. Use this option to change that cutoff.

### `--star_max_overlap_mismatch [float]`
This controls the mismatch rate allowed in the merging step mentioned above. The default is 0.01, meaning the mismatch rate between the two reads in the overlap region cannot exceed 0.01 (basically, perfect match required). Use this option to change this cutoff.

### `--star_twopass`
By default, we don't ask STAR to run two-pass alignment. Two-pass alignment is very useful for novel splice junction discovery, but this is usually not the focus of our services. Use this option to turn two-pass alignment on.

### `--percent_mapped_cutoff [float]`
By default, we filter samples by requiring at least 5% of reads uniquely mapped to genome by STAR. Those samples failing this filter will be discarded after the alignment step. Use this option to change this cutoff.

## Read Counting
### `--read_quant_method [str]`
We offer two ways of conducting read counting, 'STAR_featureCounts' and 'STAR_Salmon'. Both rely on alignments generated by STAR. STAR_featureCounts uses [featureCounts of the Subread package](http://subread.sourceforge.net/) while STAR_Salmon uses [Salmon](https://combine-lab.github.io/salmon/). It is worth noting that we are not using the alignment-free mode of Salmon. The main difference between the two methods are their quantification strategy. featureCounts relies on simple overlap between aligned reads and exons, disards multi-mapped reads. Salmon estimates the number of reads that originate from each transcript, taking into account potential biases. It uses a probabilistic model to assign reads to transcripts, even when reads may map to multiple transcripts. Default is 'STAR_featureCounts'.<br>

For the STAR_featureCounts method, there are several settings that are critical to get this counting right, and they depend on the attributes used in the annotation GTF file. If you are using a standard ENSEMBL genome, then you don't need to change any of the following options. If you are not, you need to carefully examine your GTF file and change the following option accordingly.

### `--fc_extra_attributes [str]`
By default, the pipeline uses `gene_name` as additional gene identifiers apart from ENSEMBL identifiers in the pipeline.
This behaviour can be modified by specifying `--fc_extra_attributes` when running the pipeline, which is passed on to featureCounts as an `--extraAttributes` parameter.
Note that you can also specify more than one desired value, separated by a comma:
``--fc_extra_attributes gene_id,...``

### `--fc_group_features [str]`
By default, the pipeline uses `gene_id` as the atrribute type to group features. This means reads from different exons of the same gene are summed. Use this option to change this, which is passed on to featureCounts as an `-g` parameter.

### `--fc_count_type [str]`
By default, the pipeline uses `exon` as the feature to count. This means only reads that overlap with exons will be counted. Use this option to change this to, for example, `gene`, which is passed to featureCounts as an `-t` parameter.

### `--fc_biotype_group_features [str]`
In addition to counting reads by genes, we also perform an additional counting step by biotypes of genes. This is a QC step that will give an estimate of the composition of the RNA library, which will be ploted in the report. By default, in this step, reads that have the same `gene_biotype` are grouped together and counted. Biotype information may be saved in some other attribute in some GTF files. Use this option to change this attribute. This is passed to featureCounts as an `-g` parameter.

### `--fc_biotype_count_type [str]`
This is the feature to count, equivalent of `--fc_count_type`, but for biotype analysis. Default is `exon`.

### `--skip_biotype_qc`
Sometimes, the biotype counting and QC step mentioned above may not be possible or useful, for example, some GTF files may have incomplete or no biotype attributes. Use this option to skip the biotype counting step.

### `--gene_detection_method [str]`
Use this option to choose which metrics to use for calculating the numbers of genes detected. Available: `reads`, `fpkm`, `tpm`. Default is `reads`.

## Alignment QC

### Skipping QC steps
The pipeline contains a large number of quality control steps. Sometimes, it may not be desirable to run all of them if time and compute resources are limited.
The following options make this easy:

* `--skip_qc` -                Skip **all QC steps**, apart from MultiQC
* `--skip_fastqc` -            Skip FastQC
* `--skip_rseqc` -             Skip RSeQC
* `--skip_preseq` -            Skip Preseq
* `--skip_markduplicates` -    Skip Picard MarkDuplicates (and dupRadar)
* `--skip_dupradar` -          Skip dupRadar
* `--skip_qualimap` -          Skip Qualimap
* `--skip_multiqc` -           Skip MultiQC

### `--generate_bigwig`
Use this option to generate genomic coverage bigWig files. If the protocol is stranded, the pipeline will generate two bigWig files for each sample, one for each strand. 

## Comparisons

### `--skip_deseq2`
Use this option to skip DESeq2, and gProfiler by extension. 

### `--deseq2_fdr [float]`
Use this option to specify the false discovery rate (FDR) cutoff used in DESeq2. Genes whose FDR is smaller than this will be considered significantly differentially expressed. Default is 0.05.

### `--deseq2_lfc [float]`
Use this option to specify the fold change cutoff in DESeq2. Please note that this is not a simple filter, this cutoff changes the null hypothesis in the statistical tests, therefore affects the p-values and FDRs in DESeq2. This is passed to DESeq2 as a `Log2FoldChange` parameter, therefore, please use Log2 values of intended fold change cutoff. Default is 0.

### `--comparisons`
By default, the pipeline will compare all pairwise combinations of sample groups. If this is not desirable, use this option to specify the path to a CSV file that describes which sample groups you want to comapre. The CSV file should have the following format:
```
group_1,group_2
Experiment1,Control
Experiment2,Control
```
The header must be the same as shown. All group labels should also appear in the design CSV file "group" column.

### `--gprofiler_fdr [float]`
Use this option to specify the FDR cutoff used in gProfiler. Default is 0.05.

### `--dtu_analysis`
Use this option to ask the pipeline to conduct Differential Transcript Usage(DTU) analysis, adapted from [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6178912/).

### `--dexseq_fdr [float]`
Use this option to specify the FDR cutoff used in DEXSeq in DTU analysis. Default is 0.05

### `--prop_filter_transcript_counts [float]`
This is a filtering criteria in DTU analysis. This requires a minimum proportion of samples that a transcript must have >=10 reads. Any transcripts below the threshold will be discarded before the statistical analysis. Default is 0.5.

### `--prop_filter_transcript_props [float]`
This is a filtering criteria in DTU analysis. This requires a minimum proportion of samples that a transcript must have >=10% of the reads of its gene. Any transcripts below the threshold will be discarded before the statistical analysis. Default is 0.5.

### `--prop_filter_gene_counts [float]`
This is a filtering criteria in DTU analysis. This requires a minimum proportion of samples that the gene of a transcript must have >=10 reads. Any transcripts below the threshold will be discarded before the statistical analysis. Default is 1, meaning a gene and its transcripts must have at least 10 reads in all samples to be included in DTU analysis.

### `--heatmap_group_order`
Use this option to order the Top Gene Expression Patterns Heatmap by group label instead of similarity between samples.

## Job resources

### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Resource allocation
The computing resources allocated to each type of jobs are listed in the [base.config](../conf/base.config) file. It has been tested on genomes that are similar or smaller than the human genome. If more computing resources are needed, say for a much larger custom genome, it is best to modify the `base.config` file. Remember that you would need a different JobQueue with larger instance as well.

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters. Please make sure to also set the `-w/-work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.

### `--awsregion`
The AWS region to run your job in. Default is set to `us-east-1` but can be adjusted to your needs.

## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `-work-dir`
The working directory where intermediate files will be saved.

**NB:** Single hyphen (core Nextflow option)

### `--name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML.

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file.

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--multiqc_config`
Specify a path to a custom MultiQC configuration file.

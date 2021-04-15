# Zymo-Research/nxf-rnaseq: Usage

## Table of contents

<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:true -->
* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Running full RNAseq pipeline](#running-full-rnaseq-pipeline)
  * [Running comparison only](#running-comparison-only)
  * [Running URL generation only](#running-url-generation-only)
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
  * [`--trim_nextseq [int]`](#--trim_nextseq-int)
  * [`--adapter_overlap [int]`](#--adapter_overlap-int)
  * [`--min_read_length [int]`](#--min_read_length-int)
  * [`--save_trimmed`](#--save_trimmed)
* [Alignment](#--alignment)
  * [`--save_unaligned`](#--save_unaligned)
  * [`--save_secondary_alignments`](#--save_secondary_alignments)
  * [`--star_min_overlap [int]`](#--star_min_overlap-int)
  * [`--star_max_overlap_mismatch [float]`](#--star_max_overlap_mismatch-float)
  * [`--star_twopass`](#--star_twopass)
  * [`--percent_mapped_cutoff [float]`](#--percent_mapped_cutoff-float)
* [Read Counting](#read-counting)
  * [`--fc_extra_attributes [str]`](#--fc_extra_attributes-str)
  * [`--fc_group_features [str]`](#--fc_group_features-str)
  * [`--fc_count_type [str]`](#--fc_count_type-str)
  * [`--fc_biotype_group_features [str]`](#--fc_biotype_group_features-str)
  * [`--fc_biotype_count_type [str]`](#--fc_biotype_count_type-str)
  * [`--skip_biotype_qc`](#--skip_biotype_qc)
  * [`--gene_detection_method [str]`](#--gene_detection_method-str)
* [Skipping QC steps](#skipping-qc-steps)
* [Comparisons](#comparisons)
  * [`--skip_deseq2`](#--skip_deseq2)
  * [`--deseq2_fdr [float]`](#--deseq2_fdr-float)
  * [`--deseq2_lfc [float]`](#--deseq2_lfc-float)
  * [`--comparisons`](#--comparisons)
  * [`--gprofiler_fdr [float]`](#--gprofiler_fdr-float)
* [CloudFront settings](#cloudfront-settings)
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
  * [`--kit_QC`](#--kit_QC)
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
nextflow run Zymo-Research/nxf-rnaseq \
  --genome GRCh38 \
  --design "<design CSV file on S3>" \
  -profile awsbatch \
	-work-dir "<work dir on S3>" \
	--awsqueue "<SQS ARN>" \
	--outdir "<output dir on S3>" \
	-r "1.4.0" \
	--cloudfront_private_key "<CloudFront private key file on S3>" \
	--protocol zymo_ribofree \
	--name "<study name>"
```

This command will retrieve the pipeline code from GitHub. This requires that you have setup appropriate GitHub access previlidges. If you have downloaded the pipeline code, you can subsitute `Zymo-Research/nxf-rnaseq` with `main.nf`.

### Running comparison only
Sometimes, one would need to rerun the comparison part of the pipeline without reruning sample processing, for example, with a different FDR cutoff or discarding some bad samples. These is an alternative entrypoint to the pipeline starting with combined read counts of genes, running only `DESeq2`, `gProfiler` and report. To run the pipeline via this entrypoint, the command is as follows:

```bash
nextflow run Zymo-Research/nxf-rnaseq \
  --genome GRCh38 \
  --design "<design CSV file on S3>" \
  -profile awsbatch \
	-work-dir "<work dir on S3>" \
	--awsqueue "<SQS ARN>" \
	--outdir "<output dir on S3>" \
	-r "1.4.0" \
	--cloudfront_private_key "<CloudFront private key file on S3>" \
	--name "<study name>" \
  --merged_counts "<merged counts file on S3>"
```

If you have run the full pipeline on the same data before, you can use the merged counts file generated by the pipeline, usually in the `featureCounts` folder, as the input. If not, you must provide a TSV file similar to that. Note, if you don't need pathway enrichment analysis, then `--genome` can be omitted.

### Running URL generation only
Sometimes, one would need to regenerate URLs for downloading files because the old ones have expired. To run the pipeline for this purpose only, the command is as follows:

```bash
nextflow run Zymo-Research/nxf-rnaseq \
  -profile awsbatch \
	-work-dir "<work dir on S3>" \
	--awsqueue "<SQS ARN>" \
	--outdir "<output dir on S3>" \
	-r "1.4.0" \
	--cloudfront_private_key "<CloudFront private key file on S3>" \
	--name "<study name>" \
  --file_locations "<a text file with file locations on S3>"
```

If you have run the full pipeline on the same data before, which is usually the case, you can use the file locations text file generated by the pipeline, usually in the `download_data` folder, as the input. If not, you must provide a text file with a different file location on each line.

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull Zymo-Research/nxf-rnaseq
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [Zymo-Research/nxf-rnaseq releases page](https://github.com/Zymo-Research/nxf-rnaseq/releases) and find the latest version number - numeric only (eg. `2.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Main arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. While there are multiple profiles originally in the pipeline listed below, Zymo strongly recommends using the `awsbatch` profile, and we have only tested this profile. We do not guarantee other profiles would work.

* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/rnaseq`](http://hub.docker.com/r/nfcore/rnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

When using `awsbatch` profile, one must supply [other options related to AWS batch](#aws-batch-specific-parameters), and supply the locations of [work directory](#-work-dir) and [output directory](#--outdir) on AWS S3.

### `--design`
You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--design 'path/to/data/sample_sheet.csv'
```

#### Full analysis
The `group` identifier should be identical when you have multiple replicates from the same experimental group.

The `sample` identifier should be used to name your samples. That's how they will appear in the report.

These two columns must contain only alphanumerical characters or underscores, and must start with a letter. 

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

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with AWS, the configuration is set up to use a copy of our own [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resources and several others that we have added.

### `--genome`
There are many different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh37` or `--genome GRCh38`
* Mouse
  * `--genome GRCm38`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. All you need to do is to edit the [iGenomes config file](../conf/igenomes.config).

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      star    = '<path to the star index folder>'
      fasta   = '<path to the genome fasta file>' // Used if no star index given
      gtf     = '<path to the genome gtf file>'
      bed12   = '<path to the genome bed file>' // Generated from GTF if not given
      // Other relevant settings 
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

The "GRCh38" set in the [iGenomes config file](../conf/igenomes.config) is a good example on how we at Zymo are doing this for our genome resources.

| Setting        | Description                                                                    |
|----------------|--------------------------------------------------------------------------------|
| `star`         | Path to the STAR index folder                                                  |
| `fasta`        | Path to the genome FASTA file (not required if STAR index is provided)         |
| `gtf`          | Path to the annotation GTF file                                                |
| `bed12`        | Path to the annotation BED file                                                |
| `gprofiler`    | Organism ID used in gProfiler, if applicable                                   |
| `bacteria`     | Indicates this is a bacteria, which will disallow introns in STAR alignment    |
| `ensembl_web`  | The domain of Ensembl website the organism belongs to, if applicable           |
| `rRNA_gtf`     | Path to a GTF file containing additional rRNA regions                          |
| `csi_index`    | Indicates building BAM indices using csi, applicable to very long chromosomes  |

### `--ercc_spikein`
Use this option to indicate that ERCC spike-in was added to the samples. Possible values are `1` and `2`, which indicates Mix 1 and Mix 2 of [ERCC spike-in](https://www.thermofisher.com/order/catalog/product/4456739#/4456739), respectively. When this option is activated, reads that are not aligned to the main genome are aligned again to ERCC92 reference genome and counted. The counts and concentrations of ERCC spike-in transcripts will be displayed in scatter plots in the report.

## Adapter Trimming

### `--protocol`
Use this to specify which library kit or protocol was used to make the RNAseq library. The possible options are `illumina`, `zymo_ribofree`, `zymo_3mrna`, `zymo_3mrna_nodedup`, `pico` with `illumina` being the default. We have chosen the appropriate trimming parameters such as whether to trim bps off 5' or 3' of R1 or R2 and adapter sequences for these protocols. Please refer to [this function](../libs/parse_presets.nf) for more details. You can also modify this file to add your own protocols. 

### `--trim_nextseq [int]`
Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails.

### `--adapter_overlap [int]`
Instructs Trim Galore to require at least this many bp of overlap with adapters to trim a read. Default is 1.

### `--min_read_length [int]`
Instructs Trim Galore to discard reads shorter than this. Default is 20.

### `--save_trimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this flag to copy these files when complete.

## Alignment

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

### `--percent_mapped_cutoff`
By default, we filter samples by requiring at least 5% of reads uniquely mapped to genome by STAR. Those samples failing this filter will be discarded after the alignment step. Use this option to change this cutoff.

## Read Counting
Read counting is conducted using [featureCounts of the Subread package](http://subread.sourceforge.net/). There are several settings that are critical to get this counting right, and they depend on the attributes used in the annotation GTF file. If you are using a standard ENSEMBL genome, then you don't need to change any of the following options. If you are not, you need to carefully examine your GTF file and change the following option accordingly.

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

## Skipping QC steps
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

### `--gene_detection_method [str]`
Use this option to choose which metrics to use for calculating the numbers of genes detected. Available: `reads`, `fpkm`, `tpm`. Default is `reads`.

## Comparisons

### `--skip_deseq2`
Use this option to skip DESeq2, and gProfiler by extension. 

### `--deseq2_fdr [float]`
Use this option to specify the false discovery rate (FDR) cutoff used in DESeq2. Genes whose FDR is smaller than this will be considered significantly differentially expressed. (Default is `0.05`)

### `--deseq2_lfc [float]`
Use this option to specify the fold change cutoff in DESeq2. Please note that this is not a simple filter, this cutoff changes the null hypothesis in the statistical tests, therefore affects the p-values and FDRs in DESeq2. This is passed to DESeq2 as a `Log2FoldChange` parameter, therefore, please use Log2 values of intended fold change cutoff. (Default is `0.585`, which is Log2(1.5))

### `--comparisons`
By default, the pipeline will compare all pairwise combinations of sample groups. If this is not desirable, use this option to specify the path to a CSV file that describes which sample groups you want to comapre. The CSV file should have the following format:
```
group_1,group_2
Experiment1,Control
Experiment2,Control
```
The header must be the same as shown. All group labels should also appear in the design CSV file "group" column.

### `gprofiler_fdr [float]`
Use this option to specify the FDR cutoff used in gProfiler. (Default is `0.05`)

## CloudFront settings
The pipeline provides a way for customer to download original and intermediate data, and final results through AWS CloudFront. The use of AWS CloudFront can be customized, though most of the settings have default values.

| Setting                      | Description                        | Default                                 |
|------------------------------|------------------------------------|-----------------------------------------|
| `--cloudfront_origin_path`   | The origin path for Cloudfront     | `s3://zymo-epiquest2`                   |
| `--cloudfront_domain_name`   | The domain name for Cloudfront     | `https://d2f6fa8431ldcd.cloudfront.net` |
| `--cloudfront_private_key_ID`| The private key ID for Cloudfront  | `APKAIVZC5GR5ZEN7576A`                  |
| `--cloudfront_private_key`   | The private key file path          | `false`                                 |
| `--cloudfront_link_duration` | Number of days links are active    | `60`                                    |
| `--deliver_fastqs`           | Whether to deliver original FASTQs | `false`                                 |

Only when the private key file path is provided, the file downloading function part of the pipeline will work. Please contact Zymo AWS admin for the private key file.

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

### `--kit_QC`
Use this option to provide a CSV file specifying samples for generating scatter plots for the purpose of QC of Zymo-Seq RiboFree Total RNA Library Kit. Must be used with `--protocol zymo_ribofree`. The CSV file must have the following format:
```
library_1,library_2
Dep1,Dep2
Dep1,UnDep1
```
The file header cannot be changed. All sample labels must appear in the `sample` column of the design CSV file.
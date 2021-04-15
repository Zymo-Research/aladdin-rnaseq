# RNA-Seq Nextflow Pipeline

### Introduction

This is a bioinformatics analysis pipeline used for RNA sequencing data. This pipeline was originally adpated from [nf-core/rnaseq](https://github.com/nf-core/rnaseq) pipeline version 1.4.2. Many changes have been made so the current pipeline is very different from the original version or its subsequent versions.

The workflow processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)), aligns the reads ([STAR](https://github.com/alexdobin/STAR), generates gene counts ([featureCounts](http://bioinf.wehi.edu.au/featureCounts/)), and performs extensive quality-control on the results ([RSeQC](http://rseqc.sourceforge.net/), [Qualimap](http://qualimap.bioinfo.cipf.es/), [dupRadar](https://bioconductor.org/packages/release/bioc/html/dupRadar.html), [Preseq](http://smithlabresearch.org/software/preseq/), [MultiQC](http://multiqc.info/)). See the [sample report documentation](https://github.com/Zymo-Research/service-pipeline-documentation/blob/master/docs/how_to_use_RNAseq_report.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

### Updates made by Zymo

We have added sample overview analysis and differential gene expression analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), and a pathway enrichment analysis using [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost). Other updates include improvements of biotype QC process, a section in the report for downloading data, the ability to handle ERCC spike-in, the ability to generate QC results for Zymo RiboFree kit production, other appearance/style changes to the report and various bug fixes or optimizations. 

### Quick Start

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
1. The options `-design` and `--genome` are required.
2. The options `-profile`, `-work-dir`, `--outdir`, and `--awsqueue` are required only when running pipelines on AWS Batch, but are highly recommended.
3. The option `-r` helps pin workflow to a specific release on GitHub. Note: releases prior to 1.2.0 are deprecated, releases after 2.0.0 are recommeded.
4. The option `--cloudfront_private_key` is required only when you want a file download section in the report. Note: This option is named `--cloudfrontPrivateKey` before release 2.0.0.
5. The option `--protocol` is used to identify kits/protocols used to generate the RNAseq library. This will affect trimming and strandedness parameters. Available options are `illumina`, `zymo_ribofree`, `zymo_3mrna`, `zymo_3mrna_nodedup`, and `pico`. The default values is `illumina`. Note: prior to release 2.0.0, `--protocol` is not available, one can use options such as `--zymo_ribofree` to achieve the same purpose.
6. The option `--name` will add a title to the report.
7. The design CSV file must have the following format.
```
group,sample,read_1,read_2
Control,Sample1,s3://mybucket/this_is_s1_R1.fastq.gz,s3://mybucket/this_is_s1_R2.fastq.gz
Control,Sample2,s3://mybucket/this_is_s2_R1.fastq.gz,s3://mybucket/this_is_s2_R2.fastq.gz
Experiment,Sample3,s3://mybucket/that_is_s3_R1.fastq.gz,
Experiment,Sample4,s3://mybucket/that_be_s4_R1.fastq.gz,
```
The header line must be present and cannot be changed. Sample labels and group names must contain only alphanumerical characters or underscores, and must start with a letter. Full S3 paths of R1 and R2 FASTQ files should be provided. Mixing of Paired-end and Single-end data is allowed. If you do not wish to carry out comparisons of samples, i.e., DESeq2 and g:Profiler analysis, simply leave the Group column blank (but keep the coma).<br>

There are many other options built in the pipeline to customize your run and handle specific situations, please refer to the [Usage Documentation](docs/usage.md).

### Credits
These scripts were originally written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/), part of [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden, by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).

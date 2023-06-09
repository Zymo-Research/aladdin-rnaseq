#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = OrderedDict()

regexes['RNAseq pipeline'] = ['v_pipeline.txt', r"(\S+)"]
regexes['Nextflow'] = ['v_nextflow.txt', r"(\S+)"]
regexes['FastQC'] = ['v_fastqc.txt', r"FastQC v(\S+)"]
regexes['UMI-tools'] = ['v_umi_tools.txt', r"version: (\S+)"]
regexes['Trim Galore!'] = ['v_trim_galore.txt', r"version (\S+)"]
regexes['Cutadapt'] = ['v_cutadapt.txt', r"(\S+)"]
regexes['STAR'] = ['v_star.txt', r"STAR_(\S+)"]
regexes['Samtools'] = ['v_samtools.txt', r"samtools (\S+)"]
regexes['Preseq'] = ['v_preseq.txt', r"Version: (\S+)"]
regexes['Picard MarkDuplicates'] = ['v_markduplicates.txt', r"Version:(\S+)"]
regexes['dupRadar'] = ['v_dupRadar.txt', r"(\S+)"]
regexes['RSeQC'] = ['v_rseqc.txt', r"read_duplication.py ([\d\.]+)"]
regexes['Qualimap'] = ['v_qualimap.txt', r"QualiMap v\.(\S+)"]
regexes['featureCounts'] = ['v_featurecounts.txt', r"featureCounts v(\S+)"]
regexes['DESeq2'] = ['v_DESeq2.txt', r"(\S+)"]
regexes['gProfiler'] = ['v_gProfiler.txt', r"gprofiler-official==(\S+)"]

results = OrderedDict()

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'Software Versions'
plot_type: 'html'
description: 'Software versions are collected at run time from the software output. This pipeline is adapted from <a href="https://github.com/nf-core/rnaseq" target="_blank">nf-core RNAseq pipeline</a>.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    # Only display versions of softwares that were actually used.
    if v:
        print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        if v:
            f.write("{}\t{}\n".format(k,v))
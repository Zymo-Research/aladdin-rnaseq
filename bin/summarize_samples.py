#!/usr/bin/env python

###
# This program collects read quantification results from all samples
# and summarize them into combined counts/TPM/FPKM files. 
# It also counts how many genes were detected in each sample.
# The definition of "detected gene" can be changed via options.
###

import argparse
import logging
import json
import pandas as pd
import os

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

def summarize_featurecounts(results, extra_attributes):
    """Combine featureCounts results and calculate FPKM/TPM"""

    dfs = []
    gene_lengths = None
    index_cols = ['Geneid']
    if extra_attributes:
        index_cols += extra_attributes.split(',')
    
    for f in results:
        logger.info("Processing {}".format(f))
        data = pd.read_csv(f, sep="\t", comment="#")
        data = data.set_index(index_cols)
        if gene_lengths is None:
            gene_lengths = data["Length"]
        # Counts in the last column
        data = data.iloc[:,-1]
        # Clean up the sample name
        data.name = data.name.replace(".genome.primary.bam","").replace(".dedupped.bam","")
        dfs.append(data)

    # Combine and output raw counts
    counts = pd.concat(dfs, axis=1)
    counts = counts.reindex(sorted(counts.columns),axis=1)
    counts.to_csv("merged_gene_counts.txt", index=True, sep="\t")

    # Calculate FPKM and TPM
    fpm = counts.apply(lambda x:x.div(x.sum()).mul(1e6))
    fpkm = fpm.div(gene_lengths, axis=0).mul(1e3)
    fpkm.to_csv("gene_FPKM.txt", sep="\t", index=True, float_format="%.3f")
    fpk = counts.div(gene_lengths, axis=0).mul(1e3)
    tpm = fpk.apply(lambda x:x.div(x.sum()).mul(1e6))
    tpm.to_csv("gene_TPM.txt", sep="\t", index=True, float_format="%.3f")

def summarize_salmon(results, filename, level):
    """Combine Salmon results and calculate FPKM"""

    counts = []
    fpkms = []
    tpms = []

    for f in results:
        filepath = os.path.join(f, filename)
        logger.info("Processing {}".format(filepath))
        data = pd.read_csv(filepath, sep="\t")
        data = data.set_index("Name")
        counts.append(data['NumReads'].rename(f))
        tpms.append(data['TPM'].rename(f))
        # Calculate FPKM
        fpm = data["NumReads"].div(data["NumReads"].sum()).mul(1e6)
        fpkm = fpm.div(data['EffectiveLength']).mul(1e3)
        fpkm.name = f
        fpkms.append(fpkm)

    # Output
    counts = pd.concat(counts, axis=1)
    counts = counts.reindex(sorted(counts.columns),axis=1)
    counts.to_csv("merged_{}_counts.txt".format(level), index=True, sep="\t")
    fpkms = pd.concat(fpkms, axis=1)
    fpkms = fpkms.reindex(sorted(fpkms.columns),axis=1)
    fpkms.to_csv("{}_FPKM.txt".format(level), index=True, sep="\t")
    tpms = pd.concat(tpms, axis=1)
    tpms = tpms.reindex(sorted(tpms.columns),axis=1)
    tpms.to_csv("{}_TPM.txt".format(level), index=True, sep="\t")

def count_genes(metric, cutoff, namespace, extra_attributes=None):
    """Count No. detected genes """

    # How the metric will appear in the report
    metric_formal = {'reads': 'No. exonic reads', 'fpkm':'FPKM', 'tpm':'TPM'}
    # Base info for the general stats table
    gs_dict = {
        'id': 'num_genes_detected',
        'plot_type': 'generalstats',
        'pconfig': {
            'num_genes': {
                'title': 'Numbers of genes detected',
                'description': "Numbers of genes with {} >= {:.1f}".format(metric_formal[metric], cutoff),
                'namespace': namespace,
                'min': 0,
                'format': '{:,.0f}' 
            }
        }
    }
    
    # Determine data file based on metric selected
    data_files = {'reads':'merged_gene_counts.txt', 'fpkm':'gene_FPKM.txt', 'tpm':'gene_TPM.txt'}
    data_file = data_files[metric]
    # Read input
    data = pd.read_csv(data_file, sep="\t", index_col=0)
    if extra_attributes:
        extra_cols = extra_attributes.split(',')
        data = data.drop(extra_cols, axis=1)
    # Count genes meeting cutoff
    detected = data.apply(lambda x:sum(x >= cutoff))
    # Covert to intended format
    gs_dict['data'] = detected.to_frame('num_genes').to_dict('index')
    # Write the output to files
    with open('num_genes_detected_mqc.json', 'w') as ofl:
        json.dump(gs_dict, ofl, indent=4)

def summarize_samples(results, quant_method, metric, cutoff, extra_attributes):
    """Main function to carry out all steps"""
    
    # Combined read counts
    if quant_method == 'star_featurecounts':
        summarize_featurecounts(results, extra_attributes)
        count_genes(metric, cutoff, 'featureCounts', extra_attributes)
    else:
        summarize_salmon(results, 'quant.genes.sf', 'gene')
        summarize_salmon(results, 'quant.sf', 'transcript')
        count_genes(metric, cutoff, 'Salmon', None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Summarize read quantification results and count the numbers of genes detected""")
    parser.add_argument("results", type=str, nargs="+", help="File(s) or Folder(s) with read quantification results")
    parser.add_argument("-q", "--quant_method", dest="quant_method", type=str.lower, default="star_featurecounts",
                        choices=["star_featurecounts", "star_salmon"], help="Which read quantification method was used")
    parser.add_argument("-m", "--metric", dest="metric", type=str.lower, default="reads", choices=["reads","fpkm","tpm"], help="Which metric to use to determine if a gene is detected") 
    parser.add_argument("-c", "--cutoff", dest="cutoff", type=float, default=1.0, help="Cutoff for selected metric")
    parser.add_argument("-e", "--extra_attributes", dest="extra_attributes", type=str, default=None, help="Extra attributes featureCounts was asked to collect")
    args = parser.parse_args()
    summarize_samples(args.results, args.quant_method, args.metric, args.cutoff, args.extra_attributes)

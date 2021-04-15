#!/usr/bin/env python

###
# This program counts how many genes were detected in each sample.
# The definition of "detected gene" can be changed via options.
###

import argparse
import logging
import json
import pandas as pd
from collections import defaultdict

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

def count_genes(source, metric, cutoff):
    
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
                'namespace': 'featureCounts',
                'min': 0,
                'format': '{:,.0f}' 
            }
        }
    }
    
    # Read input
    data = pd.read_csv(source, sep="\t", index_col=0)
    if metric == 'reads':
        data = data.iloc[:,1:] # Skip gene_name column
    # Count genes meeting cutoff
    detected = data.apply(lambda x:sum(x >= cutoff))
    # Covert to intended format
    gs_dict['data'] = detected.to_frame('num_genes').to_dict('index')
    # Write the output to files
    with open('num_genes_detected_mqc.json', 'w') as ofl:
        json.dump(gs_dict, ofl, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Count the numbers of genes detected""")
    parser.add_argument("source", type=str, help="Source data, either read counts, FPKM, or TPM values")
    parser.add_argument("-m", "--metric", dest='metric', type=str.lower, default='reads', choices=['reads','fpkm','tpm'], help="Which metric to use to determine if a gene is detected") 
    parser.add_argument("-c", "--cutoff", dest='cutoff', type=float, default=1.0, help="Cutoff for selected metric")
    args = parser.parse_args()
    count_genes(args.source, args.metric, args.cutoff)

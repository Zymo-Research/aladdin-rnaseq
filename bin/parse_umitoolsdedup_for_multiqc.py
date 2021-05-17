#!/usr/bin/env python

###
# This program captures the input and output read numbers from umi-tools dedup logs from each sample and adds the results to the general_stats section.
###

import argparse
import logging
import json
import os
from collections import defaultdict

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

def get_dedup_stats(dedup_logs):

    # Base info for the general stats table
    gs_dict = {
        'id': 'reads_post_dedup',
        'plot_type': 'generalstats',
        'pconfig': {
            'kept_dedupped': {
                'title': '% Dedupped',
                'description': '% Reads post dedup',
                'namespace': 'UMI-tools',
                'min': 0,
                'format': '{:.02%}' 
            }
        }
    }

    data = defaultdict(lambda: defaultdict(float))
    for f in dedup_logs:
        # Extract sample name
        sname = os.path.basename(f).split('.')[0]
        with open(f, 'r') as f:
            lines = f.readlines()
        input_reads = int(lines[-6].split(' ')[6])
        output_reads = int(lines[-5].split(' ')[7])
        data[sname]['kept_dedupped'] = output_reads/input_reads

    # Add data to the table dict
    gs_dict['data'] = data
    # Write the output to files
    with open('percent_reads_post_dedup_mqc.json', 'w') as prpd:
        json.dump(gs_dict, prpd, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Parse UMI-tools dedup results to report the number of reads output post dedup""")
    parser.add_argument("dedup_logs", type=str, nargs='+', help="UMI-tools dedup logfiles")
    args = parser.parse_args()
    get_dedup_stats(args.dedup_logs)
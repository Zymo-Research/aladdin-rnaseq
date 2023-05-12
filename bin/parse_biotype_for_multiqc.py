#!/usr/bin/env python

###
# This program aggregates biotype QC results and does the following parsing:
# 1. Only display the top n and mandatory biotypes, and collapsing the rest into "Others",
#    to limit the number categories shown in MultiQC report.
# 2. Calculate percentages of requested categories and display them in general_stats section
#    of MultiQC report.
###

import argparse
import logging
import os
import json
import re
from collections import OrderedDict, defaultdict

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

biotype_dict = {
    'id': 'biotype-counts',
    'section_name': 'Biotype Counts',
    'description': "shows reads overlapping genomic features of different biotypes, counted by <a href='http://bioinf.wehi.edu.au/featureCounts'>featureCounts</a>.",
    'plot_type': 'bargraph',
    'anchor': 'featurecounts_biotype',
    'pconfig': {
        'id': 'featureCounts_biotype_plot',
        'title': 'featureCounts: Biotypes',
        'cpswitch_counts_label': "Number of Reads",
        "cpswitch_c_active": False
    }
}

biotype_gs_dict = {
    'id': 'biotype-gs',
    'plot_type': 'generalstats',
    'pconfig': dict()
}

log_patterns = {
    'Unassigned_NoFeatures': r"Unassigned_NoFeatures\t(\d+)",
    'Unassigned_Ambiguity': r"Unassigned_Ambiguity\t(\d+)"
}

def mqc_feature_stat(counts_files, features, top_n):
    
    biotype_data = OrderedDict()
    biotype_gs_data = dict()
    # Sort counts file by name
    for bfile in sorted(counts_files):
        # Extract sample name
        sname = os.path.basename(bfile).split('.')[0].replace("_biotype", "")
        # Try to parse and read biocount file 
        fcounts = dict()
        try:
            with open(bfile, 'r') as bfl:
                for ln in bfl:
                    if ln.startswith('#'):
                        continue
                    ft, cn = ln.strip().split('\t')
                    fcounts[ft] = float(cn)
        except:
            logger.error("Trouble reading the biocount file {}".format(bfile))
            return
        # Try to read and parse the corresponding featureCounts summary file
        if os.path.isfile(bfile+'.summary'):
            try:
                with open(bfile+'.summary', 'r') as sfl:
                    txt = sfl.read()
                    for k,v in log_patterns.items():
                        m = re.search(v, txt)
                        fcounts[k] = float(m.group(1))
            except:
                logger.error("Trouble reading the featureCounts summary file {}.summary".format(bfile))
        else:
            logger.warning("featureCounts summary file {}.summary NOT found.".format(bfile))
            for k in log_patterns:
                fcounts[k] = 0
        # Try to read and parse ERCC featureCounts file
        if os.path.isfile(sname+'_ERCC.featureCounts.txt.summary'):
            try:
                with open(sname+'_ERCC.featureCounts.txt.summary', 'r') as efl:
                    txt = efl.read()
                    m = re.search(r"Assigned\t(\d+)", txt)
                    fcounts['ERCC spike-in'] = float(m.group(1))
            except:
                logger.error("Trouble reading the ERCC featureCounts summary file {}_ERCC.featureCounts.txt.summary".format(sname))
        # Calculate total
        total_count = sum(fcounts.values())
        if total_count == 0:
            logger.error("No biocounts found, exiting")
            return
        # Calculate percentage for each requested feature
        fpercent = dict()
        for ft in features:
            if ft in fcounts:
                fpercent[ft] = fcounts[ft] / total_count * 100
                # Add general stats column config if not added already
                if ft not in biotype_gs_dict['pconfig']:
                    biotype_gs_dict['pconfig'][ft] = {
                        'title': '% {}'.format(ft),
                        'namespace': 'Biotype Counts',
                        'description': '% reads overlapping {} features'.format(ft),
                        'max': 100,
                        'min': 0,
                        'scale': 'RdYlGn-rev',
                        'format': '{:.2f}%'
                    }
            else:
                logger.warning("Requested feature {} not found in biocount file {}".format(ft, bfile))
                fpercent[ft] = 'N/A'
        # Prepare the output dict for general stats
        biotype_gs_data[sname] = fpercent
        # Record the biotype counts
        biotype_data[sname] = fcounts
    # Finalize output dict for general stats    
    biotype_gs_dict['data'] = biotype_gs_data
    
    # Sort categories by total counts in all samples
    total_counts = defaultdict(float)
    for data in biotype_data.values():
        # Not including unassigned reads
        for k, v in data.items():
            if k not in log_patterns:
                total_counts[k] += v
    sorted_biotypes = sorted(total_counts, key=total_counts.get, reverse=True)
    # Select categories to keep
    keep = sorted_biotypes[:top_n]
    # Always add tracked categories
    for f in features:
        if f not in keep:
            keep.append(f)
    # Collapse the other categories into "others"
    for sname, data in biotype_data.items():
        # Use OrderedDict to keep order of categories from most to least
        d = OrderedDict(dict.fromkeys(keep, 0))
        d['others'] = 0
        for k, v in data.items():
            if k in keep or k in log_patterns or k == 'ERCC spike-in':
                d[k] = v
            else:
                d['others'] += v
        biotype_data[sname] = d
    # Add the data to the output dict
    biotype_dict['data'] = biotype_data

    # Write the output to files
    with open('biotype_mqc.json', 'w') as ofl:
        json.dump(biotype_dict, ofl, indent=4)
    with open('biotype_gs_mqc.json', 'w') as ofl:
        json.dump(biotype_gs_dict, ofl, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Merge and parse biotype QC results""")
    parser.add_argument("biocounts", type=str, nargs='+', help="File(s) with all biocounts")
    parser.add_argument("-f", "--features", dest='features', required=True, nargs='+', help="Features to display in general stats")
    parser.add_argument("-n", "--top_n_biotypes", type=int, dest='top_n', default=5, help="Top n biotypes to display in MultiQC report")
    args = parser.parse_args()
    mqc_feature_stat(args.biocounts, args.features, args.top_n)

#!/usr/bin/env python

###
# This program extracts number of reads aligned from STAR log
# and total number of read assigned from Salmon log, and generate
# MultiQC report sections
###

import argparse
import logging
import json
import re

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)
    
def parse_star_log(star_log):
    """Extract total number of reads aligned from STAR"""

    logger.info("Parsing STAR log {}".format(star_log))
    regexes = {
        "uniquely_mapped": r"Uniquely mapped reads number \|\s+(\d+)",
        "multimapped": r"Number of reads mapped to multiple loci \|\s+(\d+)"
    }
    parsed_data = {}
    with open(star_log) as fh:
        content = fh.read()
        for key, pattern in regexes.items():
            m = re.search(pattern, content)
            if m:
                parsed_data[key] = float(m.group(1))
    total = parsed_data['uniquely_mapped'] + parsed_data['multimapped']
    return total

def parse_dedup_log(dedup_log):
    """Extract total number of reads after UMI-tools deduplication"""
    
    logger.info("Parsing UMI-tools dedup log {}".format(dedup_log))
    pattern = r"Number of reads out: (\d+)"
    with open(dedup_log) as fh:
        content = fh.read()
        m = re.search(pattern, content)
        return int(m.group(1))

def parse_salmon_log(salmon_log):
    """Extract stats from Salmon quant log"""

    logger.info("Parsing Salmon quant log {}".format(salmon_log))
    regexes = {
        "total_reads": r"Total # of mapped reads :\s+(\d+)",
        "assigned": r"Counted (\d+) total reads in the equivalence classes"
    }
    parsed_data = {}
    with open(salmon_log) as fh:
        content = fh.read()
        for key, pattern in regexes.items():
            m = re.search(pattern, content, re.MULTILINE)
            if m:
                parsed_data[key] = float(m.group(1))
    return parsed_data["assigned"], parsed_data["total_reads"]

def parse_salmon_quant_stats(salmon_log, star_log, dedup_log):
    """Parse Salmon quant log and either STAR log or UMI-tools dedup for MultiQC"""

    # Parse Salmon log
    num_reads_assigned, num_reads_procssed = parse_salmon_log(salmon_log)
    # Parse STAR log
    if star_log is not None:
        sample_name = star_log.replace('Log.final.out', '')
        num_reads_aligned = parse_star_log(star_log)
    # Parse UMI-tools dedup log
    else:
        sample_name = dedup_log.replace('.dedup_LOGFILE', '')
        num_reads_aligned = parse_dedup_log(dedup_log)
    # Do some calculations
    num_reads_not_transcriptome = num_reads_aligned - num_reads_procssed
    num_reads_filtered = num_reads_procssed - num_reads_assigned
    pct_assigned = num_reads_assigned / num_reads_aligned * 100
    # Build the Multiqc input files
    section_dict = {
        'id': 'reads_assigned_STAR_Salmon',
        'section_name': 'Read quantification (STAR_Salmon)',
        'description': ("This section shows reads assigned to genes by <a href='http://combine-lab.github.io/salmon/'>Salmon</a> "
                        "from <a href='https://github.com/alexdobin/STAR'>STAR</a> alignments."),
        'plot_type': 'bargraph',
        'anchor': 'salmon_quant',
        'pconfig': {
            'id': 'salmon_quant_plot',
            'title': 'Salmon quant results from STAR alignments',
            'cpswitch_counts_label': "Numbers of Reads"
        },
        'data' : {
            sample_name: {
                'assigned to genes': num_reads_assigned,
                'not aligned to any transcript': num_reads_not_transcriptome,
                'filtered by Salmon': num_reads_filtered
            }
        }
    }
    gs_dict = {
        'id': 'reads_assigned_STAR_Salmon_gs',
        'plot_type': 'generalstats',
        'pconfig' : {
            'pct_assigned': {
                'title': '% Assigned',
                'namespace': 'Salmon quant',
                'description': '% reads assigned to genes',
                'max': 100,
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:.2f}%'
            },
            'reads_assigned': {
                'title': 'M reads Assigned',
                'namespace': 'Salmon quant',
                'description': 'No. reads assigned to genes (in millions)',
                'min': 0,
                'scale': 'PuBu'
            }
        },
        'data' : {
            sample_name: {
                'pct_assigned': pct_assigned,
                'reads_assigned': num_reads_assigned / 1000000
            }
        }
    }
    # Write the output to files
    with open("{}_salmon_quant_mqc.json".format(sample_name), 'w') as ofl:
        json.dump(section_dict, ofl, indent=4)
    with open("{}_salmon_quant_gs_mqc.json".format(sample_name), 'w') as ofl:
        json.dump(gs_dict, ofl, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Merge and parse STAR and Salmon quant logs""")
    parser.add_argument("-l", "--salmon_log", type=str, dest='salmon_log', required=True, help="File path of Salmon quant log")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--star_log", type=str, dest='star_log', help="File path of STAR log")
    group.add_argument("-d", "--dedup_log", type=str, dest='dedup_log', help="File path of UMI-tools dedup log")
    args = parser.parse_args()
    parse_salmon_quant_stats(args.salmon_log, args.star_log, args.dedup_log)
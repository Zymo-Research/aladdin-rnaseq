#!/usr/bin/env python

###
# This programs simply adds the ERCC concentration to the ERCC read counts table
###

import argparse
import pandas as pd

def add_ercc_info(counts_fn, info_fn, mixture, output_name):
    counts = pd.read_csv(counts_fn, sep="\t", index_col=0)
    info = pd.read_csv(info_fn, sep="\t", index_col=1) #ERCC IDs are in 2nd column
    conc_col_name = "concentration in Mix {} (attomoles/ul)".format(mixture)
    concentration = info[[conc_col_name,'length']]
    concentration.columns = ['concentration','length']
    counts = counts.join(concentration)
    counts.to_csv(output_name, sep="\t", index=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Add ERCC concentration to the counts matrix""")
    parser.add_argument("-c", "--counts", dest="counts", type=str, required=True, help="Read counts file")
    parser.add_argument("-i", "--info", dest="info", type=str, required=True, default="ERCC info file")
    parser.add_argument("-m", "--mixture", dest="mixture", type=str, required=True,
                        choices=['1', '2'], help="ERCC Mixure number")
    parser.add_argument("-o", "--output_name", dest="output_name", type=str, default="ERCC_counts_info.tsv",
                        help="Output file name")
    args = parser.parse_args()
    add_ercc_info(args.counts, args.info, args.mixture, args.output_name)
#!/usr/bin/env python

###
# This program calculates FPKM and TPM values
# from featureCounts read counts.
# Calculation is according to FPKM and TPM definitions.
###

import argparse
import logging
import pandas as pd

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

def calculate_fpkm_tpm(merged_counts, gene_length):

    # Read input files
    counts = pd.read_csv(merged_counts, sep="\t", index_col=0)
    counts = counts.iloc[:,1:] # Dicard gene_name column
    lengths = pd.read_csv(gene_length, sep="\t", index_col=0)
    # Calculate fpkm
    fpm = counts.apply(lambda x:x.div(x.sum()).mul(1e6))
    fpkm = fpm.div(lengths['Length'], axis=0).mul(1e3)
    fpkm.to_csv("FPKM.tsv", sep="\t", index=True, float_format='%.3f')
    # Calculate tpm
    fpk = counts.div(lengths['Length'], axis=0).mul(1e3)
    tpm = fpk.apply(lambda x:x.div(x.sum()).mul(1e6))
    tpm.to_csv("TPM.tsv", sep="\t", index=True, float_format='%.3f')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Calculate FPKM,TPM values from featureCounts merged counts""")
    parser.add_argument("merged_counts", type=str, help="Merged read counts file")
    parser.add_argument("-l", "--length", dest='gene_length', type=str, required=True, help="Gene length file") 
    args = parser.parse_args()
    calculate_fpkm_tpm(args.merged_counts, args.gene_length)
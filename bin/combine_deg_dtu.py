#!/usr/bin/env python

import argparse
import pandas as pd
import os

def combine_results(deg_results, dtu_results):
    deg_suffix = "_DESeq_results.tsv"
    dtu_suffix = "_DTU_analysis_DEXSeq_results.tsv"
    deg = dict()
    dtu = dict()
    for fn in deg_results:
        name = os.path.basename(fn)
        if deg_suffix not in name:
            raise ValueError("Unexpected naming pattern for DEG results file {}!".format(name))
        else:
            name = name.replace(deg_suffix, "")
        data = pd.read_csv(fn, sep="\t")
        data = data[["Gene ID", "Gene Name", "padj"]]
        data = data.rename(columns={"padj":"DEG_padj"})
        deg[name] = data
    for fn in dtu_results:
        name = os.path.basename(fn)
        if dtu_suffix not in name:
            raise ValueError("Unexpected naming pattern for DTU results file {}!".format(name))
        else:
            name = name.replace(dtu_suffix, "")
        data = pd.read_csv(fn, sep="\t")
        data = data[["Gene ID", "Gene Name", "gene_padj"]].drop_duplicates()
        data = data.rename(columns={"gene_padj":"DTU_padj"})
        dtu[name] = data
    for name, deg_data in deg.items():
        try:
            dtu_data = dtu[name]
            data = pd.merge(deg_data, dtu_data, how="outer")
            data = data.sort_values(["DEG_padj", "DTU_padj", "Gene ID"])
            data.to_csv("{}_combined_DEG_DTU_padj.tsv".format(name), sep="\t", index=False)
            data.to_excel("{}_combined_DEG_DTU_padj.xlsx".format(name), index=False)
        except KeyError:
            continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Combine DTU results with DEG results into one file so it's easy to check whether a gene is DEG and/or DTU.""")
    parser.add_argument("-g", "--deg_results", dest="deg_results", nargs="+", type=str, help="DEG results files")
    parser.add_argument("-t", "--dtu_results", dest="dtu_results", nargs="+", type=str, help="DTU results files")
    args = parser.parse_args()
    combine_results(args.deg_results, args.dtu_results)
#!/usr/bin/env python

import argparse
import logging
import pandas as pd
import os
from gprofiler import GProfiler

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

def run_gprofiler(deseq_result, organism, deseq_fdr, gprofiler_fdr,
                  sources=["GO:BP", "REAC", "KEGG"], domain_scope="annotated",
                  sname=None):
    
    # If sample name not given use file name
    if not sname:
        sname = os.path.splitext(os.path.basename(deseq_result))[0]
        sname = sname.replace("_DESeq_results", "")
    
    # Extract two group names
    try:
        g1, g2 = sname.split('_vs_')
    except:
        logger.error("The sample name are not in the format of '[Group1]_vs_[Group2]")
        return
    
    # Read the DESeq results
    try:
        data = pd.read_csv(deseq_result, sep="\t")
        # Extract up and down-regulated genes
        query = dict()
        up_genes = data.loc[(data["padj"]<deseq_fdr)&(data["log2FoldChange"]>0),]\
                   .sort_values("padj")["Gene ID"].tolist()
        logger.info("Found {} genes up-regulated in {}".format(len(up_genes), g1))
        if len(up_genes)>0:
            query["Higher in "+g1] = up_genes
        down_genes = data.loc[(data["padj"]<deseq_fdr)&(data["log2FoldChange"]<0),]\
                     .sort_values("padj")["Gene ID"].tolist()
        logger.info("Found {} genes up-regulated in {}".format(len(down_genes), g2))
        if len(down_genes)>0:
            query["Higher in "+g2] = down_genes
        if len(query) == 0:
            query = ["empty"]
    except:
        logger.error("Trouble parsing DESeq result file {}".format(deseq_result))
        return

    # Run gProfiler
    gp = GProfiler(return_dataframe=True)
    try:
        result = gp.profile(query=query, organism=organism, sources=sources,
                            user_threshold=gprofiler_fdr, all_results=True,
                            domain_scope=domain_scope, ordered=True)
    except:
        logger.error("Trouble running gProfiler")
        return
    
    # Parse and output results
    result = result.rename(columns={"native":"ID","query":"expression pattern"})
    result.to_csv(sname+"_gProfiler_results.tsv", sep="\t", index=False)
    result.to_excel(sname+"_gProfiler_results.xlsx", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Run gProfiler pathway enrichment analysis from DESeq2 results""")
    parser.add_argument("deseq_results", type=str, nargs='+', help="TSV file(s) with DESeq2 results")
    parser.add_argument("-o", "--organism", dest="organism", required=True, help="Organism code")
    parser.add_argument("-q", "--deseqFDR", dest="deseq_fdr", type=float, default=0.05, help="FDR used in DESeq2")
    parser.add_argument("-p", "--gprofilerFDR", dest="gprofiler_fdr", type=float, default=0.05, help="FDR to use in gProfiler")
    args = parser.parse_args()
    for f in args.deseq_results:
        logger.info("Processing {}...".format(f))
        run_gprofiler(f, args.organism, args.deseq_fdr, args.gprofiler_fdr)
#!/usr/bin/env python

###
# This program produces a transcript-to-gene mapping TSV file
# that is used in tximport from GTF file. It also collects 
# extraAttributes such as gene_name from GTF file.
###

import argparse
import logging
import re
from collections import OrderedDict

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

def tx2gene(gtf, extra_attributes, feature, outfile):
    """Make a transcript-to-gene mapping TSV file from GTF file"""

    pattern = r'([^"]+) "([^"]*)"; ?' # GTF attributes format
    mappings = OrderedDict()
    output_header = ["transcript_id", "gene_id"]
    if extra_attributes:
        extra_attributes = extra_attributes.split(',') # Parse the extra_attributes string
        output_header += extra_attributes
    
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if cols[2] == feature: # Only use selected feature type
                attributes = dict()
                for match in re.finditer(pattern, cols[8]):
                    key, value = match.group(1), match.group(2)
                    if value:
                        attributes[key] = value
                if ("gene_id" not in attributes) or ("transcript_id" not in attributes):
                    logger.error("gene_id or transcript_id not found in {}".format(line.strip()))
                # Check if there is any conflict gene_ids
                elif attributes["transcript_id"] in mappings:
                    if mappings[attributes["transcript_id"]]["gene_id"] != attributes["gene_id"]:
                        raise Exception("Conflict gene_ids found for {}".format(attributes["transcript_id"]))
                else:
                    mappings[attributes["transcript_id"]] = attributes
    
    if not mappings:
        raise Exception("No mappings found!")

    with open(outfile, "w") as fh:
        fh.write("\t".join(output_header))
        fh.write("\n")
        for transcript, attributes in mappings.items():
            fh.write("{}\t{}".format(transcript, attributes["gene_id"]))
            for a in extra_attributes:
                fh.write("\t{}".format(attributes.get(a, '')))
            fh.write("\n")            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Make a transcript-to-gene mapping file from GTF file""")
    parser.add_argument("gtf", type=str, help="GTF file")
    parser.add_argument("-e", "--extra_attributes", dest="extra_attributes", type=str, default="gene_name",
                        help="Extra attributes to include in the output file")
    parser.add_argument("-f", "--feature", dest="feature", type=str, default="transcript",
                        help="Feature type from which to extract transcript-to-gene mapping information, refers to 3rd column in GTF file")
    parser.add_argument("-o", "--outfile", dest="outfile", type=str, default="tx2gene.tsv",
                        help="Output file name")
    args = parser.parse_args()
    tx2gene(args.gtf, args.extra_attributes, args.feature, args.outfile)
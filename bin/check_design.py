#!/usr/bin/env python

import argparse
import re
import pandas as pd
import sys

def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)

def match_legal_pattern(label):
    '''Helper function to check if labels are legal'''
    legal_pattern = r"^[a-zA-Z][a-zA-Z0-9_]*$"
    if not re.match(legal_pattern, label):
        print_error("Sample/Group label contain illegal characters or does not start with letters", "Label", label)
    # MultiQC reserved strings
    illegal_patterns = ["_tophat", "ReadsPerGene", "_star_aligned", "_fastqc", "_counts", "Aligned", "_slamdunk", "_bismark", "_SummaryStatistics", \
                        "_duprate", "_vep", "ccs", "_NanoStats", r"_trimmed$", r"_val$", r"_mqc$", r"short_summary_$", r"_summary$", r"_matrix$", r"_$"]
    for p in illegal_patterns:
        if re.match(p, label):
            print_error("Sample/Group label {} contain string reserved by pipeline or MultiQC. This may cause a problem.".format(label),
                        "Reserved string", p)
    return True
    
def check_fastq_suffix(filename):
    '''Helper function to check if filenames have the correct suffix'''
    FQ_EXTENSIONS = (".fq.gz", ".fastq.gz")
    if not filename.endswith(FQ_EXTENSIONS):
        print_error("FASTQ path has invalid extension", "Path", filename)
    return True

def check_all_se_or_all_pe(group):
    '''Helper function to check if all runs of same sample are either all single-ended or all paired-ended'''
    return group['read_2'].count() == 0 or group['read_2'].count() == group['read_1'].count()

def check_design(DesignFileIn, DesignFileOut, DEseq2GroupLabels, Comparison):
    """
    This function checks that the samplesheet follows the following structure:

    sample,read_1,read_2,group
    sample1,s1_run1_R1.fastq.gz,s1_run1_R2.fastq.gz,groupA
    sample1,s1_run2_R1.fastq.gz,s1_run2_R2.fastq.gz,groupA
    sample2,s2_run1_R1.fastq.gz,,groupB
    sample3,s3_run1_R1.fastq.gz,s3_run1_R2.fastq.gz,groupB

    Use same sample labels for different sequencing runs of the same sample.
    """

    design = pd.read_csv(DesignFileIn, index_col=False)

    # Check column names. Additional columns are OK but will be ignored.
    HEADER = {"sample", "read_1", "read_2", "group"}
    assert HEADER.issubset(set(design.columns)), "Design file must contain {} columns".format(','.join(HEADER))

    # Check if all sample labels legal
    assert design["sample"].map(match_legal_pattern).all()

    # Make sure group labels are either all missing or all present
    assert design["group"].isna().all() or design["group"].notnull().all(), "Group labels missing in some samples but not others!"

    # Check if all group labels legal
    assert design["group"].map(match_legal_pattern, na_action="ignore").all()
    
    # Check FASTQ extension
    assert design["read_1"].map(check_fastq_suffix).all()
    assert design["read_2"].map(check_fastq_suffix, na_action="ignore").all()
        
    # Make sure there is no duplication of FASTQ locations
    assert design["read_1"].duplicated().any() == False, "There are duplications within Read 1 paths!"
    assert design["read_2"].dropna().duplicated().any() == False, "There are duplications within Read 2 paths!"

    # Make sure group labels are consistent for the same sample, different runs
    group_label_consistent = design.groupby("sample")["group"].apply(lambda x:len(x.unique())==1)
    if not group_label_consistent.all():
        inconsistent_group_labels = group_label_consistent[~group_label_consistent].index.tolist()
        print_error("Group labels for samples not consistent across different runs", "Samples", ",".join(inconsistent_group_labels))

    # Make sure there are no mixtures of single and paired-end data within the same sample, different runs
    sample_all_se_or_pe = design.groupby("sample").apply(check_all_se_or_all_pe)
    if not sample_all_se_or_pe.all():
        bad_samples = sample_all_se_or_pe[~sample_all_se_or_pe].index.tolist()
        print_error("Single-end and paired-end data cannot be mixed for the same sample", "Samples", ",".join(bad_samples))

    if Comparison is not None:
        comparisons = pd.read_csv(Comparison, index_col=None)
        HEADER = ['group_1', 'group_2']
        assert list(comparisons.columns) == HEADER, "Header of group comparison file must be {}!".format(','.join(HEADER))
        assert comparisons['group_1'].isin(design['group']).all(), "Some labels in group_1 column not present in design file!"
        assert comparisons['group_2'].isin(design['group']).all(), "Some labels in group_2 column not present in design file!"

    # Output the design file after all check passed
    design.to_csv(DesignFileOut, index=False)

    # Make a Group label files for DESeq2
    design[['group','sample']].drop_duplicates().to_csv(DEseq2GroupLabels, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Sanity check the RNAseq design CSV file, the group comparison file, and the scatterplot datasets file.""")
    parser.add_argument("DesignFileIn", type=str, help="Input design CSV file")
    parser.add_argument("-c", "--comparison", dest="Comparison", help="Group comparison CSV file")
    args = parser.parse_args()
    check_design(args.DesignFileIn, "checked_"+args.DesignFileIn, "DESeq2_"+args.DesignFileIn, args.Comparison)

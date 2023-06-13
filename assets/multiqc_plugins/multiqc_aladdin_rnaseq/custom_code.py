#!/usr/bin/env python
"""
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_aladdin_rnaseq_version = get_distribution("multiqc_aladdin_rnaseq").version

# Add default config options that can be overriden by user config
def plugin_before_config():
    
    # Use the aladdin template by default
    config.template = 'aladdin'
    
# Add additional config options
def plugin_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsed.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', False) is True:
        return None

    log.info("Running rnaseq MultiQC Plugins v{}".format(config.multiqc_aladdin_rnaseq_version))

    # Add to the search patterns used by modules
    if 'Trim_Galore' not in config.sp:
        config.update_dict( config.sp, { 'Trim_Galore' : { 'contents': 'This is cutadapt', 'shared': True } } )
    if 'plot_ERCC' not in config.sp:
        config.update_dict( config.sp, { 'plot_ERCC': { 'fn' : '*plot_ERCC.csv' } } )
    if 'plot_sample_distance/heatmap' not in config.sp:
        config.update_dict( config.sp, { 'plot_sample_distance/heatmap': [ {'fn' : '*sample_distance_matrix*'}, {'fn' : '*sample_similarity_matrix*'} ] } )
    if 'plot_sample_distance/pca' not in config.sp:
        config.update_dict( config.sp, { 'plot_sample_distance/pca': [ {'fn' : '*sample_pca_plot.tsv'}, {'fn' : '*sample_PCA_plot.tsv'}, {'fn' : '*sample_MDS_plot.tsv'}, {'fn' : '*sample_mds_plot.tsv'} ] } )
    if 'plot_gene_heatmap' not in config.sp:
        config.update_dict( config.sp, { 'plot_gene_heatmap': { 'fn' : '*gene_heatmap.tsv' } } )
    if 'DESeq2' not in config.sp:
        config.update_dict( config.sp, { 'DESeq2' : { 'fn' : '*DESeq_results.tsv' } } )
    if 'gProfiler' not in config.sp:
        config.update_dict( config.sp, { 'gProfiler' : { 'fn' : '*gProfiler_results.tsv' } } )
    if 'trimming_2step/1st' not in config.sp:
        config.update_dict( config.sp, { 'trimming_2step/1st': { 'fn' : '*_1st_adapter_trimming_report.txt' } } )
    if 'trimming_2step/2nd' not in config.sp:
        config.update_dict( config.sp, { 'trimming_2step/2nd': { 'fn' : '*_2nd_polyA_trimming_report.txt' } } )
    
    # Some additional filename cleaning
    config.fn_clean_exts.extend([
        '_plot_ERCC',
        '_DESeq_results',
        '_gProfiler_results',
        '_trimmed_first'
    ])

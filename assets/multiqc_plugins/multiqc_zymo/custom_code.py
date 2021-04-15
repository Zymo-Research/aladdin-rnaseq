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
config.multiqc_zymo_version = get_distribution("multiqc_zymo").version

# Add default config options that can be overriden by user config
def plugin_before_config():
    
    # Use the zymo template by default
    config.template = 'zymo'
    
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

    log.info("Running rnaseq MultiQC Plugins v{}".format(config.multiqc_zymo_version))

    # Add to the search patterns used by modules
    if 'Trim_Galore' not in config.sp:
        config.update_dict( config.sp, { 'Trim_Galore' : { 'contents': 'This is cutadapt', 'shared': True } } )
    if 'plot_ERCC' not in config.sp:
        config.update_dict( config.sp, { 'plot_ERCC': { 'fn' : '*plot_ERCC*' } } )
    if 'plot_sample_distance/heatmap' not in config.sp:
        config.update_dict( config.sp, { 'plot_sample_distance/heatmap': [ {'fn' : '*sample_distance_matrix*'}, {'fn' : '*sample_similarity_matrix*'} ] } )
    if 'plot_sample_distance/pca' not in config.sp:
        config.update_dict( config.sp, { 'plot_sample_distance/pca': [ {'fn' : '*sample_pca_plot*'}, {'fn' : '*sample_PCA_plot*'}, {'fn' : '*sample_MDS_plot*'}, {'fn' : '*sample_mds_plot*'} ] } )
    if 'plot_gene_heatmap' not in config.sp:
        config.update_dict( config.sp, { 'plot_gene_heatmap': { 'fn' : '*gene_heatmap*' } } )
    if 'DESeq2' not in config.sp:
        config.update_dict( config.sp, { 'DESeq2' : { 'fn' : '*DESeq_results*' } } )
    if 'gProfiler' not in config.sp:
        config.update_dict( config.sp, { 'gProfiler' : { 'fn' : '*gProfiler_results*' } } )
    if 'download_data' not in config.sp:
        config.update_dict( config.sp, { 'download_data' : { 'fn' : '*download_links*.json', 'shared': True } } )

    
    # Some additional filename cleaning
    config.fn_clean_exts.extend([
        '_R1',
        '_R2',
        '_plot_ERCC',
        '_DESeq_results',
        '_gProfiler_results'
    ])

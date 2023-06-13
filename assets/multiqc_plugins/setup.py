#!/usr/bin/env python
"""
Setup code for Aladdin RNAseq pipeline MultiQC plugin.
For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '1.0.0'

setup(
    name = 'multiqc_aladdin_rnaseq',
    version = version,
    description = "MultiQC plugins for rnaseq pipeline",
    packages = find_packages(),
    include_package_data = True,
    install_requires = ['multiqc==1.9'],
    entry_points = {
        'multiqc.templates.v1': [
            'aladdin = multiqc_aladdin_rnaseq.templates.aladdin'
        ],
        'multiqc.modules.v1': [
            'Trim_Galore = multiqc_aladdin_rnaseq.modules.Trim_Galore:MultiqcModule',
            'plot_ERCC = multiqc_aladdin_rnaseq.modules.plot_ERCC:MultiqcModule',
            'plot_sample_distance = multiqc_aladdin_rnaseq.modules.plot_sample_distance:MultiqcModule',
            'plot_gene_heatmap = multiqc_aladdin_rnaseq.modules.plot_gene_heatmap:MultiqcModule',
            'DESeq2 = multiqc_aladdin_rnaseq.modules.DESeq2:MultiqcModule',
            'gProfiler = multiqc_aladdin_rnaseq.modules.gProfiler:MultiqcModule',
            'trimming_2step = multiqc_aladdin_rnaseq.modules.trimming_2step:MultiqcModule'
        ],
        'multiqc.hooks.v1': [
            'before_config = multiqc_aladdin_rnaseq.custom_code:plugin_before_config',
            'execution_start = multiqc_aladdin_rnaseq.custom_code:plugin_execution_start'
        ]
    }
)

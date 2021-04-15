#!/usr/bin/env python
"""
Setup code for Zymo Research MultiQC plugin.
For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '1.0.0'

setup(
    name = 'multiqc_zymo',
    version = version,
    description = "MultiQC plugins for rnaseq pipeline",
    packages = find_packages(),
    include_package_data = True,
    install_requires = ['multiqc==1.9'],
    entry_points = {
        'multiqc.templates.v1': [
            'zymo = multiqc_zymo.templates.zymo'
        ],
        'multiqc.modules.v1': [
            'Trim_Galore = multiqc_zymo.modules.Trim_Galore:MultiqcModule',
            'plot_ERCC = multiqc_zymo.modules.plot_ERCC:MultiqcModule',
            'plot_sample_distance = multiqc_zymo.modules.plot_sample_distance:MultiqcModule',
            'plot_gene_heatmap = multiqc_zymo.modules.plot_gene_heatmap:MultiqcModule',
            'DESeq2 = multiqc_zymo.modules.DESeq2:MultiqcModule',
            'gProfiler = multiqc_zymo.modules.gProfiler:MultiqcModule',
            'download_data = multiqc_zymo.modules.download_data:MultiqcModule'
        ],
        'multiqc.hooks.v1': [
            'before_config = multiqc_zymo.custom_code:plugin_before_config',
            'execution_start = multiqc_zymo.custom_code:plugin_execution_start'
        ]
    }
)
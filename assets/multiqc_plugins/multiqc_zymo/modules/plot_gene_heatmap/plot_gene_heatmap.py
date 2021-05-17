#!/usr/bin/env python

"""MultiQC plugin module to plot a heatmap of gene expression.
Originally written for DESeq2 results, but can be modified for other heatmaps."""

from __future__ import print_function
import logging

from multiqc import config
from multiqc.plots import heatmap
from multiqc.modules.base_module import BaseMultiqcModule

import pandas as pd

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', False) is True:
            return None

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name = "Top gene expression patterns",
            target = "",
            anchor = "plot_gene_heatmap"
        )

        # Initialize the data dict
        self.plot_heatmap_data = dict()
        # Find the data files
        for f in self.find_log_files("plot_gene_heatmap", filehandles=True):
            self.add_data_source(f)
            # Parse the data file
            parsed_data = pd.read_csv(f["f"], sep="\t", index_col=0)
            if len(parsed_data):
                self.plot_heatmap_data[f["fn"]] = parsed_data
            else:
                log.debug("Could not parse heatmap data in {}".format(f["fn"]))
                raise UserWarning

        # Filter out samples matching ignored sample names
        self.plot_heatmap_data = self.ignore_samples(self.plot_heatmap_data)

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(self.plot_heatmap_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} heatmap reports".format(len(self.plot_heatmap_data)))

        # Write parsed report data to a file
        if len(self.plot_heatmap_data):
            self.write_data_file(self.plot_heatmap_data, "multiqc_plot_gene_heatmap")

        # Create heatmap
        for filename, data in self.plot_heatmap_data.items():
            pconfig = {
                "id": "gene_heatmap",
                "title": "Expression patterns of top genes",
                "xTitle": "Samples",
                "yTitle": "Genes",
                "square": False
            }
            heatmap_plot_html = heatmap.plot(data.values.tolist(), data.columns.tolist(), data.index.tolist(), pconfig)
            # Add a report section
            self.add_section(
                    description = ("Transformed read counts of top genes with highest variance, "
                                   "calculated using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). "
                                   "Read counts transformed using 'rlog' algorithm in DESeq2 and plotted in Log2 scale after centering per gene. "
                                   "A static version of this figure can be download on the results page on Aladdin platform."),
                    plot = heatmap_plot_html
            )

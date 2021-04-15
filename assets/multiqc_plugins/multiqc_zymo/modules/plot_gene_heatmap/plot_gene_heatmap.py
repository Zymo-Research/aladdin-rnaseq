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
                    "xTitle": "Samples",
                    "square": False
                    }
            if "isomirs" in filename.lower():
                pconfig["id"] = "miRNA_heatmap"
                pconfig["title"] = "Expression patterns of top miRNAs"
                pconfig["yTitle"] = "miRNAs"
                description = ("Normalized read counts of top miRNAs with highest variance, "
                               "calculated using [isomiRs](https://www.bioconductor.org/packages/release/bioc/html/isomiRs.html). "
                               "Values plotted in Log2 scale after centering per miRNA. "
                               "A static version of this figure can be "
                               " download in the [Download data section](#download_data).")
            elif "deseq2" in filename.lower():
                pconfig["id"] = "gene_heatmap"
                pconfig["title"] = "Expression patterns of top genes"
                pconfig["yTitle"] = "genes"
                description = ("Normalized read counts of top genes with highest variance, "
                               "calculated using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). "
                               "Values plotted in Log2 scale after centering per gene. "
                               "A static version of this figure can be "
                               "download in the [Download data section](#download_data).")
                
            heatmap_plot_html = heatmap.plot(data.values.tolist(), data.columns.tolist(), data.index.tolist(), pconfig)
            # Add a report section
            self.add_section(
                    description = description,
                    plot = heatmap_plot_html
                    )

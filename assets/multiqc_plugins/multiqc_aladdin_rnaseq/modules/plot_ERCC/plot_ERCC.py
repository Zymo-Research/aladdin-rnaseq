#!/usr/bin/env python

""" MultiQC plugin module to plot ERCC spike-in read counts and actual concentrations """

from __future__ import print_function
import logging

from multiqc import config
from multiqc.plots import scatter
from multiqc.modules.base_module import BaseMultiqcModule

from collections import OrderedDict, defaultdict
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
            name = "ERCC expression vs concentration",
            target = "",
            anchor = "plot_ERCC",
            info = "This section plots the known concentrations of ERCC RNA spike-in against "\
                    "their relative experssion levels (corrected by lengths)."
        )

        # Get whether to plot data points with zero counts from config
        plot_zero_count = getattr(config, "plot_ERCC_zero_count", True)

        # Initialize the data dict
        self.plot_ERCC_data = dict()
        # Find the data files
        for f in self.find_log_files("plot_ERCC", filehandles=True):
            self.add_data_source(f)
            # Parse the data file
            self.plot_ERCC_data = self.parse_ERCC_counts(f["f"], plot_zero_count)

        # Filter out samples matching ignored sample names
        self.plot_ERCC_data = self.ignore_samples(self.plot_ERCC_data)

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(self.plot_ERCC_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found ERCC spike-in data for {} samples.".format(len(self.plot_ERCC_data)))

        # Sort the data by sample name, so that tabs are consistent every time.
        self.plot_ERCC_data = OrderedDict(sorted(self.plot_ERCC_data.items()))

        # Write parsed report data to a file
        self.write_data_file(self.plot_ERCC_data, "multiqc_plot_ERCC")

        # Create scatter plot
        ylab = "Relative expression level (RPK or FPK)"
        pconfig = {
            "id": "ERCC_plot",
            "title": "Relative expression levels and concentrations of ERCC spike-in",
            "ylab": "Relative expression level (RPK or FPK)",
            "xlab": "Known concentration (attomoles/ul)",
            "xLog": True,
            "yLog": True,
            "data_labels": [{"name": sample, "ylab": ylab} for sample in self.plot_ERCC_data.keys()]
        }
        # Plot data must be a list of dicts
        plot_data = [{"ERCC transcripts":v} for v in self.plot_ERCC_data.values()]
        scatter_plot_html = scatter.plot(plot_data, pconfig)

        # Add a report section with the line plot
        helptext = "The relative expression levels is calcuated as:"
        if plot_zero_count:
            helptext += " (read count+1) / (transcript length/1000)"
        else:
            helptext += " read count/(transcript length/1000)"
        self.add_section(
            #description = "This plot shows how expression levels of ERCC spike-in transcripts correlate with their concentrations.",
            helptext = helptext,
            plot = scatter_plot_html
        )

    def parse_ERCC_counts(self, filehandle, plot_zero_count):
        ''' Parse the ERCC counts csv file '''

        data = pd.read_csv(filehandle, sep="\t", index_col=0)
        parsed_data = defaultdict(list)
        for idx, row in data.iterrows():
            for col in row.keys():
                if col not in ["concentration", "length"]:
                    if float(row[col])+int(plot_zero_count) > 0:
                        parsed_data[col].append({
                            "x": float(row["concentration"]),
                            "y": (float(row[col])+int(plot_zero_count)) * 1000 / int(row["length"]),
                            "name": idx
                        })
        return parsed_data
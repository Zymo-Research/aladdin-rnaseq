#!/usr/bin/env python

""" MultiQC plugin module to plot ERCC spike-in read counts and actual concentrations """

from __future__ import print_function
import logging

from multiqc import config
from multiqc.plots import scatter
from multiqc.modules.base_module import BaseMultiqcModule

import csv
from collections import OrderedDict

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
            parsed_data = self.parse_ERCC_counts(f["f"], plot_zero_count)
            if parsed_data:
                self.plot_ERCC_data[f["s_name"]] = parsed_data
            else:
                log.debug("Could not find any ERCC counts in {}".format(f["fn"]))
                raise UserWarning

        # Filter out samples matching ignored sample names
        self.plot_ERCC_data = self.ignore_samples(self.plot_ERCC_data)

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(self.plot_ERCC_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.plot_ERCC_data)))

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
        scatter_plot_html = scatter.plot(list(self.plot_ERCC_data.values()), pconfig)

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

        data = list(csv.DictReader(filehandle, delimiter="\t"))
        parsed_data = dict()
        transcripts = []
        for row in data:
            if float(row["count"])+int(plot_zero_count) > 0:
                transcripts.append({
                    "x": float(row["concentration"]),
                    "y": (float(row["count"])+int(plot_zero_count)) * 1000 / int(row["length"]),
                    "name": row["ID"]
                    })
        if len(transcripts):
            parsed_data["ERCC transcripts"] = transcripts

        return parsed_data
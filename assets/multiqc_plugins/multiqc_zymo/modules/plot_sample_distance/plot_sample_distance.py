#!/usr/bin/env python

"""MultiQC plugin module to plot distances/similarities between samples.
Originally written for RNA-Seq, but should apply to most NGS applications.
Given a distance matrix and/or a PCA data matrix, make heatmap and/or PCA plot."""

from __future__ import print_function
import logging

from multiqc import config
from multiqc.plots import scatter,heatmap
from multiqc.modules.base_module import BaseMultiqcModule

import pandas as pd
from collections import defaultdict

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', False) is True:
            return None

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name = "Distances/similarities between samples",
            target = "",
            anchor = "plot_sample_distance",
            info = "This section plots the distances or similarities between samples in the form of "
                    "heatmap, PCA, and/or MDS plots."
        )

        # Initialize the data dict
        self.plot_heatmap_data = dict()
        # Find the data files
        for f in self.find_log_files("plot_sample_distance/heatmap", filehandles=True):
            self.add_data_source(f)
            # Parse the data file
            parsed_data = pd.read_csv(f["f"], sep="\t", index_col=0)
            if len(parsed_data):
                self.plot_heatmap_data[f["fn"]] = parsed_data
            else:
                log.debug("Could not parse sample distance/similarity matrix data in {}".format(f["fn"]))
                raise UserWarning

        self.plot_pca_data = dict()
        for f in self.find_log_files("plot_sample_distance/pca", filehandles=True):
            parsed_data = self.parse_pca(f["f"])
            if parsed_data:
                self.plot_pca_data[f["fn"]] = parsed_data
            else:
                log.debug("Could not parse MDS/PCA plot data in {}".format(f["fn"]))
                raise UserWarning

        # Filter out samples matching ignored sample names
        self.plot_heatmap_data = self.ignore_samples(self.plot_heatmap_data)
        self.plot_pca_data = self.ignore_samples(self.plot_pca_data)

        # Nothing found - raise a UserWarning to tell MultiQC
        if (len(self.plot_heatmap_data) == 0) and (len(self.plot_pca_data) == 0):
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} heatmap reports".format(len(self.plot_heatmap_data)))
        log.info("Found {} PCA reports".format(len(self.plot_pca_data)))

        # Write parsed report data to a file
        if len(self.plot_heatmap_data):
            self.write_data_file(self.plot_heatmap_data, "multiqc_plot_sample_distance_heatmap")
        if len(self.plot_pca_data):
            self.write_data_file(self.plot_pca_data, "multiqc_plot_sample_distance_PCA")

        # Create heatmap
        for filename, data in self.plot_heatmap_data.items():
            pconfig = {
                    "xTitle": "Samples",
                    "yTitle": "Samples"
                    }
            if "sample_distance_matrix" in filename.lower():
                pconfig["id"] = "sample_distance_matrix"
                pconfig["title"] = "Sample distance matrix"
                section_name = "Distance matrix of samples"
                section_anchor = "sample_distance_heatmap"
                section_description = ("The distances (Euclidean distance) between samples "
                                       "are visualized here in the form of heatmap. Smaller "
                                       "values indicate higher similarity between samples.")
            elif "sample_similarity_matrix" in filename.lower():
                pconfig["id"] = "sample_similarity_matrix"
                pconfig["title"] = "Sample similarity matrix"
                section_name = "Similarity matrix of samples"
                section_anchor = "sample_similarity_heatmap"
                section_description = ("The similarities (Pearson correlation coefficient) "
                                       "between samples are visualized here in the form of "
                                       "heatmap. Larger values indicate higher similarity "
                                       "between samples.")
            if "isomirs" in filename.lower():
                section_description += (" The similarities were calculated using Log2 values "
                                        "of normalized read counts of all mature miRNAs using "
                                        "[isomiRs](https://www.bioconductor.org/packages/release/bioc/html/isomiRs.html).")
            elif "deseq" in filename.lower():
                section_description += (" The similarities were calculated using normalized and "
                                        "'rlog' transformed read counts of all genes using "
                                        "[DESeq2](https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html).") 
            
            heatmap_plot_html = heatmap.plot(data.values.tolist(), data.columns.tolist(), data.index.tolist(), pconfig)
            # Add a report section
            self.add_section(
                    name = section_name,
                    anchor = section_anchor,
                    description = section_description,
                    plot = heatmap_plot_html
                    )

        # Create PCA/MDS plot
        for filename, data in self.plot_pca_data.items():
            pconfig = {
                    "xlab": data["xlab"],
                    "ylab": data["ylab"],
                    "colors": data["colors"],
                    "marker_line_colour": data["colors"],
                    "square": True
                    }
            if 'pca_plot' in filename.lower():
                pconfig["id"] = "sample_pca_plot"
                pconfig["title"] = "PCA plot of samples"
                section_name = "Principal component analysis of samples"
                section_anchor = "sample_PCA_plot"
                section_description = ("[Principal component analysis](https://en.wikipedia.org/wiki/Principal_component_analysis) "
                    "was conducted to visualize the distance/similarity between samples.")
            elif 'mds_plot' in filename.lower():
                pconfig["id"] = "sample_mds_plot"
                pconfig["title"] = "MDS plot of samples"
                section_name = "Multidimensional scaling analysis of samples"
                section_anchor = "sample_MDS_plot"
                section_description = ("[Multidimensional scaling](https://en.wikipedia.org/wiki/Multidimensional_scaling) "
                    "was conducted to visualize the distance/similarity between samples.")
            if "isomirs" in filename.lower():
                section_description += " Top 500 miRNAs with highest variance among samples were used to make this plot."
            elif "deseq" in filename.lower():
                section_description += " Top 500 genes with highest variance among samples were used to make this plot."
            scatter_plot_html = scatter.plot(data["data"], pconfig)
            self.add_section(
                    name = section_name,
                    anchor = section_anchor,
                    description = section_description,
                    plot = scatter_plot_html
                    )

    def parse_pca(self, filehandle):
        ''' Parse the PCA/MDS plot data file '''

        parsed_data = dict()
        data = pd.read_csv(filehandle, sep="\t", index_col=0)
        if "group" in data.columns.str.lower():
            group_col = "group"
        elif "condition" in data.columns.str.lower():
            group_col = "condition"
        else:
            log.debug("'condition' or 'group' column NOT FOUND in PCA/MDS plot data file.")
            raise UserWarning
            return None
        # PCA axis labels and data are stored in the first two columns
        parsed_data["xlab"] = data.columns[0]
        parsed_data["ylab"] = data.columns[1]
        # Default colors
        default_colors = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9',
                          '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1']
        parsed_data["colors"] = dict()
        # If there is no group information, return original data
        if data[group_col].isna().all():
            data = data.iloc[:,[0,1]]
            data.columns = ['x','y']
            parsed_data["data"] = data.to_dict('index')
        else:
        # If there is group information, group data points by group
            datapoints = defaultdict(list)
            for index, row in data.iterrows():
                datapoints[row[group_col]].append({
                    "x": row[0],
                    "y": row[1],
                    "name": index
                    })
            parsed_data["data"] = datapoints
            for idx, group in enumerate(datapoints.keys()):
                cidx = idx
                while cidx >= len(default_colors):
                    cidx -= len(default_colors)
                parsed_data["colors"][group] = default_colors[cidx]
        # return results
        return parsed_data

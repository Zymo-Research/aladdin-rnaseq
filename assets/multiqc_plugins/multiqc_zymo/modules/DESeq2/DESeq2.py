#!/usr/bin/env python

"""MultiQC plugin module to process DESeq2 differential gene expression results.
It outputs a table tallying the numbers of up- and down-regulated genes and 
generates MA-plots for all combinations of conditions."""

from __future__ import print_function
import logging

from multiqc import config
from multiqc.plots import scatter,table
from multiqc.modules.base_module import BaseMultiqcModule

import pandas as pd
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
            name = "Differential gene expression",
            target = "DESeq2",
            href = "https://bioconductor.org/packages/release/bioc/html/DESeq2.html",
            anchor = "DESeq2",
            info = (" calculates the expression levels of genes and conducts "
            "statistical analysis of differential gene expression.")
        )

        # Get the pvalue cutoff from config
        alpha = getattr(config, "DESeq2_alpha", 0.1)
        link_prefix = getattr(config, "gene_link_prefix", None)
        if link_prefix is not None:
            link_prefix = '<a href=\"https://{}/id/'.format(link_prefix)

        # Initialize the data dict
        self.deseq_results = dict()
        gene_label = "Gene Name"
        # Find the data files
        for f in self.find_log_files("DESeq2", filehandles=True):
            self.add_data_source(f)
            # Parse the data file
            parsed_data = pd.read_csv(f["f"], sep="\t")
            # Sanity check
            cols = ["baseMean", "log2FoldChange", "padj", "lfcShrink"]
            for col in cols:
                if col not in parsed_data.columns:
                    log.debug("{} column missing in DESeq2 results.".format(col))
                    parsed_data = None
                    break
            if "Gene Name" not in parsed_data.columns:
                gene_label = "Gene ID"
            if parsed_data is not None:
                self.deseq_results[f["s_name"]] = parsed_data
            else:
                log.debug("Could not parse DESeq2 results in {}".format(f["fn"]))
                raise UserWarning

        # Filter out samples matching ignored sample names
        self.deseq_results = self.ignore_samples(self.deseq_results)

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(self.deseq_results) ==  0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} DESeq2 reports".format(len(self.deseq_results)))

        # Sort the data by sample name, so that tabs are consistent every time.
        self.deseq_results = OrderedDict(sorted(self.deseq_results.items()))

        # Write parsed report data to a file
        if len(self.deseq_results):
            self.write_data_file(self.deseq_results, "multiqc_plot_DESeq2_results")

        # Create table
        table_data = OrderedDict()
        # Count no. of genes in each category
        for sample_name, data in self.deseq_results.items():
            table_data[sample_name] = {
                    "filtered": data["padj"].isna().sum(),
                    "not_deg": (data["padj"]>=alpha).sum(),
                    "up_deg": ((data["padj"]<alpha) & (data["log2FoldChange"]>0)).sum(),
                    "down_deg": ((data["padj"]<alpha) & (data["log2FoldChange"]<0)).sum()
                    }
        # Set table headers
        headers = OrderedDict()
        headers["up_deg"] = {
                "title": "Higher in condition 1",
                "description": "Numbers of differentially expressed genes with higher expression in condition 1"
                }
        headers["down_deg"] = {
                "title": "Higher in condition 2",
                "description": "Numbers of differentially expressed genes with higher expression in condition 2"
                }
        headers["not_deg"] = {
                "title": "Not differentially expressed",
                "description": "Numbers of genes that were not differentially expressed between conditions 1 and 2"
                }
        headers["filtered"] = {
                "title": "Did not pass filter",
                "description": "Numbers of genes that failed to pass filter in DESeq2 because low expression level or outliers"
                }
        table_config = {
                "id": "Differential_Expression_Table",
                "col1_header": "Comparison(cond.1_cond.2)",
                "format": "{:,.0f}", # No decimal
                "sortRows": False
                }
        table_plot_html = table.plot(table_data, headers, table_config)
        self.add_section(
                name = "Summary table of gene differential expression",
                anchor = "differential_expression_table",
                description = ("General statistics of differentially expressed "
                               "genes in pairwise comparisons. Genes with adjusted "
                               "p-values smaller than {} were considered "
                               "differentially expressed.").format(alpha),
                plot = table_plot_html
                )
        
        # Create scatter plots
        xy_plot_data = OrderedDict()
        for sample_name, data in self.deseq_results.items():
            cond1, cond2 = sample_name.split("_vs_")
            cols = [cond1, cond2, gene_label]
            deg = data.loc[data["padj"]<alpha, cols]
            # Only display the first 1000 genes
            deg = deg.iloc[:1000,]
            non_deg = data.loc[data["padj"]>=alpha, cols]
            # Downsample the non-DEGs so display won't be a problem
            if len(non_deg) > 1000:
                non_deg = non_deg.sample(1000)
            # Rename the columns for plotting
            col_rename = {
                    cond1: "x",
                    cond2: "y",
                    gene_label: "name"
                    }
            deg = deg.rename(columns = col_rename)
            non_deg = non_deg.rename(columns = col_rename)
            # Convert DataFrame to list of dicts
            plot_data = {}
            plot_data["DEG"] = deg.to_dict("records")
            plot_data["nonDEG"] = non_deg.to_dict("records")
            xy_plot_data[sample_name] = plot_data
        # Set the figure config
        pconfig = {
                "id": "Scatter_plot",
                "xlab": "Normalized mean counts in {}".format(cond1),
                "ylab": "Nomarlized mean counts in {}".format(cond2),
                "xLog": True,
                "yLog": True,
                "colors": { "DEG":"#f45b5b", "nonDEG":"#434348" },
                "marker_line_colour": { "DEG":"#f45b5b", "nonDEG":"#434348" },
                "square": True,
                "marker_size": 3,
                "data_labels": [{"name": sample} for sample in xy_plot_data.keys()]
                }
        scatter_plot_html = scatter.plot(list(xy_plot_data.values()), pconfig)
        self.add_section(
                name = "Scatter plot",
                anchor = "scatter_plot",
                description = ("Scatter plot is a simple and straightforward way"
                " to visualize differential gene expression results."
                " Expression levels of genes in one condition are shown on X-axis"
                " while those in other are shown on Y-axis.<br>"
                "The scatter plots here include differentially expressed genes"
                " (up to the first 1000 genes) and up to 1000 randomly selected"
                " non-differentially expressed genes. You can download the"
                " scatter plots with all genes in the [Download data section](#download_data).<br>"
                "Red dots represent differentially expressed genes (adjusted p-values<{}). "
                "Grey dots represent non-differentially expressed genes. "
                "Count data transformation was carried out using the 'rlog' method in DESeq2").format(alpha),
                plot = scatter_plot_html
                )
                
        # Create MA plots
        ma_plot_data = OrderedDict()
        cols = ["baseMean", "lfcShrink", gene_label]
        for sample_name, data in self.deseq_results.items():
            # Separate DEGs and non-DEGs
            deg = data.loc[data["padj"]<alpha, cols]
            # Only display the first 1000 genes
            deg = deg.iloc[:1000,]
            non_deg = data.loc[data["padj"]>=alpha, cols]
            # Downsample the non-DEGs so display won't be a problem
            if len(non_deg) > 1000:
                non_deg = non_deg.sample(1000)
            # Rename the columns for plotting
            col_rename = {
                    "baseMean": "x",
                    "lfcShrink": "y",
                    gene_label: "name"
                    }
            deg = deg.rename(columns = col_rename)
            non_deg = non_deg.rename(columns = col_rename)
            # Convert DataFrame to list of dicts
            plot_data = {}
            plot_data["DEG"] = deg.to_dict("records")
            plot_data["nonDEG"] = non_deg.to_dict("records")
            ma_plot_data[sample_name] = plot_data
        # Set the figure config
        pconfig = {
                "id": "MA_plot",
                "xlab": "Normalized mean read counts",
                "ylab": "Log2(Fold Change)",
                "xLog": True,
                "colors": { "DEG":"#f45b5b", "nonDEG":"#434348" },
                "marker_line_colour": { "DEG":"#f45b5b", "nonDEG":"#434348" },
                "square": True,
                "marker_size": 3,
                "data_labels": [{"name": sample} for sample in ma_plot_data.keys()]
                }
        scatter_plot_html = scatter.plot(list(ma_plot_data.values()), pconfig)
        self.add_section(
                name = "MA plot",
                anchor = "ma_plot",
                description = ("[MA plot](https://en.wikipedia.org/wiki/MA_plot)"
                " is a type of visualization of differential gene expression results"
                " often used in publications. Expression levels are shown on X-axis"
                " while log2 of fold changes are shown on Y-axis.<br>"
                "The MA plots here include differentially expressed genes (up to"
                " the first 1000 genes) and up to 1000 randomly selected non-differentially"
                " expressed genes. You can download the MA plots with all genes"
                " in the [Download data section](#download_data).<br>"
                "Red dots represent differentially expressed genes (adjusted p-values<{}). "
                "Grey dots represent non-differentially expressed genes. "
                "Shrinkage of effect size was carried out using the 'normal' method in DESeq2").format(alpha),
                plot = scatter_plot_html
                )
        
        # Create snapshot tables for each comparison
        headers = OrderedDict()
        headers["gene_link"] = {
                "title": "Gene name",
                "description": "Gene name"
                }
        headers["baseMean"] = {
                "title": "Mean counts",
                "description": "Mean read counts in both conditions",
                "format": "{:,.0f}"
                }
        headers["log2FoldChange"] = {
                "title": "Log2 Fold change",
                "description": "Log2 Fold change",
                "format": "{:,.2f}"
                }
        headers["padj"] = {
                "title": "False discovery rate",
                "description": "False discovery rate",
                "format": "{:,.1e}",
                "min": 0,
                "max": 1
                }
        table_config = {
                "id": "top_deg_table",
                "col1_header": "Rank",
                "sortRows": False,
                "scale": False
                }
        for sample_name, data in self.deseq_results.items():
            cond1, cond2 = sample_name.split("_vs_")
            # Sort by FDR and get top 50 genes
            data = data.sort_values("padj")
            genes = data.copy().head(50)
            # Add link to the gene names
            if link_prefix is not None:
                genes["gene_link"] = genes.apply(lambda x:link_prefix+x["Gene ID"]+'" target="_blank">'+x[gene_label]+"</a>", axis=1)
            else:
                genes["gene_link"] = genes[gene_label]
            genes = genes[["gene_link", "baseMean", "log2FoldChange", "padj"]]
            genes = genes.to_dict("records")
            # Make the top gene set table
            table_data = OrderedDict()
            for idx, row in enumerate(genes):
                table_data[str(idx+1)] = row
            table_config["id"] = "Top_DEGs_"+sample_name
            table_plot_html = table.plot(table_data, headers, table_config)
            self.add_section(
                    name = "Top differentially expressed genes in comparison {} vs. {}".format(cond1, cond2),
                    anchor = "top_degs_"+sample_name,
                    description = ("Top 50 differentially expressed genes, "
                                   "ranked by FDR, in comparison {} vs. {}. "
                                   "Full DESeq2 results can be downloaded in "
                                   "the [Download data section](#download_data).<br>"
                                   "Genes with positive Log2 fold changes have "
                                   "higher expression in {}. Genes with negative "
                                   "Log2 fold changes have higher expression in "
                                   "{}.").format(cond1, cond2, cond1, cond2),
                    plot = table_plot_html
                    )
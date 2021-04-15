#!/usr/bin/env python

"""MultiQC plugin module to process gProfiler gene set enrichment results.
It outputs a table tallying the numbers of up- and down-regulated gene sets and 
generates tables of top 20 gene sets in each comparison."""

from __future__ import print_function
import logging

from multiqc import config
from multiqc.plots import table
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
            name = "Gene set enrichment analysis",
            target = "g:Profiler",
            href = "https://biit.cs.ut.ee/gprofiler/gost",
            anchor = "gProfiler",
            info = (" performs functional enrichment analysis also known as "
                    "gene set enrichment analysis on input gene list. "
                    "It maps genes to known functional information sources "
                    "and detects statistically significantly enriched terms.")
        )

        # Get the pvalue cutoff from config
        alpha = getattr(config, "gProfiler_alpha", 0.05)

        # Initialize the data dict
        self.gProfiler_results = dict()
        # Find the data files
        for f in self.find_log_files("gProfiler", filehandles=True):
            self.add_data_source(f)
            # Parse the data file
            parsed_data = pd.read_csv(f["f"], sep="\t")
            # Sanity check
            cols = ["ID", "name", "p_value", "expression pattern"]
            for col in cols:
                if col not in parsed_data.columns:
                    log.debug("{} column missing in gProfiler results.".format(col))
                    parsed_data = None
                    break
            if parsed_data is not None:
                self.gProfiler_results[f["s_name"]] = parsed_data[cols]
            else:
                log.debug("Could not parse gProfiler results in {}".format(f["fn"]))
                raise UserWarning

        # Filter out samples matching ignored sample names
        self.gProfiler_results = self.ignore_samples(self.gProfiler_results)

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(self.gProfiler_results) ==  0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} gProfiler reports".format(len(self.gProfiler_results)))

        # Sort the data by sample name, so that tabs are consistent every time.
        self.gProfiler_results = OrderedDict(sorted(self.gProfiler_results.items()))

        # Write parsed report data to a file
        if len(self.gProfiler_results):
            self.write_data_file(self.gProfiler_results, "multiqc_gProfiler_results")

        # Create summary table
        table_data = OrderedDict()
        # Count no. of genes in each category
        for sample_name, data in self.gProfiler_results.items():
            cond1, cond2 = sample_name.split("_vs_")
            cond1 = "Higher in " + cond1
            cond2 = "Higher in " + cond2
            table_data[sample_name] = {
                    "not_enriched": (data["p_value"]>=alpha).sum(),
                    "up_gs": ((data["p_value"]<alpha) & (data["expression pattern"]==cond1)).sum(),
                    "down_gs": ((data["p_value"]<alpha) & (data["expression pattern"]==cond2)).sum()
                    }
        # Set table headers
        headers = OrderedDict()
        headers["up_gs"] = {
                "title": "Higher in condition 1",
                "description": "Numbers of gene sets enriched with higher expression in condition 1"
                }
        headers["down_gs"] = {
                "title": "Higher in condition 2",
                "description": "Numbers of gene sets enriched with higher expression in condition 2"
                }
        headers["not_enriched"] = {
                "title": "Not enriched",
                "description": ("Numbers of gene sets that were not enriched "
                                "with higher expression in either conditions")
                }
        table_config = {
                "id": "gene_set_enrichment_table",
                "col1_header": "Comparison(cond.1_cond.2)",
                "format": "{:,.0f}", # No decimal
                "sortRows": False
                }
        table_plot_html = table.plot(table_data, headers, table_config)
        self.add_section(
                name = "Summary table of gene set enrichment analysis",
                anchor = "gProfiler_summary_table",
                description = ("General statistics of gene set enrichment "
                               "analysis in pairwise comparisons. Gene sets "
                               "with false discovery rate smaller than {} "
                               "were considered enriched.").format(alpha),
                plot = table_plot_html
                )
        
        # Create snapshot tables for each comparison
        # Reset the table headers
        headers = OrderedDict()
        headers["name"] = {
                "title": "Gene set name",
                "description": "Gene set name",
                }
        headers["ID"] = {
                "title": "Gene set category",
                "description": "The source of gene sets"
                }
        headers["p_value"] = {
                "title": "Adjusted p-value",
                "description": "Adjusted p-value"
                }
        headers["expression pattern"] = {
                "title": "Expression pattern",
                "description": "Expression pattern of the gene set"
                }
        table_config = {
                "id": "top_gene_sets_table",
                "col1_header": "Rank",
                "format": "{:,.3f}",
                "sortRows": False,
                "scale": False
                }
        for sample_name, data in self.gProfiler_results.items():
            cond1, cond2 = sample_name.split("_vs_")
            if len(data)>0:
                # Get top 30 gene sets
                gene_set = data.copy().head(30)
                # Add link to the gene set names
                gene_set["name"] = gene_set[["ID", "name"]].apply(self.find_address, axis=1)
                gene_set["ID"] = gene_set["ID"].map(self.parse_gene_set_category)
                gene_set = gene_set.to_dict("records")
                # Make the top gene set table
                table_data = OrderedDict()
                for idx, row in enumerate(gene_set):
                    table_data[str(idx+1)] = row
                table_config["id"] = "Top_Gene_Sets_"+sample_name
                table_plot_html = table.plot(table_data, headers, table_config)
                self.add_section(
                        name = "Top enriched gene sets in comparison {} vs. {}".format(cond1, cond2),
                        anchor = "top_gene_sets_"+sample_name,
                        description = ("Top 30 gene sets, ranked by p-value, in"
                                       " comparison {} vs. {}. Full g:Profiler "
                                       "results can be downloaded in the "
                                       "[Download data section](#download_data)."
                                      ).format(cond1, cond2),
                        plot = table_plot_html
                        )
            else:
                self.add_section(
                        name = "Top enriched gene sets in comparison {} vs. {}".format(cond1, cond2),
                        anchor = "top_gene_sets_"+sample_name,
                        description = ("g:Profiler analysis was not carried out"
                                       " because no gene belonging to the"
                                       " tested gene set(s) were differentially"
                                       " expressed in comparison {} vs. {}").format(cond1, cond2)
                        )                    
    
    def find_address(self, row):
        gene_set_id, name = row
        if gene_set_id.startswith("GO:"):
            address = "http://amigo.geneontology.org/amigo/term/"+gene_set_id
        elif gene_set_id.startswith("REAC:"):
            gene_set_id = gene_set_id.replace("REAC:", "")
            address = "https://reactome.org/content/detail/"+gene_set_id
        elif gene_set_id.startswith("KEGG:"):
            gene_set_id = gene_set_id.replace("KEGG:", "")
            address = "https://www.genome.jp/dbget-bin/www_bget?pathway:map"+gene_set_id
        return "<a href='{}' target='_blank'>{}</a>".format(address, name)
    
    def parse_gene_set_category(self, gene_set_id):
        if gene_set_id.startswith("GO:"):
            return "GO:Biological Process"
        elif gene_set_id.startswith("REAC:"):
            return "Reactome Pathway"
        elif gene_set_id.startswith("KEGG:"):
            return "KEGG Pathway"

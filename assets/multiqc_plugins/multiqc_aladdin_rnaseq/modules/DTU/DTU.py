#!/usr/bin/env python

"""MultiQC module to display differential transcript usage (DTU) results.
It plots a scatterplot showing which genes have DTU or DGE or both.
It also plots a table of the most significant DTU genes and transcripts.
"""

from __future__ import print_function
import logging

from multiqc import config
from multiqc.plots import scatter,table
from multiqc.modules.base_module import BaseMultiqcModule

import pandas as pd
import numpy as np
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
            name = "Differential transcript usage",
            target = "",
            anchor = "DTU",
            info = ("This section display the results of differential transcript usage (DTU) analysis. "
            "Methods of this analysis were adapted from <a href=\"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6178912\" target=\"_blank\">this paper</a>. "
            "The analysis consists of two stages with their own adjusted p-values: "
            "a screening stage that asks if a gene has DTU; and a confirmation stage that checks transcripts for DTU if their gene passed screening stage. "
            "Please see the paper for details and code.")
        )

        # Get the pvalue cutoff from config
        deseq_alpha = getattr(config, "DESeq2_alpha", 0.05)
        dexseq_alpha = getattr(config, "DEXSeq_alpha", 0.05)
        link_prefix = getattr(config, "gene_link_prefix", None)
        if link_prefix is not None:
            link_prefix = '<a href=\"https://{}/id/'.format(link_prefix)
        
        # Find and load DESeq2 and DEXSeq results
        self.deseq_results = dict()
        self.dexseq_results = dict()
        for f in self.find_log_files('DTU/DESeq2', filehandles=True):
            self.parse_deseq_results(f)
        for f in self.find_log_files('DTU/DEXSeq', filehandles=True):
            self.parse_dexseq_results(f)
        if self.dexseq_results is None:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        # Write parsed report data to a file
        self.write_data_file(self.dexseq_results, 'multiqc_DTU')

        # Summary table for DTU
        if len(self.dexseq_results):
            self.dtu_summary(dexseq_alpha)

        # Scatterplot of DEG and DTU genes
        if len(self.deseq_results) and len(self.dexseq_results):
            self.gene_scatterplot(-np.log10(deseq_alpha), -np.log10(dexseq_alpha))

        # Table of top DTU genes
        self.top_dtu_table(dexseq_alpha, link_prefix)

    def parse_deseq_results(self, f):
        """Parse DESeq2 results to get DEG genes and its adjusted p-values"""
        self.add_data_source(f)
        parsed_data = pd.read_csv(f["f"], sep="\t")
        try:
            if 'Gene Name' in parsed_data.columns:
                self.deseq_results[f["s_name"]] = parsed_data[['Gene ID', 'Gene Name', 'padj']]
            else:
                self.deseq_results[f["s_name"]] = parsed_data[['Gene ID', 'padj']]
        except:
            log.error("Failed to parse DESeq2 results in {}".format(f['fn']))

    def parse_dexseq_results(self, f):
        """Parse DEXSeq results"""
        self.add_data_source(f)
        parsed_data = pd.read_csv(f["f"], sep="\t")
        cols = ['Gene ID', 'Transcript ID', 'screening_padj','confirmation_padj']
        for col in cols:
            if col not in parsed_data.columns:
                log.error("Failed to parse DEXSeq results in {} because missing {} column".format(f['fn'], col))
        self.dexseq_results[f["s_name"]] = parsed_data

    def dtu_summary(self, dexseq_alpha):
        """Make a summary table for DTU analysis"""
        table_data = OrderedDict()
        for sample_name, dexseq_results in self.dexseq_results.items():
            dtu_gene_padj = dexseq_results[['Gene ID', 'screening_padj']].drop_duplicates()
            table_data[sample_name] = {
                "not_dtu": (dtu_gene_padj["screening_padj"]>=dexseq_alpha).sum(),
                "dtu": (dtu_gene_padj["screening_padj"]<dexseq_alpha).sum(),
            }
        # Configure the table
        headers = OrderedDict()
        headers["dtu"] = {
            "title": "Genes involved in DTU",
            "description": "Numbers of genes involved in differential transcript usage"
        }
        headers["not_dtu"] = {
            "title": "Genes not involved in DTU",
            "description": "Numbers of genes not involved in differential transcript usage"
        }
        table_config = {
            "id": "Differential_Transcript_Usage_Table",
            "col1_header": "Comparison",
            "format": "{:,.0f}", # No decimal
            "sortRows": False
        }
        table_plot_html = table.plot(table_data, headers, table_config)
        self.add_section(
            name = "Summary table of differential transcript usage",
            anchor = "dtu_summary_table",
            description = ("Summary of genes involved in differential transcript usage (DTU) or not in pairwise comparisons. "
                           "Not all genes are included in this table. Pre-filtering of transcripts were carried out to remove "
                           "transcripts and genes with low read counts or low proportion of expression within the gene. "
                           "Statistical analysis were not carried out on those genes/transcripts. "
                           "Please see [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6178912/) for filtering criteria."),
            plot = table_plot_html
        )

    def gene_scatterplot(self, deseq_alpha, dexseq_alpha):
        """Make a scatterplot of log scaled adjusted p-values for DEG and DTU analysis for each gene"""
        # Function to transform to -log10 scale, with capping and jitter
        def transform(x, min_val, jitter_stdev):
            if x<=min_val:
                return -np.log10(min_val) + np.random.normal(0, jitter_stdev)
            else:
                return -np.log10(x)
        scatterplot_data = OrderedDict()
        # Collect DESeq2 and DEXSeq data, combine and parse
        for sample_name, dexseq_results in self.dexseq_results.items():
            group1, group2 = sample_name.split("_vs_")
            if sample_name in self.deseq_results:
                deseq_results = self.deseq_results[sample_name]
            elif "{}_vs_{}".format(group2, group1) in self.deseq_results:
                deseq_results = self.deseq_results["{}_vs_{}".format(group2, group1)]
            else:
                log.debug("DEXSeq comparison {} has no corresponding DESeq2 results. Skipping...".format(sample_name))
                continue
            # Combine both DGE and DTU results
            dtu_gene_padj = dexseq_results[['Gene ID', 'screening_padj']].drop_duplicates()
            combined_padjs = pd.merge(deseq_results, dtu_gene_padj, on='Gene ID', how='outer').fillna(1)
            combined_padjs['padj'] = combined_padjs['padj'].apply(transform, min_val=1e-20, jitter_stdev=0.25)
            combined_padjs['screening_padj'] = combined_padjs['screening_padj'].apply(transform, min_val=1e-20, jitter_stdev=0.25)
            # Rename columns for plotting
            if 'Gene Name' in combined_padjs.columns:
                combined_padjs = combined_padjs.drop('Gene ID', axis=1)
                combined_padjs = combined_padjs.rename(columns={'Gene Name':'name', 'padj':'x', 'screening_padj':'y'})
            else:
                combined_padjs = combined_padjs.rename(columns={'Gene ID':'name', 'padj':'x', 'screening_padj':'y'})
            # Separate categories, only show 1000 genes for each
            dtu_only = combined_padjs.loc[(combined_padjs['y']>dexseq_alpha)&(combined_padjs['x']<=deseq_alpha), ].sort_values('y',ascending=False).iloc[:500,]
            dge_only = combined_padjs.loc[(combined_padjs['y']<=dexseq_alpha)&(combined_padjs['x']>deseq_alpha), ].sort_values('x',ascending=False).iloc[:500,]
            both = combined_padjs.loc[(combined_padjs['y']>dexseq_alpha)&(combined_padjs['x']>deseq_alpha), ].sort_values(['y','x'],ascending=False).iloc[:500,]
            neither = combined_padjs.loc[(combined_padjs['y']<=dexseq_alpha)&(combined_padjs['x']<=deseq_alpha), ]
            if len(neither) > 100:
                neither = neither.sample(100)
            # Prepare plot data
            plot_data = dict()
            plot_data['DTU'] = dtu_only.to_dict('records')
            plot_data['DGE'] = dge_only.to_dict('records')
            plot_data['Both'] = both.to_dict('records')
            plot_data['Neither'] = neither.to_dict('records')
            scatterplot_data[sample_name] = plot_data
        # Make the scatterplot
        pconfig = {
            "id" : "DGE_DTU_genes_scatterplot",
            "xlab" : "-log10 adjusted p-value for DGE",
            "ylab" : "-log10 adjusted p-value for DTU",
            "colors" : { 'DTU':'#7cb5ec', 'DGE':'#90ed7d', 'Both':'#f45b5b', 'Neither':'#434348' },
            "marker_line_colour" : { 'DTU':'#7cb5ec', 'DGE':'#90ed7d', 'Both':'#f45b5b', 'Neither':'#434348' },
            "marker_size" : 3,
            "square": True,
            "data_labels" : [{'name':k, "xlab":"-log10 adjusted p-value for DGE", "ylab":"-log10 adjusted p-value for DTU"} for k in scatterplot_data.keys()]
        }
        scatter_plot_html = scatter.plot(list(scatterplot_data.values()), pconfig)
        self.add_section(
            name = "DGE and DTU of genes",
            anchor = "dge_dtu_genes_scatterplot",
            description = ('This plot summarizes the results of both differential gene expression(DGE) '
            'and differential transcript usage(DTU) in each pairwise comparison. Genes in blue are genes involved only in DTU, '
            'i.e. the total expression level of a gene did not change, but proportions of transcripts changed. '
            'Green are those involved only in DGE, i.e. the proportions of transcripts did not change but there were '
            'an increase or decrease on the gene level. Red are those involved in both. Grey are those involved in neither.'
            'Genes plotted for DGE/DTU/both are capped at the top 500 genes and 100 random selected genes involved in neither are plotted. '
            'X,Y axes plot -log10 adjusted p-values of DGE and DTU analysis, respectively. '
            'Plotted p-values are capped at around 1e-20. Please see the results files for actual p-values.'),
            plot = scatter_plot_html
        )
        
    def top_dtu_table(self, dexseq_alpha, link_prefix):
        """Make a talbe listing the top DTU genes and transcripts"""
        # Configure the table
        headers = OrderedDict()
        headers["gene_link"] = {
                "title": "Gene name",
                "description": "Gene name",
                "scale": False
                }
        headers["transcript_link"] = {
                "title": "Transcript ID",
                "description": "Transcript ID",
                "scale": False
                }
        headers["screening_padj"] = {
                "title": "Screening adjusted p-value",
                "description": "Adjusted p-values for the screening stage, i.e. whether a gene is involved in DTU",
                "format": "{:,.1e}",
                "scale": False
                }
        headers["confirmation_padj"] = {
                "title": "Confirmation adjusted p-value",
                "description": "Adjusted p-values for the confirmation stage, i.e. whether a transcript is involved in DTU",
                "format": "{:,.1e}",
                "scale": False
                }
        headers['Log2foldchange'] = {
                "title": "Log2 fold change",
                "description": "Log2 fold change of expression levels of this transcript between the two groups",
                "format": "{:,.2f}",
                "min": -5,
                "max": 5,
                "bars_zero_centrepoint": True,
                "scale": "RdYlBu"
                }
        table_config = {
                "id": "Top_DTU_table",
                "col1_header": "Rank",
                "sortRows": False
                }
        # Make one table for each comparison
        for sample_name, dexseq_results in self.dexseq_results.items():
            group1, group2 = sample_name.split("_vs_")
            # Make a subset of up to top 30 DTU genes
            top30_genes = dexseq_results.drop_duplicates(subset=['Gene ID']).head(30)['Gene ID']
            plot_data = dexseq_results.loc[(dexseq_results['Gene ID'].isin(top30_genes)) & (dexseq_results['screening_padj']<dexseq_alpha), ]
            if len(plot_data):
                # Create a link for the gene and transcript
                if 'Gene Name' in plot_data.columns:
                    gene_label = 'Gene Name'
                else:
                    gene_label = 'Gene ID'
                if link_prefix is not None:
                    plot_data['gene_link'] = plot_data.apply(lambda x:link_prefix+x['Gene ID']+'" target="_blank">'+x[gene_label]+"</a>", axis=1)
                    plot_data['transcript_link'] = plot_data.apply(lambda x:link_prefix+x['Transcript ID']+'" target="_blank">'+x['Transcript ID']+"</a>", axis=1)
                else:
                    plot_data['gene_link'] = plot_data[gene_label]
                    plot_data['transcript_link'] = plot_data['Transcript ID']
                plot_data['Log2foldchange'] = plot_data["log2fold_{}_{}".format(group1,group2)]
                plot_data = plot_data.to_dict('records')
                idx_gene = 0
                idx_transcript = 1
                gene_id = ''
                table_data = OrderedDict()
                for row in plot_data:
                    if row['Gene ID'] == gene_id:
                        idx_transcript += 1
                    else:
                        gene_id = row['Gene ID']
                        idx_gene += 1
                        idx_transcript = 1
                    table_data["{}.{}".format(idx_gene,idx_transcript)] = row
                table_plot_html = table.plot(table_data, headers, table_config)
                self.add_section(
                    name = "Top genes involved in DTU in comparison {} vs. {}".format(group1, group2),
                    anchor = "top_dtu_genes_"+sample_name,
                    description = ("Top genes and their transcripts involved in DTU, in comparison {} vs. {}. "
                    "Up to the top 30 genes that passed screening stage are shown here. "
                    "Adjusted p-values for screening stage and confirmation stages for all transcripts of the same gene are shown. "
                    "Transcripts with positive Log2 fold changes have higher expression in {}. "
                    "Transcripts with negative Log2 fold changes have higher expression in {}. "
                    "Please refer to the results file for full list of genes and transcripts.").format(group1, group2, group1, group2),
                    plot = table_plot_html
                    )
            else:
                self.add_section(
                    name = "Top genes involved in DTU in comparison {} vs. {}".format(group1, group2),
                    anchor = "top_dtu_genes",
                    description = ("No genes were found to be involved in differential transcript usage. "
                    "Please refer to the results file for full list of genes and their p-values.")
                )

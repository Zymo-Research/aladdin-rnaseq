#!/usr/bin/env python

"""
MultiQC module to parse output from the 2-step trimming (adapter+polyA) used for 3'mRNAseq kit.
Code modified from MultiQC v1.9 cutadapt module.
The 2-step process runs cutadapt twice, but we want to combine the stats of those two steps to simplify the report for customers.
"""

from __future__ import print_function
from collections import OrderedDict, defaultdict
import logging
import re

from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger('multiqc')

class MultiqcModule(BaseMultiqcModule):
    """
    2-step trimming report Module Class, combining two cutadapt outputs into one report section.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Read Triming',
            anchor='trimming',
            target="",
            info="This section summarizes the results of read trimming including both adapter and poly-A triming."
        )

        # Find and load any Trim Galore reports
        self.cutadapt_data = dict()
        self.cutadapt_length_counts = dict()
        self.cutadapt_length_exp = dict()
        self.cutadapt_length_obsexp = dict()

        for f in self.find_log_files('trimming_2step/1st', filehandles=True):
            self.parse_cutadapt_logs(f, '1st')

        for f in self.find_log_files('trimming_2step/2nd', filehandles=True):
            self.parse_cutadapt_logs(f, '2nd')

        # Filter to strip out ignored sample names
        self.cutadapt_data = self.ignore_samples(self.cutadapt_data)

        if len(self.cutadapt_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.cutadapt_data)))

        # Calculate a few extra numbers of our own
        for s_name, d in self.cutadapt_data.items():
            if '1st' in d and '2nd' in d:
                # Percent trimmed
                self.cutadapt_data[s_name]['percent_trimmed'] = float(d['1st']['bp_processed'] - d['2nd']['bp_written']) / d['1st']['bp_processed'] * 100
                # Percent reads pass filter
                # SE reads only for current 2-step trimming
                self.cutadapt_data[s_name]['r_written'] = d['2nd']['r_written']
                self.cutadapt_data[s_name]['r_too_short'] = d['1st'].get('r_too_short', 0) + d['2nd'].get('r_too_short', 0)
                self.cutadapt_data[s_name]['r_too_long'] = d['1st'].get('r_too_long', 0) + d['2nd'].get('r_too_long', 0)
                self.cutadapt_data[s_name]['r_too_many_N'] = d['1st'].get('r_too_many_N', 0) + d['2nd'].get('r_too_many_N', 0)
                self.cutadapt_data[s_name]['percent_pass_filter'] = float(self.cutadapt_data[s_name]['r_written']) / d['1st']['r_processed'] * 100
            else:
                log.error("Either 1st or 2nd step trimming log missing for {}".format(s_name))
                raise UserWarning

        # Write parsed report data to a file
        self.write_data_file(self.cutadapt_data, 'multiqc_cutadapt')

        # Basic Stats Table
        self.cutadapt_general_stats_table()

        # Bar plot with number of reads trimmed
        self.cutadapt_filtered_barplot()

        # Trimming Length Profiles
        self.cutadapt_length_trimmed_plot()


    def parse_cutadapt_logs(self, f, step):
        """ Go through log file looking for cutadapt output """
        fh = f['f']
        # Removed support for cutadapt v1.6 or older.
        regexes = {
            "bp_processed": "Total basepairs processed:\s*([\d,]+) bp",
            "bp_written": "Total written \(filtered\):\s*([\d,]+) bp",
            "quality_trimmed": "Quality-trimmed:\s*([\d,]+) bp",
            "r_processed": "Total reads processed:\s*([\d,]+)",
            "pairs_processed": "Total read pairs processed:\s*([\d,]+)",
            "r_with_adapters": "Reads with adapters:\s*([\d,]+)",
            "r1_with_adapters": "Read 1 with adapter:\s*([\d,]+)",
            "r2_with_adapters": "Read 2 with adapter:\s*([\d,]+)",
            "r_too_short": "Reads that were too short:\s*([\d,]+)",
            "pairs_too_short": "Pairs that were too short:\s*([\d,]+)",
            "r_too_long": "Reads that were too long:\s*([\d,]+)",
            "pairs_too_long": "Pairs that were too long:\s*([\d,]+)",
            "r_too_many_N": "Reads with too many N:\s*([\d,]+)",
            "pairs_too_many_N": "Pairs with too many N:\s*([\d,]+)",
            "r_written": "Reads written \(passing filters\):\s*([\d,]+)",
            "pairs_written": "Pairs written \(passing filters\):\s*([\d,]+)",
        }
        s_name = None
        log_section = None
        for l in fh:
            # New log starting
            # Get sample name from end of command line params
            if l.startswith('Command line parameters'):
                f_name = l.strip().split()[-1]
                s_name = self.clean_s_name(f_name, f['root'])
                # Clear data if sample already exists, otherwise add the dict for this step
                if s_name in self.cutadapt_data:
                    if step in self.cutadapt_data[s_name]:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                        self.cutadapt_data[s_name] = { step : defaultdict(int) }
                    else:
                        self.cutadapt_data[s_name][step] = defaultdict(int)
                else:
                    self.cutadapt_data[s_name] = { step : defaultdict(int) }

            if s_name is not None:
                # Search regexes for overview stats and sum R1 and R2 stats
                for k, r in regexes.items():
                    match = re.search(r, l)
                    if match:
                        self.cutadapt_data[s_name][step][k] += int(match.group(1).replace(',', ''))

                # Starting a new section
                if '===' in l:
                    log_section = l.strip().strip('=').strip()

                # Histogram showing lengths trimmed, only for the 1st step
                if step == '1st':
                    if 'length' in l and 'count' in l and 'expect' in l:
                        plot_sname = s_name
                        if log_section is not None:
                            plot_sname = '{} - {}'.format(plot_sname, log_section)
                        self.cutadapt_length_counts[plot_sname] = dict()
                        self.cutadapt_length_exp[plot_sname] = dict()
                        self.cutadapt_length_obsexp[plot_sname] = dict()
                        # Nested loop to read this section while the regex matches
                        for l in fh:
                            r_seqs = re.search("^(\d+)\s+(\d+)\s+([\d\.]+)", l)
                            if r_seqs:
                                a_len = int(r_seqs.group(1))
                                self.cutadapt_length_counts[plot_sname][a_len] = int(r_seqs.group(2))
                                self.cutadapt_length_exp[plot_sname][a_len] = float(r_seqs.group(3))
                                if float(r_seqs.group(3)) > 0:
                                    self.cutadapt_length_obsexp[plot_sname][a_len] = float(r_seqs.group(2)) / float(r_seqs.group(3))
                                else:
                                    # Cheating, I know. Infinity is difficult to plot.
                                    self.cutadapt_length_obsexp[plot_sname][a_len] = float(r_seqs.group(2))
                            else:
                                break

        if s_name is not None:
            self.add_data_source(f, s_name)

    def cutadapt_general_stats_table(self):
        """ Take the parsed stats from the Trim Galore report and add it to the
        basic stats table at the top of the report """

        headers = {}
        headers['percent_trimmed'] = {
            'title': '% BP Trimmed',
            'description': '% Total Base Pairs trimmed',
            'namespace': 'cutadapt',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['percent_pass_filter'] = {
            'title': '% Reads PF',
            'description': '% Reads pass various length and quality filters',
            'namespace': 'cutadapt',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn'
        }
        self.general_stats_addcols(self.cutadapt_data, headers)

    def cutadapt_filtered_barplot(self):
        """ Bar plot showing proportion of reads trimmed """

        pconfig = {
            'id': 'cutadapt_filtered_reads_plot',
            'title': 'cutadpat: Filtered Reads',
            'ylab': 'Counts'
        }

        # Just SE data for now
        cats = OrderedDict()
        cats['r_written'] = { 'name': 'Reads passing filters' }
        cats['r_too_short'] = { 'name': 'Reads that were too short' }
        cats['r_too_long'] = { 'name': 'Reads that were too long' }
        cats['r_too_many_N'] = { 'name': 'Reads with too many N' }

        self.add_section(
            name = 'Filtered Reads',
            anchor = 'cutadapt_filtered_reads',
            description = 'This plot shows the number of reads retained or removed by cutadpat after adapter and poly-A trimming.',
            plot = bargraph.plot(self.cutadapt_data, cats, pconfig)
        )

    def cutadapt_length_trimmed_plot (self):
        """ Generate the trimming length plot """

        pconfig = {
            'id': 'cutadapt_trimmed_sequences_plot',
            'title': 'cutadapt: Lengths of Trimmed Sequences',
            'ylab': 'Counts',
            'xlab': 'Length Trimmed (bp)',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': [
                {'name': 'Counts', 'ylab': 'Count'},
                {'name': 'Obs/Exp', 'ylab': 'Observed / Expected'}
            ]
        }

        self.add_section(
            name = 'Trimmed Sequence Lengths',
            anchor = 'cutadapt_trimmed_sequences',
            description = ('This plot shows the number of reads with certain lengths of adapter trimmed. '
                           'Quality trimmed and hard trimmed sequences are not included. '
                           'Poly-A trimming is not included, either.'),
            helptext = '''
            Obs/Exp shows the raw counts divided by the number expected due to sequencing errors.
            A defined peak may be related to adapter length.
            See the [cutadapt documentation](http://cutadapt.readthedocs.org/en/latest/guide.html#how-to-read-the-report)
            for more information on how these numbers are generated.
            ''',
            plot = linegraph.plot([self.cutadapt_length_counts, self.cutadapt_length_obsexp], pconfig)
        )

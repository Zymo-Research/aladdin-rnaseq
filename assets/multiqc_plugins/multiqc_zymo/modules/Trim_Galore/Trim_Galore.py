#!/usr/bin/env python

"""
MultiQC module to parse output from Trim Galore! 
Code modified from MultiQC v1.9 cutadapt module.
Although cutadapt module can understand Trim Galore! logs, it misses key information like reads discarded for being too short and others.
This module aims to address that and also plots R1 and R2 trimming stats separately.
"""

from __future__ import print_function
from collections import OrderedDict, defaultdict
import logging
import re

from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Trim Galore module class, parses its trimming_report.txt files.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Trim Galore', anchor='trim_galore',
        href='https://github.com/FelixKrueger/TrimGalore',
        info="is a wrapper around Cutadapt and FastQC to consistently apply adpater"\
         " and quality trimming to FastQ files.")

        # Find and load any Trim Galore reports
        self.trim_galore_data = dict()
        self.trim_galore_length_counts = dict()
        self.trim_galore_length_exp = dict()
        self.trim_galore_length_obsexp = dict()

        for f in self.find_log_files('Trim_Galore', filehandles=True):
            self.parse_trim_galore_logs(f)

        # Filter to strip out ignored sample names
        self.trim_galore_data = self.ignore_samples(self.trim_galore_data)

        if len(self.trim_galore_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.trim_galore_data)))

        # Calculate a few extra numbers of our own
        for s_name, d in self.trim_galore_data.items():
            # Percent trimmed
            if 'bp_processed' in d and 'bp_written' in d:
                self.trim_galore_data[s_name]['percent_trimmed'] = float(d['bp_processed'] - d['bp_written']) / d['bp_processed'] * 100
            # Percent reads pass filter
            if 'r_processed' in d:
                # PE reads
                if 'R1' in d and 'R2' in d:
                    self.trim_galore_data[s_name]['pairs_written'] = d['r_processed'] / 2 - d.get('pairs_too_short',0) - d.get('pairs_too_many_N',0)
                    self.trim_galore_data[s_name]['percent_pass_filter'] = float(self.trim_galore_data[s_name]['pairs_written']) / d['r_processed'] * 200
                # SE reads
                else:
                    self.trim_galore_data[s_name]['r_written'] = d['r_processed'] - d.get('r_too_short',0) - d.get('r_too_long',0) - d.get('r_too_many_N',0)
                    self.trim_galore_data[s_name]['percent_pass_filter'] = float(self.trim_galore_data[s_name]['r_written']) / d['r_processed'] * 100

        # Write parsed report data to a file
        self.write_data_file(self.trim_galore_data, 'multiqc_trim_galore')

        # Basic Stats Table
        self.trim_galore_general_stats_table()

        # Bar plot with number of reads trimmed
        self.trim_galore_filtered_barplot()

        # Trimming Length Profiles
        self.trim_galore_length_trimmed_plot()


    def parse_trim_galore_logs(self, f):
        """ Go through log file looking for Trim Galore output """
        fh = f['f']
        # Trim galore always run cutadapt separately for R1 and R2, so stats for R1 and R2 are irrelavant.
        # Trim galore has a few lines appended to the end of cutadapt outputs, stating no. reads removed.
        # Removed support for cutadapt v1.6 or older.
        regexes = {
            'bp_processed': "Total basepairs processed:\s*([\d,]+) bp",
            'bp_written': "Total written \(filtered\):\s*([\d,]+) bp",
            'bp_quality_trimmed': "Quality-trimmed:\s*([\d,]+) bp",
            'r_processed': "Total reads processed:\s*([\d,]+)",
            'r_with_adapters': "Reads with adapters:\s*([\d,]+)",
            'r_too_short': "Sequences removed because they became shorter than the length cutoff of \d+ bp:\s*(\d+)",
            'r_too_long': "Sequences removed because after trimming they were longer than the maximum length cutoff of \d+ bp:\s*(\d+)",
            'r_too_many_N': "Sequences removed because they contained more Ns than the cutoff of \d+:\s*(\d+)",
            'pairs_too_short': "Number of sequence pairs removed because at least one read was shorter than the length cutoff \(\d+ bp\):\s*(\d+)",
            'pairs_too_many_N': "Number of sequence pairs removed because at least one read contained more N\(s\) than the specified limit of \d+:\s*(\d+)",
        }
        s_name = None
        log_section = None
        read_number = 'R1'
        for l in fh:
            # New log starting
            # Get sample name from end of command line params
            if l.startswith('Command line parameters'):
                f_name = l.strip().split()[-1]
                # Figure out whether R1 or R2
                if '_R1.fastq' in f_name:
                    read_number = 'R1'
                elif '_R2.fastq' in f_name:
                    read_number = 'R2'
                s_name = self.clean_s_name(f_name, f['root'])
                # Clear data if sample already exists and same R1/R2 has been found before
                if s_name in self.trim_galore_data:
                    if read_number in self.trim_galore_data[s_name]:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                        self.trim_galore_data[s_name] = defaultdict(int)
                else:
                    self.trim_galore_data[s_name] = defaultdict(int)
                # Record that data for R1/R2 has been found
                self.trim_galore_data[s_name][read_number] = 1

            if s_name is not None:
                # Search regexes for overview stats and sum R1 and R2 stats
                for k, r in regexes.items():
                    match = re.search(r, l)
                    if match:
                        self.trim_galore_data[s_name][k] += int(match.group(1).replace(',', ''))

                # Starting a new section
                if '===' in l:
                    log_section = l.strip().strip('=').strip()

                # Histogram showing lengths trimmed
                if 'length' in l and 'count' in l and 'expect' in l:
                    plot_sname = "{}_{}".format(s_name, read_number)
                    if log_section is not None:
                        plot_sname = '{} - {}'.format(plot_sname, log_section)
                    self.trim_galore_length_counts[plot_sname] = dict()
                    self.trim_galore_length_exp[plot_sname] = dict()
                    self.trim_galore_length_obsexp[plot_sname] = dict()
                    # Nested loop to read this section while the regex matches
                    for l in fh:
                        r_seqs = re.search("^(\d+)\s+(\d+)\s+([\d\.]+)", l)
                        if r_seqs:
                            a_len = int(r_seqs.group(1))
                            self.trim_galore_length_counts[plot_sname][a_len] = int(r_seqs.group(2))
                            self.trim_galore_length_exp[plot_sname][a_len] = float(r_seqs.group(3))
                            if float(r_seqs.group(3)) > 0:
                                self.trim_galore_length_obsexp[plot_sname][a_len] = float(r_seqs.group(2)) / float(r_seqs.group(3))
                            else:
                                # Cheating, I know. Infinity is difficult to plot.
                                self.trim_galore_length_obsexp[plot_sname][a_len] = float(r_seqs.group(2))
                        else:
                            break

        if s_name is not None:
            self.add_data_source(f, s_name)

    def trim_galore_general_stats_table(self):
        """ Take the parsed stats from the Trim Galore report and add it to the
        basic stats table at the top of the report """

        headers = {}
        headers['percent_trimmed'] = {
            'title': '% BP Trimmed',
            'description': '% Total Base Pairs trimmed',
            'namespace': 'Trim_Galore',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['percent_pass_filter'] = {
            'title': '% Reads PF',
            'description': '% Reads pass various length and quality filters',
            'namespace': 'Trim_Galore',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn'
        }
        self.general_stats_addcols(self.trim_galore_data, headers)

    def trim_galore_filtered_barplot(self):
        """ Bar plot showing proportion of reads trimmed """

        pconfig = {
            'id': 'trim_galore_filtered_reads_plot',
            'title': 'Trim Galore: Filtered Reads',
            'ylab': 'Counts'
        }

        # We just use all categories. If a report is generated with a mixture
        # of SE and PE data then this means quite a lot of categories.
        # Usually, only a single data type is used though - in that case
        # any categories with 0 across all samples will be ignored.
        cats = OrderedDict()
        cats['pairs_written'] = { 'name': 'Pairs passing filters' }
        cats['r_written'] = { 'name': 'Reads passing filters' }
        cats['pairs_too_short'] = { 'name': 'Pairs that were too short' }
        cats['r_too_short'] = { 'name': 'Reads that were too short' }
        cats['pairs_too_long'] = { 'name': 'Pairs that were too long' }
        cats['r_too_long'] = { 'name': 'Reads that were too long' }
        cats['pairs_too_many_N'] = { 'name': 'Pairs with too many N' }
        cats['r_too_many_N'] = { 'name': 'Reads with too many N' }

        self.add_section(
            name = 'Filtered Reads',
            anchor = 'trim_galore_filtered_reads',
            description = 'This plot shows the number of reads (SE) / pairs (PE) removed by Trim Galore.',
            plot = bargraph.plot(self.trim_galore_data, cats, pconfig)
        )

    def trim_galore_length_trimmed_plot (self):
        """ Generate the trimming length plot """

        pconfig = {
            'id': 'trim_galore_trimmed_sequences_plot',
            'title': 'Trim Galore: Lengths of Trimmed Sequences',
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
            anchor = 'trim_galore_trimmed_sequences',
            description = ('This plot shows the number of reads with certain lengths of adapter trimmed. '
                           'Quality trimmed and hard trimmed sequences are not included.'),
            helptext = '''
            Obs/Exp shows the raw counts divided by the number expected due to sequencing errors.
            A defined peak may be related to adapter length.
            See the [cutadapt documentation](http://cutadapt.readthedocs.org/en/latest/guide.html#how-to-read-the-report)
            for more information on how these numbers are generated.
            ''',
            plot = linegraph.plot([self.trim_galore_length_counts, self.trim_galore_length_obsexp], pconfig)
        )

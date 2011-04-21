import sys
from matplotlib import pyplot as plt
import math
import os
import ConfigParser
from optparse import OptionParser
from IPython.Debugger import Tracer
debug = Tracer()
# Provide:
# ts_ids 
# if no ts_ids are provided, plot the 10 most significantly distinct 3utrs.

# Program will read UTR_SETTINGS to get the datasets from which the ts_ids
# should be plotted from.


def get_settings(settings_file, here):
    """Read settings file for parameters needed to run the script."""

    conf = ConfigParser.ConfigParser()
    conf.read(settings_file)

    if 'PLOTTING' not in conf.sections():
        print('Section "PLOTTING" not found in {0}'.format(settings_file))
        sys.exit()

    dsets = [dset for (dset, val) in dict(conf.items('PLOTTING')).items() if val]
    # Now you have the name of the dsets. Get the paths of the summary and
    # coverage files for the dsets, based on their pre-determined paths.

    # coverage files are found in here/temp_files/covered_dset
    # summary files are found in here/output/utr_dset
    coverage = dict((dset, os.path.join(here, 'temp_files', 'covered_'+dset)) for
                dset in dsets)

    summary = dict((dset, os.path.join(here, 'output', 'utr_'+dset)) for dset
                   in dsets)

    # Check that all files exist
    for needed_file in coverage.values() + summary.values():
        if not os.path.isfile(needed_file):
            print('File:\n{0}\nnot found. Exiting...'.format(needed_file))
            sys.exit()

    return (summary, coverage)

def summary_reader(summary):

    # Make a relationship between dset and index for future use

    # Prepare a dictionary of its_id with relative lenghts of compartments
    cuml_99_dict = {}

    # Get the ids of the transcripts whose extensions do not overlap exons in
    # the gencode annotation. (it's a set() by default)
    exons = '/users/rg/jskancke/phdproject/3UTR/gencode5/just_exons/just_exons.bed'

    for (dset, summary_path) in summary.items():
        # Create a dict for each dset
        cuml_99_dict[dset] = {}
        for line in open(summary_path, 'rb'):
            # Skip lines that don't start with chr (it's fast enough)
            if line[:3] != 'chr':
                continue
            #(chrm, beg, end, ts_ID, strand, cuml99, int_mean_99,\
             #int_mean_99, int_mean_annotation, int_mean_annotation,\
             #PAS_type, pas_distance, rpkm) = line.split()
            cuml_99_dict[dset][line.split()[3]] = line.split()

    return cuml_99_dict

def parse_arguments():
    """Set up the argument parsing. Return """

    desc = """%prog plots the read-coverage of 3'UTRs. A single transcript ID
    can be supplied with the -i option, and a file with a list of transcript
    IDs can be supplied with the -f option. If no options are defined, the top
    5 most significant 3UTRs will be plotted.
    For %prog to work, make sure that the datasets whose read-coverage should
    be plotted are listed under the 'PLOTTING' section in the file
    UTR_COVERAGE. The names of the datasets should be the same as the names
    given in the 'DATASETS' section, also in UTR_COVERAGE."""

    how_to = "usage: %prog [option] [argument]"

    parser = OptionParser(usage = how_to, description = desc)

    parser.add_option("-i", "--transcript",
                      default = False,
                      dest = 'ts_id_single',
                      action = 'store',
                      help = "transcript ID whose 3'UTR will be plotted")

    parser.add_option("-f", "--transcriptFile",
                      default = False,
                      dest = 'ts_id_file',
                      action = 'store',
                      help = "file with one or more transcript IDs for plotting")

    (options, args) = parser.parse_args()

    # if both options were given, complain.
    if options.ts_id_single and options.ts_id_file:
        parser.error("Two options received. Use either -i or -f but not both.")

    return (options.ts_id_single, options.ts_id_file)

def top_ts_id_targets(summary_dict):
    """Go through summary dict and get five TS that have a within-length
    relative difference of at least 0.3. Select the five based on rpkm values
    (for all conditions)."""

    # Make a new dictionary that holds the relative_lenghts 
    conditions = summary_dict.keys()
    rel_len_keepers = {}

    # For each condition, add to cumul99 value if rpkm is high enough
    for condition in conditions:
        for ts_id, info in summary_dict[condition].items():
            (chrm, beg, end, ts_ID, strand, cuml99, int_mean_99, int_mean_99,
             int_mean_annotation, int_mean_annotation, PAS_type, pas_distance,\
             rpkm) = info
            rpkm = float(rpkm)
            # Add only those with rpkm over 2
            if rpkm > 2:
                if ts_id in rel_len_keepers:
                    rel_len_keepers[ts_id].append(float(cuml99))
                else:
                    rel_len_keepers[ts_id] = [float(cuml99)]

    # Keep only the transcripts that 1) have rel_lens for all conditions and 2)
    # have a minimum distance of 0.4 btween two of the rel_lens.
    nr_conditions = len(conditions)
    restricted_len_dict = {}
    for ts_id, rel_lens in rel_len_keepers.items():
        if len(rel_lens) < nr_conditions:
            continue
        else:
            # Calculate distances between all rrel_dist
            parw_dist = [abs(val1-val2) for val1 in rel_lens for val2 in
                         rel_lens]
            sig_dist = ['yes' for val in parw_dist if val > 0.4]
            if 'yes' in sig_dist:
                restricted_len_dict[ts_id] = rel_lens

    # Get a list of ts_ids with average rpkm
    ts_rpkm = []
    for ts_id in restricted_len_dict:
        avrg_rpkm = sum(float(summary_dict[condition][ts_id][-1]) for condition in
                        conditions)/len(conditions)
        ts_rpkm.append((avrg_rpkm, ts_id))

    ts_rpkm.sort(reverse=True)

    # Return top five
    return [tup[1] for tup in ts_rpkm[:5]]


def coverage_reader(coverage, summary_dict, ts_id_single, ts_id_file):
    # If ts_id_single, use only this one.
    if ts_id_single:
        key_ids = [ts_id_single]

    # If ts_id_file, read the ts_ids in the file.
    elif ts_id_file:
        key_ids = [line.strip() for line in open(ts_id_file, 'rb')]

    # If neither are supplied, print the top 5 most conspicuous targets
    if not ts_id_single or ts_id_file:
        key_ids = top_ts_id_targets(summary_dict)

    exp_nr = len(summary_dict.keys())
    covr_dict = dict((key, [[] for bla in range(exp_nr)]) for key in key_ids)

    # Get the coverage from the coverage-files
    for (indx, (dset, cvrg_path)) in enumerate(coverage.items()):
        for line in open(cvrg_path):
            (chrm, beg, end, ts_id, d, strnd, rel_pos, covrg) = line.split()
            if ts_id in key_ids:
                covr_dict[ts_id][indx].append(int(covrg))

    return covr_dict, key_ids

def plot_3utr_ends(ts_id, covr_dict, strand, exp_nr, dsets, figdir):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    tseries = [moving_average(covr_dict[ts_id][i], width=0) for i in range(exp_nr)]

    if strand == '-':
        tseries[:] = [list(reversed(tser)) for tser in tseries]
    for indx, ser in enumerate(tseries):
        ax.plot(ser, label=dsets[indx], linewidth=4)

    # TODO turn transcript ID to gene ID.
    ax.set_title('Transcript ID: ' + ts_id)
    ax.set_xlabel('Nucleotide from 3UTR start')
    ax.set_ylabel('Nucleotide coverage depth')
    ax.legend(loc=1, ncol=1, shadow=True)
    max_covr = max(sum(covr_dict[ts_id], []))
    ax.set_ylim(0, max_covr + math.floor(max_covr*0.3))
    filename = ts_id
    fig.savefig(os.path.join(figdir, filename + '.png'), format='png', dpi=300)
    #fig.show() # you only see figures from within ipython. must save figures to file.

def moving_average(dataseries, width):
    """ Divide the dataseries into blocks of size width"""

    # If no extension is called for, return input untouched.
    if width == 0:
        return dataseries

    newseries = []
    for pos, val in enumerate(dataseries):
        newseries.append(sum(dataseries[pos-width:pos+width])/width*2)

    return newseries

def utr_plotter(coverage, covr_dict, summary_dict, key_ids, imagedir):
    # Get the names of the dsets
    dsets = coverage.keys()
    exp_nr = len(coverage.keys())

    # Print the coverages in plots. Take care of strand.
    for ts_id in key_ids:
        strand = summary_dict[summary_dict.keys()[0]][ts_id][4]
        plot_3utr_ends(ts_id, covr_dict, strand, exp_nr, dsets, imagedir)


def main():
    here = os.path.dirname(os.path.realpath(__file__))
    settings_file = os.path.join(here, 'UTR_SETTINGS')
    imagedir = os.path.join(here, 'images')
    if not os.path.exists(imagedir):
        os.makedirs(imagedir)

    # Parse the arguments. Get either ts_ID, ts_ID_file, or neither.
    (ts_id_single, ts_id_file) = parse_arguments()

    # Get the paths to the summary and coverage files
    (summary, coverage) = get_settings(settings_file, here)

    # Get dict[dset][ts_ids] = line.split()
    summary_dict = summary_reader(summary)

    # Get coverage from coverage files
    covr_dict, key_ids = coverage_reader(coverage, summary_dict, ts_id_single,
                                         ts_id_file)
    # Plot the utrs
    utr_plotter(coverage, covr_dict, summary_dict, key_ids, imagedir)


if __name__ == '__main__':
    main()

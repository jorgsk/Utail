"""
Read the utr_settings file, download the split-map data, and convert this into
bed format. (Check if you have already downloaded and converted -- don't do it
again).
"""
########################################
# only get the debug function if run from Ipython #
def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

if run_from_ipython():
    from IPython.Debugger import Tracer
    debug = Tracer()
else:
    def debug():
        print ('Warning: debugging mark present.')
        pass
########################################

import ConfigParser
import os
from subprocess import Popen, PIPE

here = os.path.dirname(os.path.realpath(__file__))
savedir = os.path.join(here, 'splitmapped_reads')
settings_file = os.path.join(here, 'UTR_SETTINGS')

# create the directory if it doesn't exist
if not os.path.isdir(savedir):
    os.makedirs(savedir)

conf = ConfigParser.ConfigParser()
conf.optionxform = str
conf.read(settings_file)

datasets = dict((dset, files.split(':')) for dset,
                files in conf.items('DATASETS'))

# Go through all the items in 'datsets'. Pop the directories from the list.
# They are likely to be shortcuts.
for (dset, dpaths) in datasets.items():
    for pa in dpaths:
        if os.path.isdir(pa):
            datasets.pop(dset)

gtfpaths = {}

import glob
for (dset, dpaths) in datasets.items():
    root_path = os.path.dirname(os.path.dirname(dpaths[0]))
    splitdir = os.path.join(root_path, 'splitmapping')
    gtfpaths[dset] = glob.glob(splitdir+'/*single.unique.gtf.gz*')

bedpaths = {}

# call zcat on the files and merge them to bed
for dset, gtf_paths in gtfpaths.items():

    outfile = os.path.join(savedir, dset+'_splitmapped_reads.bed')

    for path in gtf_paths:
        cmd = ['zcat', path]
        p = Popen(cmd, stdout=PIPE)

        for line in p.stdout:

            # split output
            (chrm, d, d, beg, end, d, strand) = line.split()[:7]

            # bed format
            bed_output = '\t'.join([chrm, beg, end, '0', '0', strand]) + '\n'

            # write to std out
            outfile.write(bed_output)

    bedpaths[dset] = outfile

# call mergeBed on the dset files (use annotation_parser version)
from annotation_parser import merge_output

stranded = True
bedpaths = merge_output(bedpaths, stranded)



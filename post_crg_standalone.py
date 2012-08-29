"""
The purpose of this is only to output the polyA sites.
"""
from __future__ import division
print('Loading modules ...\n')
import os
import sys
import shutil
import ConfigParser
from multiprocessing import Pool
from multiprocessing import cpu_count
from operator import attrgetter

import re

# only get the debug function if run from Ipython
def run_from_ipython():
    try:
        __IPYTHON__active
        return True
    except NameError:
        return False

if run_from_ipython():
    from IPython.Debugger import Tracer
    #from IPython.core.debugger import Tracer
    debug = Tracer()
else:
    def debug(): pass

from subprocess import Popen
from subprocess import PIPE
import time
import math

# Your own imports
import annotation_parser as genome

class Settings(object):
    """
    An instance of this class holds all the settings parameters obtained from
    the UTR_SETTINGS file. Useful for passing to downstream code.
    """

    def __init__(self, datasets, read_limit, max_cores, chr1, gem_index):

        self.datasets = datasets
        # these two might be included later
        self.polyA_annotation_path = 0
        self.transcript_annotation_path = 0

        self.read_limit = read_limit
        self.max_cores = max_cores
        self.chr1 = chr1
        self.gem_index = gem_index

    def DEBUGGING(self):
        """
        For modifying the settings from UTR_FINDER -- only to be used in when
        debugging!
        """

        self.chr1 = True
        #self.chr1 = False
        #self.read_limit = False
        self.read_limit = 1000000 # less than 10000 no reads map
        self.max_cores = 3

def verify_access(f):
    """
    Try to access file 'f'. If not successful, which means that the file either
    doesn't exist or you don't have access to it, abort program.
    """
    try:
        open(f, 'rb')
    except:
        print('Could not access {0}.\nCheck if file exists, and if it exists '\
              'check permissions'.format(f))
        sys.exit()

def get_dataset_paths(dataset_filepaths):
    """
    Make a mapping from dsetId to paths and description
    """

    datasets = {}
    for line in open(dataset_filepaths, 'rb'):
        if line.startswith('#'):
            continue

        exp_id, description, paths = line.split('\t')
        paths = paths.rstrip() # remove end of line character

        subset = {'description': description, 'paths': paths.split(';')}

        datasets[exp_id] = subset

    return datasets

def read_settings(settings_file):
    """
    Read the settings and get all the settings parameters. These should be used
    to create a 'Settings' object.
    """

    conf = ConfigParser.SafeConfigParser()
    conf.optionxform = str
    conf.read(settings_file)

    expected_fields = ['DATASETS', 'ANNOTATION', 'CPU_CORES', 'RESTRICT_READS',
                       'CHROMOSOME1', 'GEM_MAPPER']

    missing = set(conf.sections()) - set(expected_fields)

    if len(missing) == 0:
        pass
    else:
        print('The following options sections are missing: {}'.format(missing))
        sys.exit()

    dataset_filepaths = conf.items('DATASETS')[0][1]
    datasets = get_dataset_paths(dataset_filepaths)

    # check if the datset files are actually there...
    #for dset, files in datasets.items():
        #[verify_access(f) for f in files]

    # cpu cores
    try:
        max_cores = conf.getint('CPU_CORES', 'max_cores')
    except ValueError:
        max_cores = cpu_count()-1
        if max_cores < 1:
            max_cores = 1

    # restrict number of reads from source
    try:
        read_limit = conf.getint('RESTRICT_READS', 'restrict_reads')
    except ValueError:
        read_limit = conf.getboolean('RESTRICT_READS', 'restrict_reads')

    # restrict to chromosome 1
    chr1 = conf.getboolean('CHROMOSOME1', 'only_chr1')

    # Get the gem_index_file
    gem_index = conf.get('GEM_MAPPER', 'gem_mapper_index')

    # you need the index file for the gem-mapper
    verify_access(gem_index+'.blf') # Check if is a file

    return(datasets, read_limit, max_cores, chr1, gem_index)

def make_directories(here, dirnames, DEBUGGING):
    """
    For each name in dirnames, return a list of paths to newly created
    directories in folder 'here'. Don't overwrite folders if they exist.
    """
    outdirs = []

    if DEBUGGING:
        # if DEBUGGING, make a DEBUGGING output directory
        debug_dir = os.path.join(here, 'DEBUGGING')
        if not os.path.exists(debug_dir):
            os.makedirs(debug_dir)

        # modify 'here' to poin to the debugging dir
        here = debug_dir


    for dirname in dirnames:

        this_dir = os.path.join(here, dirname)

        if not os.path.exists(this_dir):
            os.makedirs(this_dir)

        outdirs.append(this_dir)

    return outdirs

def run_dataset_map():
    """
    Create a mapping from labExpId to bam file locations under
    ENCODE/gingerascarrie...

    For this you require the files:
        1) Dataset.info.txt, which maps datasets to folders
        2) hg19_RNA_dashboard_files.txt.crg, maps datsets to labExpId

    You have copied both these files to your local drive.

    UPDATE: the bam files from David's pipeline do not contain the unmapped
    reads. wtflord. sticking to gem. I still need to make the mapping though.

    1) Make the dataset -> folder mapping K532_WC_2 will be the second
    replicate (1 being the first replicate ..) of cell, compartment.

    2) make the labExpId -> datset mapping
    """
    from glob import glob

    ging_folder = '/users/rg/dgonzalez/Projects/ENCODE/GingerasCarrie'

    # 1) dataset -> folder mapping
    dset_info = 'post_project_files/Dataset.info.txt'
    dset_path = {}

    for line in open(dset_info, 'rb'):
        set_code, cell_line, compartment, longshort, rpl = line.split()

        # get only the long poly(A) samples
        if longshort != 'LONGPOLYA':
            continue

        # only get gingeras datsets
        if set_code[:4] != 'Ging':
            continue

        # get rid of two pesky datasets
        if ('BROAD' in cell_line) or ('SKNSHRA' in cell_line):
            continue

        # don't get spikein or poly(A) minus
        if ('spike' in set_code) or ('Minus' in set_code):
            continue

        #if 'BROAD' in line:
            #debug()

        dataset = '_'.join([cell_line, compartment.lower(), rpl])

        folder = set_code[4:]

        dset_path[dataset] = glob(ging_folder+'/'+folder+'/SAM/*.map.gz')

    # 2) dataset -> labExpId
    labExpId_info = 'post_project_files/all_Gingeras_samples.tsv'
    dset_labExpId = {}

    headers = 0

    for line in open(labExpId_info, 'rb'):

        # skip comments
        if line.startswith('#'):
            headers = line.split()
            headers[0] = headers[0][1:]
            continue

        # make some rudementary filtes
        if ('RnaSeq' not in line) or ('longPolyA' not in line):
            continue

        # make a dictionary of the line
        line_info = dict(zip(headers, line.split()))

        # re-filter just to be sure some filters
        if line_info['dataType'] != 'RnaSeq':
            continue

        if line_info['rnaExtract'] != 'longPolyA':
            continue

        # obtain cell line, compartment, and repliate info
        cell_line = line_info['cell']

        # if the cell is hela or h1, change the wording to make it compatible
        # with the other textfile :(
        if cell_line == 'HeLa-S3':
            cell_line = 'HELAS3'

        elif cell_line == 'H1-hESC':
            cell_line = 'H1HESC'

        elif cell_line == 'MCF-7':
            cell_line = 'MCF7'

        elif cell_line == 'HepG2':
            cell_line = 'HEPG2'

        compartment =  line_info['localization']

        # not all entries have a replicate
        if 'replicate' in line_info:
            rpl = line_info['replicate']
        else:
            rpl = '1'

        dataset = '_'.join([cell_line, compartment, rpl])

        dset_labExpId[dataset] = line_info['labExpId']

    outp_file = 'post_project_files/id_descr_paths.tsv'
    outp_handle = open(outp_file, 'wb')

    outp_handle.write('#labExpID\tDescription\tPaths\n')

    # write to file 
    for dset, paths in dset_path.items():
        labId = dset_labExpId[dset]
        outp_handle.write('\t'.join([labId, dset, ';'.join(paths)]) + '\n')

    outp_handle.close()

    return outp_file

def main():
    """
    This method is called if script is run as __main__.
    """

    # obtain a file that maps from labID to paths and previous dset-info-keys
    # XXX not part of the pipeline; just here for convenience
    #map_file = run_dataset_map()

    # Set debugging mode. This affects the setting that are in the debugging
    # function (called below). It also affects the 'temp' and 'output'
    # directories, respectively.
    DEBUGGING = True #
    #DEBUGGING = False

    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    # Make directories needed by downstream code
    # If debugging, make debugging output
    dirnames = ['temp_files', 'source_bedfiles', 'polyA_output']
    (tempdir, beddir, output_dir) = make_directories(here, dirnames, DEBUGGING)

    # Location of settings file
    settings_file = os.path.join(here, 'polyA_settings.txt')

    # Create the settings object from the settings file
    settings = Settings(*read_settings(settings_file))

    # This option should be set only in case of debugging. It makes sure you
    # just run chromosome 1 and only extract a tiny fraction of the total reads.
    if DEBUGGING:
        settings.DEBUGGING()

    # The program reads a lot of information from the annotation. The annotation
    # object will hold this information (file-paths and datastructures).
    #print('Reading settings ...\n')
    #annotation = Annotation(settings.annotation_path,
                            #settings.annotation_format,
                            #settings.annotated_polyA_sites)

    # Create a pool of processes; one dataset will take up one process.
    my_pool = Pool(processes = settings.max_cores)
    results = []

    # Apply all datasets to the pool
    t1 = time.time()

    # exp_id and dset_reads are as given in UTR_SETTINGS
    for exp_id, dset_info in settings.datasets.items():

        # The arguments needed for the pipeline
        arguments = (exp_id, dset_info, tempdir, output_dir, settings,
                     DEBUGGING, here)

        ###### FOR DEBUGGING ######
        akk = pipeline(*arguments)
        ##########################
        debug()

        #result = my_pool.apply_async(pipeline, arguments)
        #results.append(result)

    my_pool.close()
    my_pool.join()

    # we don't return anything, but get results anyway
    [result.get() for result in results]

    # Print the total elapsed time
    print('Total elapsed time: {0}\n'.format(time.time()-t1))
    ###################################################################

def pipeline(exp_id, dset_info, tempdir, output_dir, settings, DEBUGGING,
             here):
    """
    Get reads, get polyA reads, cluster polyA reads, get coverage, combine it in
    a 3UTr object, do calculations on the object attributes, write calculation
    to output files ... this is where it all happens: the PIPEline.
    """

    # Give new names to some parameters to shorten the code
    read_limit = settings.read_limit

    just_dset = '_'.join(exp_id.split('_')[:-1])
    polyA_path = os.path.join(tempdir, 'polyA_reads_'+just_dset+'.fa')

    if DEBUGGING:
        cache_dir = os.path.join(os.path.join(here, 'DEBUGGING'),
                                 'polyA_cache')
    else:
        cache_dir = os.path.join(here, 'polyA_output')

    if not os.path.isdir(cache_dir):
        os.makedirs(cache_dir)

    print('Obtaining reads from mapping for {0} ...\n'.format(exp_id))
    all_reads = get_bed_reads(dset_info['paths'], exp_id, read_limit, tempdir,
                              polyA_path)

    # Extract info from reads
    (unmapped_reads, total_reads, total_mapped_reads, acount, tcount) = all_reads

    # PolyA pipeline: remove low-quality reads, remap, and -> .bed-format:

    # 1) Process reads by removing those with low-quality, and removing the
    #    leading Ts and/OR trailing As.
    processed_reads, avrg_read_len = process_reads(unmapped_reads)

    # 2) Map the surviving reads to the genome and return unique ones
    print('Remapping poly(A)-reads for {0} + noise filtering...'.format(exp_id))
    polyA_bed_path = map_reads(processed_reads, avrg_read_len, settings)

    # 3) Save the polyA_bed_path to the cache dir
    polydone_path = os.path.splitext(os.path.split(polyA_path)[1])[0]
    polyA_cached_path = os.path.join(cache_dir, polydone_path+\
                                     '_processed_mapped.bed')

    shutil.copyfile(polyA_bed_path, polyA_cached_path)

    debug()

def process_reads(pA_reads_path):
    """
    Remove reads that are too short or have a poor nucleotide composition.
    """
    processed_reads = os.path.splitext(pA_reads_path)[0]+'_processed.fas'
    outfile = open(processed_reads, 'wb')

    # Go through the file two lines at a time. If the next line does not begin
    # with '>', append line to last entry (it's a split fasta-file).

    length_sum = 0
    tot_reads = 0

    for line in open(pA_reads_path, 'rb'):
        # add lines until there are 2 entries in linepair
        (seq, at, tail_info) = line.split()
        seqlen = len(seq)
        if seqlen > 25:
            As = seq.count('A')
            Ts = seq.count('T')
            # only save if A/T frequencies are not abnormal
            if (As/seqlen < 0.70) and (Ts/seqlen < 0.70):
                length_sum += seqlen
                tot_reads += 1
                # trim the title
                outfile.write('>{0}\t{1}\n'.format(at, tail_info)+seq+'\n')

    if tot_reads > 0:
        avrg_len = length_sum/float(tot_reads)
    else:
        avrg_len = 0

    outfile.close()

    return processed_reads, avrg_len

def map_reads(processed_reads, avrg_read_len, settings):
    """
    Map the processed reads using gem-mapper. Use the average read length to
    determine the number of mismatches for the mapper according to the following
    scheme where X is the average read length:
        * if X < 50, then 1 mismatch
        * if 50 < X < 100, then 2 mismatchs
        * if 100 < X, then 3 mismatches
    Regardless, only accept uniquely mapping reads.
    """

    # File path of the mapped reads
    mapped_reads = os.path.splitext(processed_reads)[0]+'_mapped'

    # Naming the final output
    polybed_path = os.path.splitext(processed_reads)[0] + '_mapped.bed'

    # How many mismatches depend on read length
    if avrg_read_len < 50:
        mismatch_nr = 1
    elif 50 < avrg_read_len < 100:
        mismatch_nr = 2
    elif 100 < avrg_read_len:
        mismatch_nr = 3

    ### mapping trimmed reads
    command = "gem-mapper -I {0} -i {1} -o {2} -q ignore -m {3}"\
            .format(settings.gem_index, processed_reads, mapped_reads, mismatch_nr)

    p = Popen(command.split())
    p.wait()

    # Accept mismatches according to average read length
    acceptables = {1: set(('1:0', '0:1')), 2: set(('1:0:0', '0:1:0', '0:0:1')),
                   3: set(('1:0:0:0', '0:1:0:0', '0:0:1:0', '0:0:0:1'))}

    acceptable = acceptables[mismatch_nr]
    getstrand = {'R':'-', 'F':'+'}
    start_re = re.compile('[0-9]*')

    reads_file = open(polybed_path, 'wb')

    # count the number of noisy reads and total reads
    noisecount = 0
    allcount = 0

    for line in open(mapped_reads + '.0.map', 'rb'):
        (at, tail_info, seq, mapinfo, position) = line.split('\t')

        # Acceptable reads and poly(A) reads are mutually exclusive.
        if mapinfo in acceptable:
            allcount += 1
            # Get chromosome, strand, and beg
            (chrom, rest) = position.split(':')
            strand = getstrand[rest[0]]
            beg = start_re.match(rest[1:]).group()

            ## Don't write if this was a noisy read
            if read_is_noise(chrom, strand, beg, seq, at, tail_info, settings):
                noisecount += 1
                continue

            reads_file.write('\t'.join([chrom, beg, str(int(beg)+len(seq)), '.',
                                      '.', strand]) + '\n')
    # close file
    reads_file.close()

    # Write to logfile
    if allcount > 0:
        vals = (noisecount, allcount, noisecount/float(allcount))

        noiseinf = '\nNoise reads: {0}\nTotal reads: {1}\nNoise ratio: {2:.2f}\n'\
               .format(*vals)

        noiselog = open('NOISELOG.LOG', 'ab')
        noiselog.write('-'*80+'\n')
        noiselog.write(polybed_path)
        noiselog.write(noiseinf)
        noiselog.write('-'*80+'\n')
        noiselog.close()

    return polybed_path

def read_is_noise(chrm, strand, beg, seq, at, tail_info, settings):
    """
    Get the original sequence. If the A/T count on the geneomic sequence is
    within 30% of the A/T count on the genome, ignore the read. This corresponds
    to 4:1, 5:1, 6:1, 7:2, 8:2, 9:2, 10:3, 11:3, 12:3, 14:, 4
    40 % would be
    to 4:1, 5:2, 6:2, 7:2, 8:4, 9:3, 10:4, 11:4, 12:5, 14:, 5
    of course, this is a simply solution. In fact, the longer the tail, the less
    likely it is that the reads follow each other. But this will have to do.

    However, keep in mind that you accept all kinds of reads, irrespective of
    quality measure. To do this properly, you need to take qualities into
    account.
    """
    frac = 0.3
    # get the sequences! note: the sequence fetcher automatically
    # reverse-transribes regions on the - strand. In other words, it
    # always returns sequences in the 5->3 direction.
    # However
    seq_len = len(seq)
    beg = int(beg)
    end = beg+seq_len
    # for getting the nucleotide counts fast
    nuc_dicts = dict(b.split('=') for b in tail_info.split(':'))
    strip_len = sum([int(val) for val in nuc_dicts.values()])

    ## it came from the - strand
    if (at == 'T' and strand == '+') or (at == 'A' and strand == '-'):
        # get beg-strip_len:beg
        seq_dict = {'seq': (chrm, beg-strip_len-1, beg, strand)}
        d = genome.get_seqs(seq_dict, settings.hgfasta_path)

    ## it came from the + strand
    if (at == 'T' and strand == '-') or (at == 'A' and strand == '+'):
        # get beg+seq_len:beg+seq_len+strip_len
        seq_dict = {'seq': (chrm, end, end+strip_len+1, strand)}
        d = genome.get_seqs(seq_dict, settings.hgfasta_path)

    genome_count = d['seq'].count(at)
    tail_count = int(nuc_dicts[at])

    winsize = int(math.floor(tail_count*frac))
    window = range(tail_count-winsize, tail_count) +\
            range(tail_count, tail_count+winsize)

    # return as true noise 
    if genome_count in window:
        return True
    # return as not-noise
    else:
        return False

    pass

def get_bed_reads(dset_paths, exp_id, read_limit, tempdir, polyA_path):
    """
    Get reads from file. Determine file-type. If gem, extract from gem and
    convert to bed. If .bed, concatenate the bedfiles and convert them to the
    desired internal format. Works like a wrapper around zcat.

    Write unmapped reads with leading poly(T) or tailing poly(A) to .bed file.
    The read is written to the bed-file as stripped of polyA or polyT stretch.
    Either 5 contiguous A/Ts or 6 A/Ts in the last/first 7 nucleotides must be
    present for the read to be considered as a putative poly(A) read.
    """

    # Keep track on the number of reads you get for the sake of RPKM
    total_reads = 0
    total_mapped_reads = 0
    acount = 0
    tcount = 0

    # File object for writing to
    polyA_file = open(polyA_path, 'wb')

    # Accept up to two mismatches. Make as set for speedup.
    acceptable = set(('1:0:0', '0:1:0', '0:0:1'))

    # nucleotides
    nucleotides = set(['G', 'A', 'T', 'C'])

    # Run zcat with -f to act as noram cat if the gem-file is not compressed
    cmd = ['zcat', '-f'] + dset_paths
    f = Popen(cmd, stdout=PIPE)

    # You need an overall method for determining tails.
    # Two pieces of information
    # 1) minimum length of tail
    # 2) nr of non-A/T nucleotides ratio 1:3, 1:4, 1:5?
    # 4) Commence with a minimum UTR length of 4 and 1:5 ratio
    # 5) let these be optional in the final version of Utail
    # 7) ratio denominator must b
    # 8) also accept and iteratively trim 2/10 2/12!

    min_length = 5
    ratio_denominator = 6
    rd = ratio_denominator

    min_A = 'A'*min_length
    min_T = 'T'*min_length

    # Make regular expression for getting read-start
    #start_re = re.compile('[0-9]*')
    trail_A = re.compile('A{{{0},}}'.format(min_length))
    lead_T = re.compile('T{{{0},}}'.format(min_length))
    non_A = re.compile('[^A]')
    non_T = re.compile('[^T]')

    for map_line in f.stdout:
        #
        try:
            (ID, seq, quality, mapinfo, position) = map_line.split('\t')
        except ValueError:
            print 'read error from line {0} in {1} '.format(map_line, dset_paths)
            continue

        total_reads +=1

        # If short, only get up to 'limit' of reads
        if read_limit and total_reads > read_limit:
            break

        # count mapped reads
        if mapinfo in acceptable:
            total_mapped_reads +=1

        # Look for unmapped reads
        if mapinfo[:5] == '0:0:0':
            if seq[:2] == 'NN':
                seq = seq[2:]

            # If more than two bases are ambigious, discard.
            if seq.count('N') > 2:
                continue

            # Check for putative poly(T)-head. Remove tail and write to file.
            if (seq[:min_length] == min_T) or (seq[:rd].count('T') >=rd-1)\
               or (seq[:2*rd].count('T') >=2*rd-2):

                t_strip = strip_tailT(seq, lead_T, non_T, min_length, rd,
                                      min_T)
                striplen = len(t_strip)
                if striplen > 25:
                    seqlen = len(seq)
                    tail = seq[:seqlen-striplen]
                    tail_rep = ':'.join([l+'='+str(tail.count(l)) for
                                         l in nucleotides])
                    polyA_file.write(t_strip + ' T '+tail_rep+'\n')
                    tcount +=1

            # Check for putative poly(A)-tail. Remove tail and write to file.
            if (seq[-min_length:] == min_A) or (seq[-rd:].count('A') >=rd-1)\
               or (seq[-2*rd:].count('A') >=2*rd-2):

                a_strip = strip_tailA(seq, trail_A, non_A, min_length, rd,
                                      min_A)
                striplen = len(a_strip)
                if striplen > 25:
                    seqlen = len(seq)
                    tail = seq[-(seqlen-striplen):]
                    tail_rep = ':'.join([l+'='+str(tail.count(l)) for
                                         l in nucleotides])
                    polyA_file.write(a_strip + ' A '+tail_rep+'\n')
                    acount += 1

    polyA_file.close()

    print('\n{0}: Nr. of total reads: {1}'.format(exp_id, total_reads))
    print('{0}: Nr. of readsReads with poly(A)-tail: {1}'.format(exp_id, acount))
    print('{0}: Nr. of readsReads with poly(T)-tail: {1}\n'.format(exp_id, tcount))

    return (polyA_path, total_reads, total_mapped_reads, acount, tcount)

def strip_tailT(seq, lead_T, non_T, min_length, rd, min_T):
    """
    Strip the polyT tail of a read iteratively.
    """

    if seq[:min_length] == min_T:
        # Remove all leading Ts
        seq = strip_tailT(lead_T.sub('', seq, count=1), lead_T, non_T,
                          min_length, rd, min_T)

    if seq[:rd].count('T') >=rd-1:
        # First remove the non-T character
        # Then remove all leading Ts
        seq = strip_tailT(lead_T.sub('', non_T.sub('', seq, count=1), count=1),
                           lead_T, non_T, min_length, rd, min_T)

    if seq[:2*rd].count('T') >=2*rd-2:
        # First remove the non-T characters
        seq = strip_tailT(lead_T.sub('', non_T.sub('', seq, count=2), count=1),
                           lead_T, non_T, min_length, rd, min_T)
        # Then remove all leading Ts
    return seq

def strip_tailA(seq, trail_A, non_A, min_length, rd, min_A):
    """
    Strip the polyA tail of a read iteratively. Strip things ending in a stretch
    of As and stript 1/min_length and 2/min_length*2
    """

    if seq[-min_length:] == min_A:
        # Remove all trailing As on reverse sequence; then reverse
        # again. rexexp only works from the beginning
        seq = strip_tailA(trail_A.sub('', seq[::-1], count=1)[::-1], trail_A,
                           non_A, min_length, rd, min_A)

    if seq[-rd:].count('A') >=rd-1:
        # First remove the non-A character
        # Then remove all trailing As
        seq = strip_tailA(trail_A.sub('', non_A.sub('', seq[::-1], count=1),
                           count=1)[::-1], trail_A, non_A, min_length, rd, min_A)

    if seq[-2*rd:].count('A') >=2*rd-2:
        # First remove the non-A characters
        # Then remove all trailing As
        seq = strip_tailA(trail_A.sub('', non_A.sub('', seq[::-1], count=2),
                           count=1)[::-1], trail_A, non_A, min_length, rd, min_A)

    return seq

if __name__ == '__main__':
    main()

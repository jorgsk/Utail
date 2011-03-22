"""
As a first approach, calculate the UTR-lengths of transcript that are annotated
to have only 1 utr.

Do this for different cellular compartments.

See if there is a difference in the different compartments.
"""
from __future__ import division
import os
import sys
import shutil
import numpy as np
import cPickle
import ConfigParser
from multiprocessing import Pool
from multiprocessing import cpu_count
from pprint import pprint

import re

def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

# only get the debug functio prop
if run_from_ipython():
    from IPython.Debugger import Tracer
    debug = Tracer()
else:
    def debug(): pass

from subprocess import Popen
from subprocess import PIPE
import time
import math
import matplotlib.pyplot as plt

# Your own imports
import annotation_analysis_progress as genome

class Settings(object):
    """Store the settings obtained from file"""
    def __init__(self, datasets, annotation_path, utrfile_provided, read_limit,
                 max_cores, chr1, hgfasta_path, polyA, min_utrlen, extendby,
                 del_reads):

        self.datasets = datasets
        self.annotation_path = annotation_path
        self.utrfile_provided = utrfile_provided
        self.read_limit = read_limit
        self.max_cores = max_cores
        self.chr1 = chr1
        self.hgfasta_path = hgfasta_path
        self.polyA = polyA
        self.min_utrlen = min_utrlen
        self.extendby = int(extendby)
        self.del_reads = del_reads

class Annotation(object):
    """A convenience class that holds the files and data structures relevant to
    the annotation that has been supplied."""

    def __init__(self, annotation_path):
        """Initiate as empty strings and fill in as they are being created."""
        self.path = annotation_path
        self.utr_bed_path = ''
        self.utr_bed_dict = ''
        self.annotated_polyA_path = ''
        self.annotated_polyA_dict = ''
        self.utr_seq_dict = ''

class FinalOutput(object):
    """This class takes care of writing output to file. Avoids passing a million
    parameters here and there."""

    def __init__(self, chrm, beg, end, strand, ts_ID, rpkm, extendby,
                 first_covr, sequence, polyA_cluster):

        # variables you have to initialize
        self.chrm = chrm
        # Note that beg and end might be subject to extensions! See the
        # nonextended versions further down
        self.beg = beg
        self.end = end
        self.ts_ID = ts_ID
        self.strand = strand
        self.rpkm = rpkm
        self.extendby = extendby # how far has this annotation been extended
        self.sequence = sequence
        self.polyA_cluster = polyA_cluster

        # variables that change depending on if reads map to this 3UTR
        self.has_PAS = 'NA'
        self.pas_pos ='NA'
        self.ext_mean_99 = 'NA'
        self.int_mean_99 = 'NA'
        self.ext_mean_annot = 'NA'
        self.int_mean_annot = 'NA'
        self.cuml_rel_size = 'NA'
        self.polyA_cumul_rel_pos = 'NA'
        self.polyA_support = 0

        # the coverage vectors
        self.covr_vector = [first_covr]
        self.cumul_covr = [first_covr]

        # Get the non-extended end as well!
        if extendby > 0:
            if strand == '+':
                self.end_nonextended = self.end - extendby
            if strand == '-':
                self.beg_nonextended = self.beg + extendby

    def header_dict(self):

        return dict((('chrm', self.chrm), ('beg', self.frmt(self.beg)),
                    ('end', self.frmt(self.end)),
                    ('ts_ID', self.frmt(self.ts_ID)),
                    ('strand', self.strand),
                    ('cuml_rel_size', self.frmt(self.cuml_rel_size)),
                    ('3utr_extended_by', self.frmt(self.extendby)),
                    ('mean_int_covrg_cuml_point', self.frmt(self.int_mean_99)),
                    ('mean_ext_covrg_cuml_point', self.frmt(self.ext_mean_99)),
                    ('mean_int_covrg_annotation', self.frmt(self.int_mean_annot)),
                    ('mean_ext_covrg_annotation', self.frmt(self.ext_mean_annot)),
                    ('polyA_support_for_cuml_point',self.frmt(self.polyA_support)),
                    ('cuml_point_for_polyA_site',self.frmt(self.polyA_cumul_rel_pos)),
                    ('PAS_type', self.has_PAS),
                    ('PAS_distance_from_cuml_point', self.frmt(self.pas_pos)),
                    ('RPKM', self.frmt(self.rpkm))
                     ))

    def frmt(self, element):
        """Return float objects with four decimals. Return all other objects as
        they were"""

        if type(element) is float:
            return format(element, '.4f')
        if type(element) is int:
            return str(element)
        else:
            return element

    def header_order(self):
        return """
        chrm
        beg
        end
        3utr_extended_by
        ts_ID
        strand
        cuml_rel_size
        mean_int_covrg_cuml_point
        mean_ext_covrg_cuml_point
        mean_int_covrg_annotation
        mean_ext_covrg_annotation
        polyA_support_for_cuml_point
        cuml_point_for_polyA_site
        PAS_type
        PAS_distance_from_cuml_point
        RPKM
        """.split()


    def is_empty(self):
        """Determine if non-extended coverage vector is empty"""

        if self.strand == '-':
            covr_slice = iter(self.covr_vector[self.extendby:])
        elif self.strand == '+':
            covr_slice = iter(self.covr_vector[:-self.extendby])

        # Check if covr_vector is non-empty in non-extended 3UTR
        for val in covr_slice:
            if val > 0:
                return False

        return True

    def calculate_output(self, pas_list, utrs, transl):
        # Don't do anything if coverage vector is empty
        if self.is_empty():
            return

        # Call different strand-specific functions
        if self.strand == '-':
            # calculate the cumulative values
            self.cumul_minus()
            # calculate has_polyA values
            self.polyA_minus(pas_list, utrs)

        if self.strand == '+':
            # calculate the cumulative values
            self.cumul_plus()
            # calculate has_polyA values
            self.polyA_plus(pas_list, utrs)

        # calculate the PAS and pas distance
        self.set_pas(pas_list, utrs)

    def cumul_minus(self):
        """Note that rel_pos is relative to the un-extended 3UTR. Thus, a
        rel_pos value of 3 here would be the same as self.beg + extension + 3."""
        covr_vector = self.covr_vector
        extendby = self.extendby
        cumul_covr = self.cumul_covr
        # get the cumulated coverage resulting from the extension
        ext_cumul = self.cumul_covr[extendby-1]
        # subtract this from extension
        if ext_cumul > 0:
            cumul_covr = [val-ext_cumul for val in cumul_covr[extendby:]]
        else:
            cumul_covr = cumul_covr[extendby:]

        # Get normalized cuml-coverage of non-extended 3UTR
        covr_sum = sum(covr_vector[extendby:])
        self.norm_cuml = [1-val/covr_sum for val in cumul_covr]

        # Test for a special case where only last entry has value
        if covr_sum == covr_vector[-1]:
            rel_pos = 1
            cuml_rel_size = rel_pos/float(self.end-self.beg)

        # Get the point where 99.5% of reads have landed
        for ind, el in enumerate(self.norm_cuml):
            if el < 0.995:
                rel_pos = ind
                length = float(self.end-self.beg)
                self.cuml_rel_size = (length-rel_pos)/length
                break

        # Save relative position with the object
        self.rel_pos = rel_pos

        # rel_pos according to extended 3utr
        ext_rel_pos = rel_pos + extendby

        # Then calculate the mean coverage on both sides of this.
        # Note for ext_mean_99: ext_rel_pos - extendby = ind
        self.ext_mean_99 = sum(covr_vector[ind:ext_rel_pos])/extendby
        self.int_mean_99 = sum(covr_vector[ext_rel_pos:ext_rel_pos +extendby])/extendby

        # Get the mean values 'extendby' around the annotated end too
        self.ext_mean_annot = sum(covr_vector[:extendby])/extendby
        self.int_mean_annot = sum(covr_vector[extendby: 2*extendby])/extendby

    def cumul_plus(self):
        covr_vector = self.covr_vector
        extendby = self.extendby
        # Get normalized cuml-coverage of un-extended 3UTR
        cumul_covr = self.cumul_covr[:-extendby]
        covr_sum = sum(covr_vector[:-extendby])
        self.norm_cuml = [val/covr_sum for val in cumul_covr]

        # Test special case where only first entry has value
        if covr_sum == covr_vector[0]:
            rel_pos = 1
            cuml_rel_size = rel_pos/float((self.end-self.beg))

        # Get the point where 99.5% of reads have landed
        for ind, el in enumerate(reversed(self.norm_cuml)):
            if el < 0.995:
                length = self.end-self.beg
                rel_pos = length - ind
                self.cuml_rel_size = rel_pos/float(length)
                break

        # Save relative position with the object
        # This relative position is 
        self.rel_pos = rel_pos

        # rel_pos according to extended 3utr
        ext_rel_pos = rel_pos - extendby

        # The calculate the mean coverage on both sides of this.
        # Note for ext_mean_99: ext_rel_pos + extendby = rel_pos
        self.ext_mean_99 = sum(covr_vector[ext_rel_pos:rel_pos])/extendby
        self.int_mean_99 = sum(covr_vector[ext_rel_pos-extendby:ext_rel_pos])/extendby

        # Get the mean values extendby around the annotated end too
        self.ext_mean_annot = sum(covr_vector[-extendby:])/extendby
        self.int_mean_annot = sum(covr_vector[-2*extendby:-extendby])/extendby


    def set_pas(self, pas_list, utrs):

        ts_ID = self.ts_ID

        # rel pos is different depending on strand because we are dealing with a
        # reverse-complemented sequence that has been extended
        if self.strand == '+':
            rel_pos = self.rel_pos
        if self.strand == '-':
            rel_pos = len(self.sequence) - self.extendby - self.rel_pos

        rel_seq = self.sequence[rel_pos-40:rel_pos+10]
        for pas in pas_list:
            try:
                (self.has_PAS, pas_indx) = (pas, rel_seq.index(pas))
                self.pas_pos = pas_indx - 40
                break
            except ValueError:
                (self.has_PAS, self.pas_pos) = ('NA', 'NA')


    def polyA_plus(self, pas_list, utrs):
        # Look for polyA-read support -- and get the cumulative
        # percentage at which the majority of the polyAreads are
        # found.

        # Get the absolute end-position of the 99.5%
        end_pos = self.beg + self.rel_pos
        # Check if this position has polyA reads within 50 nt
        self.polyA_cumul_rel_pos = 0
        for pAsite in self.polyA_cluster[0]:
            if (end_pos > pAsite-75) and (end_pos < pAsite+75):
                # get cumulative value at the polyA site
                rel_polyA = pAsite - self.beg
                # If the polyA is outside the end, report it as the last
                # value.
                if rel_polyA > self.end_nonextended - self.beg:
                    self.polyA_cumul_rel_pos = self.norm_cuml[-1]
                else:
                    self.polyA_cumul_rel_pos = self.norm_cuml[rel_polyA]

                # A general marker for polyA close to end_pos
                self.polyA_support = 1

                break

    def polyA_minus(self, pas_list, utrs):

        # Get the absolute end position of the 99.5%
        # self.rel_pos is relative the non-extended utr.
        end_pos = self.beg_nonextended + self.rel_pos
        # Check if this position has polyA reads within 50 nt
        for pAsite in self.polyA_cluster[0]:
            if (end_pos > pAsite-75) and (end_pos < pAsite+75):
                # Relative polyA site -- considering non-extended UTR!
                rel_polyA = pAsite - self.beg_nonextended
                # If the polyA is in the extension, report it as the last value
                if rel_polyA < 0:
                    self.polyA_cumul_rel_pos = self.norm_cuml[0]
                else:
                    self.polyA_cumul_rel_pos = self.norm_cuml[rel_polyA]

                # A general marker for polyA close to end_pos
                self.polyA_support = 1

                break

    def write_output(self, outobject):
        """Format the output as desired, then save"""

        # Create a list of formatted output
        output_order = self.header_order()
        output_dict = self.header_dict()

        output = [output_dict[hedr] for hedr in output_order]

        outobject.write('\t'.join(output) + '\n')


def zcat_wrapper(bed_reads, read_limit, out_path, polyA, polyA_path):
    """Function run by processes. Get only uniquely mapped reads with up to two
    mismatches. Convert gem-data into bed-format. If polyA is true, get poly(A)
    reads as well"""

    # File objects for writing to
    out_file = open(out_path, 'wb')
    polyA_file = open(polyA_path, 'wb')

    # Accept up to two mismatches. Make as set for speedup.
    acceptable = set(('1:0:0', '0:1:0', '0:0:1'))
    # Make regular expression for getting read-start
    start_re = re.compile('[0-9]*')
    trail_A = re.compile('A{5,}') # must be used on reversed string
    lead_T = re.compile('T{5,}')
    non_A = re.compile('[^A]')
    non_T = re.compile('[^T]')

    # A dictionary of strands
    getstrand = {'R':'-', 'F':'+'}

    # Keep track on the number of reads you get for the sake of RPKM
    total_reads = 0
    # Run zcat with -f to act as noram cat if the gem-file is not compressed
    cmd = ['zcat', '-f'] + bed_reads
    f = Popen(cmd, stdout=PIPE, stderr=PIPE)
    if read_limit:
        count = 0
    for map_line in f.stdout:
        (ID, seq, quality, mapinfo, position) = map_line.split('\t')

        if read_limit: # If short, only get up to 'limit' of reads
            count +=1
            if count > read_limit:
                break

        # Acceptable and poly(A) reads are mutually exclusive.
        if mapinfo in acceptable:
            # Get chromosome
            chrom, rest = position.split(':')
            # Get the strand
            strand = getstrand[rest[0]]
            # Get read beg
            beg = start_re.match(rest[1:]).group()

            # Write to file
            out_file.write('\t'.join([chrom, beg, str(int(beg)+len(seq)), '.',
                                      '.', strand]) + '\n')
            total_reads = total_reads + 1

        # When looking for poly(A) reads, filter by non-uniquely mapped reads
        if polyA == True:
            if mapinfo[:5] == '0:0:0':
                if seq.startswith('NN'):
                    seq = seq[2:]

                # If more than two ambigious, discard.
                if seq.count('N') > 2:
                    continue

                # Check for putative poly(A)-tail. Remove tail and write to file.
                if (seq[-5:] == 'AAAAA') or (seq[-7:].count('A') >=6):
                    polyA_file.write(strip_tailA(seq, trail_A, non_A) + '\n')

                if (seq[:5] == 'TTTTT') or (seq[:7].count('T') >=6):
                    polyA_file.write(strip_tailT(seq, lead_T, non_T) + '\n')

    out_file.close()
    polyA_file.close()

    return total_reads

def strip_tailA(seq, trail_A, non_A):

    if seq[-5:] == 'AAAAA':
        # Remove all trailing As on reverse sequence; then reverse
        # again
        seq = strip_tailA(trail_A.sub('', seq[::-1], count=1)[::-1], trail_A,
                           non_A)

    if seq[-7:].count('A') >=6:
        # First remove the non-A character
        # Then remove all trailing As
        seq = strip_tailA(trail_A.sub('', non_A.sub('', seq[::-1], count=1),
                           count=1)[::-1], trail_A, non_A)

    return seq

def strip_tailT(seq, lead_T, non_T):

    if seq[:5] == 'TTTTT':
        # Remove all leading Ts
        seq = strip_tailT(lead_T.sub('', seq, count=1), lead_T, non_T)

    if seq[:7].count('T') >=6:
        # First remove the non-T character
        # Then remove all leading Ts
        seq = strip_tailT(lead_T.sub('', non_T.sub('', seq, count=1), count=1),
                           lead_T, non_T)

    return seq

def coverage_wrapper(dset_id, filtered_reads, utrfile_path, options):
    """ Do intersections of a_files and b_files. The output will be in terms of
    the lines in the b_file and how much each line is covered by lines in the
    a_file.
    """

    # The gold standard for renaming paths :)
    (dirpath, basename) = os.path.split(filtered_reads)
    out_path = os.path.join(dirpath, 'covered_'+dset_id)
    outfile = open(out_path, 'wb')

    cmd = ['coverageBed', options, '-a', filtered_reads, '-b', utrfile_path]

    f = Popen(cmd, shell = False, stdout=outfile)
    f.wait()

    outfile.close()

    return out_path

def get_pas_list(inverted=False):
    # The known PAS in order of common appearance
    known_pas = ['AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA',
                 'AATACA', 'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA',
                 'AATAGA']
    if inverted:
        transl = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
        return [''.join([transl[letr] for letr in reversed(pas)]) for pas in
                  known_pas]
    else:
        return known_pas


def ends_wrapper(dset_id, coverage, utrs, utr_seqs, rpkm, extendby,
                 polyA_cluster):
    """Putting together all the info on the 3UTRs and writing it to file"""

    # Renaming paths for the output file
    (dirpath, basename) = os.path.split(coverage)
    out_path = os.path.join(dirpath, 'utr_'+dset_id)
    outfile = open(out_path, 'wb')

    pas_list = get_pas_list(inverted = False)
    transl = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}

    # First get line 1 as something to start from
    read_from = open(coverage, 'rb')
    line1 = read_from.next()
    (chrm, beg, end, ts_ID, d, strand, rel_pos, covr) = line1.split()

    # Create an output-instance
    this_utr = FinalOutput(chrm, int(beg), int(end), strand, ts_ID, rpkm[ts_ID],
                           extendby, int(covr), utr_seqs[ts_ID],
                           polyA_cluster[ts_ID])

    # Write the header of the file
    outfile.write('\t'.join(this_utr.header_order()) + '\n')

    # Assert that the file information is the same as you started with
    assert utrs[ts_ID] == (chrm, int(beg), int(end), strand), 'Mismatch'

     # staring reading from line nr 2
    for line in read_from:
        (chrm, beg, end, ts_ID, d, strand, pos, covr) = line.split('\t')

        # If ts_ID remains the same, we're on the same 3UTR
        if ts_ID == this_utr.ts_ID:
            # Note: we also add coverage from extension; this is corrected for later
            covr = int(covr)
            this_utr.covr_vector.append(covr) # add coverage and continue
            this_utr.cumul_covr.append(this_utr.cumul_covr[-1] + covr)

        # Ooops, this is the first line of a new UTR
        else:
            # Calculate output values like cumulative coverage, polyA sites, etc..
            this_utr.calculate_output(pas_list, utrs, transl)
            # Save output to file
            this_utr.write_output(outfile)

            # Update to the new utr and start the loop from scratch
            this_utr = FinalOutput(chrm, int(beg), int(end), strand, ts_ID,
                                   rpkm[ts_ID], extendby, int(covr),
                                   utr_seqs[ts_ID], polyA_cluster[ts_ID])

            # Assert that the next utr has correct info
            assert utrs[ts_ID] == (chrm, int(beg), int(end), strand), 'Mismatch'

    outfile.close()

    return out_path

def verify_access(f):
    try:
        open(f, 'rb')
    except:
        print('Could not access {0}'.format(f))
        sys.exit()

def read_settings(settings_file):
    """Read settings file for parameters needed to run the script."""

    conf = ConfigParser.ConfigParser()
    conf.read(settings_file)

    expected_fields = ['DATASETS', 'ANNOTATION', 'CPU_CORES', 'RESTRICT_READS',
                       'CHROMOSOME1', 'SUPPLIED_3UTR_BEDFILE', 'HG_FASTA',
                       'POLYA_READS', 'MIN_3UTR_LENGTH', 'EXTEND', 'PLOTTING',
                       'DELETE_READS']

    missing = set(conf.sections()) - set(expected_fields)
    if len(missing) == 0:
        pass
    else:
        print('The following options sections are missing: {}'.format(missing))

    # datasets and annotation
    datasets = dict((dset, files.split(':')) for dset, files in conf.items('DATASETS'))
    annotation = conf.get('ANNOTATION', 'annotation')
    # check if the files are actually there...
    for dset, files in datasets.items():
        [verify_access(f) for f in files]

    # set minimum length of 3 utr
    try:
        min_utrlen = conf.getint('MIN_3UTR_LENGTH', '3utrlen')
    except ValueError:
        min_utrlen = 100

    # cpu cores
    try:
        max_cores = conf.getint('CPU_CORES', 'max_cores')
    except ValueError:
        max_cores = cpu_count()-1

    # restrict number of reads from source
    try:
        read_limit = conf.getint('RESTRICT_READS', 'restrict_reads')
    except ValueError:
        read_limit = conf.getboolean('RESTRICT_READS', 'restrict_reads')

    # supplied 3utr bedfile
    try:
        utr_bedfile_path = conf.getboolean('SUPPLIED_3UTR_BEDFILE',
                                           'utr_bed_path')
    except ValueError:
        utr_bedfile_path = conf.get('SUPPLIED_3UTR_BEDFILE', 'utr_bed_path')
        verify_access(utr_bedfile_path) # Check if is a file

    # restrict to chromosome 1
    chr1 = conf.getboolean('CHROMOSOME1', 'chr1')

    # human genome fasta file
    hg_fasta = conf.get('HG_FASTA', 'hg_fasta')

    # by how much should the 3utrs be extended?
    extendby = conf.get('EXTEND', 'extend_by')

    # should you delete reads?
    del_reads = conf.getboolean('DELETE_READS', 'delete_reads')

    # poly(A) reads -- get them / dont get them / supply them
    try:
        polyA = conf.getboolean('POLYA_READS', 'polya')
    except ValueError:
        polyA = conf.get('POLYA_READS', 'polya')
        verify_access(polyA) # Check if is a file

    return(datasets, annotation, utr_bedfile_path, read_limit, max_cores, chr1,
          hg_fasta, polyA, min_utrlen, extendby, del_reads)


def get_chromosome1(annotation, beddir):
    (basepath, extension) = os.path.splitext(annotation)
    (dirpath, basename) = os.path.split(basepath)
    out_path = os.path.join(beddir, basename) + '_chr1' + extension
    # If the file already exists, don't make it again
    if os.path.isfile(out_path):
        return out_path

    # Make the directory if it doesn't exist
    if not os.path.exists(beddir):
        os.makedirs(beddir)

    # The file to write to
    out_file = open(out_path, 'wb')

    t1 = time.time()
    print('Extracting chromosome 1 from the annotation ...')
    for line in open(annotation, 'rb'):
        if line.split('\t')[0] == 'chr1':
            out_file.write(line)
    print('\tTime taken to get chromsome 1: {0}\n'.format(time.time()-t1))

    return out_path

def get_polyA_sites_path(settings, beddir):
    """Get polyA_sites_file path from annotation via genome module"""

    # If options are not set, make them a 0-string define options
    if settings.chr1:
        chrm = '_chr1'
    else:
        chrm = ''

    if settings.utrfile_provided:
        user_provided = '_user_provided'
    else:
        user_provided = ''

    # utr path
    polyA_sites_filename = 'polyA_sites' + chrm + user_provided + '.bed'
    polyA_site_bed_path = os.path.join(beddir, polyA_sites_filename)

    ## If the file already exists, don't make it again
    if os.path.isfile(polyA_site_bed_path):
        return polyA_site_bed_path

    # If a utrfile is provided, never re-make from annotation
    if settings.utrfile_provided:
        polyA_site_bed_path = shape_provided_bed(polyA_site_bed_path, settings)

    # If utrfile is not provided, get it yourself from a provided annotation
    if not settings.utrfile_provided:
        if settings.chr1:
            # path or dict?
            annotation_path = get_chr1_annotation(settings.annotation_path, beddir)

        t1 = time.time()
        print('polyA sites bedfile not found. Generating from annotation ...')
        reload(genome)
        genome.get_polyA_sites_bed(annotation_path, polyA_site_bed_path, settings)

        print('\tTime taken to generate polyA-bedfile: {0}\n'\
              .format(time.time()-t1))

    return polyA_site_bed_path

def get_utr_path(settings, beddir):
    """Get 3utr-file path from annotation via genome module"""

    # If options are not set, make them a 0-string define options
    if settings.extendby:
        ext = '_extendby_'+str(settings.extendby)
    else:
        ext = ''

    if settings.chr1:
        chrm = '_chr1'
    else:
        chrm = ''

    if settings.utrfile_provided:
        user_provided = '_user_provided'
    else:
        user_provided = ''

    # utr path
    utr_filename = '3utrs' + chrm + ext + user_provided + '.bed'
    utr_bed_path = os.path.join(beddir, utr_filename)

    # If the file already exists, don't make it again
    if os.path.isfile(utr_bed_path):
        return utr_bed_path

    # If a utrfile is provided, never re-make from annotation
    if settings.utrfile_provided:
        utr_bed_path = shape_provided_bed(utr_bed_path, settings)

    # If utrfile is not provided, get it yourself from a provided annotation
    if not settings.utrfile_provided:
        if settings.chr1:
            # path or dict?
            annotation_path = get_chr1_annotation(settings.annotation, beddir)

        t1 = time.time()
        print('3UTR-bedfile not found. Generating from annotation ...')
        genome.get_3utr_bed(annotation_path, utr_bed_path, settings)

        print('\tTime taken to generate 3UTR-bedfile: {0}\n'\
              .format(time.time()-t1))

    return utr_bed_path

def get_chr1_annotation(annotation, beddir):
    """Split the original (GENCODE) annotation into one with only chrm1 """

    (name, suffix) = os.path.splitext(os.path.basename(annotation))
    filename = name + '_chr1' + suffix
    outpath = os.path.join(beddir, filename)

    # If the file already exists, don't make it again
    if os.path.isfile(outpath):
        return outpath


    t1 = time.time()
    print('Separating chr1 from the annotation ...')
    outhandle = open(outpath, 'wb')

    for line in open(annotation, 'rd'):
        if line[:5] == 'chr1\t':
            outhandle.write(line)

    outhandle.close()

    print('\tTime taken to separate chr1: {0}\n'.format(time.time()-t1))

    return outpath

def shape_provided_bed(utr_bed_path, settings):
    """Go through provided bedfile and shape it according to settings. Save in
    local directory of this script"""

    outfile = open(utr_bed_path, 'wb')

    for line in open(settings.utrfile_provided, 'rb'):
        (chrm, beg, end, name, val, strand) = line.split()[:6]

        end = int(end)
        beg = int(beg)

        # Skip lines that are not chrm1 if option is set
        if settings.chr1:
            if chrm not in ['chr1', 'Chr1']:
                continue

        # Skip short utrs
        if end - beg < settings.utrlen:
            continue

        # Extend by the nr of nucleotides specified
        if settings.extendby:
            if strand == '+':
                end = end + settings.extendby
            if strand == '-':
                beg = beg - settings.extendby

        outfile.write('\t'.join([chrm, str(beg), str(end), name, val,
                                 strand])+'\n')
    outfile.close()

    return utr_bed_path

def save_output(final_dict, output_dir):
    """Copy files in out_dict to the output folder"""

    for ID, final_path in final_dict.items():
        out_path = os.path.join(output_dir, os.path.basename(final_path))
        shutil.copyfile(final_path, out_path)

def get_rpkm(reads, utrfile_path, total_reads, utrs, extendby):
    """Run coverageBed of reads on provided 3UTR and get RPKM. Correct for any
    extension made to the 3UTRs"""
    # Information needed for rpkm:
    # TOTAL reads
    # Length of UTR
    # # of reads landing in the UTR.
    rpkm = {}
    # If the bed-file has been extended, we need to unextend it.
    if extendby:
        temp_bed_path = os.path.join(os.path.dirname(utrfile_path), 'temp_bed')
        temp_handle = open(temp_bed_path, 'wb')
        for line in open(utrfile_path, 'rb'):

            (chrm, beg, end, ts_id, d, strand) = line.split()

            # Reverse the extensions so you get correct RPKM!
            if strand == '+':
                end = int(end) - extendby
            if strand == '-':
                beg = int(beg) + extendby

            temp_handle.write('\t'.join([chrm, str(beg), str(end), ts_id, d,
                                         strand]) + '\n')
        temp_handle.close()

        utrfile_path = temp_bed_path

    p = Popen(['coverageBed', '-a', reads, '-b', utrfile_path], shell=False,
                   stdout=PIPE)

    for line in p.stdout:
        (chrm, bl, bl, ts_id, d, d, reads_covering) = line.split('\t')[:7]
        utr_length = utrs[ts_id][2] - utrs[ts_id][1]
        rpkm[ts_id] = ((10**9)*int(reads_covering))/(total_reads*utr_length)

    return rpkm

def process_reads(pA_reads_path):
    """ Remove reads that are too short or have a poor nucleotide composition.
    Re
    """
    processed_reads = os.path.splitext(pA_reads_path)[0]+'_processed.fas'
    outfile = open(processed_reads, 'wb')

    # Go through the file two lines at a time. If the next line does not begin
    # with '>', append line to last entry (it's a split fasta-file).

    for line in open(pA_reads_path, 'rb'):
        # add lines until there are 2 entries in linepair
        seq = line.rstrip()
        seqlen = len(seq)
        if seqlen > 25:
            As = seq.count('A')
            Ts = seq.count('T')
            # only save if A/T frequencies are not abnormal
            if (As/seqlen < 0.45) and (Ts/seqlen < 0.45):
                # trim the title
                outfile.write('>read\n'+line)

    outfile.close()

    return processed_reads

def map_reads(processed_reads):
    """ Map the processed reads using gem-mapper """
    ## Andrea's GEM index of hg19
    g_ind='/users/rg/atanzer/DATA/GEM_indices/Genomes/H.sapiens.genome.hg19.main'
    mapped_reads = os.path.splitext(processed_reads)[0]+'_mapped'

    # Naming the final output
    base_dir = os.path.dirname(os.path.split(mapped_reads)[0])
    polybed_path = os.path.splitext(processed_reads)[0] + '_mapped.bed'

    # mapping trimmed reads
    command = "gem-mapper -I {0} -i {1} -o {2} -q ignore -m 1"\
            .format(g_ind, processed_reads, mapped_reads)

    p = Popen(command.split())
    p.wait()

    # Accept up to one mismatch. Make as set for speedup.
    acceptable = set(('1:0', '0:1'))
    getstrand = {'R':'-', 'F':'+'}
    start_re = re.compile('[0-9]*')

    reads_file = open(polybed_path, 'wb')

    for line in open(mapped_reads + '.0.map', 'rb'):
        (ID, seq, mapinfo, position) = line.split('\t')

        # Acceptable and poly(A) reads are mutually exclusive.
        if mapinfo in acceptable:
            # Get chromosome, strand, and beg
            (chrom, rest) = position.split(':')
            strand = getstrand[rest[0]]
            beg = start_re.match(rest[1:]).group()

            # Write to file in .bed format
            reads_file.write('\t'.join([chrom, beg, str(int(beg)+len(seq)), '.',
                                      '0', strand]) + '\n')

    return polybed_path

def remap_polyAs(pA_reads_path):
    """Get fasta poly(A) reads, remap them, and return the uniquely mapped.
    However, check if the poly(A) reads have been obtained for this dataset
    before. In that case, just return the path, and use the previously obtained
    results."""

    # Check if a ...unique.bed file already exists. In that case, return.

    # Process reads by removing those with low-quality, removing the leading Ts
    # and/OR trailing As.
    processed_reads = process_reads(pA_reads_path)

    # Map the surviving reads to the genome and return unique ones
    unique_reads = map_reads(processed_reads)

    return unique_reads

def get_bed_reads(dset_reads, dset_id, read_limit, tempdir, polyA):
    """ Determine datatype according to file-suffix and act accordingly. """

    # Path of .bed output
    out_path = os.path.join(tempdir, 'reads_'+dset_id+'.bed')
    # Define the path of the polyA file
    polyA_path = os.path.join(tempdir, 'polyA_reads_'+dset_id+'.fa')

    # allowed suffixes:
    # [gem, map, gz]
    # [gem, map]
    # [bed]
    # [bed.gz]
    ok_sufx = ['gem', 'map', 'gz', 'bed']

    nr_files = len(dset_reads)
    # more robust way of building a suffix. Start with a '.' separated list of
    # the file name. Proceed backwards, adding to reverese_suffix if the entry
    # is in the allowed group.
    dotsplit = os.path.basename(dset_reads[0]).split('.')
    suflist = [sufpart for sufpart in reversed(dotsplit) if sufpart in ok_sufx]

    suffix = '.'.join(reversed(suflist))

    # If in gem-format, go through the file with zcat -f
    if suffix in ['gem.map.gz', 'gem.map']:
        print('Obtaining reads from mapping for {0} ...\n'.format(dset_id))
        # Get only the uniquely mapped reads (up to 2 mismatches)
        total_reads = zcat_wrapper(dset_reads, read_limit, out_path, polyA,
                                   polyA_path)

    # If in bed-format, also run zcat -f on all the files to make one big
    # bedfile. How to restrict reads in this case?? No restriction?
    elif suffix in ['bed.gz', 'bed']:
        total_reads = concat_bedfiles(dset_reads, out_path, polyA, polyA_path)

    else:
        print('Non-valid suffix: {0}. Allowed suffixes are .gem.map.gz,\
              .gem.map, and .gem.map.gz'.format(suffix))
        sys.exit()

    return (out_path, polyA_path, total_reads)

def concat_bedfiles(dset_reads, out_path, polyA, polyA_path):

    # File objects for writing to
    out_file = open(out_path, 'wb')
    polyA_file = open(polyA_path, 'wb')

    cmd = ['zcat', '-f'] + dset_reads
    f = Popen(cmd, stdout=out_file, stderr=PIPE)
    f.wait()

    # If polyA is a path to a bedfile (or bedfiles) concatenate these too
    if type(polyA) == str:
        cmd = ['zcat', '-f'] + polyA
        f = Popen(cmd, stdout=polyA_file, stderr=PIPE)
        f.wait()

    # Return the number of line numbers (= number of reads)
    return sum(1 for line in open(out_file, 'rb'))

def get_polyAutr(polyAbed, utrfile_path):
    """ Intersect poly(A) reads with 3UTR file. Return dictionary where each
    3UTR is listed with its poly(A) sites."""
    # A dictionary to hold the ts_ID -> poly(A)-reads relation
    utr_polyAs = {}

    cmd = ['intersectBed', '-wb', '-a', polyAbed, '-b', utrfile_path]

    # Run the above command -- outside the shell -- and loop through output
    f = Popen(cmd, stdout=PIPE)
    for line in f.stdout:
        (polyA, ts) = (line.split()[:6], line.split()[6:])
        ts_id = ts[3]
        if not ts_id in utr_polyAs:
            utr_polyAs[ts_id] = [tuple(polyA)]
        else:
            utr_polyAs[ts_id].append(tuple(polyA))

    return utr_polyAs

def cluster_loop(ends):
    clustsum = 0
    clustcount = 0
    this_cluster = []
    clusters = []
    for indx, val in enumerate(ends):
        ival = int(val)

        clustsum = clustsum + ival
        clustcount += 1
        mean = clustsum/clustcount

        # If dist between new entry and cluster mean is < 20, keep in cluster
        if abs(ival - mean) < 20:
            this_cluster.append(ival)
        else: # If not, start a new cluster, and save the old one
            clusters.append(this_cluster)
            clustsum = ival
            clustcount = 1
            this_cluster = [ival]

    # Append the last cluster
    clusters.append(this_cluster)
    # Get only the clusters with length more than one
    res_clus = [clus for clus in clusters if len(clus) > 1]
    # Get the mean of the clusters with length greater than one
    mixed_cluster = [int(math.floor(sum(clus)/len(clus))) for clus in res_clus]

    # Return the clusters and their average
    return (mixed_cluster, res_clus)

def cluster_polyAs(utr_polyAs, utrs):
    """Cluster the poly(A) reads for each ts_id. Choose the pTTS site depending
    on the strand of the ts."""
    plus_values = {'+': [], '-':[]}
    minus_values = {'+': [], '-':[]}

    polyA_cluster = {}

    for ts_id, polyAs in utr_polyAs.iteritems():
        ts_info = utrs[ts_id]
        real_strand = ts_info[3]

        # Getting statistics for which strand the polyA reads map to
        nr_polyA_reads = len(polyAs)
        plus_ratio = sum(1 for tup in polyAs if tup[5] == '+')/nr_polyA_reads
        minus_ratio = sum(1 for tup in polyAs if tup[5] == '-')/nr_polyA_reads

        if real_strand == '+':
            plus_values['+'].append(plus_ratio)
            plus_values['-'].append(minus_ratio)
        if real_strand == '-':
            minus_values['+'].append(plus_ratio)
            minus_values['-'].append(minus_ratio)

        # For + strands, choose beg as polyA_site; choose end for - strands.
        if real_strand == '+':
            #ends = [tup[1] for tup in polyAs]
            ends = sorted([tup[1] for tup in polyAs if tup[5] == '-'])
        if real_strand == '-':
            #ends = [tup[2] for tup in polyAs]
            ends = sorted([tup[2] for tup in polyAs if tup[5] == '+'])

        # Getting the actual clusters
        # TODO you can get statistics from these clusters to make a
        # distance-distribution figure.
        # This figure shows that your method is sound because it looks like
        # all the other polyA-distribution figures.
        polyA_cluster[ts_id] = cluster_loop(ends)

    # For those ts_id that don't have a cluster, give them an empty list; ad hoc
    for ts_id in utrs:
        if ts_id not in polyA_cluster:
            polyA_cluster[ts_id] = [[],[]]

    # Statistics on the ratio of + and - mapped reads for the genes that are in
    # the positive or negative strand. RESULT: for paired end reads, the polyA
    # reads map to the OPPOSITE strand 90% of the times.

    plus_avrg = {'+': sum(plus_values['+'])/len(plus_values['+']),
                 '-': sum(plus_values['-'])/len(plus_values['-'])}

    minus_avrg = {'+': sum(minus_values['+'])/len(minus_values['+']),
                 '-': sum(minus_values['-'])/len(minus_values['-'])}

    return polyA_cluster


def pipeline(dset_id, dset_reads, tempdir, output_dir, utr_seqs, settings,
             annotation):
    """Get reads, get polyA reads, cluster polyA reads, get coverage, do
    calculations, write to output... this is where it all happens"""

    t0 = time.time()
    # Set some settings for brievity
    read_limit = settings.read_limit
    polyA = settings.polyA
    # Get the reads in bed format. And the polyA reads if this option is set.
    # As well get the total number of reads for getting RPKM later.
    (bed_reads, polyA_bed, total_nr_reads) = get_bed_reads(dset_reads, dset_id,
                                                      read_limit, tempdir,
                                                      polyA)

    # If polyA is a string, it is a string to polyA sequences in bed-format
    if type(polyA) == str:
        polyA_bed_path = polyA
        utr_polyAs = get_polyAutr(polyA_bed_path, annotation.utrfile_path)

    elif polyA == True:
        # Do the polyA: extract, trim, remap, and -> .bed-format.
        # However, if a poly(A).bed file exists for this dset, don't redo the
        # job.
        print('Processing poly(A) reads for {0}...'.format(dset_id))
        polyA_bed_path = remap_polyAs(polyA_bed)

        # If polyA, find where the polyAs overlap with the 3UTR in the first place
        # run overlap-bed or something similar
        # To verify the correctness of your polyA reads, to a quick check on the
        # genome browser. If you have no overlaps, then something is wrong.
        utr_polyAs = get_polyAutr(polyA_bed_path, settings.utrfile_path)

    # Cluster the poly(A) reads for each ts_id
    polyA_cluster = cluster_polyAs(utr_polyAs, annotation.utrs)

    #debug()
    # Get the RPKM by running intersect-bed of the reads on the 3utr
    print('Obtaining RPKM for {0} ...\n'.format(dset_id))
    rpkm = get_rpkm(bed_reads, annotation.utrfile_path, total_nr_reads,
                    annotation.utrs, settings.extendby)

    # Cover the annotated 3UTR with the reads
    print('Getting read-coverage for the 3UTRs for {0} ...\n'.format(dset_id))
    options = '-d'
    coverage = coverage_wrapper(dset_id, bed_reads, annotation.utrfile_path, options)
    # Add coverage file_locations to dict

    print('Using coverage to get 3UTR ends for {0} ...\n'.format(dset_id))
    # Get transcript 3utr endings as determined by read coverage
    utr_ends = ends_wrapper(dset_id, coverage, annotation.utrs, utr_seqs, rpkm,
                            settings.extendby, polyA_cluster)

    print('Total time for {0}: {1}\n'.format(dset_id, time.time() - t0))

    # If requested, delete the reads in .bed format after use
    if settings.del_reads:
        os.remove(bed_reads)

    # Return a list with the salient file paths of output files
    return {dset_id: {'coverage': coverage, 'utr_ends': utr_ends}}

def remove_directory(some_dir):
    try:
        shutil.rmtree(some_dir)
    except OSError:
        print('Unable to remove {0}'.format(some_dir))


def output_analyzer(final_output, utrs, utrs_path):
    """Present graphically the different lengths of the different tissues"""

    ##########
    # Restrict final output according to cell_line
    filtered_output = dict((key, final_output[key]) for key in final_output if
                        'K562' in key)
    #filtered_output = dict((key, final_output[key]) for key in final_output if
                        #'HeLa' in key)
    #filtered_output = dict((key, final_output[key]) for key in final_output if
                        #'GM12878' in key)
    #filtered_output = final_output

    # Do basic statistics and make a box-plot
    statistics_boxplot(filtered_output, utrs, utrs_path)

def rel_len_variation(measure_dict, dset_indx):
    # Output: histograms of differences

    strict_rel_len = {}
    without_zeros = 0
    for its_id, el_sizes in measure_dict.items():
        if 0 in el_sizes:
            continue
        else:
            without_zeros += 1
            # Calculate distances between all rel_dist
            parw_dist = [abs(val1-val2) for val1 in el_sizes for val2 in
                         el_sizes]
            sig_dist = ['yes' for val in parw_dist if val > 0.3]
            if 'yes' in sig_dist:
                strict_rel_len[its_id] = el_sizes

    pprint(without_zeros)
    pprint(len(strict_rel_len))


def statistics_boxplot(filtered_output, utrs, utrs_path):
    """Do basic statistics and make boxplots """

    # Make a relationship between dset and index for future use
    dset_indx = dict((dset, indx) for indx, dset in enumerate(filtered_output))
    exp_nr = len(filtered_output)

    # Prepare a dictionary of its_id with relative lenghts of compartments
    cuml_99_dict = dict((ts_id, [0 for val in range(exp_nr)]) for ts_id in utrs)

    # Get the ids of the transcripts whose extensions do not overlap exons in
    # the gencode annotation. (it's a set() by default)
    exons = '/users/rg/jskancke/phdproject/3UTR/gencode5/just_exons/just_exons.bed'

    rpkm_threshold = 1
    for dset, path in filtered_output.items():
        for line in open(path, 'rb'):
            # Skip lines that don't start with chr (it's fast enough)
            if line[:3] != 'chr':
                continue
            (chrm, beg, end, ts_ID, strand, cuml99, int_mean_99, ext_mean_99,
             int_mean_annot, ext_mean_annot, PAS_type,\
             pas_distance, rpkm) = line.split()
            rpkm = float(rpkm)
            # Get 3utrs with rpkm above a set threshold
            if rpkm > rpkm_threshold:
                cuml_99_dict[ts_ID][dset_indx[dset]] = float(cuml99)

    # Get how many of the transcripts are 0 or > 0.999 in rel_len.
    # count 0, > 99
    one_or_zero = [[0,0], [0,0], [0,0]]
    zeros = [0,0,0,0]

    for ts_ID, rel_values in cuml_99_dict.items():
        if rel_values.count(0) == 3:
            zeros[3] += 1
        if rel_values.count(0) == 2:
            zeros[2] += 1
        if rel_values.count(0) == 1:
            zeros[1] += 1
        if rel_values.count(0) == 0:
            zeros[0] += 1

        for ind, val in enumerate(rel_values):
            if val == 0:
                one_or_zero[ind][0] += 1
            if val > 0.95:
                one_or_zero[ind][1] += 1

    print('One or zero: {0}'.format(one_or_zero))
    print('Zeros: {0}').format(zeros)

    # Find these things [0.3, 0.5, 0.8] as opposed to [0.3, 0.3, 0.3] and [0.8,
    # 0.8, 0.8]
    print("Total transcripts: {0}".format(len(utrs)))

    #rel_len_variation(cuml_99_dict, dset_indx)

    # Do some statistics and boxpolot for your measure of choice
    measure_dict = cuml_99_dict
    # Boxplot needs[[3,4,....,4], [4,3,...3]]
    # Average normalized_len for each condition
    for_boxplot = [[] for val in range(exp_nr)]
    avr_len = dict((dset, 0) for dset in dset_indx)

    # For boxplot, get only the utrs with expression in all tissues
    for dset, dset_ind in dset_indx.items():
        for_boxplot[dset_ind] = [val[dset_ind] for (key, val) in
                                 measure_dict.items() if 0 not in val]

        dset_vals = [val[dset_ind] for (key, val) in measure_dict.items()]
        avr_len[dset] = (np.mean(dset_vals), np.std(dset_vals))

    print('For boxplot: {0}'.format(len(for_boxplot[0])))

    # Split the utrs into size categories, since you are doing relative size
    # increases.
    ranges = [(0,100), (101, 500), (501, 1000), (1001, 1500), (1501, 2000),
              (2501, 3000), (3001, 3500), (3501, 4000), (4001, 20000)]

    # Make boxplot
    boxplot(for_boxplot, dset_indx.keys())


def boxplot(arrays, dsets):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.boxplot(arrays)
    ax.set_xticklabels(dsets)
    ax.set_ylabel('Relative 3UTR length', size=20)
    ax.set_title(dsets[0].split('_')[0])
    ax.set_ylim(0,1.2)
    fig.show()


def make_utr_plots(final_output, utrs, coverage):
    """ """
    # Choose a ts_ID from final_output where there is a clear difference between
    # the different tissues (one has 0.4, one has 0.6, one has 0.8 rel_len)
    # from the final_output files. Also make sure they have reasonably high
    # RPKMs. Then look through the coverage files for these ts_ids, and
    # calculate the average coverage per 10 reads or so and produce a plot of
    # this.

    # Select the transcripts that have rpkm over 2
    exp_nr = len(final_output)
    rel_len_dict = dict((ts_id, [0]*exp_nr) for ts_id in utrs)
    index = 0
    count = 0
    for dset, path in final_output.items():
        # Skip the first line starting with '#'
        for line in open(path, 'rb'):
            # Skip lines that don't start with chr
            if line[:3] != 'chr':
                continue
            (chrm, beg, end, ts_ID, last_utr_read, strand, rel_size, mean_cvr,
             std_cvr, PAS_type, pas_distance, rpkm) = line.split()
            rpkm = float(rpkm)
            if rpkm > 2:
                count = count+1
                rel_len_dict[ts_ID][index] = float(rel_size)
        index = index + 1

    print count

    # Keep only the transcripts that 1) have no 0 rel_lens AND 2) have a minimum
    # distance of 0.4 btween two of the rel_lens.
    strict_rel_len = {}
    for its_id, el_sizes in rel_len_dict.items():
        if 0 in el_sizes:
            continue
        else:
            # Calculate distances between all rel_dist
            parw_dist = [abs(val1-val2) for val1 in el_sizes for val2 in
                         el_sizes]
            sig_dist = ['yes' for val in parw_dist if val > 0.4]
            if 'yes' in sig_dist:
                strict_rel_len[its_id] = el_sizes

    print len(strict_rel_len)

    # If more than 5 subjects, make a random sample of the ones in the set
    if len(strict_rel_len) > 1:
        strict_rel_len = dict((rankey, strict_rel_len[rankey]) for rankey in
                              strict_rel_len.keys()[:3])

    # For each element in strict_rel_len, do a graphical presentation of the
    # read distribution
    # Method: go through each coverage file and grep the coverage into a list.
    # Make averages over that list. Print the output.
    key_ids = strict_rel_len.keys()
    covrg_dict = dict((key, [[] for bla in range(exp_nr)]) for key in key_ids)
    indx=0
    for dset, cvrg_path in coverage.items():
        for line in open(cvrg_path):
            (chrm, beg, end, ts_id, d, strnd, rel_pos, covrg) = line.split()
            if ts_id in key_ids:
                covrg_dict[ts_id][indx].append(int(covrg))
        indx = indx+1

    for ts_id, covregs in covrg_dict.items():
        plot_3utr_ends(ts_id, covregs, utrs[ts_id][3])

def plot_3utr_ends(ts_id, coverages, strand):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    new_sries0 = moving_average(coverages[0])
    new_sries1 = moving_average(coverages[1])
    ax.plot(moving_average(coverages[0]))
    ax.plot(moving_average(coverages[1]))
    fig.show()
    pass

def moving_average(dataseries):
    # Divide the dataseries into blocks of size 10nt.
    newseries = []
    for pos, val in enumerate(dataseries):
        newseries.append(sum(dataseries[pos-3:pos+3])/6)

    return newseries

    # Divide the dataseries into averages of a fixed length, so that all utrs
    # can be compared.
    pass

def get_utrdict(utrfile_path):
    utr_dict = {}
    for line in open(utrfile_path, 'rb'):
        (chrm, beg, end, name, value, strand) = line.split()
        utr_dict[name] = (chrm, int(beg), int(end), strand)

    return utr_dict

def extend_utrs(utrs, extendby, extended_path):
    """Take a .bedfile and extend the 3UTRs. Save only the extensions."""

    extended_handle = open(extended_path, 'wb')

    for line in open(utrs, 'rb'):
        (chrm, beg, end, ts_id, d, strnd) = line.split()
        end, beg = int(end), int(beg)
        length = end-beg
        extension = int(math.ceil(length*extendby))
        if strnd == '+':
            end = end + extension
            beg = end
        if strnd == '-':
            beg = beg - extension
            end = beg

        extended_handle.write('\t'.join([chrm, str(beg), str(end), ts_id, d,
                                         strnd]) + '\n')
    extended_handle.close()

def intersect(filea, fileb):
    """ Intersect filea with fileb. Report back only the entries in filea that
    did NOT intersect with anything in fileb (the '-v' parameter)"""

    cmd = ['intersectBed', '-v', '-a', filea, '-b', fileb]

    # Run the above command
    f = Popen(cmd, shell=False, stdout=PIPE)

    # From the output, get the ts_ID (column 4) that had no overlap
    return set([line.split()[3] for line in f.stdout])


def get_no_overlaps(exons, utrs):

    here = os.path.dirname(os.path.realpath(__file__))
    phdproj = '/users/rg/jskancke/phdproject/3UTR/'

    # These utrs are filtered: they are the longest utr of the genes
    utrs = phdproj+'Characterizer/Work_in_progress/source_bedfiles/3utrs.bed'

    extendby = 0.5
    # make a .bed file with just the extensions of the 3utr.bed file
    extended_path = os.path.join(os.path.dirname(utrs), 'utr_extensions.bed')

    # If the file doesn't alraedy exist, make it
    if os.path.isfile(extended_path):
        if os.path.getsize(extended_path) == 0:
            extend_utrs(utrs, extendby, extended_path)
    else:
        extend_utrs(utrs, extendby, extended_path)

    filea = extended_path
    fileb = exons
    # Get only those elements in fileA that do NOT intersect elements in fileB
    no_overlaps = intersect(filea, fileb)

    return no_overlaps

def get_polyA_dict(utr_path, polyA_path):
    polyA_dict = {}

    cmd = ['intersectBed', '-s', '-wb', '-a', polyA_path , '-b', utr_path]

    # Run the above command -- outside the shell -- and loop through output
    f = Popen(cmd, stdout=PIPE)
    for line in f.stdout:

        (chrm, beg, end, d, d, strnd, d, d, d, ts_id, d, d) = line.split()

        if ts_id in polyA_dict:
            polyA_dict[ts_id].append((chrm, beg, end, strnd))
        else:
            polyA_dict[ts_id] = [(chrm, beg, end, strnd)]

    return polyA_dict

def main():

    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    # For storing temporary files
    tempdir = os.path.join(here, 'temp_files')
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    # For storing the 3utr_annotation-bedfile if none is provided
    beddir = os.path.join(here, 'source_bedfiles')
    if not os.path.exists(beddir):
        os.makedirs(beddir)

    # Output directory
    output_dir = os.path.join(here, 'output')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Location of settings file
    settings_file = os.path.join(here, 'UTR_SETTINGS')

    # Get the necessary variables from the settings file
    settings = Settings(*read_settings(settings_file))

    # Whether to re-simulate or to get results from pickle
    #simulate = False
    simulate = True

    DEBUGGING = True
    #DEBUGGING = False
    if DEBUGGING:
        settings.chr1 = True
        settings.read_limit = 100000
        #read_limit = False
        settings.max_cores = 3
        #polyA = False
        #polyA = True
        settings.polyA = '/users/rg/jskancke/phdproject/3UTR/Characterizer/'\
                'Work_in_progress/temp_files/polyA_reads_k562_whole_cell_'\
                'processed_mapped_in_3utr.bed'

    #proj_dir = '__'.join(datasets.keys()+['Chr1'+str(chr1)]+['Reads'+str(read_limit)])
    #if not os.path.exists(proj_dir):
        #os.makedirs(proj_dir)

    # the Annotation instance that holds file-paths and datastructures obtained
    # from the annotation
    annotation = Annotation(settings.annotation_path)

    annotation.utrfile_path = get_utr_path(settings, beddir)
    # Get dictionary with utr-info
    print('Making 3UTR data structures ...\n')
    annotation.utrs = get_utrdict(annotation.utrfile_path)

    # file path to annotated polyA sites as obtained form annotation
    annotation.polyA_sites_path = get_polyA_sites_path(settings, beddir)
    # Get dictionary with polyA sites by intersecting with utr-bed
    print('Making polyA_sites data structures ...\n')
    annotation.polyA_sites_dict = get_polyA_dict(annotation.utrfile_path,
                                                annotation.polyA_sites_path)

    # Pickle the final results. Initiate the pickle object.
    pickled_final = os.path.join(output_dir, 'pickled_result_paths')

    ##################################################################
    if simulate:
        #Get dictionary of fasta-sequences from 'utrs'
        tx = time.time()
        fasta_bed = os.path.join(tempdir, 'bed_for_fasta.bed')
        print('Fetching the sequences of the annotated 3UTRs ...')
        utr_seqs = genome.get_seqs(annotation.utrs, settings.hgfasta_path)
        print('\tTime to get sequences: {0}\n'.format(time.time() - tx))

        # Create a pool of processes; one dataset will take up one process.
        my_pool = Pool(processes = settings.max_cores)
        results = []

        # Apply all datasets to the pool
        t1 = time.time()

        for dset_id, dset_reads in settings.datasets.items():

            arguments = (dset_id, dset_reads, tempdir, output_dir, utr_seqs,
                         settings, annotation)

            ###### WORK IN PROGRESS
            #akk = pipeline(dset_id, dset_reads, tempdir, output_dir, utr_seqs,
                           #settings, annotation)

            #debug()

            result = my_pool.apply_async(pipeline, arguments)
            results.append(result)

        # Wait for all procsses to finish
        my_pool.close()
        my_pool.join()

        # Get the paths from the final output
        dsets = settings.datasets.keys()
        outp = [result.get() for result in results]
        print('Total elapsed time: {0}\n'.format(time.time()-t1))

        # create output dictionaries
        coverage, final_output = {}, {}

        # Fill the dicts with paths
        for line in outp:
            for key, value in line.items():
                coverage[key] = value['coverage']
                final_output[key] = value['utr_ends']

        # Put these dicts in a dict again and pickle for future use
        dumpme = {'coverage': coverage, 'final_output': final_output}
        cPickle.dump(dumpme, open(pickled_final, 'wb'))

        # Copy output from temp-dir do output-dir
        save_output(final_output, output_dir)
    ##################################################################

    if not simulate:

        file_path_dict = cPickle.load(open(pickled_final, 'rb'))

        final_output = file_path_dict['final_output']
        coverage = file_path_dict['coverage']

    # Present output graphically
    #output_analyzer(final_output, utrs, utrfile_path)

    # Make some plots that compare utr_ends between tissues
    #make_utr_plots(final_output, utrs, coverage)

    # Delete temporary files
    #remove_directory(tempdir)

if __name__ == '__main__':
    main()

# You're going about this the wrong way. Here's how to get a good marker for
# 3UTR end.
# 1) Get all the annotated polyA sites (determined by 3UTR end) from gencode.
# 2) Get the self.polyA_cluster
# 3) Get the self.coverage diff on both sides of polyA_cluster sites
# 4) Get the cumul value and relative_pos value for each polyA site
# 5) For each gene, output the above variables to a file
# 5) Use that information to get a better cumul value cutoff for when a gene is
# finished

# NOTE from Friday: the GENCODE poly(A) reads are good and plentiful. You can
# rely on them,
# how many 'novel' poly(A) sites can we reliably verify?
# "poly(A) reads as a marker for RNA 3' ends"
# You can download the SVM and polyDB tracks from the browser. You can check of
# your putative transcript ends fall on these, close to these, or not close to
# these. Based on the coverage around the site, you can make a conclusion about
# the 3'UTR end at this point..

# The first thing to do with your reads is to go through some rounds of
# filtering.

# IDEA: when you merge the 3utrs together, can you keep the annotated 3utr ends?
# This list needs to be checked against he poly(As) you discover. How many are
# novel? How many are predicted?

# TODO for meeting with Roderic:
# 1) Point 1 from below
# 2) Pictures from the genome-browser with the poly(A) reads

    # For future:
    # 1) whole-cell, nucleus, cytosol 3utr lengths. Don't you always expect the
    # whole-cell to have the same length as the longest length of nucleus and
    # cytocol? Do some numbers on this.
    # 2) Can we identify some 3utrs that are longer than the annotation?
    # 3) A list of candidate genes that have different length (AND where the
    # result from 1) holds true? AND they have the highest rpkms?)
    # 4) You can evaluate the accuracy of the polyA reads by checking the
    # coverage before and after them.

# TODO AFTER meeting:
# Implement the python solution for getting sequences fast
# 1) Put the temp poly(A) reads in a different folder. Thus you can re-obtain
# them easily, because it takes so much time (and RAM) to remap them.

#* ways to visualize the 3'utr and the PAS motif; similar to what I've done with PERL; if python does not have a similar library we can implement this in a perl script;

"""
Read rna-seq read-files, a human genome annotation, a fasta file of the
entire human genome, and an inpute file UTR_SETTINGS. Based on this, write from
one to three output files, depending on the settings, that give information
about the usage of 3UTRs in the rna-seq reads.

The core of the program is a pipeline where rna-seq reads are extracted and
treated in several steps. The steps are

1) Extract the reads from a gem-mapping

1.1) If requested, also extract the polyA reads
1.2) Treat the polyA reads: trimming of poly(A) tail + remapping to the genome
1.3) Intersect the remapped polyA reads with the 3UTR
1.4) Cluster the polyA reads that map to the 3UTR (stochastic cleavage)

2) Get the RPKM of the 3UTRs: run intersectBed on the 3UTR bedfile and the readfile

3) Get read coverage of the 3UTRs: run coverageBed -d on the 3UTRs and the reads

4) Go through the coverage file one line at a time;
4.1) Each line is one base; each 3UTR exon has from 200-X0000 lines
4.2) While on the same utr, keep adding coverage data
4.3) When a line with a different UTR is encountered, pause, and calculate
output for the previous 3UTR. If however the previous 3UTR was part of a
multi-exon.......... begin here.

... TODO YOU ARE DOCUMENTING AND YOU ARE WORRIED THAT THE FULL RUN IS
CRASHING WHEN USING ALL READS BUT NOT WHEN USING 10000. NEED TO DEBUG THAT
BETTER.
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
from operator import attrgetter

import re

# only get the debug function if run from Ipython
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
    def debug(): pass

from subprocess import Popen
from subprocess import PIPE
import time
import math
import matplotlib.pyplot as plt

# Your own imports
import annotation_analysis_progress as genome

class Settings(object):
    """Store the settings obtained from the settings file"""
    def __init__(self, datasets, annotation_path, utrfile_provided, read_limit,
                 max_cores, chr1, hgfasta_path, polyA, min_utrlen, extendby,
                 del_reads, cumul_tuning):

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
        self.cumul_tuning = cumul_tuning

    def DEBUGGING(self):

        self.chr1 = True
        #self.read_limit = False
        self.read_limit = 100000
        #self.read_limit = 1000
        self.max_cores = 3
        #self.polyA = True
        #self.polyA = '/users/rg/jskancke/phdproject/3UTR/the_project/polyA_file/'\
                #'polyA_reads_k562_whole_cell_processed_mapped_in_3utr.bed'
        self.polyA = '/users/rg/jskancke/phdproject/3UTR/the_project/polyA_file/'\
                'polyA_reads_hela-s3_cytoplasm_processed_mapped_in_3utr.bed'

class Annotation(object):
    """A convenience class that holds the files and data structures relevant to
    the annotation that has been supplied."""

    def __init__(self, annotation_path):
        """Initiate as empty strings and fill in as they are being created."""
        self.path = annotation_path
        self.utrfile_path = ''
        self.utr_exons = ''
        self.a_polyA_sites_path = ''
        self.a_polyA_sites_dict = ''

    def get_utrdict(self):
        utr_dict = {}
        for line in open(self.utrfile_path, 'rb'):
            (chrm, beg, end, name, value, strand) = line.split()
            utr_dict[name] = (chrm, int(beg), int(end), strand)

        return utr_dict

    def get_polyA_dict(self):

        polyA_dict = {}

        polyA_path = self.a_polyA_sites_path
        utr_path = self.utrfile_path

        cmd = ['intersectBed', '-s', '-wb', '-a', polyA_path, '-b', utr_path]

        # Run the above command -- outside the shell -- and loop through output
        f = Popen(cmd, stdout=PIPE)
        for line in f.stdout:

            (chrm, beg, end, d, d, strnd, d, d, d, utr_id, d, d) = line.split()

            if utr_id in polyA_dict:
                polyA_dict[utr_id].append((chrm, beg, end, strnd))
            else:
                polyA_dict[utr_id] = [(chrm, beg, end, strnd)]

        # For the multi-exon utrs, add [] (they will probably not intersect with
        # annotated polyA sites.
        for utr_id in self.utr_exons.iterkeys():
            if utr_id not in polyA_dict:
                polyA_dict[utr_id] = []

        return polyA_dict

class UTR(object):
    """This class takes care of writing output to file. Avoids passing a million
    parameters here and there."""

    def __init__(self, chrm, beg, end, strand, val, utr_ID, rpkm, extendby,
                 first_covr, sequence, polyA_reads, a_polyA_sites):

        # variables you have to initialize
        self.chrm = chrm
        # Note that beg and end might be subject to extensions! See the
        # nonextended versions further down
        self.beg_ext = int(beg)
        self.end_ext = int(end)
        self.val = int(val)
        self.utr_ID = utr_ID
        self.strand = strand
        self.rpkm = rpkm
        self.extendby = extendby # how far has this annotation been extended
        self.sequence = sequence
        self.polyA_reads = polyA_reads # the polyA reads
        self.a_polyA_sites = [int(site[1]) for site in a_polyA_sites] # annotated sites

        # Assume not a multi exon
        if self.val == 1:
            self.multi_exon = False
        elif self.val > 1:
            self.multi_exon = True
        else:
            print('Exon number must be an integer larger than 1')
            sys.exit()

        self.exon_nr =int(utr_ID[-1])

        # If multi exon, update values
        self.tot_exnr = self.val

        # the coverage vectors
        self.covr_vector = [int(first_covr)]
        self.cumul_covr = [int(first_covr)]

        # Get the non-extended begs and ends as well!
        if extendby > 0:
            if strand == '+':
                self.end_nonext = self.end_ext - extendby
                self.beg_nonext = self.beg_ext
            if strand == '-':
                self.end_nonext = self.end_ext
                self.beg_nonext = self.beg_ext + extendby

        # Get the UNEXTENDED length
        self.length_ext = self.end_ext - self.beg_ext
        self.length_nonext = self.end_nonext - self.beg_nonext

        # if multi exon, make a beg-end list (extended and nonextended) for the
        # utrs
        if self.multi_exon:
            self.begends_ext = [(self.beg_ext, self.end_ex_extt)]
            self.begends_nonext = [(self.beg_nonext, self.end_nonext)]

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

    def update_utr(self, new_exon):
        """Update the first exon in the utr (which by sorting is the one with
        the smallest utr.beg) with the next exon."""

        # Update RPKM with weights on the lenghts of the two exons
        self.rpkm = ((self.rpkm/self.length) +
                     new_exon.rpkm/new_exon.length)/(self.length + new_exon.length)
        # Update length
        self.length += new_exon.length

        # Update begends lists to know where all exons begin and end
        self.begends.append((new_exon.beg, new_exon.end))
        if self.strand == '+':
            self.begends_nonextended.append((new_exon.beg, new_exon.end_nonextended))
        if self.strand == '-':
            self.begends_nonextended.append((new_exon.beg_nonextended, new_exon.end))

        # Update polyA_reads if any
        if new_exon.polyA_reads[0] != []:
            # If only one new polyA read
            if len(new_exon.polyA_reads[0]) == 1:
                self.polyA_reads[0].append(new_exon.polyA_reads[0][0])
                self.polyA_reads[1].append(new_exon.polyA_reads[1][0])

            else:
                # If more than one cluster, be smart about it
                # This code possibly also works for == 1, but too lazy to check
                avrges = new_exon.polyA_reads[0]
                supprt = new_exon.polyA_reads[1]
                for (pA_site, pA_support) in zip(avrges, supprt):
                    self.polyA_reads[0].append(pA_site)
                    self.polyA_reads[1].append(pA_support)

        # Update annotated polyA sites
        for site in new_exon.a_polyA_sites:
            self.a_polyA_sites.append(site)

        # Update the sequence, depending on strand
        # If + strand, new.exon is downstream self.exon
        if self.strand == '+':
            self.sequence = self.sequence + new_exon.sequence

        # If - strand, new.exon is upstream self.exon
        if self.strand == '-':
            self.sequence = new_exon.sequence + self.sequence

def frmt(self, element):
    """Return float objects with four decimals. Return all other objects as
    they were"""

    if type(element) is float:
        return format(element, '.4f')
    if type(element) is int:
        return str(element)
    else:
        return element

class FullLength(object):
    """Class for writing the UTR end file"""

    def __init__(self, utr_ID):
        self.utr_ID = utr_ID

        # variables that change depending on if reads map to 3UTR
        self.has_PAS = 'NA'
        self.pas_pos ='NA'
        self.ext_mean_99 = 'NA'
        self.int_mean_99 = 'NA'
        self.ext_mean_annot = 'NA'
        self.int_mean_annot = 'NA'
        self.cuml_rel_size = 'NA'
        self.polyA_support = 'NA'
        self.rel_pos = 'NA'

    def header_dict(self, this_utr):
        """Return  """
        return dict((('chrm', this_utr.chrm), ('beg', self.frmt(this_utr.beg_nonext)),
                    ('end', self.frmt(this_utr.end_nonext)),
                    ('utr_ID', self.frmt(this_utr.utr_ID[:-2])),
                    ('strand', this_utr.strand),
                    ('3utr_extended_by', self.frmt(this_utr.extendby)),
                    ('epsilon_coord', self.frmt(self.rel_pos)),
                    ('epsilon_rel_size', self.frmt(self.cuml_rel_size)),
                    ('epsilon_mean_int_covrg', self.frmt(self.int_mean_99)),
                    ('epsilon_mean_ext_covrg', self.frmt(self.ext_mean_99)),
                    ('annotation_mean_int_covrg', self.frmt(self.int_mean_annot)),
                    ('annotation_mean_ext_covrg', self.frmt(self.ext_mean_annot)),
                    ('epsilon_polyA_support', self.frmt(self.polyA_support)),
                    ('epsilon_PAS_type', self.has_PAS),
                    ('epsilon_PAS_distance', self.frmt(self.pas_pos)),
                    ('3utr_RPKM', self.frmt(this_utr.rpkm))
                     ))

    def frmt(self, element):
        """Return float objects with four decimals. Return all other objects as
        they came in"""

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
        strand
        utr_ID
        epsilon_coord
        epsilon_rel_size
        epsilon_mean_int_covrg
        epsilon_mean_ext_covrg
        annotation_mean_int_covrg
        annotation_mean_ext_covrg
        epsilon_polyA_support
        epsilon_PAS_type
        epsilon_PAS_distance
        3utr_RPKM
        """.split()

    def write_header(self, outfile):
        outfile.write('\t'.join(self.header_order()) + '\n')

    def calculate_output(self, pas_list, this_utr):
        # Don't do anything if coverage vector is empty
        if this_utr.is_empty():
            return

        # Call different strand-specific functions
        if this_utr.strand == '-':
            # calculate the cumulative values
            self.cumul_minus(this_utr)

        if this_utr.strand == '+':
            # calculate the cumulative values
            self.cumul_plus(this_utr)

        # see if the end position has polyA support (if polyA reads are given)
        if this_utr.polyA_reads != []:
            self.get_polyA_support(this_utr)

        # calculate the PAS and pas distance for 'length'
        self.get_pas(pas_list, this_utr)


    def cumul_minus(self, this_utr):
        """Note that rel_pos is relative to the un-extended 3UTR. Thus, a
        rel_pos value of 3 here would be the same as self.beg + extension + 3."""
        covr_vector = this_utr.covr_vector
        extendby = this_utr.extendby
        cumul_covr = this_utr.cumul_covr
        # NOTE the coverage arrays are in terms of the non-extended one

        # get the cumulated coverage resulting from the extension
        ext_cumul = cumul_covr[extendby-1]
        # subtract this from extension
        if ext_cumul > 0:
            cumul_covr = [val-ext_cumul for val in cumul_covr[extendby:]]
        else:
            cumul_covr = cumul_covr[extendby:]

        # Get normalized cuml-coverage of non-extended 3UTR. save to this_utr
        covr_sum = sum(covr_vector[extendby:])
        this_utr.norm_cuml = [1-val/covr_sum for val in cumul_covr]

        # Test for a special case where only the last entry has coverage
        if covr_sum == covr_vector[-1]:
            rel_pos = 1
            self.cuml_rel_size = rel_pos/float(this_utr.length_nonext)

        # Get the non_extended utr-relative position where 99.5% of the reads
        # have landed
        for ind, el in enumerate(this_utr.norm_cuml):
            if el < 0.995:
                rel_pos = ind
                length_nonext = this_utr.length_nonext
                self.cuml_rel_size = (length_nonext-rel_pos)/length_nonext
                break

        # Save relative cumulative position with the this utr for later usage
        this_utr.rel_pos = rel_pos
        #Note that rel_pos is relative to the un-extended 3UTR.

        # rel_pos according to extended 3utr
        ext_rel_pos = rel_pos + extendby

        # Then calculate the mean coverage on both sides of this.
        # Note for ext_mean_99: ext_rel_pos - extendby = ind
        self.ext_mean_99 = sum(covr_vector[ind:ext_rel_pos])/extendby
        self.int_mean_99 = sum(covr_vector[ext_rel_pos:ext_rel_pos +extendby])/extendby

        # Get the mean values 'extendby' around the annotated end too
        self.ext_mean_annot = sum(covr_vector[:extendby])/extendby
        self.int_mean_annot = sum(covr_vector[extendby: 2*extendby])/extendby

    def cumul_plus(self, this_utr):
        """Get cumulative read coverage"""
        covr_vector = this_utr.covr_vector
        extendby = this_utr.extendby

        # Get normalized cuml-coverage of UN-extended 3UTR
        cumul_covr = this_utr.cumul_covr[:-extendby] # non-extended cumul_covr
        covr_sum = sum(covr_vector[:-extendby]) # non-extended sum(covr_vector)
        this_utr.norm_cuml = [val/covr_sum for val in cumul_covr]

        # Test special case where only first entry has value
        if covr_sum == covr_vector[0]:
            rel_pos = 1
            self.cuml_rel_size = rel_pos/float((this_utr.length_nonext))

        # Get the point where 99.5% of reads have landed
        for ind, el in enumerate(reversed(this_utr.norm_cuml)):
            if el < 0.995:
                length_nonext = this_utr.length_nonext
                rel_pos = length_nonext - ind
                self.cuml_rel_size = rel_pos/float(length_nonext)
                break

        # Save relative position (relative to extended 3utr) with the object
        this_utr.rel_pos = rel_pos

        # rel_pos according to extended 3utr
        ext_rel_pos = rel_pos - extendby

        # Calculate the mean coverage on both sides of this.
        # Note for ext_mean_99: ext_rel_pos + extendby = rel_pos
        self.ext_mean_99 = sum(covr_vector[ext_rel_pos:rel_pos])/extendby
        self.int_mean_99 = sum(covr_vector[ext_rel_pos-extendby:ext_rel_pos])/extendby

        # Get the mean values extendby around the annotated end too
        self.ext_mean_annot = sum(covr_vector[-extendby:])/extendby
        self.int_mean_annot = sum(covr_vector[-2*extendby:-extendby])/extendby


    def get_pas(self, pas_list, this_utr):

        utr_ID = this_utr.utr_ID

        # rel pos is different depending on strand because we are dealing with a
        # reverse-complemented sequence that has been extended
        if this_utr.strand == '+':
            rel_pos = this_utr.rel_pos
        if this_utr.strand == '-':
            rel_pos = len(this_utr.sequence) - this_utr.extendby - this_utr.rel_pos

        rel_seq = this_utr.sequence[rel_pos-40:rel_pos]
        for pas in pas_list:
            try:
                (self.has_PAS, pas_indx) = (pas, rel_seq.index(pas))
                self.pas_pos = pas_indx - 40
                break
            except ValueError:
                (self.has_PAS, self.pas_pos) = ('NA', 'NA')


    def get_polyA_support(self, this_utr):
        """ Look for polyA-read support -- and get the cumulative
         percentage at which the majority of the polyAreads are
         found.
         """

        # Get the absolute end-position of the 99.5% 
        end_pos = this_utr.beg_nonext + this_utr.rel_pos

        # Check if this position has polyA reads within 50 nt
        for pAsite in this_utr.polyA_reads[0]:
            if pAsite-40 < end_pos < pAsite+40:
                # A general marker for polyA close to end_pos
                self.polyA_support = 1
                break

    def write_output(self, outobject, this_utr):
        """Format the output as desired, then save"""

        # Create a list of formatted output
        output_order = self.header_order()
        output_dict = self.header_dict(this_utr)

        output = [output_dict[hedr] for hedr in output_order]

        outobject.write('\t'.join(output) + '\n')

class PolyAReads(object):

    def __init__(self, utr_ID):
        self.utr_ID = utr_ID

        # variables to be printed
        self.polyRead_sites = 'NA'
        self.rel_polyRead_sites = 'NA'
        self.read_coverage_change = 'NA'
        self.annotation_support = 'NA'
        self.PAS_list = 'NA'

    def write_header(self, outfile):
        outfile.write('\t'.join(self.header_order()) + '\n')

    def frmt(self, element):
        """Return float objects with four decimals. Return all other objects as
        they were"""

        if type(element) is float:
            return format(element, '.4f')
        if type(element) is int:
            return str(element)
        else:
            return element

    def header_dict(self, this_utr, polA_nr, pAcoord, nr_supp_pA, covR, covL,
                    annotpA_dist, nearbyPAS, PAS_dist):

        return dict((
                    ('utr_ID', this_utr.utr_ID[:-2]),
                    ('polyA_number', self.frmt(polA_nr)),
                    ('strand', this_utr.strand),
                    ('polyA_coordinate', self.frmt(pAcoord)),
                    ('number_supporting_reads', self.frmt(nr_supp_pA)),
                    ('coverage_50nt_downstream', self.frmt(covR)),
                    ('coverage_50nt_upstream', self.frmt(covL)),
                    ('annotated_polyA_distance', self.frmt(annotpA_dist)),
                    ('nearby_PAS', nearbyPAS),
                    ('PAS_distance', self.frmt(PAS_dist)),
                    ))

    def header_order(self):
        return """
        utr_ID
        polyA_number
        strand
        polyA_coordinate
        number_supporting_reads
        coverage_50nt_downstream
        coverage_50nt_upstream
        annotated_polyA_distance
        nearby_PAS
        PAS_distance
        """.split()

    def write_output(self, outobject, this_utr):
        """Now there is 0 or > 1 output for each UTR. Need to loop through the
        number of pA sites and write output for each loop"""

        # Don't write for utrs without polyA-reads
        if self.polyRead_sites == 'NA':
            return

        # The column-order in which the output should be printed
        output_order = self.header_order()

        # The output relies on that the calculations have all been 'in correct
        # order' in calculated_output()
        for indx, site in enumerate(self.polyRead_sites):
            polAnr = indx + 1
            pAcord = site
            nr_supp_pA = len(this_utr.polyA_reads[1][indx])
            (covL, covR) = self.read_coverage_change[indx]
            annotpA_dist = self.annotation_support[indx]
            # One site might have several PAS of several creeds
            (nearbyPAS, PAS_dist) = self.select_PAS(indx)

            # Get the output dictionary with updated values
            output_dict = self.header_dict(this_utr, polAnr, pAcord, nr_supp_pA,
                                           covR, covL, annotpA_dist, nearbyPAS,
                                           PAS_dist)

            output = [output_dict[hedr] for hedr in output_order]

            outobject.write('\t'.join(output) + '\n')


    def select_PAS(self, index):
        if self.PAS_list[index] == []:
            return ('NA', 'NA')
        # If not, choose the best PAS
        if len(self.PAS_list[index]) == 1:
            return self.PAS_list[index][0]
        else:
            # Go through pas-list from best to worst. take the first you get
            paslist = get_pas_list()
            for pas in paslist:
                for this_pas, pasdist in self.PAS_list[index]:
                    if pas == this_pas:
                        return (this_pas, pasdist)

    def calculate_output(self, pas_list, this_utr):
        """
        For each reported polyA site as determined by poly(A) reads, check if
        1) The coverage changes 50nt on both sides of the site
        2) If there is an annotated polyA site nearby
        3) If there is a PAS nearby
        4) The relative length of the 3UTR at this site
        5) The normalized cumulative read value at this site

        And if the user supplies the information:

        6) SVM nearby?
        7) polyADB nearby?
        """

        # Don't do anything if coverage vector is empty
        if this_utr.is_empty():
            return

        self.polyRead_sites = this_utr.polyA_reads[0]
        # If there are no polyA reads in the UTR, pass on
        if self.polyRead_sites == []:
            return
        # Get the relative to extended utr(!) location of the polyA sites.
        self.rel_polyRead_sites = [pos-this_utr.beg_ext for pos in
                                       self.polyRead_sites]

        ############# TEST #################
        # for a test, print the distance from the end of the polyA reads as well
        # as their abundance
        #polyRead_sites_count = [len(sites) for sites in this_utr.polyA_reads[1]]

        #if self.strand  == '+':
            #end_dist = [self.end_nonextended - pos for pos in polyRead_sites]
            #print (self.strand, sorted(zip(end_dist, polyRead_sites_count)))

        #if self.strand == '-':
            #end_dist = [pos - self.beg_nonextended for pos in polyRead_sites]
            #print (self.strand, sorted(zip(end_dist, polyRead_sites_count)))
        ###################################

        #1) Get coverage on both sides of polyA read
        # This variable should be printed relative to polyA_read count divided
        # by the rpkm
        rccp = [(sum(this_utr.covr_vector[point:point-50])/50,
                  sum(this_utr.covr_vector[point:point+50])/50) for point in
                 self.rel_polyRead_sites]

        self.read_coverage_change = rccp

        #2) Is there an annotated polyA site nearby?
        # Report if there is one within +/- 40 nt and report the distance. Also
        # report if no distance is found.
        self.annotation_support = self.read_annotation_support(this_utr)
        # Result: a lot of support for the final site; less support for 'within'
        # sites

        #3) If there is a PAS nearby? This is secondary information. Look 40 nt
        #downstream and look for PAS. The SVM is a better indicator.
        self.PAS_list = self.get_PAS(pas_list, this_utr)

    def get_PAS(self, pas_list, this_utr):
        """Go through the -40 from the polyA read average. Collect any PAS you
        find."""
        pA_read_PAS = []
        for rpoint in self.rel_polyRead_sites:
            rel_seq = this_utr.sequence[rpoint-40:rpoint]
            found = False
            thispoint = []
            for pas in pas_list:
                try:
                    (PAS_type, pas_indx) = (pas, rel_seq.index(pas))
                    pas_rel_pos = 40 - pas_indx # AAT-FROMHERE-AAA
                    thispoint.append((PAS_type, pas_rel_pos))
                    found = True
                    # is this biased by finding the first? let's say there are
                    # two TATAAA -- now you only find the most distant one
                except ValueError:
                    (PAS_type, pas_rel_pos) = ('NA', 'NA')
            pA_read_PAS.append(thispoint)

        return pA_read_PAS

    def read_annotation_support(self, this_utr):
        """Check if the read polyA points have annotated points near to them"""
        supp = []
        for rpoint in self.polyRead_sites:
            found = False
            for apoint in this_utr.a_polyA_sites:
                if rpoint-40 < apoint < rpoint+40:
                    found = True
                    found_point = apoint
                    found_distance = rpoint-apoint
                    break
            if found:
                supp.append(found_distance)
            if not found:
                supp.append('NA')

        return supp

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
    f = Popen(cmd, stdout=PIPE)
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
            out_file.write('\t'.join([chrom, beg, str(int(beg)+len(seq)), '0',
                                      '0', strand]) + '\n')
            total_reads = total_reads + 1

        # When looking for poly(A) reads, filter by non-uniquely mapped reads
        if polyA == True:
            if mapinfo[:5] == '0:0:0':
                if seq[:2] == 'NN':
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

    f = Popen(cmd, stdout=outfile)
    f.wait()

    outfile.close()

    return out_path

def get_pas_list():
    # The known PAS in order of common appearance
    return set(['AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA',
            'AATACA', 'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA',
            'AATAGA'])

def join_multiexon_utr(multi_exon_utr):
    # Sort utr-exon instances according to utr-beg
    multi_exon_utr.sort(key = attrgetter('beg'))
    # Select the first exon in the utr as the 'main' utr. Add information from
    # the other utrs to this utr.
    main_utr = multi_exon_utr[0]
    for new_exon in multi_exon_utr[1:]:
        main_utr.update_utr(new_exon)

    return main_utr

def output_writer(dset_id, coverage, annotation, utr_seqs, rpkm, extendby,
                 polyA_reads, settings):
    """Putting together all the info on the 3UTRs and writing to files. Write
    one file mainly about the length of the 3UTR, and write another file about
    the polyA sites found in the 3UTR."""

    a_polyA_sites_dict = annotation.a_polyA_sites_dict
    utr_exons = annotation.utr_exons

    (dirpath, basename) = os.path.split(coverage)

    # Paths and file objecutr for the two output files (length and one polyA)
    length_outpath = os.path.join(dirpath, 'utr_'+dset_id)
    polyA_outpath = os.path.join(dirpath, 'polyA_'+dset_id)
    length_outfile = open(length_outpath, 'wb')
    polyA_outfile = open(polyA_outpath, 'wb')

    # list of PAS hexamers
    pas_list = get_pas_list()

    # Create the multi-exon dictionary where unfinshed multi-exon utrs will be
    # stored. When all exons of a multi exons have been added, they will be
    # combined and written to file.
    multi_exons = {}

    # First get line 1 as something to start from
    read_from = open(coverage, 'rb')

    line1 = read_from.next()
    (chrm, beg, end, utr_ID, val, strand, rel_pos, covr) = line1.split()

    # Create a UTR-instance
    this_utr = UTR(chrm, beg, end, strand, val, utr_ID, rpkm[utr_ID], extendby,
                   covr, utr_seqs[utr_ID], polyA_reads[utr_ID],
                   a_polyA_sites_dict[utr_ID])

    # Create instances for writing to two output files
    length_output = FullLength(utr_ID)
    pAread_output = PolyAReads(utr_ID)

    # Write the headers of the length and polyA output files
    length_output.write_header(length_outfile)
    pAread_output.write_header(polyA_outfile)

    # If tuning the cumulative length, open a file for this
    if settings.cumul_tuning:
        basedir = os.path.split(dirpath)[0]
        outfile = 'cumul_' + dset_id + '.stat'
        tuning_handle = open(os.path.join(basedir, 'output', outfile), 'wb')
        #write header for the polyA-cumulative stats
        tuning_handle.write('\t'.join(['utr_id', 'epsilon_relative',
                                       'pA_to_cumul_dist', 'pA_cumul',
                                       'covr_minus', 'covr_plus', 'rpkm',
                                       'utr-length', 'strand']) + '\n')

    # Assert that the file information is the same as you started with
    assert utr_exons[utr_ID] == (chrm, int(beg), int(end), strand), 'Mismatch'

    # staring reading from line nr 2
    for line in read_from:
        (chrm, beg, end, utr_ID, val, strand, pos, covr) = line.split('\t')

        # Check if you stay on the same UTR-exon in the coverage file
        if utr_ID == this_utr.utr_ID:

            this_utr.covr_vector.append(int(covr))
            this_utr.cumul_covr.append(this_utr.cumul_covr[-1] + int(covr))

        # if not, this is the first line of a new UTR
        else:

            calculate_and_write = True

            # check if the previous utr was a multi-exon utr
            if this_utr.tot_exnr > 1:
                # If so, check if it is complete
                utr_name = this_utr.utr_ID[:-2]

                if utr_name in multi_exons:
                    multi_exons[utr_name].append(this_utr)
                else:
                    multi_exons[utr_name] = [this_utr]

                if len(multi_exons[utr_name]) == this_utr.tot_exnr:
                    # If complete, join the utrs into a this_utr object
                    this_utr = join_multiexon_utr(multi_exons[utr_name])

                    # Create instances for writing to file
                    length_output = FullLength(utr_ID)
                    pAread_output = PolyAReads(utr_ID)

                else:
                    # if incomplete, continue without writing anything
                    calculate_and_write = False

            if calculate_and_write:

                # Calculate output values like 99.5% length, cumulative coverage, etc
                length_output.calculate_output(pas_list, this_utr)
                pAread_output.calculate_output(pas_list, this_utr)

                # Save output to files
                length_output.write_output(length_outfile, this_utr)
                pAread_output.write_output(polyA_outfile, this_utr)

                # If tuning, calculate the tuning variables and write to file
                if settings.cumul_tuning:
                    calc_write_tuning(tuning_handle, length_output, this_utr)

            # Update to the new utr and start the loop from scratch
            this_utr = UTR(chrm, beg, end, strand, val, utr_ID, rpkm[utr_ID],
                           extendby, covr, utr_seqs[utr_ID],
                           polyA_reads[utr_ID], a_polyA_sites_dict[utr_ID])

            # If not a multi exon, make output instances for writing to file
            if not this_utr.multi_exon:
                length_output = FullLength(utr_ID)
                pAread_output = PolyAReads(utr_ID)

            # Assert that the next utr has correct info
            assert utr_exons[utr_ID] == (chrm, int(beg), int(end), strand), 'Mismatch'

    # Close opened files
    length_outfile.close()
    polyA_outfile.close()

    if settings.cumul_tuning:
        tuning_handle.close()

    return (polyA_outpath, length_outpath)

def calc_write_tuning(tuning_handle, length_output, this_utr):
    """Write output important for tuning internal parameters"""

    if this_utr.is_empty():
        return

    # Get the absolute end-position of the 99.5% 
    # (this_utr.rel_pos is relative to the non-extended UTR)
    end_pos = this_utr.beg_nonext + this_utr.rel_pos

    # Get the cumulative coverage and +/- 50nt coverage of the closest pA site
    write_output = False
    for pAsite in this_utr.polyA_reads[0]:
        if pAsite-50 < end_pos < pAsite+50:
            # Get the absolute distance from the polyA site
            pA_dist = pAsite-end_pos

            # Get the relative-to-exon position of the polyA site 
            rel_pA_pos = pAsite - this_utr.beg_nonext

            # Test this by visualizing in the browser.
            # Upload the polyA reads and scrutinize the coverage vector

            # Get the coverage 50 nt on both sides of the polyA site
            if this_utr.strand == '+':
                d_stream_covr = sum(this_utr.covr_vector[rel_pA_pos-50:rel_pA_pos])/50
                u_stream_covr = sum(this_utr.covr_vector[rel_pA_pos:rel_pA_pos+50])/50

            if this_utr.strand == '-':
                d_stream_covr = sum(this_utr.covr_vector[rel_pA_pos:rel_pA_pos+50])/50
                u_stream_covr = sum(this_utr.covr_vector[rel_pA_pos-50:rel_pA_pos])/50

            # 
            debug()
            # RESULT: end-pos (thereby rel-pos) seems wrong! :S it's 
            # Rel_pos is according to the EXTENDED 3utr!
            # RESULT: polyA clusters are in correct place.
            # TODO fig rel_pos! I want it to be relative to the annotation!

            # Get the cumulative coverage of the polyA site
            # However, watch out for reads that land outside the non-extended
            # region

            # For now I simply assign them to the last value of the region.. but
            # it's not satisfactory!

            if (rel_pA_pos < 0) or (rel_pA_pos >= this_utr.length):
                if this_utr.strand == '+':
                    cumul_pA = this_utr.norm_cuml[-1]
                if this_utr.strand == '-':
                    cumul_pA = this_utr.norm_cuml[0]
            else:
                cumul_pA = this_utr.norm_cuml[rel_pA_pos]

            write_output = True

            break

    if write_output:
        d_stream_covr = str(d_stream_covr)
        u_stream_covr = str(u_stream_covr)
        pA_dist = str(pA_dist)
        cumul_pA = str(cumul_pA)
        rpkm = str(this_utr.rpkm)
        length = str(this_utr.length)
        default_pos = str(this_utr.rel_pos)
        strand = str(this_utr.strand)
        tuning_handle.write('\t'.join([this_utr.utr_ID[:-2], default_pos,
                                       pA_dist, cumul_pA, d_stream_covr,
                                       u_stream_covr, rpkm, length, strand]) +
                            '\n')

    # Write utr_id, distances, cumulative coverage, rpkm, and length of utr


def verify_access(f):
    try:
        open(f, 'rb')
    except:
        print('Could not access {0}.\nCheck if file exists, and if it exists '\
              'check permissions'.format(f))
        sys.exit()

def read_settings(settings_file):
    """Read settings file for parameters needed to run the script."""

    conf = ConfigParser.ConfigParser()
    conf.read(settings_file)

    expected_fields = ['DATASETS', 'ANNOTATION', 'CPU_CORES', 'RESTRICT_READS',
                       'CHROMOSOME1', 'SUPPLIED_3UTR_BEDFILE', 'HG_FASTA',
                       'POLYA_READS', 'MIN_3UTR_LENGTH', 'EXTEND', 'PLOTTING',
                       'DELETE_READS', 'UTR_LENGTH_TUNING']

    missing = set(conf.sections()) - set(expected_fields)
    if len(missing) == 0:
        pass
    else:
        print('The following options sections are missing: {}'.format(missing))
        sys.exit()

    # datasets and annotation
    datasets = dict((dset, files.split(':')) for dset, files in conf.items('DATASETS'))
    annotation = conf.get('ANNOTATION', 'annotation')
    # check if the files are actually there...
    datasets.pop('carrie')
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

    # tuning of the 99.5% value
    cuml_tuning = conf.getboolean('UTR_LENGTH_TUNING', 'tuning')

    # if polyA is false, you can't do tuning
    if cuml_tuning and not polyA:
        print('To tune the cumulative value, you need to enable polyA reads')
        sys.exit()

    return(datasets, annotation, utr_bedfile_path, read_limit, max_cores, chr1,
          hg_fasta, polyA, min_utrlen, extendby, del_reads, cuml_tuning)


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

def get_a_polyA_sites_path(settings, beddir):
    """Get a_polyA_sites_file path from annotation via genome module"""

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
    a_polyA_sites_filename = 'a_polyA_sites' + chrm + user_provided + '.bed'
    polyA_site_bed_path = os.path.join(beddir, a_polyA_sites_filename)

    ## If the file already exists, don't make it again
    if os.path.isfile(polyA_site_bed_path):
        return polyA_site_bed_path

    # If a utrfile is provided, never re-make from annotation
    if settings.utrfile_provided:
        polyA_site_bed_path = shape_provided_bed(polyA_site_bed_path, settings)

    # If utrfile is not provided, get it yourself from a provided annotation
    if not settings.utrfile_provided:
        if settings.chr1:
            settings.annotation_path = get_chr1_annotation(settings, beddir)

        t1 = time.time()
        print('Annotated polyA sites-file not found. Generating from annotation ...')
        genome.get_a_polyA_sites_bed(settings, polyA_site_bed_path)

        print('\tTime taken to generate polyA-bedfile: {0}\n'\
              .format(time.time()-t1))

    return polyA_site_bed_path

def get_utr_path(settings, beddir):
    """Get 3utr-file path from annotation via genome module. In the genome
    module, a bedfile is made from the annotation."""

    # Put together the name of the output file. The name depends on the options
    # that were used to get it.
    # If options are not set, make them a 0-string
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

    # If a utrfile is provided, never re-make from annotation, but shape it to
    # your needs? This is a bit unfinished.
    if settings.utrfile_provided:
        utr_bed_path = shape_provided_bed(utr_bed_path, settings)

    # If utrfile is not provided, get it yourself from a provided annotation
    if not settings.utrfile_provided:
        if settings.chr1:
            settings.annotation_path_chr1 = get_chr1_annotation(settings, beddir)

        t1 = time.time()
        print('3UTR-bedfile not found. Generating from annotation ...')
        #reload(genome)
        genome.get_3utr_bed_all_exons(settings, utr_bed_path)

        print('\tTime taken to generate 3UTR-bedfile: {0}\n'\
              .format(time.time()-t1))

    return utr_bed_path

def get_chr1_annotation(settings, beddir):
    """Split the original (GENCODE) annotation into one with only chrm1 """

    (name, suffix) = os.path.splitext(os.path.basename(settings.annotation_path))
    filename = name + '_chr1' + suffix
    outpath = os.path.join(beddir, filename)

    # If the file already exists, don't make it again
    if os.path.isfile(outpath):
        return outpath

    t1 = time.time()
    print('Separating chr1 from the annotation ...')
    outhandle = open(outpath, 'wb')

    for line in open(settings.annotation_path, 'rd'):
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

def get_rpkm(reads, utrfile_path, total_reads, utrs, extendby, dset_id):
    """Run coverageBed of reads on provided 3UTR and get RPKM. Correct for any
    extension made to the 3UTRs"""
    # Information needed for rpkm:
    # TOTAL reads
    # Length of UTR
    # # of reads landing in the UTR.
    rpkm = {}
    # If the bed-file has been extended, we need to unextend it.
    if extendby:
        temp_bed_path = os.path.join(os.path.dirname(utrfile_path), dset_id+
                                     '_temp_bed')
        temp_handle = open(temp_bed_path, 'wb')
        for line in open(utrfile_path, 'rb'):

            (chrm, beg, end, utr_id, val, strand) = line.split()

            # Reverse the extensions so you get correct RPKM!
            # Only extend those that were extended!
            # consider exNr_TotExNr
            # The 1_4 should be unextended for the - strand
            # the 4_4 should be unextended for the + strand
            # the 1_1 should be unextended for both strands

            ex_nr = int(utr_id.split('_')[-1])
            tot_exnr = int(val)

            if ex_nr == tot_exnr == 1:
                if strand == '+':
                    end = int(end) - extendby
                if strand == '-':
                    beg = int(beg) + extendby

            if ex_nr == tot_exnr > 1:
                if strand == '+':
                    end = int(end) - extendby

            if (tot_exnr > 1) and (ex_nr == 1):
                if strand == '-':
                    beg = int(beg) + extendby

            temp_handle.write('\t'.join([chrm, str(beg), str(end), utr_id, val,
                                         strand]) + '\n')
        temp_handle.close()

        utrfile_path = temp_bed_path

    p = Popen(['coverageBed', '-a', reads, '-b', utrfile_path], stdout=PIPE)

    for line in p.stdout:
        (chrm, bl, bl, utr_id, d, d, reads_covering) = line.split('\t')[:7]
        utr_length = utrs[utr_id][2] - utrs[utr_id][1]
        rpkm[utr_id] = ((10**9)*int(reads_covering))/(total_reads*utr_length)

    # Delete the tempfile
    if extendby:
        os.remove(temp_bed_path)

    return rpkm

def process_reads(pA_reads_path):
    """ Remove reads that are too short or have a poor nucleotide composition.
    Re
    """
    processed_reads = os.path.splitext(pA_reads_path)[0]+'_processed.fas'
    outfile = open(processed_reads, 'wb')

    # Go through the file two lines at a time. If the next line does not begin
    # with '>', append line to last entry (it's a split fasta-file).

    length_sum = 0
    tot_reads = 0

    for line in open(pA_reads_path, 'rb'):
        # add lines until there are 2 entries in linepair
        seq = line.rstrip()
        seqlen = len(seq)
        if seqlen > 25:
            As = seq.count('A')
            Ts = seq.count('T')
            # only save if A/T frequencies are not abnormal
            if (As/seqlen < 0.50) and (Ts/seqlen < 0.50):
                length_sum += seqlen
                tot_reads += 1
                # trim the title
                outfile.write('>read\n'+line)

    if tot_reads > 0:
        avrg_len = length_sum/float(tot_reads)
    else:
        avrg_len = 0

    outfile.close()

    return processed_reads, avrg_len

def map_reads(processed_reads, avrg_read_len):
    """ Map the processed reads using gem-mapper. Use the average read length to
    determine the number of mismatches for the mapper. Only accept uniquely
    mapping reads."""

    ## Andrea's GEM index of hg19
    g_ind='/users/rg/atanzer/DATA/GEM_indices/Genomes/H.sapiens.genome.hg19.main'
    mapped_reads = os.path.splitext(processed_reads)[0]+'_mapped'

    # Naming the final output
    base_dir = os.path.dirname(os.path.split(mapped_reads)[0])
    polybed_path = os.path.splitext(processed_reads)[0] + '_mapped.bed'

    # How many mismatches depend on read length
    if avrg_read_len < 50:
        mismatch_nr = 1
    elif 50 < avrg_read_len < 100:
        mismatch_nr = 2
    elif 100 < avrg_read_len:
        mismatch_nr = 3

    # mapping trimmed reads
    command = "gem-mapper -I {0} -i {1} -o {2} -q ignore -m {3}"\
            .format(g_ind, processed_reads, mapped_reads, mismatch_nr)

    p = Popen(command.split())
    p.wait()

    # Accept mismatches according to average read length
    acceptables = {1: set(('1:0', '0:1')), 2: set(('1:0:0', '0:1:0', '0:0:1')),
                   3: set(('1:0:0:0', '0:1:0:0', '0:0:1:0', '0:0:0:1'))}

    acceptable = acceptables[mismatch_nr]
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
    processed_reads, avrg_read_len = process_reads(pA_reads_path)

    # Map the surviving reads to the genome and return unique ones
    unique_reads = map_reads(processed_reads, avrg_read_len)

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
    # Building a suffix: Start with a '.' separated list of
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
    f = Popen(cmd, stdout=out_file)
    f.wait()

    # If polyA is a path to a bedfile (or bedfiles) concatenate these too
    if type(polyA) == str:
        cmd = ['zcat', '-f'] + polyA
        f = Popen(cmd, stdout=polyA_file)
        f.wait()

    # Return the number of line numbers (= number of reads)
    return sum(1 for line in open(out_file, 'rb'))

def get_polyA_utr(polyAbed, utrfile_path):
    """ Intersect poly(A) reads with 3UTR file. Return dictionary where each
    3UTR is listed with its poly(A) sites."""
    # A dictionary to hold the utr_id -> poly(A)-reads relation
    utr_polyAs = {}

    cmd = ['intersectBed', '-wb', '-a', polyAbed, '-b', utrfile_path]

    # Run the above command -- outside the shell -- and loop through output
    f = Popen(cmd, stdout=PIPE)
    for line in f.stdout:
        (polyA, utr) = (line.split()[:6], line.split()[6:])
        utr_id = utr[3]
        if not utr_id in utr_polyAs:
            utr_polyAs[utr_id] = [tuple(polyA)]
        else:
            utr_polyAs[utr_id].append(tuple(polyA))

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
    """Cluster the poly(A) reads for each gene. Choose the pTTS site depending
    on the strand of the ts. Return dictionary with (utr_id, clustered pA sites)"""

    # Test if there are no polyA reads: return empty lists 
    if utr_polyAs == {}:
        polyA_reads = dict((utr_id, [[],[]]) for utr_id in utrs)
        return polyA_reads

    plus_values = {'+': [], '-':[]}
    minus_values = {'+': [], '-':[]}

    polyA_reads = {}

    for utr_id, polyAs in utr_polyAs.iteritems():
        utr_info = utrs[utr_id]
        real_strand = utr_info[3]

        # Getting statistics for which strand the polyA reads map to
        # ONLY UNCOMMENT WHEN RUNNING WITH ALL READS
        #nr_polyA_reads = len(polyAs)
        #plus_ratio = sum(1 for tup in polyAs if tup[5] == '+')/nr_polyA_reads
        #minus_ratio = sum(1 for tup in polyAs if tup[5] == '-')/nr_polyA_reads

        #if real_strand == '+':
            #plus_values['+'].append(plus_ratio)
            #plus_values['-'].append(minus_ratio)
        #if real_strand == '-':
            #minus_values['+'].append(plus_ratio)
            #minus_values['-'].append(minus_ratio)

        # For + strands, choose beg as polyA_site; choose end for - strands.
        # Select only the reads with opposite strand! Make use of bias :)
        # TODO do statistics on how many in which strand; count the number
        # depending on direction so you can justify removing those with only one
        # read. This will help you screening out fake poly(A)s for non-stranded
        # data as well :)

        if real_strand == '+':
            ends = sorted([tup[2] for tup in polyAs if tup[5] == '-'])
        if real_strand == '-':
            ends = sorted([tup[1] for tup in polyAs if tup[5] == '+'])

        # Getting the actual clusters
        polyA_reads[utr_id] = cluster_loop(ends)

    # For those utr_id that don't have a cluster, give them an empty list; ad hoc
    for utr_id in utrs:
        if utr_id not in polyA_reads:
            polyA_reads[utr_id] = [[],[]]

    # Statistics on the ratio of + and - mapped reads for the genes that are in
    # the positive or negative strand. RESULT: for paired end reads, the polyA
    # reads map to the OPPOSITE strand 90% of the times.

    #plus_avrg = {'+': sum(plus_values['+'])/len(plus_values['+']),
                 #'-': sum(plus_values['-'])/len(plus_values['-'])}

    #minus_avrg = {'+': sum(minus_values['+'])/len(minus_values['+']),
                 #'-': sum(minus_values['-'])/len(minus_values['-'])}

    return polyA_reads


def pipeline(dset_id, dset_reads, tempdir, output_dir, utr_seqs, settings,
             annotation):
    """Get reads, get polyA reads, cluster polyA reads, get coverage, do
    calculations, write to output... this is where it all happens"""

    # Give new names to some parameters to shorten the code
    utrfile_path = annotation.utrfile_path
    extendby = settings.extendby
    read_limit = settings.read_limit
    polyA = settings.polyA

    t0 = time.time()

    # Get the normal reads (in bed format). Get the polyA reads if this option is set.
    # As well get the total number of reads for calculating the RPKM later.
    (bed_reads, polyA_bed, total_nr_reads) = get_bed_reads(dset_reads, dset_id,
                                                      read_limit, tempdir, polyA)

    # If polyA is a string, it is a path of a bedfile with to polyA sequences
    if type(polyA) == str:
        polyA_bed_path = polyA
        utr_polyAs = get_polyA_utr(polyA_bed_path, utrfile_path)

    # If polyA is True, trim the extracte polyA reads, remap them, and save the
    # uniquely mappign ones to bed format
    elif polyA == True:
        # Do the polyA: extract, trim, remap, and -> .bed-format.
        print('Processing poly(A) reads for {0}...'.format(dset_id))
        polyA_bed_path = remap_polyAs(polyA_bed)

        # Get a dictionary for each utr_id with its overlapping polyA reads
        utr_polyAs = get_polyA_utr(polyA_bed_path, utrfile_path)

    # Cluster the poly(A) reads for each utr_id
    if polyA:
        polyA_reads = cluster_polyAs(utr_polyAs, annotation.utr_exons)

    # Get the RPKM
    print('Obtaining RPKM for {0} ...\n'.format(dset_id))
    rpkm = get_rpkm(bed_reads, utrfile_path, total_nr_reads, annotation.utr_exons,
                    extendby, dset_id)

    # Cover the 3UTRs with the reads
    print('Getting read-coverage for the 3UTRs for {0} ...\n'.format(dset_id))
    options = '-d'
    coverage_path = coverage_wrapper(dset_id, bed_reads, utrfile_path, options)

    # Iterate through the coverage files and write to output files
    print('Writing output files... {0} ...\n'.format(dset_id))
    output_files = output_writer(dset_id, coverage_path, annotation, utr_seqs,
                                 rpkm, extendby, polyA_reads, settings)
    print('Total time for {0}: {1}\n'.format(dset_id, time.time() - t0))

    # If requested, delete the read-files after use
    if settings.del_reads:
        os.remove(bed_reads)

    # Get the paths to the two output files explicitly
    (polyA_output, length_output) = output_files

    # Return a dictionary with the file paths of the output files
    return {dset_id: {'coverage': coverage_path, 'length': length_output,
                      'polyA':polyA_output}}

def remove_directory(some_dir):
    try:
        shutil.rmtree(some_dir)
    except OSError:
        print('Unable to remove {0}'.format(some_dir))


def output_analyzer(final_output, utrs, utrs_path):
    """Present graphically the different lengths of the different tissues"""

    ##########
    # Restrict final output according to cell_line
    #filtered_output = dict((key, final_output[key]) for key in final_output if
                        #'K562' in key)
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
    for iutr_id, el_sizes in measure_dict.items():
        if 0 in el_sizes:
            continue
        else:
            without_zeros += 1
            # Calculate distances between all rel_dist
            parw_dist = [abs(val1-val2) for val1 in el_sizes for val2 in
                         el_sizes]
            sig_dist = ['yes' for val in parw_dist if val > 0.3]
            if 'yes' in sig_dist:
                strict_rel_len[iutr_id] = el_sizes

    pprint(without_zeros)
    pprint(len(strict_rel_len))


def statistics_boxplot(filtered_output, utrs, utrs_path):
    """Do basic statistics and make boxplots """

    # Make a relationship between dset and index for future use
    dset_indx = dict((dset, indx) for indx, dset in enumerate(filtered_output))
    exp_nr = len(filtered_output)

    # Prepare a dictionary of iutr_id with relative lenghts of compartments
    cuml_99_dict = dict((utr_id, [0 for val in range(exp_nr)]) for utr_id in utrs)

    # Get the ids of the transcripts whose extensions do not overlap exons in
    # the gencode annotation. (it's a set() by default)
    exons = '/users/rg/jskancke/phdproject/3UTR/gencode5/just_exons/just_exons.bed'

    rpkm_threshold = 1
    for dset, path in filtered_output.items():
        for line in open(path, 'rb'):
            # Skip lines that don't start with chr (it's fast enough)
            if line[:3] != 'chr':
                continue
            (chrm, beg, end, utr_id, strand, cuml99, int_mean_99, ext_mean_99,
             int_mean_annot, ext_mean_annot, PAS_type,\
             pas_distance, rpkm) = line.split()
            rpkm = float(rpkm)
            # Get 3utrs with rpkm above a set threshold
            if rpkm > rpkm_threshold:
                cuml_99_dict[utr_id][dset_indx[dset]] = float(cuml99)

    # Get how many of the transcripts are 0 or > 0.999 in rel_len.
    # count 0, > 99
    one_or_zero = [[0,0], [0,0], [0,0]]
    zeros = [0,0,0,0]

    for utr_id, rel_values in cuml_99_dict.items():
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
    # Choose a utr_id from final_output where there is a clear difference between
    # the different tissues (one has 0.4, one has 0.6, one has 0.8 rel_len)
    # from the final_output files. Also make sure they have reasonably high
    # RPKMs. Then look through the coverage files for these utr_ids, and
    # calculate the average coverage per 10 reads or so and produce a plot of
    # this.

    # Select the transcripts that have rpkm over 2
    exp_nr = len(final_output)
    rel_len_dict = dict((utr_id, [0]*exp_nr) for utr_id in utrs)
    index = 0
    count = 0
    for dset, path in final_output.items():
        # Skip the first line starting with '#'
        for line in open(path, 'rb'):
            # Skip lines that don't start with chr
            if line[:3] != 'chr':
                continue
            (chrm, beg, end, utr_id, last_utr_read, strand, rel_size, mean_cvr,
             std_cvr, PAS_type, pas_distance, rpkm) = line.split()
            rpkm = float(rpkm)
            if rpkm > 2:
                count = count+1
                rel_len_dict[utr_id][index] = float(rel_size)
        index = index + 1

    print count

    # Keep only the transcripts that 1) have no 0 rel_lens AND 2) have a minimum
    # distance of 0.4 btween two of the rel_lens.
    strict_rel_len = {}
    for iutr_id, el_sizes in rel_len_dict.items():
        if 0 in el_sizes:
            continue
        else:
            # Calculate distances between all rel_dist
            parw_dist = [abs(val1-val2) for val1 in el_sizes for val2 in
                         el_sizes]
            sig_dist = ['yes' for val in parw_dist if val > 0.4]
            if 'yes' in sig_dist:
                strict_rel_len[iutr_id] = el_sizes

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
            (chrm, beg, end, utr_id, d, strnd, rel_pos, covrg) = line.split()
            if utr_id in key_ids:
                covrg_dict[utr_id][indx].append(int(covrg))
        indx = indx+1

    for utr_id, covregs in covrg_dict.items():
        plot_3utr_ends(utr_id, covregs, utrs[utr_id][3])

def plot_3utr_ends(utr_id, coverages, strand):
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

def extend_utrs(utrs, extendby, extended_path):
    """Take a .bedfile and extend the 3UTRs. Save only the extensions."""

    extended_handle = open(extended_path, 'wb')

    for line in open(utrs, 'rb'):
        (chrm, beg, end, utr_id, d, strnd) = line.split()
        end, beg = int(end), int(beg)
        length = end-beg
        extension = int(math.ceil(length*extendby))
        if strnd == '+':
            end = end + extension
            beg = end
        if strnd == '-':
            beg = beg - extension
            end = beg

        extended_handle.write('\t'.join([chrm, str(beg), str(end), utr_id, d,
                                         strnd]) + '\n')
    extended_handle.close()

def intersect(filea, fileb):
    """ Intersect filea with fileb. Report back only the entries in filea that
    did NOT intersect with anything in fileb (the '-v' parameter)"""

    cmd = ['intersectBed', '-v', '-a', filea, '-b', fileb]

    # Run the above command
    f = Popen(cmd, stdout=PIPE)

    # From the output, get the utr_id (column 4) that had no overlap
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

def polyA_analysis(final_outp_polyA):
    myfile = final_outp_polyA['k562_whole_cell']
    for line in myfile:
        debug()

    pass

def make_directories(here, dirnames):
    """For each name in dirnames, return a list of paths to newly created
    directories in folder 'here'."""
    outdirs = []

    for dirname in dirnames:

        this_dir = os.path.join(here, dirname)

        if not os.path.exists(this_dir):
            os.makedirs(this_dir)

        outdirs.append(this_dir)

    return outdirs

def main():

    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    # Make directories needed by downstream code
    dirnames = ['temp_files', 'source_bedfiles', 'output']
    (tempdir, beddir, output_dir) = make_directories(here, dirnames)

    # Location of settings file
    settings_file = os.path.join(here, 'UTR_SETTINGS')
    # Get the necessary variables from the settings file and create the settings
    # object. This object will be sent around the program, so that settings are
    # always acessable.
    settings = Settings(*read_settings(settings_file))

    # When a simulation is over, the paths to the output files are pickled. When
    # simulate is False, the program will not read rna-seq files but will
    # instead try to get the output files from the last simulation.
    #simulate = False
    simulate = True

    # This option should be set only in case of debugging. It makes sure you
    # just run chromosome 1 and only extract a tiny fraction of the total reads.
    DEBUGGING = True
    #DEBUGGING = False
    if DEBUGGING:
        settings.DEBUGGING()

    # The program reads a lot of information from the annotation. The annotation
    # object will hold this information (file-paths and datastructures).
    annotation = Annotation(settings.annotation_path)

    # Check if 3UTRfile has been made or provided; if not, get it from annotation
    annotation.utrfile_path = get_utr_path(settings, beddir)

    # Get dictionary with (chrm, beg, end, strand) values for each 3utr-exon key
    print('Making 3UTR data structures ...\n')
    annotation.utr_exons = annotation.get_utrdict()

    # You extract all annotated polyA sites into a bedfile: a_polyA_sites_path
    annotation.a_polyA_sites_path = get_a_polyA_sites_path(settings, beddir)

    # You intersect the annotated polyA files with the 3UTR bedfile you got from
    # the annotation. Put the intersected annotated polyA sites in a dictionary.
    print('Making a_polyA_sites data structures ...\n')
    annotation.a_polyA_sites_dict = annotation.get_polyA_dict()

    # Pickle the final results. Initiate the pickle object.
    pickled_final = os.path.join(output_dir, 'pickled_result_paths')

    ##################################################################
    if simulate:
        # For all 3UTR exons, get the genomic sequence. The sequence is needed
        # to look for PAS sites.
        tx = time.time() # time the operation
        fasta_bed = os.path.join(tempdir, 'bed_for_fasta.bed')
        print('Fetching the sequences of the annotated 3UTRs ...')
        utr_seqs = genome.get_seqs(annotation.utr_exons, settings.hgfasta_path)
        print('\tTime to get sequences: {0}\n'.format(time.time() - tx))

        # Create a pool of processes; one dataset will take up one process.
        my_pool = Pool(processes = settings.max_cores)
        results = []

        # Apply all datasets to the pool
        t1 = time.time()

        # dset_id and dset_reads are as given in UTR_SETTINGS
        for dset_id, dset_reads in settings.datasets.items():

            # The arguments needed for the pipeline
            arguments = (dset_id, dset_reads, tempdir, output_dir, utr_seqs,
                         settings, annotation)

            ###### WORK IN PROGRESS
            akk = pipeline(dset_id, dset_reads, tempdir, output_dir, utr_seqs,
                           settings, annotation)
            debug()


            result = my_pool.apply_async(pipeline, arguments)
            results.append(result)

        #debug()
        # Wait for all procsses to finish
        my_pool.close()
        my_pool.join()

        dsets = settings.datasets.keys()
        # Get the paths from the final output
        outp = [result.get() for result in results]
        print('Total elapsed time: {0}\n'.format(time.time()-t1))

        # create output dictionaries
        coverage, final_outp_length, final_outp_polyA = {}, {}, {}

        # Fill the dicts with paths
        for line in outp:
            for key, value in line.items():
                coverage[key] = value['coverage']
                final_outp_polyA[key] = value['polyA']
                final_outp_length[key] = value['length']

        # Put these dicts in a dict again and pickle for future use
        dumpme = {'coverage': coverage, 'final_output_polyA':
                  final_outp_polyA, 'final_output_length': final_outp_length}

        cPickle.dump(dumpme, open(pickled_final, 'wb'))

        # Copy output from temp-dir do output-dir
        save_output(final_outp_polyA, output_dir)
        save_output(final_outp_length, output_dir)
    ##################################################################
    # NOTE TO SELF: everything starting from there should be in a separate
    # script: utr_output_analysis.py or smth similar. It's only been put here
    # for conveniece for checking the output during work.

    if not simulate:

        file_path_dict = cPickle.load(open(pickled_final, 'rb'))

        final_outp_polyA = file_path_dict['final_output_polyA']
        final_outp_length = file_path_dict['final_output_length']
        coverage = file_path_dict['coverage']

    # Present output graphically
    #output_analyzer(final_output, utrs, utrfile_path)
    #polyA_analysis(final_outp_polyA)

    # Make some plots that compare utr_ends between tissues
    #make_utr_plots(final_output, utrs, coverage)

    # Delete temporary files
    #remove_directory(tempdir)

if __name__ == '__main__':
    main()

# TODO: make wig-file from the coverage-file
# RESULT: you need to save it on your own server, accesable from ftp/http
# TODO: quality check the output to .stats files for coverage before/after
#       + don't do minus/plus: do downstream and upstream
# TODO: write the genomic rel_pos in the length file
# TODO: generate the figures.
# TODO: don't leave out polyA clusters from the 'wrong' strand. filter that
# TODO: verify the polyA reads clusters with polyAdb + other results
# TODO: go through the entire program and improve documentation.
# TODO: verify the 3UTR length with annotation + smth else

# PAPER IDEAS:
    # 1) Reproduce everything they have done with ESTs. This shows that detailed
    # sequencing now gives the same results as all the ESTs combined! :)
    # 2) Compare with the nature paper with 36bp reads. Even if they have many
    # more unmapped reads, the read length is critical to get good information.
    # Run the whole pipelien with _JUST_THE_EXTENSIOSN! you can make an
    # extension.bed file and mark it with the gene you extend from. extend until
    # you meet another annotated object, and at most a certain distance.
    # 3) Consult Hagen about what he would expect from the different
    # compartments
    # 4) The postoc's speech was that methylation downstream a poly(A) site can
    # contribute to polyadenylation, thuogh Hagen injected that it could
    # contribute to the slowing and subsequent release of the polymerase
    # 5) Don't remove the polyA clusters of different direction. A few of them
    # are genuinly from the other strand! Remove instread just those that have a
    # support of 1 read.
    # 6) Give the absolute number of poly(A) reads in the different
    # compartments. Seems like there are a lot more poly(A) reads in the
    # cytoplasm than in the nucleus


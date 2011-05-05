"""
Read rna-seq read-files, a human genome annotation, a fasta file of the
entire human genome, and an inpute file UTR_SETTINGS. Based on this, write from
one to three output files, depending on the settings, that give information
about the usage of 3UTRs in the rna-seq reads.

The core of the program is a pipeline where rna-seq reads are extracted and
treated in several steps. The steps are

1) Extract the reads from a gem-mapping

    If run with **polyA = true** in the settings file
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

Dependencies:
    * pyFasta
    * bedGraphToBigWig
    * bedTools
    * numpy (can easily be removed -- just mean and std used)
    * matplotlib (optional -- for plotting)
    * python 2.6.4 or greater

"""

from __future__ import division
print('Loading modules ...\n')
import os
import sys
import shutil
import cPickle
import ConfigParser
from multiprocessing import Pool
from multiprocessing import cpu_count
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

# Your own imports
import annotation_parser as genome

class Settings(object):
    """
    Convenience class.One instance is created from this class: it holds all the
    settings parameters obtained from the UTR_SETTINGS file.
    """
    def __init__(self, datasets, annotation_path, utrfile_provided, read_limit,
                 max_cores, chr1, hgfasta_path, polyA, min_utrlen, extendby,
                 cumul_tuning, bigwig, bigwig_datasets, bigwig_savedir,
                 bigwig_url):

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
        self.cumul_tuning = cumul_tuning
        self.bigwig = bigwig
        self.bigwig_datasets = bigwig_datasets
        self.bigwig_savedir = bigwig_savedir
        self.bigwig_url = bigwig_url

        # for debugging only
        self.bed_reads = False

    def DEBUGGING(self):
        """
        Modify settings for debugging ONLY!
        """

        self.chr1 = True
        #self.chr1 = False
        #self.read_limit = False
        #self.read_limit = 10000000
        self.read_limit = 100000
        self.max_cores = 3
        #self.polyA = True
        #self.polyA = False
        self.polyA = '/users/rg/jskancke/phdproject/3UTR/the_project/temp_files'\
                '/polyA_reads_HeLa-S3_Whole_Cell_processed_mapped.bed'
        #self.bed_reads = False
        #self.bed_reads = '/users/rg/jskancke/phdproject/3UTR/the_project/temp_files'\
                #'/reads_K562_Whole_Cell.bed'

class Annotation(object):
    """
    A convenience class. Refers to the files and data structures relevant to
    the annotation that has been supplied.
    """

    def __init__(self, annotation_path):
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
    """
    Class for UTR-objects. Holds all the relevant information about a UTR-object
    that has been calculated in the pipeline. Attributes are either
    variables that are written to output files, or variables that assist in
    calculating the output variables.

    The UTR object is created from individual exon objects. Often, there is only
    one exon for the UTR object. In these cases, the beg_ext/non_ext variables
    are clear. However, when a UTR object consists of two or more exons, the
    internal exons have not been extended, so for them the beg/end_ext and
    beg/end_nonext attributes have the same value.
    """

    def __init__(self, chrm, beg, end, strand, val, utr_ID, rpkm, extendby,
                 first_covr, sequence, polyA_reads, a_polyA_sites):

        # variables you have to initialize
        self.chrm = chrm

        # ASSUME that the exon has been extended. This might not be the case. 
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
        self.rel_pos = 'NA' # it might not be updated

        # The length of the extended UTR
        self.length_ext = self.end_ext - self.beg_ext

        # Assume not a multi exon
        if self.val == 1:
            self.multi_exon = False
        elif self.val > 1:
            self.multi_exon = True

        # The exon_nr for this exon in the utr.
        self.this_exnr = int(utr_ID.split('_')[-1])

        # Total number of exons in UTR
        self.tot_exnr = self.val

        # initiate the coverage vectors with the first 'coverage'
        self.covr_vector = [int(first_covr)]
        self.cumul_covr = [int(first_covr)]

        # If this exon has been extended, get the non-extended begs and ends. If
        # the exon has not been extended, set them as the same value as the
        # beg.ext. This is confusing, but I'm too tired to sort it out.

        # Start off by assuming no extensions: set non_extended to the same as
        # extended
        self.end_nonext = self.end_ext
        self.beg_nonext = self.beg_ext
        self.length_nonext = self.length_ext

        # TODO you are in the process of making sure that the non_extended is
        # only performed on exons that have in fact been extended.
        # If extended, update these values
        if extendby > 0:

            # This UTR has only one exon
            one_exon = (self.this_exnr == self.tot_exnr == 1)

            # This is the 3' exon on the + strand
            fin_ex_plus = (strand == '+') and (self.this_exnr == self.tot_exnr > 1)

            # This is the 3' exon on the - strand
            fin_ex_min = (strand == '-') and (self.tot_exnr > 1) and (self.this_exnr == 1)

            if one_exon:
                if strand == '+':
                    self.end_nonext = self.end_ext - extendby
                if strand == '-':
                    self.beg_nonext = self.beg_ext + extendby

            if fin_ex_plus:
                self.end_nonext = self.end_ext - extendby

            if fin_ex_min:
                self.beg_nonext = self.beg_ext + extendby

            # if either of the cases were encountered, update the unextended
            # length
            if one_exon or fin_ex_plus or fin_ex_min:
                self.length_nonext = self.end_nonext - self.beg_nonext

        # if multi exon, make a beg-end list (extended and nonext) for the
        # utrs
        if self.multi_exon:
            self.begends_ext = [(self.beg_ext, self.end_ext)]
            if extendby > 0:
                self.begends_nonext = [(self.beg_nonext, self.end_nonext)]

    def __repr__(self):
        return self.utr_ID[-8:]

    def __str__(self):
        return "\nChrm\t{0}\nBeg\t{1}\nEnd\t{2}\nStrand\t{3}\n"\
                .format(self.chrm, self.beg_nonext, self.end_nonext, self.strand)


    def is_empty(self):
        """
        Determine if non-extended coverage vector is empty, i.e, has no
        read-coverage.
        """

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
        """
        For multi-exon 3UTRs. Updates the attributes of the UTR-object with the
        next exon with the next exon. The first exon's ID will be the one
        written to file.
        """

        # Update RPKM with weights on the lenghts of the two exons
        self.rpkm = ((self.rpkm/self.length_nonext) +
                     new_exon.rpkm/new_exon.length_nonext)/\
                (self.length_nonext + new_exon.length_nonext)
        # Update length
        self.length_nonext += new_exon.length_nonext

        # Update begends lists to know where all exons begin and end
        self.begends_ext.append((new_exon.beg_ext, new_exon.end_ext))
        self.begends_nonext.append((new_exon.beg_nonext, new_exon.end_nonext))

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

def frmt(element):
    """
    Return floats as strings with four decimals.
    Return ints as strings.
    Return all other objects as they came in.
    """

    if type(element) is float:
        return format(element, '.4f')
    if type(element) is int:
        return str(element)
    else:
        return element

class FullLength(object):
    """
    Class for writing the 'length' file. Calculates the output based on
    variables from a UTR instance. Also writes the header of the 'length' output
    file.
    """

    def __init__(self, utr_ID):
        self.utr_ID = utr_ID

        self.PAS_dist = 'NA'
        self.has_PAS = 'NA'
        self.eps_dstream_coverage = 'NA'
        self.eps_ustream_coverage = 'NA'
        self.annot_dstream_coverage = 'NA'
        self.annot_ustream_coverage = 'NA'
        self.cuml_rel_size = 'NA'

        self.epsilon = 0.999

    def header_dict(self, this_utr):
        """
        Return a dictionary that maps header-entries to the UTR or FullLength
        instances that are being written to file. This ensures that each column
        contains the data that corresponds to the colum header.
        """
        return dict((
                    ('chrm', this_utr.chrm),
                    ('beg', frmt(this_utr.beg_nonext)),
                    ('end', frmt(this_utr.end_nonext)),
                    ('utr_ID', frmt(this_utr.utr_ID[:-2])),
                    ('strand', this_utr.strand),
                    ('3utr_extended_by', frmt(this_utr.extendby)),
                    ('epsilon_coord', frmt(this_utr.rel_pos)),
                    ('epsilon_rel_size', frmt(self.cuml_rel_size)),
                    ('epsilon_downstream_covrg', frmt(self.eps_dstream_coverage)),
                    ('epsilon_upstream_covrg', frmt(self.eps_ustream_coverage)),
                    ('annotation_downstream_covrg',
                     frmt(self.annot_dstream_coverage)),
                    ('annotation_upstream_covrg',
                     frmt(self.annot_ustream_coverage)),
                    ('epsilon_PAS_type', self.has_PAS),
                    ('epsilon_PAS_distance', frmt(self.PAS_dist)),
                    ('3utr_RPKM', frmt(this_utr.rpkm))
                     ))

    def header_order(self):
        """
        The order in which colums appear in the 'length' output file.
        """
        return """
        chrm
        beg
        end
        3utr_extended_by
        strand
        utr_ID
        epsilon_coord
        epsilon_rel_size
        epsilon_downstream_covrg
        epsilon_upstream_covrg
        annotation_downstream_covrg
        annotation_upstream_covrg
        epsilon_PAS_type
        epsilon_PAS_distance
        3utr_RPKM
        """.split()

    def write_header(self, outfile):
        """
        Write the header of the 'length' output file
        """
        outfile.write('\t'.join(self.header_order()) + '\n')

    def calculate_output(self, pas_patterns, this_utr):
        """
        Calculate output variables if *this_utr* is nonempty:
            1) Get cumulative coverage using the *epsilon* parameter
            2) See if there is a PAS close to the estimated end, and if so,
            return the distance to that PAS
        """
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

        # calculate the PAS and pas distance for 'length'
        self.get_pas(pas_patterns, this_utr)


    def cumul_minus(self, this_utr):
        """
        Calculate the position of cumulative coverage (e.g 0.98) relative to the
        length of the annotated 3UTR. This is refered to as the *epsilon*
        position. Calculate the coverage on both sides of this position.
        """
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
            if el < self.epsilon:
                rel_pos = ind
                length_nonext = this_utr.length_nonext
                self.cuml_rel_size = (length_nonext-rel_pos)/length_nonext
                break

        # Save relative cumulative position with the this utr for later usage
        this_utr.rel_pos = rel_pos
        #Note that rel_pos is relative to the un-extended 3UTR.

        # rel_pos according to extended 3utr
        er_pos = rel_pos + extendby

        # Calculate the mean coverage on both sides of rel_pos
        self.eps_dstream_coverage = sum(covr_vector[er_pos-50:er_pos])/float(50)
        self.eps_ustream_coverage = sum(covr_vector[er_pos:er_pos+50])/float(50)

        # Get the mean values 'extendby' around the annotated end too
        self.annot_ustream_coverage = sum(covr_vector[extendby: 2*extendby])/extendby
        self.annot_dstream_coverage = sum(covr_vector[:extendby])/extendby


    def cumul_plus(self, this_utr):
        """
        See cumul_minus
        """
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

        # Get the point where epsilon percent of reads have landed (e.g. 99.5)
        for ind, el in enumerate(reversed(this_utr.norm_cuml)):
            if el < self.epsilon:
                length_nonext = this_utr.length_nonext
                rel_pos = length_nonext - ind
                self.cuml_rel_size = rel_pos/float(length_nonext)
                break

        # Save relative position (relative to extended 3utr) with the object
        this_utr.rel_pos = rel_pos

        # Calculate the mean coverage on both sides of this.
        # Note for ext_mean_99: ext_rel_pos + extendby = rel_pos
        # ext_mean_99 -> upstream_coverage
        self.eps_dstream_coverage = sum(covr_vector[rel_pos:rel_pos+50])/float(50)
        self.eps_ustream_coverage = sum(covr_vector[rel_pos-50:rel_pos])/float(50)

        # Check for the occasions when rel_pos -50 is less than 0
        if rel_pos < 50 and rel_pos != 0:
            self.eps_ustream_coverage = sum(covr_vector[:rel_pos])/float(rel_pos)

        if rel_pos == 0:
            self.eps_ustream_coverage = 0

        # Get the mean values extendby around the annotated end too
        self.annot_dstream_coverage = sum(covr_vector[-2*extendby:-extendby])/extendby
        self.annot_ustream_coverage = sum(covr_vector[-extendby:])/extendby

    def get_pas(self, pas_patterns, this_utr):
        """
        Return all close-by PAS and their distances.
        """

        utr_ID = this_utr.utr_ID

        if this_utr.strand == '+':
            rel_pos = this_utr.rel_pos
            # For negative strand, take the length of the sequence mins rel_pos.
            # However, since the sequence has been extended, and rel_pos is from
            # the non-extended sequence, subtract the extension.
        if this_utr.strand == '-':
            rel_pos = len(this_utr.sequence) - this_utr.extendby - this_utr.rel_pos

        temp_PAS = []
        pas_seq = this_utr.sequence[rel_pos-40:rel_pos]
        for pas_exp in pas_patterns:
            temp_PAS.append([(m.group(), 40-m.end()) for m in
                            pas_exp.finditer(pas_seq)])

        if sum(temp_PAS, []) != []:
            has_PAS = sum(temp_PAS, [])
            if len(has_PAS) == 1:
                (has_PAS, PAS_dist) = ([has_PAS[0][0]], [has_PAS[0][1]])
            else:
                (has_PAS, PAS_dist) = zip(*has_PAS)

            self.has_PAS = ' '.join([str(pas) for pas in has_PAS])
            self.PAS_dist = ' '.join([str(dist) for dist in PAS_dist])


    def write_output(self, outobject, this_utr):
        """
        Write out the output as determined in 'header_order()'
        """

        # Get output dict
        output_dict = self.header_dict(this_utr)

        # Write output that corresponds to header in self.header_order
        output = [output_dict[hedr] for hedr in self.header_order()]

        outobject.write('\t'.join(output) + '\n')

class PolyAReads(object):
    """
    Returns object for writing to 'polyA' file. Contains method for writing
    header of 'polyA' file.
    """

    def __init__(self, utr_ID):
        self.utr_ID = utr_ID

        # variables to be printed
        self.polyRead_sites = 'NA'
        self.rel_polyRead_sites = 'NA'
        self.ds_covrg = 'NA'
        self.us_covrg = 'NA'
        self.annotation_support = 'NA'
        self.all_PAS = 'NA'

    def write_header(self, outfile):
        """
        Write the header of the 'polyA' output file
        """
        outfile.write('\t'.join(self.header_order()) + '\n')

    def header_dict(self, this_utr, polA_nr, pAcoord, nr_supp_pA, ds_covr,
                    us_covr, annotpA_dist, nearbyPAS, PAS_dist):
        """
        See the equivalent method for the 'FullLength' class.
        """

        return dict((
                    ('chrm', this_utr.chrm),
                    ('beg', frmt(this_utr.beg_nonext)),
                    ('end', frmt(this_utr.end_nonext)),
                    ('utr_ID', this_utr.utr_ID[:-2]),
                    ('polyA_number', frmt(polA_nr)),
                    ('strand', this_utr.strand),
                    ('polyA_coordinate', frmt(pAcoord)),
                    ('number_supporting_reads', frmt(nr_supp_pA)),
                    ('coverage_50nt_downstream', frmt(ds_covr)),
                    ('coverage_50nt_upstream', frmt(us_covr)),
                    ('annotated_polyA_distance', frmt(annotpA_dist)),
                    ('nearby_PAS', nearbyPAS),
                    ('PAS_distance', frmt(PAS_dist)),
                    ('3utr_RPKM', frmt(this_utr.rpkm))
                    ))

    def header_order(self):
        """
        See the equivalent method for the 'FullLength' class.
        """
        return """
        chrm
        beg
        end
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
        3utr_RPKM
        """.split()

    def write_output(self, outobject, this_utr):
        """
        Write the output in the order determined in 'header_order()'. For each
        UTR, there might be several polyA clusters. Each cluster gets a line in
        the output file.
        """

        # Don't write for utrs without polyA-reads
        if self.polyRead_sites == 'NA':
            return

        # The column-order in which the output should be printed
        output_order = self.header_order()

        # Write one output line for each polyA cluster
        for indx, site in enumerate(self.polyRead_sites):
            polAnr = indx + 1
            rel_polyA_site = self.rel_polyRead_sites[indx]
            nr_supp_pA = len(this_utr.polyA_reads[1][indx])
            ds_covrg = self.ds_covrg[indx]
            us_covrg = self.us_covrg[indx]
            annotpA_dist = self.annotation_support[indx]

            # One site might have several PAS sites
            (nearbyPAS, PAS_dist) = self.all_PAS[rel_polyA_site]

            # If found, Turn the PAS sites and their dists into strings
            if (nearbyPAS, PAS_dist) != ('NA', 'NA'):
                nearbyPAS = ' '.join([str(pas) for pas in nearbyPAS])
                PAS_dist = ' '.join([str(dist) for dist in PAS_dist])

            # Get the output dictionary with updated values
            output_dict = self.header_dict(this_utr, polAnr, site, nr_supp_pA,
                                           ds_covrg, us_covrg, annotpA_dist,
                                           nearbyPAS, PAS_dist)

            output = [output_dict[hedr] for hedr in output_order]

            outobject.write('\t'.join(output) + '\n')


    def calculate_output(self, pas_patterns, this_utr):
        """
        Calculate the following output:
            # The relative position of the cluster to non-extended UTR beg
            # Average coverage 50nt on both sides of the cluster
            # The distance to an annotated polyA site, if distance is less than
              40 nt.
            # The distance and type of PAS, if closer than 40 nt.
        """

        # Don't do anything if coverage vector is empty
        if this_utr.is_empty():
            return

        self.polyRead_sites = this_utr.polyA_reads[0]
        # If there are no polyA reads in the UTR, pass on
        if self.polyRead_sites == []:
            return
        # Get the relative to EXTENDED utr(!) location of the polyA sites,
        # because it's from the extended utr we take the coverage values.
        self.rel_polyRead_sites = [pos-this_utr.beg_ext for pos in
                                       self.polyRead_sites]

        ############# TEST #################
        # for a test, print the distance from the end of the polyA reads as well
        # as their abundance
        #polyRead_sites_count = [len(sites) for sites in this_utr.polyA_reads[1]]

        #if self.strand  == '+':
            #end_dist = [self.end_nonext - pos for pos in polyRead_sites]
            #print (self.strand, sorted(zip(end_dist, polyRead_sites_count)))

        #if self.strand == '-':
            #end_dist = [pos - self.beg_nonext for pos in polyRead_sites]
            #print (self.strand, sorted(zip(end_dist, polyRead_sites_count)))
        ###################################

        #1) Get coverage on both sides of polyA read
        # This variable should be printed relative to polyA_read count divided
        # by the rpkm

        # The question: are the sits in rel_polyRead_sites relative to the
        # extended or the non-extended utr?
        #if this_utr.utr_ID.startswith('ENSG00000142632_1'):
            #debug()
        #ENSG00000142632_1
        #and
        #ENSG00000097021_1

        # bug source: the 'relative' positions in rel_polyRead_sites can be
        # negative... which messes up list indexing completely of course.

        # Determine coverage relative to the vector
        left_covrg = []
        right_covrg = []

        for p in self.rel_polyRead_sites:
            right_covrg.append(sum(this_utr.covr_vector[p:p+50])/float(50))

            # Left-covrg depends on if you are close to the beg of cover vector
            if p < 50 and p != 0:
                left_covrg.append(sum(this_utr.covr_vector[:p])/float(p))
            elif p == 0:
                left_covrg.append(0)
            else:
                left_covrg.append(sum(this_utr.covr_vector[p-50:p])/float(50))


        # Find upstream/downstream from left/right depending on strand
        if this_utr.strand == '+':
            self.ds_covrg = right_covrg
            self.us_covrg = left_covrg

        if this_utr.strand == '-':
            self.ds_covrg = left_covrg
            self.us_covrg = right_covrg

        #2) Is there an annotated polyA site nearby?
        # Report if there is one within +/- 40 nt and report the distance. Also
        # report if no distance is found.
        self.annotation_support = self.read_annotation_support(this_utr)

        #3) If there is a PAS nearby? all_PAS hold all the PAS and their
        #distances for all the polyA read clusters in the 3UTR
        self.all_PAS = self.get_pas(pas_patterns, this_utr)

    def get_pas(self, pas_patterns, this_utr):
        """
        Go through the -40 from the polyA read average. Collect PAS and distance
        as you find them.
        """

        all_PAS = {}
        for rpoint in self.rel_polyRead_sites:
            pas_seq = this_utr.sequence[rpoint-40:rpoint]

            # special case if the polyA read is is early in the sequence
            if rpoint < 40:
                pas_seq = this_utr.sequence[:rpoint]

            temp_PAS = []
            for pas_exp in pas_patterns:
                temp_PAS.append([(m.group(), m.start()) for m in
                                pas_exp.finditer(pas_seq)])

            if sum(temp_PAS, []) != []:
                has_PAS = sum(temp_PAS, [])
                if len(has_PAS) == 1:
                    (has_PAS, PAS_dist) = ([has_PAS[0][0]], [has_PAS[0][1]])
                else:
                    (has_PAS, PAS_dist) = zip(*has_PAS)

                all_PAS[rpoint] = (has_PAS, PAS_dist)

            # In case no PAS are found, return 'NA'
            else:
                all_PAS[rpoint] = ('NA', 'NA')

        return all_PAS

    def read_annotation_support(self, this_utr):
        """
        Check if the polyA cluster has an annotated polyA site nearby.
        """
        supp = []
        for rpoint in self.polyRead_sites:
            found = False
            for apoint in this_utr.a_polyA_sites:
                if rpoint-40 < apoint < rpoint+40:
                    found = True
                    found_point = apoint
                    if this_utr.strand == '+':
                        found_distance = rpoint-apoint + 1
                    else:
                        found_distance = rpoint-apoint
                    # RESULT: for +, distance is one too little
                    # Otherwise, it seems good.
                    break
            if found:
                supp.append(found_distance)
            if not found:
                supp.append('NA')

        return supp

def zcat_wrapper(bed_reads, read_limit, out_path, polyA, polyA_path):
    """
    Wrapper around zcat. Call on gem-mapped reads. Write uniquely mapped
    reads (up to 2 mismatches) to .bed file.

    If polyA parameter was passed as True, write unmapped reads with leading
    poly(T) or tailing poly(A) to .bed file. The read is written to the bed-file
    as stripped of polyA or polyT stretch. Either 5 contiguous A/Ts or 6 A/Ts in
    the last/first 7 nucleotides must be present for the read to be considered
    as a putative poly(A) read.
    """

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
    """
    Strip the polyA tail of a read iteratively.
    """

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
    """
    Strip the polyT tail of a read iteratively.
    """

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
    """
    Wrapper around coverageBed. Calcuates coverage of the reads over the
    3UTR-regions, or whatever region was sent in to the program. If the region
    is large, this takes a lot of time, and the resulting file is huge.
    """

    # Rename path to 'covered_...'
    (dirpath, basename) = os.path.split(filtered_reads)
    out_path = os.path.join(dirpath, 'covered_'+dset_id)
    outfile = open(out_path, 'wb')

    cmd = ['coverageBed', options, '-a', filtered_reads, '-b', utrfile_path]

    f = Popen(cmd, stdout=outfile)
    f.wait()

    outfile.close()

    return out_path

def join_multiexon_utr(multi_exon_utr):
    """
    For multi-exon UTRs. When all exons in an UTR has been accounted for, start
    with the first exon (determined by sorting), and call
    'first_exon.update_utr(new_exon)' for all the new exons.
    """
    # Sort utr-exon instances according to utr-beg
    multi_exon_utr.sort(key = attrgetter('beg_nonext'))
    # Select the first exon in the utr as the 'main' utr. Add information from
    # the other utrs to this utr.
    main_utr = multi_exon_utr[0]
    for new_exon in multi_exon_utr[1:]:
        main_utr.update_utr(new_exon)

    return main_utr

def output_writer(dset_id, coverage, annotation, utr_seqs, rpkm, extendby,
                 polyA_reads, settings):
    """
    Putting together all the info on the 3UTRs and writing to files. Write
    one file mainly about the length of the 3UTR, and write another file about
    the polyA sites found in the 3UTR.
    """

    a_polyA_sites_dict = annotation.a_polyA_sites_dict
    utr_exons = annotation.utr_exons

    (dirpath, basename) = os.path.split(coverage)

    # Paths and file objecutr for the two output files (length and one polyA)
    length_outpath = os.path.join(dirpath, 'length_'+dset_id)
    polyA_outpath = os.path.join(dirpath, 'polyA_'+dset_id)
    length_outfile = open(length_outpath, 'wb')
    polyA_outfile = open(polyA_outpath, 'wb')

    # list of PAS hexamers
    PAS_sites = ['AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA',
                 'AATACA', 'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA',
                 'AATAGA']

    # compiled regular expressions
    pas_patterns = [re.compile(pas) for pas in PAS_sites]

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
                   covr, utr_seqs[utr_ID], polyA_reads[utr_ID]['other_strand'],
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
                length_output.calculate_output(pas_patterns, this_utr)
                pAread_output.calculate_output(pas_patterns, this_utr)

                # Save output to files
                length_output.write_output(length_outfile, this_utr)
                pAread_output.write_output(polyA_outfile, this_utr)

                # If tuning, calculate the tuning variables and write to file
                if settings.cumul_tuning:
                    calc_write_tuning(tuning_handle, length_output, this_utr)

            # Update to the new utr and start the loop from scratch
            this_utr = UTR(chrm, beg, end, strand, val, utr_ID, rpkm[utr_ID],
                           extendby, covr, utr_seqs[utr_ID],
                           polyA_reads[utr_ID]['other_strand'],
                           a_polyA_sites_dict[utr_ID])

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
    """
    Write parameters for tuning them. Used to determine the epsilon value for
    finding the cut-off where an UTR ends, and also for the before/after
    coverage ration of the polyA clusters.
    """

    if this_utr.is_empty():
        return

    # Get the absolute end-position of the 99.5%, relative to the non-extended UTR
    end_pos = this_utr.beg_nonext + this_utr.rel_pos

    # Get the cumulative coverage and +/- 50nt coverage of the closest pA site
    write_output = False
    close_sites = []

    for pAsite in this_utr.polyA_reads[0]:
        if pAsite-50 < end_pos < pAsite+50:
            close_sites.append(pAsite)
            write_output = True

    if write_output:
        # There could be 3 clusters within 50 nt of the end_pos. You need to
        # choose the most downstream cluster.
        close_sites.sort()
        if this_utr.strand == '+':
            pAsite = close_sites[-1]
        if this_utr.strand == '-':
            pAsite = close_sites[0]

        # Get the absolute distance from the polyA site
        pA_dist = pAsite-end_pos

        # Get the relative-to-non-extended position of the polyA site 
        rel_pA_pos = pAsite - this_utr.beg_nonext

        # Test this by visualizing in the browser.
        # Upload the polyA reads and scrutinize the coverage vector

        covr_vec = this_utr.covr_vector

        # Get the coverage 50 nt on both sides of the polyA site
        # No worried about extension -- it is in the end.
        if this_utr.strand == '+':
            d_stream_covr = sum(covr_vec[rel_pA_pos-50:rel_pA_pos])/float(50)
            u_stream_covr = sum(covr_vec[rel_pA_pos:rel_pA_pos+50])/float(50)

        # if negative strand and extended, become relative to the extended one
        if this_utr.strand == '-':
            if this_utr.extendby:
                rel_ex_pA_pos = rel_pA_pos + this_utr.extendby
            d_stream_covr = sum(covr_vec[rel_ex_pA_pos:rel_ex_pA_pos+50])/float(50)
            u_stream_covr = sum(covr_vec[rel_ex_pA_pos-50:rel_ex_pA_pos])/float(50)

        # Get the cumulative coverage of the polyA site
        # However, watch out for reads that land outside the non-extended
        # region
        # For now I simply assign them to the last value of the region.. but
        # it's not satisfactory!

        if (rel_pA_pos < 0) or (rel_pA_pos >= this_utr.length_nonext):
            if this_utr.strand == '+':
                cumul_pA = this_utr.norm_cuml[-1]
            if this_utr.strand == '-':
                cumul_pA = this_utr.norm_cuml[0]
        else:
            cumul_pA = this_utr.norm_cuml[rel_pA_pos]

        # TODO VERIFY CORRECTNESS OF COVERAGES and positions

        d_stream_covr = str(d_stream_covr)
        u_stream_covr = str(u_stream_covr)
        pA_dist = str(pA_dist)
        cumul_pA = str(cumul_pA)
        rpkm = str(this_utr.rpkm)
        length = str(this_utr.length_nonext)
        default_pos = str(this_utr.rel_pos)
        strand = str(this_utr.strand)
        tuning_handle.write('\t'.join([this_utr.utr_ID[:-2], default_pos,
                                       pA_dist, cumul_pA, d_stream_covr,
                                       u_stream_covr, rpkm, length, strand]) +
                            '\n')

    # Write utr_id, distances, cumulative coverage, rpkm, and length of utr


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

def read_settings(settings_file):
    """
    Read the settings and get all the settings parameters. These should be used
    to create a 'Settings' object.
    """

    conf = ConfigParser.ConfigParser()
    conf.optionxform = str
    conf.read(settings_file)

    expected_fields = ['DATASETS', 'ANNOTATION', 'CPU_CORES', 'RESTRICT_READS',
                       'CHROMOSOME1', 'SUPPLIED_3UTR_BEDFILE', 'HG_FASTA',
                       'POLYA_READS', 'MIN_3UTR_LENGTH', 'EXTEND', 'PLOTTING',
                       'UTR_LENGTH_TUNING', 'BIGWIG']

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

    # bigWig or not?
    bigwig = conf.getboolean('BIGWIG', 'bigwig')

    bigwig_datasets = []
    bigwig_savedir = ''
    bigwig_url = ''
    if bigwig:
        bigwig_datasets = conf.get('BIGWIG', 'datasets').split(':')
        bigwig_savedir = conf.get('BIGWIG', 'save_dir')
        bigwig_url = conf.get('BIGWIG', 'url')

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
          hg_fasta, polyA, min_utrlen, extendby, cuml_tuning, bigwig,
           bigwig_datasets, bigwig_savedir, bigwig_url)

def get_a_polyA_sites_path(settings, beddir):
    """
    Call the *annotation_parser* module to obtain the locations of all annotated
    utr-ends, in the program called 'a_polyA_sites'.
    """

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
    """
    Get 3utr.bed-file path from annotation via annotation_parser module.
    """

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
        genome.get_3utr_bed_all_exons(settings, utr_bed_path)

        print('\tTime taken to generate 3UTR-bedfile: {0}\n'\
              .format(time.time()-t1))

    return utr_bed_path

def get_chr1_annotation(settings, beddir):
    """
    From the annotation, extract only entries for chromosome 1 and save to a
    separate file. This file is needed for speed-runs of the pipeline, while
    still giving relavant results.
    """

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
    """
    Go through provided bedfile and shape it to confirm with internal standards
    in the program. Save in the bed-file directory.
    """

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
    """
    Copy files in out_dict from the temp_files to the output folder. These are
    the key output files that you don't want to delete by accident.
    """

    for ID, final_path in final_dict.items():
        out_path = os.path.join(output_dir, os.path.basename(final_path))
        shutil.copyfile(final_path, out_path)

def get_rpkm(reads, utrfile_path, total_reads, utrs, extendby, dset_id):
    """
    Run coverageBed of reads on provided 3UTR and get RPKM. The RPKM is for the
    un-extended 3UTR. You need: total number of reads, length of the UTR, and
    the number of reads landing in this 3UTR. The rpkm measure has its problems
    when it comes to comparison between samples, but it's easy to calculate, and
    it tells you about expression within one experiment.
    """
    # Information needed for rpkm:
    # TOTAL reads
    # Length of UTR
    # # of reads landing in the UTR.
    rpkm = {}
    # If the bed-file has been extended, we need to unextend it. HOWEVER, we
    # should ONLY un-extend the 3' exons of the UTR. Internal exons have not
    # been extended and should not be touched.
    if extendby:
        temp_bed_path = os.path.join(os.path.dirname(utrfile_path), dset_id+
                                     '_temp_bed')
        temp_handle = open(temp_bed_path, 'wb')
        for line in open(utrfile_path, 'rb'):

            (chrm, beg, end, utr_id, val, strand) = line.split()

            # Reverse the extensions so you get correct RPKM!
            # Only extend those that were extended in the first place!
            #ENSG00000078369_1_2	2	-
            # This means UTR number 1 from this gene, exon number 2 from this
            # UTR, and there are 2 utrs in total.

            this_ex_nr = int(utr_id.split('_')[-1])
            tot_exnr = int(val)

            if this_ex_nr == tot_exnr == 1:
                if strand == '+':
                    end = int(end) - extendby
                if strand == '-':
                    beg = int(beg) + extendby

            if this_ex_nr == tot_exnr > 1:
                if strand == '+':
                    end = int(end) - extendby

            if (tot_exnr > 1) and (this_ex_nr == 1):
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
    """
    Map the processed reads using gem-mapper. Use the average read length to
    determine the number of mismatches for the mapper according to the following
    scheme where X is the average read length:
        * if X < 50, then 1 mismatch
        * if 50 < X < 100, then 2 mismatchs
        * if 100 < X, then 3 mismatches
    Regardless, only accept uniquely mapping reads.
    """

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

        # Acceptable reads and poly(A) reads are mutually exclusive.
        if mapinfo in acceptable:
            # Get chromosome, strand, and beg
            (chrom, rest) = position.split(':')
            strand = getstrand[rest[0]]
            beg = start_re.match(rest[1:]).group()

            # Write to file in .bed format
            reads_file.write('\t'.join([chrom, beg, str(int(beg)+len(seq)), '0',
                                      '0', strand]) + '\n')

    return polybed_path

def get_bed_reads(dset_reads, dset_id, read_limit, tempdir, polyA):
    """
    Get reads from file. Determine file-type. If gem, extract from gem and
    convert to bed. If .bed, concatenate the bedfiles and convert them to the
    desired internal format.
    """

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
    """
    Throw zcat on the bedfiles. With the -f option, it will concatenate both
    compressed and noncompressed files.
    """

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
    """
    Call intersectBed on the poly(A) reads and on the 3UTR file. Return
    the surviving poly(A) reads in a dictionary, where each 3UTR (key) points
    to all its polyA sites.
    """
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
    """
    Cluster the polyA reads together. It is assumed that the list 'ends' comes
    in as sorted. Go through the list; the first site is a cluster; if the
    second site is within 20nt, it also becomes a cluster, and the new cluster
    site is the average of the two sites. If the next site were not withing 20
    nt, then this previous site is kept as a cluster, and the new site is the
    new cluster; and so on.
    """
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

    ## Get only the clusters with length more than one
    #res_clus = [clus for clus in clusters if len(clus) > 1]

    # Get the mean of the clusters
    averages_cluster = []

    for clus in clusters:
        # skip empty clusters
        if clus == []:
            continue
        averages_cluster.append(int(math.floor(sum(clus)/len(clus))))

    # Return the clusters and their average
    return (averages_cluster, clusters)

def cluster_polyAs(utr_polyAs, utrs, polyA):
    """
    For each UTR, save the clustered poly(A)-reads from both strands. For
    double-stranded reads, it is known that the poly(A) reads map to the
    opposite side of the strand it came from. This information is helpful in
    estimating the false positive rate.
    """

    # Ff there are no polyA reads or polyA is false: return empty lists 
    if utr_polyAs == {} or polyA == False:
        polyA_reads = dict((utr_id,{'this_strand':[[],[]], 'other_strand':[[], []]})
                           for utr_id in utrs)
        return polyA_reads

    plus_values = {'+': [], '-':[]}
    minus_values = {'+': [], '-':[]}

    polyA_reads = {}

    for utr_id, polyAs in utr_polyAs.iteritems():
        utr_info = utrs[utr_id]
        real_strand = utr_info[3]

        if real_strand == '+':
            other_strand_ends = sorted([tup[2] for tup in polyAs if tup[5] == '-'])
            this_strand_ends = sorted([tup[1] for tup in polyAs if tup[5] == '+'])
        if real_strand == '-':
            other_strand_ends = sorted([tup[1] for tup in polyAs if tup[5] == '+'])
            this_strand_ends = sorted([tup[2] for tup in polyAs if tup[5] == '-'])

        # Getting the actual clusters
        polyA_reads[utr_id] = {'this_strand': cluster_loop(this_strand_ends),
                               'other_strand': cluster_loop(other_strand_ends)}

    # For those utr_id that don't have a cluster, give them an empty list; ad hoc
    for utr_id in utrs:
        if utr_id not in polyA_reads:
            polyA_reads[utr_id] = {'this_strand':[[],[]], 'other_strand':[[],[]]}

    return polyA_reads


def pipeline(dset_id, dset_reads, tempdir, output_dir, utr_seqs, settings,
             annotation, DEBUGGING):
    """
    Get reads, get polyA reads, cluster polyA reads, get coverage, combine it in
    a 3UTr object, do calculations on the object attributes, write calculation
    to output files ... this is where it all happens: the PIPEline.
    """

    # Give new names to some parameters to shorten the code
    utrfile_path = annotation.utrfile_path
    extendby = settings.extendby
    read_limit = settings.read_limit
    polyA = settings.polyA
    # Define utr_polyAs in case polyA is false
    utr_polyAs = {}

    t0 = time.time()

    # IF DEBUGGING AND READS IN BED FORMAT ARE SUPPLIED
    # If both the bedfiles of the reads and the bedfiles of the polyA files are
    # submitted for debugging purposes, don't get reads anew
    if DEBUGGING and (type(settings.bed_reads) == str) and (type(polyA) == str):
        bed_reads = settings.bed_reads
        total_nr_reads = sum(1 for line in open(bed_reads, 'rb'))
        p_polyA_bed = polyA

    # Get the normal reads (in bed format). Get the polyA reads if this option is set.
    # As well get the total number of reads for calculating the RPKM later.
    # p_polyA_bed stands for putative polyA reads
    else:
        (bed_reads, p_polyA_bed, total_nr_reads) = get_bed_reads(dset_reads, dset_id,
                                                          read_limit, tempdir, polyA)

    # If polyA is a string, it is a path of a bedfile with to polyA sequences
    if type(polyA) == str:
        polyA_bed_path = polyA
        # Get a dictionary for each utr_id with its overlapping polyA reads
        utr_polyAs = get_polyA_utr(polyA_bed_path, utrfile_path)

    # If polyA is True, trim the extracte polyA reads, remap them, and save the
    # uniquely mappign ones to bed format
    elif polyA == True:
        # PolyA pipeline: remove low-quality reads, remap, and -> .bed-format:

        # 1) Process reads by removing those with low-quality, and removing the
        #    leading Ts and/OR trailing As.
        print('Processing poly(A) reads for {0}...'.format(dset_id))
        processed_reads, avrg_read_len = process_reads(p_polyA_bed)

        # 2) Map the surviving reads to the genome and return unique ones
        polyA_bed_path = map_reads(processed_reads, avrg_read_len)

        # Get a dictionary for each utr_id with its overlapping polyA reads
        utr_polyAs = get_polyA_utr(polyA_bed_path, utrfile_path)

    # Cluster the poly(A) reads for each utr_id. If polyA is false or no polyA
    # reads were found, return a placeholder, since this variable is assumed to
    # exist in downstream code.
    polyA_reads = cluster_polyAs(utr_polyAs, annotation.utr_exons, polyA)

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

    # Get the paths to the two output files explicitly
    (polyA_output, length_output) = output_files

    # Return a dictionary with the file paths of the output files
    return {dset_id: {'coverage': coverage_path, 'length': length_output,
                      'polyA':polyA_output}}

def make_directories(here, dirnames):
    """
    For each name in dirnames, return a list of paths to newly created
    directories in folder 'here'. Don't overwrite folders if they exist.
    """
    outdirs = []

    for dirname in dirnames:

        this_dir = os.path.join(here, dirname)

        if not os.path.exists(this_dir):
            os.makedirs(this_dir)

        outdirs.append(this_dir)

    return outdirs

# Make a bedfile for the 'length' UTRs.
def parse_length(in_length, out_length, out_header):

    in_handle = open(in_length, 'rb')
    orig_head = in_handle.next() # remove the header before reading

    out_handle = open(out_length, 'wb')
    # print the header of the output file
    out_handle.write(out_header+'\n')

    for line in in_handle:
        (chrm, beg, end, extendby, strand, ID, eps_coord) = line.split()[:7]
        # Skip the utrs without coverage
        if eps_coord == 'NA':
            continue

        # Calculate the new 'end' (strand dependent)
        eps_coord = int(eps_coord)
        if strand == '+':
            end = str(int(beg) + eps_coord)
        if strand == '-':
            beg = str(int(beg) + eps_coord)

        out_handle.write('\t'.join([chrm, beg, end, ID, '0', strand])+'\n')

    out_handle.close()
    in_handle.close()

def parse_clusters(in_clusters, out_clusters, out_header):
    """
    Parse the polyA cluster file from output and produce a bedfile of all the
    polyA clusters. The name is the utr_id and the score is the nr of supporting
    reads.
    """
    infile = open(in_clusters, 'rb')
    in_header = infile.next()

    outfile = open(out_clusters, 'wb')
    outfile.write(out_header + '\n')
    for line in infile:
        (chrm, utr_beg, utr_end, utr_id, clr_nr, strand, cl_pos,
         supp_reads) = line.split()[:8]

        cl_beg = cl_pos
        cl_end = str(int(cl_beg)+1)
        outfile.write('\t'.join([chrm, cl_beg, cl_end, utr_id, supp_reads,
                                 strand]) + '\n')

    infile.close()
    outfile.close()


def make_bigwigs(settings, annotation, here):
    """
    Make bigwig files of the polyA reads and the normal reads of the datasets
    specified under [BIGWIG] in the UTR_SETTINGS file.

    Turn the read coverage and polyA read .bedfiles into bigWig files.
    These files can in turn be viewed on the UCSC genome browser.

    Finally print out a USCS custom track line for the bigWig file.

    Steps:

        1) Intersect with 3UTR (or whatever original bedfile was used)
        2) Do .bed -> .bedgraph
        3) Do .bedgraph -> bigWig
        5) Print USCS line

    As well, get a bedfile for the relative lengths:
        1) Parse lengths_datset in output
        2) Save to output-dir
    """

    # Define shortcut variables
    dsets = settings.bigwig_datasets # the datsets we process
    utrfile_path = annotation.utrfile_path
    savedir = settings.bigwig_savedir
    url = settings.bigwig_url

    # hg19 and bed->bigwig files
    hg19 = '/users/rg/jskancke/phdproject/3UTR/the_project/source_bedfiles/hg19'
    bedGtoBigWig = '/users/rg/jskancke/programs/other/bedGraphToBigWig'

    for_wig = []
    for_length = []
    for_polyA = []

    print('Checking if files are available...\n')
    for dset in dsets:
        polyAfile = 'polyA_reads_'+ dset +'_processed_mapped.bed'
        readfile = 'reads_'+ dset +'.bed'
        lengthfile = 'length_' + dset
        clusterfile = 'polyA_' + dset

        polyA_read_path = os.path.join(here, 'temp_files', polyAfile)
        readpath = os.path.join(here, 'temp_files', readfile)
        lengthpath = os.path.join(here, 'output', lengthfile)
        cluster_polyA_path = os.path.join(here, 'output', clusterfile)

        # Check if the files are there; if so add them to list
        for bedfile in [polyA_read_path, readpath]:
            try:
                open(bedfile, 'rb')
                for_wig.append(bedfile)
            except:
                print('Not found or no access:\n{0}\nSkipping file...\n'\
                      .format(bedfile))

        # check length-file
        try:
            open(lengthpath, 'rb')
            for_length.append(lengthpath)
        except:
            print('Not found or no access:\n{0}\nSkipping file...\n'\
                  .format(lengthpath))

        # check polyA-cluster file
        try:
            open(cluster_polyA_path, 'rb')
            for_polyA.append(cluster_polyA_path)
        except:
            print('Not found or no access:\n{0}\nSkipping file...\n'\
                  .format(cluster_polyA_path))

    print('Done checking.\n')

    ######################## BedTools reads -> BigWig ##################
    # do the bedTools work to make BigWig files for the genome browser
    for dset in for_wig:

        (dirname, filename) = os.path.split(dset)

        # Adapt file names
        if filename.startswith('polyA'):
            shortname = '_'.join(filename.split('_')[:-2])
            co = '255,0,0' # red color
        if filename.startswith('reads'):
            shortname = os.path.splitext(filename)[0]
            co = '0,0,255' # blue color

        bedG_path = os.path.join(savedir, shortname + '.bedGraph')
        bigW_path = os.path.join(savedir, shortname + '.bigWig')

        dset_sort = dset+'_sorted'
        # Don't sort file if sorted file exists. DANGEROUS BUT FASTER.
        if not os.path.isfile(dset_sort):
            print('Sorting {0} ...'.format(dset))
            sort_cmd = ['sort', '-k', '1,1', dset]
            e = Popen(sort_cmd, bufsize=-1, stdout = open(dset_sort, 'wb'))
            e.wait()

        if os.path.isfile(bigW_path):
            print('')
            print('Found: {0}\nDelete it if you want to remake BigWig from source.'\
                  .format(bigW_path))
            print('')

        # Only make bigwig files if they don't exist. Delete to remake.
        else:
            cmd_intersect = ['intersectBed', '-a', dset_sort, '-b', utrfile_path]
            cmd_bedGraph = ['genomeCoverageBed', '-bg', '-i', 'stdin', '-g', hg19]
            cmd_bedGtoBW = [bedGtoBigWig, bedG_path, hg19, bigW_path]

            f = Popen(cmd_intersect, stdout = PIPE)
            g = Popen(cmd_bedGraph, stdin = f.stdout, stdout = open(bedG_path, 'wb'))

            print('Running intersectBed + genomeCoverageBed on {0} ...'\
                  .format(dset_sort))

            g.wait() # wait for bedGraph to finish
            h = Popen(cmd_bedGtoBW)

            print('Running bedGraphToBigWig on {0} ...'.format(bedG_path))
            h.wait() # wait for bigWig to finish

        UCSC = 'track type=bigWig name="{0}" description="{0}" bigDataUrl={1} '\
        'color={2} visibility=2'\
                .format(shortname, os.path.join(url, shortname + '.bigWig'), co)

        print('Provide this USCS custom track line:\n')
        print UCSC
        print('')

    ############### Parse Length output ##########################
    for dset in for_length:

        (dirname, filename) = os.path.split(dset)

        shortname = 'epsilon_' + filename

        length_path = os.path.join(savedir, shortname)

        # Header of the bedfile
        header = 'track type=bed name="{0}" description="{0}" color=0,255,0'\
                .format(shortname, os.path.join(url, length_path))
        print('Parsing the file with epsilon-lengths ...')

        # parse the output length-file to make a bedfile for how long the 3UTRs
        # are given the read coverage
        parse_length(dset, length_path, header)

        print('')
        print('Upload this file to the genome browser:\n')
        print length_path
        print('')

        # Get the upstream/downstream stuff as well
        length_ud = 'udstream_' + shortname
        length_ud_path = os.path.join(savedir, length_ud+'.bedGraph')
        header = 'track type=bed name="{0}" description="{0}" color=0,255,0'\
                .format(length_ud)
        length_udstream_covr(dset, length_ud_path, header)

        # Get the location of the PAS and their distances.
        length_PAS = 'PAS_distance_' + shortname
        length_PAS_path = os.path.join(savedir, length_PAS+'.bedGraph')
        header = 'track type=bed name="{0}" description="{0}" color=158,35,135'\
                .format(length_PAS)

        pas_dist_bed_length(dset, length_PAS_path, header)


    ############### Parse the polyA output ########################
    for dset in for_polyA:
        (dirname, filename) = os.path.split(dset)

        shortname = 'clusters_' + filename

        cluster_path = os.path.join(savedir, shortname)

        header = 'track type=bed name="{0}" description="{0}" color=0,255,255'\
                .format(shortname, os.path.join(url, cluster_path))
        print('Parsing the file with polyA-clusters ...')

        parse_clusters(dset, cluster_path, header)

        print('')
        print('Upload this file to the genome browser:\n')
        print cluster_path
        print('')

        # FOR DEBUGGING #

        # Get the upstream, downstraem coverage 
        cluster_ud = 'udstream_' + shortname
        cluster_ud_path = os.path.join(savedir, cluster_ud+'.bedGraph')
        header = 'track type=bed name="{0}" description="{0}" color=0,255,255'\
                .format(cluster_ud)
        cluster_udstream_covr(dset, cluster_ud_path, header)

        # Get the location of the PAS and their distances.
        cluster_PAS = 'PAS_distance_' + shortname
        cluster_PAS_path = os.path.join(savedir, cluster_PAS+'.bedGraph')
        header = 'track type=bed name="{0}" description="{0}" color=128,35,175'\
                .format(cluster_PAS)

        pas_dist_bed_clusters(dset, cluster_PAS_path, header)

    # While debugging: print out tracks for the coverage upstream and downstream
    # of the length and polyA output respectively.

def length_udstream_covr(in_length, out_length, out_header):
    infile = open(in_length, 'rb')
    orig_head = infile.next() # remove the header before reading

    outfile = open(out_length, 'wb')
    # print the header of the output file
    outfile.write(out_header+'\n')

    for line in infile:
        (chrm, beg, end, ext, strand, ID, eps_coord, eps_rz, eps_ds_covr,
         eps_us_covr) = line.split()[:10]
        # Skip the utrs without coverage
        if eps_coord == 'NA':
            continue
        eps_coord = int(beg) + int(eps_coord)
        ds_covr = float(eps_ds_covr)
        us_covr = float(eps_us_covr)
        # Write 50 lines upstream/downstream (depending on strand) in bedgraph
        # format for easy viewing on the browser

        # Strand differene: what is ustream dn what is dstream
        if strand == '+':
            beg = eps_coord - 50
            end = eps_coord
            outfile.write('\t'.join([chrm, str(beg), str(end), str(us_covr)])+'\n')

            beg = eps_coord
            end = eps_coord + 50
            outfile.write('\t'.join([chrm, str(beg), str(end), str(ds_covr)])+'\n')

        if strand == '-':
            beg = eps_coord - 50
            end = eps_coord
            outfile.write('\t'.join([chrm, str(beg), str(end), str(ds_covr)])+'\n')

            beg = eps_coord
            end = eps_coord + 50
            outfile.write('\t'.join([chrm, str(beg), str(end), str(us_covr)])+'\n')

    outfile.close()
    infile.close()

def cluster_udstream_covr(in_cluster, out_cluster, out_header):
    infile = open(in_cluster, 'rb')
    orig_head = infile.next() # remove the header before reading

    outfile = open(out_cluster, 'wb')
    # print the header of the output file
    outfile.write(out_header+'\n')

    for line in infile:
        (chrm, beg, end, ID, pA_nr, strand, cl_coord, supp_reads, ds_covr,
         us_covr) = line.split()[:10]

        # YOU ARE WRITING THE POLYA CLUSTER DOWNSTREAM UPSTREAM
        ds_covr = float(ds_covr)
        us_covr = float(us_covr)
        cl_coord = int(cl_coord)

        # Write 50 lines upstream/downstream (depending on strand) in bedgraph
        # format for easy viewing on the browser

        # Strand differene: what is ustream dn what is dstream
        if strand == '+':
            beg = cl_coord - 50
            end = cl_coord
            outfile.write('\t'.join([chrm, str(beg), str(end), str(us_covr)])+'\n')

            beg = cl_coord
            end = cl_coord + 50
            outfile.write('\t'.join([chrm, str(beg), str(end), str(ds_covr)])+'\n')

        if strand == '-':
            beg = cl_coord - 50
            end = cl_coord
            outfile.write('\t'.join([chrm, str(beg), str(end), str(ds_covr)])+'\n')

            beg = cl_coord
            end = cl_coord + 50
            outfile.write('\t'.join([chrm, str(beg), str(end), str(us_covr)])+'\n')

    outfile.close()
    infile.close()


def pas_dist_bed_clusters(in_PAS, out_PAS, header):
    """
    Write bedfile with positions and distances of PASes to their clusters so you
    can check the result in the genome browser.
    """
    infile = open(in_PAS, 'rb')
    in_header = infile.next()

    outfile = open(out_PAS, 'wb')
    outfile.write(header + '\n')

    for line in infile:
        (chrm, beg, end, d, d, strand, cl_coord) = line.split('\t')[:7]
        (cl_PAS, cl_PAS_dist, rpkm) = line.split('\t')[-3:]
        cl_PAS = cl_PAS.split(' ')
        cl_PAS_dist = cl_PAS_dist.split(' ')

        for (pas, dist) in zip(cl_PAS, cl_PAS_dist):
            if pas != 'NA':
                # get the beg and end of the PAS
                if strand == '+':
                    pas_beg = str(int(cl_coord) - int(dist)-7)
                    pas_end = str(int(cl_coord) - int(dist)-1)

                if strand == '-':
                    pas_beg = str(int(cl_coord) + int(dist)-1)
                    pas_end = str(int(cl_coord) + int(dist)+5)

                pas_and_dist = pas+'_'+dist

                outfile.write('\t'.join([chrm, pas_beg, pas_end, pas_and_dist])+'\n')

    outfile.close()

def pas_dist_bed_length(in_PAS, out_PAS, header):
    """
    Write bedfile with positions and distances of PASes to their clusters so you
    can check the result in the genome browser.
    """
    infile = open(in_PAS, 'rb')
    in_header = infile.next()

    outfile = open(out_PAS, 'wb')
    outfile.write(header + '\n')

    for line in infile:
        (chrm, beg, end, d, strand, d, eps_coord) = line.split('\t')[:7]
        (eps_PAS, eps_PAS_dist, rpkm) = line.split('\t')[-3:]
        eps_PAS = eps_PAS.split(' ')
        eps_PAS_dist = eps_PAS_dist.split(' ')

        for (pas, dist) in zip(eps_PAS, eps_PAS_dist):
            if pas != 'NA':
                # get the beg and end of the PAS
                if strand == '+':
                    pas_beg = str(int(eps_coord)+int(beg) - int(dist)-7)
                    pas_end = str(int(eps_coord)+int(beg) - int(dist)-1)

                if strand == '-':
                    pas_beg = str(int(eps_coord)+int(beg) + int(dist)-1)
                    pas_end = str(int(eps_coord)+int(beg) + int(dist)+5)

                pas_and_dist = pas+'_'+dist

                outfile.write('\t'.join([chrm, pas_beg, pas_end, pas_and_dist])+'\n')

    outfile.close()

    pass

def main():
    """
    This method is called if script is run as __main__.
    """

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
    simulate = False
    #simulate = True
    #settings.bigwig = False

    # This option should be set only in case of debugging. It makes sure you
    # just run chromosome 1 and only extract a tiny fraction of the total reads.
    DEBUGGING = True
    #DEBUGGING = False
    if DEBUGGING:
        settings.DEBUGGING()

    #settings.polyA = False

    # The program reads a lot of information from the annotation. The annotation
    # object will hold this information (file-paths and datastructures).
    print('Reading settings ...\n')
    annotation = Annotation(settings.annotation_path)

    # Check if 3UTRfile has been made or provided; if not, get it from annotation
    annotation.utrfile_path = get_utr_path(settings, beddir)

    # Get dictionary with (chrm, beg, end, strand) values for each 3utr-exon key
    # NOTE: for the 3'terminal exons, the beg/end is the extended value. For
    # internal exons, the beg/end are not extended.
    print('Making 3UTR data structures ...\n')
    annotation.utr_exons = annotation.get_utrdict()

    # You extract all annotated polyA sites into a bedfile: a_polyA_sites_path
    annotation.a_polyA_sites_path = get_a_polyA_sites_path(settings, beddir)

    # You intersect the annotated polyA files with the 3UTR bedfile you got from
    # the annotation. Put the intersected annotated polyA sites in a dictionary.
    print('Making data structures for annotation poly(A) sites ...\n')
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
                         settings, annotation, DEBUGGING)

            ###### WORK IN PROGRESS
            akk = pipeline(dset_id, dset_reads, tempdir, output_dir, utr_seqs,
                           settings, annotation, DEBUGGING)

            #result = my_pool.apply_async(pipeline, arguments)
            #results.append(result)

        # Wait for all procsses to finish
        debug()
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

    #if not simulate:

        #file_path_dict = cPickle.load(open(pickled_final, 'rb'))

        #final_outp_polyA = file_path_dict['final_output_polyA']
        #final_outp_length = file_path_dict['final_output_length']
        #coverage = file_path_dict['coverage']

    # if set, make bigwig files
    #if settings.bigwig:
        #make_bigwigs(settings, annotation, here)

if __name__ == '__main__':
    main()

# XXX NOTE TO SELF XXX
# Before the holiday you updated a lot of docstrings and you started with the
# documentation  (see doc folder). You had a meeting with Roderic: see message
# below.

# TODO: the path to the gem-mapper index... if people are going to use this,
# they need to provide the path in the SETTINGS file. If you do this yourself,
# it will be easier to switch to an index file on your own hard disk, makign
# loading the index into memory much faster.
# TODO: go through the entire program and improve documentation.
# TODO: write external documentation file for the program
# TODO: generate the figures. one button!
# TODO: the two below should be part of the final analysis script; one button!
# TODO: verify the polyA reads clusters with polyAdb + other results
# TODO: verify the 3UTR length with annotation + smth else
# TODO: make a function that returns all regions (5UTR, intron, CDS)
# unoverlapping other regions.

# TODO Roderic ideas:
    # 1) Run without annotation: you can create bed-file from sections around the
    # polyA reads. This would be awsome to try out.
    # 2) Must be able to run without outputting to 'length' file. In this way,
    # the program can be run to output info about length from an annotation,
    # polyA reads from annotation or not, or just polyA reads without an
    # annotation. This will require some re-writing, but not awfully much.

# TODO: Roderic conclusion:
    # next time he wants some basic statistics. The program seems ready, so it's
    # time to do some numbers on the output.

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
    # 5) Do statistics on the number of times polyA reads fall in the opposite
    # strand (for all clusters and for clusters with more than 1 read)
    # 6) Give the absolute number of poly(A) reads in the different
    # compartments. Seems like there are a lot more poly(A) reads in the
    # cytoplasm than in the nucleus


"""
Script for displaying and summarizing the results from utr_finder.py.
"""

import os
import ConfigParser
import sys

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
        print ('Warning: debugging mark present. This is a bug.')
        pass
########################################

class Length(object):
    """
    For UTR objects from the 'length' output file in the 'output' directory.
    """

    def __init__(self, input_line):

        # Read all the parameters from line
        (chrm, beg, end, utr_extended_by, strand, utr_ID, epsilon_coord,
        epsilon_rel_size, epsilon_downstream_covrg, epsilon_upstream_covrg,
        annot_downstream_covrg, annot_upstream_covrg, epsilon_PAS_type,
        epsilon_PAS_distance, utr_RPKM) = input_line.split('\t')

        self.chrm = chrm
        self.beg = str_to_intfloat(beg)
        self.end = str_to_intfloat(end)
        self.extended_by = str_to_intfloat(utr_extended_by)
        self.strand = strand
        self.utr_ID = utr_ID
        self.eps_coord = str_to_intfloat(epsilon_coord)
        self.eps_rel_size = str_to_intfloat(epsilon_rel_size)
        self.eps_downstream_covrg = str_to_intfloat(epsilon_downstream_covrg)
        self.eps_upstream_covrg = str_to_intfloat(epsilon_upstream_covrg)
        self.annot_downstream_covrg = str_to_intfloat(annot_downstream_covrg)
        self.annot_upstream_covrg = str_to_intfloat(annot_upstream_covrg)

        # PAS type and PAS distance are space-delimited
        PAS_type = epsilon_PAS_type.split(' ')
        self.eps_PAS_type = [str_to_intfloat(pas) for pas in PAS_type]

        PAS_distance = epsilon_PAS_distance.split(' ')
        self.eps_PAS_distance = [str_to_intfloat(dist) for dist in PAS_distance]

        self.RPKM = str_to_intfloat(utr_RPKM.strip())

    def __repr__(self):
        return self.utr_ID[-8:]

    def __str__(self):
        return "\nChrm\t{0}\nBeg\t{1}\nEnd\t{2}\nStrand\t{3}\n"\
                .format(self.chrm, self.beg, self.end, self.strand)

class PolyA(object):
    """
    For polyA objects from the 'polyA' output file in the 'output' directory
    """

    def __init__(self, input_line):
        #
        (chrm, beg, end, utr_ID, polyA_number, strand, polyA_coordinate,
         number_supporting_reads, coverage_50nt_downstream,
         coverage_50nt_upstream, annotated_polyA_distance, nearby_PAS,
         PAS_distance) = input_line.split('\t')

        self.chrm = chrm
        self.beg = str_to_intfloat(beg)
        self.end = str_to_intfloat(end)
        self.utr_ID = utr_ID
        self.cluster_nr = str_to_intfloat(polyA_number)
        self.strand = strand
        self.polyA_coordinate = str_to_intfloat(polyA_coordinate)
        self.nr_support_reads = str_to_intfloat(number_supporting_reads)
        self.coverage_50nt_downstream = str_to_intfloat(coverage_50nt_downstream)
        self.coverage_50nt_upstream = str_to_intfloat(coverage_50nt_upstream)
        self.annotated_polyA_distance = str_to_intfloat(annotated_polyA_distance)

        # PAS type and distance are space-delimited
        PAS_type = nearby_PAS.split(' ')
        self.nearby_PAS = [str_to_intfloat(pas) for pas in PAS_type]

        PAS_distance = PAS_distance.split(' ')
        self.PAS_distance = [str_to_intfloat(dist) for dist in PAS_distance]

    def __repr__(self):
        return self.utr_ID[-8:]

    def __str__(self):
        return "\nChrm\t{0}\nBeg\t{1}\nEnd\t{2}\nStrand\t{3}\n#\t{4}\n"\
                .format(self.chrm, self.beg, self.end, self.strand, self.cluster_nr)

class Settings(object):
    """
    Convenience class.One instance is created from this class: it the pertinent
    settings parameters obtained from the UTR_SETTINGS file.
    """
    def __init__(self, settings_file, savedir, outputdir, here):

        conf = ConfigParser.ConfigParser()
        conf.optionxform = str
        conf.read(settings_file)

        self.datasets = conf.get('PLOTTING', 'datasets').split(':')
        self.savedir = savedir
        self.outputdir = outputdir
        self.here = here

    # Return the paths of the length files
    def length_files(self):
        return dict((d, os.path.join(self.here, self.outputdir,'length_'+d))
                    for d in self.datasets)

    # Return the paths of the polyA files
    def polyA_files(self):
        return dict((d, os.path.join(self.here, self.outputdir,'polyA_' + d))
                    for d in self.datasets)

    # Return the paths of the epsilon files
    def epsilon_files(self):
        return dict((d, os.path.join(self.here, self.outputdir,'cumul_'+d+'.stat'))
                    for d in self.datasets)

def get_lengths(settings):
    """
    Return a dictionary of class
    """
    length_files = settings.length_files()

    # Check if all length files exist or that you have access
    [verify_access(f) for f in length_files.values()]

    l_dict = {}
    for dset, f in length_files.items():
        file_obj = open(f, 'rb')
        header = file_obj.next()
        # For each dset, make a dictionary:[UTR_ID]-> Length object
        # line.split()[5] is the utr_id column
        l_dict[dset] = dict((line.split()[5], Length(line)) for line in file_obj)

    return l_dict


def get_polyAs(settings):
    """
    Return a dictionary of class
    """
    length_files = settings.polyA_files()

    # Check if all length files exist or that you have access
    [verify_access(f) for f in length_files.values()]

    l_dict = {}
    for dset, f in length_files.items():
        file_obj = open(f, 'rb')
        header = file_obj.next()
        l_dict[dset] = {}
        # For each dset, make a dictionary:[UTR_ID]-> Length object
        # line.split()[5] is the utr_id column
        for line in file_obj:
            identifier = '_'.join(line.split()[3:5])
            l_dict[dset][identifier] = PolyA(line)

    return l_dict


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

def str_to_intfloat(element):
    """
    Try to return int; if not try to return float; if not return unchanged.
    """

    try:
        return int(element)
    except ValueError:
        pass

    try:
        return float(element)
    except ValueError:
        pass

    return element

def utr_length_comparison(settings):
    """
    * compare 3UTR length in general
    * compare 3UTR length UTR-to-UTR
    """
    # 1) Box plot of 3UTR lengths in all utrs

    # 2) Box plot of 

    pass

def polyadenylation_comparison(settings):
    """
    * compare 3UTR polyadenylation in general
    * compare 3UTR polyadenylation UTR-to-UTR
    """

    pass

def reads(settings):
    """
    * number of reads for each compartment
    * number of putative poly(A) reads for each compartment
    * ratio of the two above
    * number of putative poly(A) reads that map to clusters
    * number of clusters with different sizes
    * number of clusters in 'wrong' direction, as a function of the numbers of
    * poly(A) reads in that cluster
    """

    pass

def UTR_processing(settings):
    """
    * number of annotated 3UTRs (and number of genes)
        * multi and single exon
    * number of surviving 3UTRs from clustering
        * multi and single exon
    * number of survivers from overlap-screening
        * multi and single exon
    """

    pass

def other_methods_comparison(settings):
    """
    * number of matching polyA sites in polA_db
    * number of matching polyA sites in recent publications
    * number of novel polyA sites, per tissue/cell_line
    """

    pass

def data_annotation_correspondence(settings):
    """
    * number of utrs found compared to total
    * number of utrs that seem longer than annotation
    """

    pass

def rouge_polyA_sites(settings):
    """
    * polyA clusters in the different annotation regions
    * comparison with polyA_db and recent publications.
    * "extended 3UTRs": non-overlapping 3'end downstream regions with read
        coverage and poly(A) clusters
    * novel polyA sites in annotated 3UTRs: compared to polyA_db and
        publications
    * presence of PAS and PAS type for novel sites
    """

    pass

def classic_polyA_stats(settings):
    """
    * distances from polyA cluster to PAS site
    * PAS variant distribution
    * variance within polyA clusters
    * degree of usage of early and late polyA sites
    """

    pass

def epsilon_evaluation(settings):
    """
    Statistics for how epsilon is determined.
    """

    pass


def main():
    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    # Directory paths for figures and where the output lies
    (savedir, outputdir) = [os.path.join(here, d) for d in ('figures', 'output')]

    # Read UTR_SETTINGS
    settings = Settings(os.path.join(here, 'UTR_SETTINGS'), savedir, outputdir, here)

    # Get the length-datasets
    length_dict = get_lengths(settings)
    # Get the polyA-datasets
    polyA_dict = get_polyAs(settings)

    utr_length_comparison(settings, length_dict)

    polyadenylation_comparison(settings)

    reads(settings)

    data_annotation_correspondence(settings)

    UTR_processing(settings)

    other_methods_comparison(settings)

    rouge_polyA_sites(settings)

if __name__ == '__main__':
    main()

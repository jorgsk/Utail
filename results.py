"""
Script for displaying and summarizing the results from utr_finder.py.
"""

import os
import ConfigParser

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
        return [os.path.join(self.here, self.outputdir,'length_' + d) for d
                in self.datasets]

    # Return the paths of the polyA files
    def polyA_files(self):
        return [os.path.join(self.here, self.outputdir,'polyA_' + d) for d
                in self.datasets]

    # Return the paths of the epsilon files
    def epsilon_files(self):
        return [os.path.join(self.here, self.outputdir,'cumul_'+d+'.stat') for d
                in self.datasets]

def utr_length_comparison(settings):
    """
    * compare 3UTR length in general
    * compare 3UTR length UTR-to-UTR
    """
    # 1) Box plot of 3UTR lengths in all utrs

    # Data challenge: I need each 3UTR in each dataset. Double dictionary,
    # anyone? Or a class? mydict[dataset_id][3utr_id](tuple.of.things)

    # Alternatively the mydict[datset_id] as dict, with 3utr_id being class
    # objects instead of tuples; they would have all the fields as found in the
    # output file. (ps: I'm going with this one :))

    # 2) Box plot of 
    pass

def polyadenylation_comparison(settings):
    """
    * conpare 3UTR polyadenylation in general
    * conpare 3UTR polyadenylation UTR-to-UTR
    """

    #
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
    * comparison to polyA_db and recent publications.
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
    debug()

    # Get the length-datasets
    length_dict = get_lengths(settings)

    utr_length_comparison(settings)

    polyadenylation_comparison(settings)

    reads(settings)

    data_annotation_correspondence(settings)

    UTR_processing(settings)

    other_methods_comparison(settings)

    rouge_polyA_sites(settings)

if __name__ == '__main__':
    main()

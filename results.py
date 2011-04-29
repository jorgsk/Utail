"""
Script for displaying and summarizing the results from utr_finder.py.
"""

import os
import ConfigParser
import sys
import matplotlib.pyplot as plt

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

class PolyaCluster(object):
    """
    Class that contains helper functions to work on datasets.
    """

    def __init__(self, dset_list, dset_name):
        self.pAclusters = dset_list
        self.name = dset_name

    def get_supporting_reads(self):
        return [cluster.nr_support_reads for cluster in self.pAclusters]

    def get_coverage_downstream(self):
        return [cluster.coverage_50nt_downstream for cluster in self.pAclusters]

    def get_coverage_upstream(self):
        return [cluster.coverage_50nt_upstream for cluster in self.pAclusters]

    def get_annotated_distance(self):
        return [cluster.annotated_polyA_distance for cluster in self.pAclusters]

    def get_rpkm(self):
        return [cluster.rpkm for cluster in self.pAclusters]


class LengthDataset(object):
    """
    Class that contains helper functions to work on datasets.
    """

    def __init__(self, dset_list, dset_name):
        self.utrs = dset_list
        self.name = dset_name

    def get_eps_length(self, minRPKM, IDs=False):
        """
        Return a list of epsilon-lengths for this datsets. Specify a minimum
        RPKM for the UTR in question. Optionally, specify a list or set of
        IDs that should be included.
        """

        eps_lengths = []

        # If a list of utrIDs have been supplied, include only those in that list
        if IDs:
            for utr in self.utrs:
                if (utr.eps_rel_size != 'NA') and (utr.RPKM > minRPKM):
                    if utr.ID in IDs:
                        eps_lengths.append(utr.eps_rel_size)

        # Include all utr objects
        else:
            for utr in self.utrs:
                if utr.eps_rel_size != 'NA' and utr.RPKM > minRPKM:
                    eps_lengths.append(utr.eps_rel_size)

        return eps_lengths

    def expressed_IDs(self):
        """
        Return list (or set) of IDs of utrs that are expressed in this dataset. The
        criterion for 'expressed' is that eps_rel_size is not 'NA'.
        """
        return [utr.ID for utr in self.utrs if utr.eps_rel_size != 'NA']


class UTR(object):
    """
    For UTR objects from the 'length' output file in the 'output' directory.
    """

    def __init__(self, input_line):

        # Read all the parameters from line
        (chrm, beg, end, utr_extended_by, strand, ID, epsilon_coord,
        epsilon_rel_size, epsilon_downstream_covrg, epsilon_upstream_covrg,
        annot_downstream_covrg, annot_upstream_covrg, epsilon_PAS_type,
        epsilon_PAS_distance, utr_RPKM) = input_line.split('\t')

        self.chrm = chrm
        self.beg = str_to_intfloat(beg)
        self.end = str_to_intfloat(end)
        self.extended_by = str_to_intfloat(utr_extended_by)
        self.strand = strand
        self.ID = ID
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
        return self.ID[-8:]

    def __str__(self):
        return "\nChrm\t{0}\nBeg\t{1}\nEnd\t{2}\nStrand\t{3}\n"\
                .format(self.chrm, self.beg, self.end, self.strand)

class Cluster(object):
    """
    For polyA cluster objects from the 'polyA' output file in the 'output' directory
    """

    def __init__(self, input_line):
        #
        (chrm, beg, end, ID, polyA_number, strand, polyA_coordinate,
         number_supporting_reads, coverage_50nt_downstream,
         coverage_50nt_upstream, annotated_polyA_distance, nearby_PAS,
         PAS_distance, rpkm) = input_line.split('\t')

        self.chrm = chrm
        self.beg = str_to_intfloat(beg)
        self.end = str_to_intfloat(end)
        self.ID = ID
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

        self.rpkm = str_to_intfloat(rpkm.strip())

    def __repr__(self):
        return self.ID[-8:]+'_'+str(self.cluster_nr)

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

class Plotter(object):
    """
    Collection of plot-methods
    """

    def boxplot(self, arrays, dsets, ylim):
        """
        A basic box plot
        """

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.boxplot(arrays)
        ax.set_xticklabels(dsets)
        #ax.set_ylabel('Some label', size=20)
        #ax.set_title('Some title')
        ax.set_ylim(*ylim)
        fig.show()

    def scatterplot(self, dset1, dset2, label1, label2, title, xlim, ylim):
        """
        A basic scatter plot
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(dset1, dset2)
        ax.set_xlabel = label1
        ax.set_ylabel = label2
        ax.set_title(title)
        ax.set_xlim = xlim
        ax.set_ylim = ylim
        fig.show()

    def triangleplot_scatter(self, datas, titles):
        """
        A scatterplot of multiple variables. Plots all the variables in 'datas'
        aginst each other. To remove redundancy, only the upper triangle of the
        matrix is shown.
        """
        var_nr = len(datas)
        plots = range(var_nr)

        fig = plt.figure()
        axis_nr = 0

        max_plot =  var_nr*var_nr

        # Determine which should have xlabel and which should have ylabel
        xlables = range(max_plot-var_nr+1, max_plot+1)
        ylables = range(1,max_plot,var_nr)

        remove_indices = []
        for v in range(2, var_nr+1):
            mymax = v*v
            remove_indices += range(v,mymax,var_nr)

        for indx_1 in plots:
            for indx_2 in plots:
                axis_nr += 1

                ax = fig.add_subplot(var_nr, var_nr, axis_nr)
                if axis_nr not in remove_indices:
                    ax.scatter(datas[indx_2], datas[indx_1])

                # Set only labels where they should be
                if axis_nr in xlables:
                    ax.set_xlabel(titles[indx_2], size=15)
                if axis_nr in ylables:
                    ax.set_ylabel(titles[indx_1], size=15)

                # Remove axis when not needed
                if axis_nr not in xlables:
                    ax.set_xticks([])
                if axis_nr not in ylables:
                    ax.set_yticks([])

        fig.show()
        debug()


def get_lengths(settings):
    """
    Return a list of LengthDataset instances. Each LengthDataset instance is
    instanciated with a list of UTR objects and the name of the datset.
    """
    length_files = settings.length_files()

    # Check if all length files exist or that you have access
    [verify_access(f) for f in length_files.values()]

    lengths = []
    for dset_name, f in length_files.items():
        file_obj = open(f, 'rb')
        header = file_obj.next()
        utr_list = [UTR(line) for line in file_obj]

        lengths.append(LengthDataset(utr_list, dset_name))

    return lengths


def get_pAclusters(settings):
    """
    Return a list of polyA cluster instances. Each polyA cluster instance is
    instanciated with a list of Cluster objects.
    """
    length_files = settings.polyA_files()

    # Check if all length files exist or that you have access
    [verify_access(f) for f in length_files.values()]

    pAclusters = []
    for dset_name, f in length_files.items():
        file_obj = open(f, 'rb')
        header = file_obj.next()
        cluster_list = [Cluster(line) for line in file_obj]

        pAclusters.append(PolyaCluster(cluster_list, dset_name))

    return pAclusters


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

def utr_length_comparison(settings, lengths):
    """
    * compare 3UTR length in general (with boxplot)
    * compare 3UTR length UTR-to-UTR (with boxplot)
    """
    # Get the plotter object
    p = Plotter()

    # 1) Box plot of 3UTR lengths in all utrs
    # get dset names and an array of arrays [[], [],...,[]], of lengths
    names = [dset.name for dset in lengths]
    eps_lengths = [dset.get_eps_length(minRPKM=2) for dset in lengths]

    ylim = (0, 1.2)
    #p.boxplot(eps_lengths, names, ylim)

    # 2) Box plot of 3UTR lengths for utrs expressed in all datasets.
    # Get the IDs of utrs that are expressed in each datset
    utrID_sets = [set(dset.expressed_IDs()) for dset in lengths] # list of sets

    # Do the intersection of these sets to get only the ones that are expressed
    # in all datsets
    common = set.intersection(*utrID_sets)

    # Get the lengths of only the commonly expressed utrs:
    eps_lengths_c = [dset.get_eps_length(minRPKM=2, IDs=common) for dset in lengths]
    p.boxplot(eps_lengths_c, names, ylim)

def polyadenylation_comparison(settings, polyAs):
    """
    * compare 3UTR polyadenylation in general
    * compare 3UTR polyadenylation UTR-to-UTR
    """
    # How to compare polyAdenylation between compartments?
    p = Plotter()

    # What is a good framework to call correlations from?
    # Is every 
    #
    # 0) Control: what is the correlation between rpkm and # of supporting
    # polyA-reads (first, second, third, etc)? First: scatterplots.
    for dset in polyAs:

        # Get array of data for triangle-scatterplot
        arrays = []
        arrays.append(dset.get_supporting_reads())
        arrays.append(dset.get_coverage_downstream())
        arrays.append(dset.get_coverage_upstream())
        #arrays.append(dset.get_annotated_distance())
        arrays.append(dset.get_rpkm())

        # Titles for the above variables
        #titles = ['Nr Supporting Reads', 'Downstream Coverage',
                  #'Upstream Coverage', 'Distance to annotated TTS', 'RPKM']
        titles = ['Nr Supporting Reads', 'Downstream Coverage',
                  'Upstream Coverage', 'RPKM']

        p.triangleplot_scatter(arrays, titles)


    def get_annotated_distance(self):
        return [cluster.annotated_polyA_distance for cluster in self.pAclusters]
    debug()
    # 0) Control: what is the correlation between # supporting reads; change in
    # coverage; PAS-type; PAS-distance; rpkm; and total reads for experiment?
    #
    # 1) Usage of multiple polyadenylation sites in each UTR
    # 2) When multiple sites, what is the read frequency of the sites?
    # 3) Compare utr-by-utr: is there a difference in polyA usage?

    pass

def reads(settings):
    """
    * number of reads for each compartment
    * number of putative poly(A) reads for each compartment
    * ratio of the two above
    * number of putative poly(A) reads that map to pAclusters
    * number of pAclusters with different sizes
    * number of pAclusters in 'wrong' direction, as a function of the numbers of
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
    * polyA pAclusters in the different annotation regions
    * comparison with polyA_db and recent publications.
    * "extended 3UTRs": non-overlapping 3'end downstream regions with read
        coverage and poly(A) pAclusters
    * novel polyA sites in annotated 3UTRs: compared to polyA_db and
        publications
    * presence of PAS and PAS type for novel sites
    """

    pass

def classic_polyA_stats(settings):
    """
    * distances from polyA cluster to PAS site
    * PAS variant distribution
    * variance within polyA pAclusters
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
    lengths = get_lengths(settings)
    # Get the polyA-datasets (polyA cluster datasets)
    polyAs = get_pAclusters(settings)

    #utr_length_comparison(settings, lengths)

    polyadenylation_comparison(settings, polyAs)

    reads(settings)

    data_annotation_correspondence(settings)

    UTR_processing(settings)

    other_methods_comparison(settings)

    rouge_polyA_sites(settings)

if __name__ == '__main__':
    main()

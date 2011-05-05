"""
Script for displaying and summarizing the results from utr_finder.py.
"""

from __future__ import division
import os
import ConfigParser
import sys
import matplotlib.pyplot as plt
from operator import attrgetter
import math

import numpy as np

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

# Horrible, global variables

first_pas = 'AATAAA'
second_pas = 'ATTAAA'
top_pas = set([first_pas, second_pas])
lower_pas = set(['AGTAAA', 'TATAAA', 'CATAAA', 'GATAAA', 'AATATA', 'AATACA',
                 'AATAGA', 'AAAAAG', 'ACUAAA'])

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
        return [cluster.dstream_covrg for cluster in self.pAclusters]

    def get_coverage_upstream(self):
        return [cluster.ustream_covrg for cluster in self.pAclusters]

    def get_annotated_distance(self):
        return [cluster.annotated_polyA_distance for cluster in self.pAclusters]

    def get_rpkm(self):
        return [cluster.rpkm for cluster in self.pAclusters]


class UTRDataset(object):
    """
    Class that contains helper functions to work on datasets.
    """

    def __init__(self, dset_dict, dset_name):
        self.utrs = dset_dict
        self.name = dset_name

    def __repr__(self):
        return self.name

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

        # The UTR potentially has polyA read clusters. They will be added to
        # this list. The number of clusters will also be added.
        self.clusters = []
        self.cluster_nr = 0

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
         number_supporting_reads, dstream_covrg, ustream_covrg,
         annotated_polyA_distance, nearby_PAS, PAS_distance,
         rpkm) = input_line.split('\t')

        self.cluster_nr = str_to_intfloat(polyA_number)
        self.polyA_coordinate = str_to_intfloat(polyA_coordinate)
        self.nr_support_reads = str_to_intfloat(number_supporting_reads)
        self.dstream_covrg = str_to_intfloat(dstream_covrg)
        self.ustream_covrg = str_to_intfloat(ustream_covrg)
        self.annotated_polyA_distance = str_to_intfloat(annotated_polyA_distance)

        # PAS type and distance are space-delimited
        PAS_type = nearby_PAS.split(' ')
        self.nearby_PAS = [str_to_intfloat(pas) for pas in PAS_type]

        PAS_distance = PAS_distance.split(' ')
        self.PAS_distance = [str_to_intfloat(dist) for dist in PAS_distance]

    def __repr__(self):
        return str(self.polyA_coordinate)

    def __str__(self):
        return "\nCl nr\t{0}\nCoord\t{1}\nRead nr\t{2}\n"\
                .format(self.cluster_nr,self.polyA_coordinate,self.nr_support_reads)

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
        ylables = range(1, max_plot, var_nr)

        remove_indices = []
        for v in range(2, var_nr+1):
            mymax = v*v
            remove_indices += range(v, mymax, var_nr)

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

    def last_three_clustersites(self, clus_list):
        """
        The input are the last X clusters (counting from the end -- 'first'
        is the last cluster in the utr). You want to make a sublot with two
        rows, each with boxplots: first row are the coverage ratio of the first,
        second, and third cluster, and the second row the read coverage.

        Clus list comes in order: element 0 is the 3'most element, and the -1
        element is the 5' most.

        The data has been normalized. For each gene, the feature (ud_ratio or
        coverage count), are normalized to the largest value. After
        normalization, all values are subtracted by 1; in this way, a negative
        value indicates that downstream is larger than upstream. Further, a cap
        is set on the values before normalization. Otherwise, since we divide by
        numbers close to zero, values can be very high. A ratio cap is set to 10.
        """

        fig = plt.figure()

        cluster_nr = len(clus_list)
        ratios = [dic['ud_ratio'] for dic in clus_list]
        supports = [dic['support'] for dic in clus_list]
        plot_array = ratios + supports # 1 to 6 traverse through

        max_plotnr = len(plot_array)
        row_nr = 2
        col_nr = max_plotnr/2

        for plotnr in range(1, max_plotnr+1):
            array = plot_array[plotnr-1]

            # mean and std
            median = np.median(array)
            std = np.std(array)
            mean = np.mean(array)
            ax = fig.add_subplot(2, col_nr, plotnr)

            median = format(median, '.2f')
            std = format(std, '.2f')
            mean = format(mean, '.2f')
            lbl = 'median: '+median+'\nmean: '+mean+'\nstd: '+std

            ax.boxplot(array)

            # Plot text right onto the image
            ax.text(0.55,0.9,lbl, size=13)

            # Set y limits depending on if log(ratio) or read count
            if plotnr in range(1, cluster_nr+1):
                ax.set_ylim(-3.2, 13.2)

            if plotnr in range(cluster_nr+1, max_plotnr+1):
                ax.set_ylim(-1,60)

            if plotnr == 1:
                ax.set_ylabel('Log2-ratio of upstream/downstream coverage',size=20)

            if plotnr == max_plotnr - cluster_nr+1:
                ax.set_ylabel('Poly(A)-read count', size=20)

            if plotnr in range(cluster_nr+1, max_plotnr+1):
                ax.set_xlabel("Poly(A) cluster nr. {0} from 3' end"\
                              .format(plotnr-cluster_nr))

        plt.show()


def get_utrs(settings):
    """
    Return a list of UTR instances. Each UTR instance is
    instanciated with a list of UTR objects and the name of the datset.

    1) Read through the length file, create UTR objects for each line
    2) If present, read through the polyA file, updating the UTR object created
    in step 1

    """
    # Do everything for the length files
    length_files = settings.length_files()
    cluster_files = settings.polyA_files()

    # Check if all length files exist or that you have access
    [verify_access(f) for f in length_files.values()]

    all_utrs = []
    for dset_name in settings.datasets:

        lengthfile = open(length_files[dset_name], 'rb')
        clusterfile = open(cluster_files[dset_name], 'rb')

        # Get headers
        l_header = lengthfile.next()
        c_header = clusterfile.next()

        # Create the utr objects
        utr_dict = dict((line.split()[5], UTR(line)) for line in lengthfile)

        # Add cluster objects to the utr objects (line[3] is UTR_ID
        for line in clusterfile:
            utr_dict[line.split()[3]].clusters.append(Cluster(line))

        # Update the UTRs with some new info about their clusters
        for (utr_id, utr_obj) in utr_dict.iteritems():
            utr_dict[utr_id].cluster_nr = len(utr_obj.clusters)

        all_utrs.append(UTRDataset(utr_dict, dset_name))

    return all_utrs


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

def utr_length_comparison(settings, utrs):
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

def output_control(settings, dsets):
    """
    Control: what is the correlation between # supporting reads; change in
    coverage; PAS-type; PAS-distance; and rpkm?
    """
    p = Plotter()

    ## TRIANGLE PLOT ##

    #for dset in dsets:

        #nr_supp_reads = []
        #covrg_down = []
        #covrg_up = []
        #rpkm = []

        #for (utr_id, utr) in dset.utrs.iteritems():

            ## For all polyA clusters get
            ## 1) NR of supporting reads
            ## 2) Coverage downstream
            ## 3) Coverage upstream
            ## 4) RPKM

            #for cls in utr.clusters:
                #nr_supp_reads.append(cls.nr_support_reads)
                #covrg_down.append(cls.dstream_covrg)
                #covrg_up.append(cls.ustream_covrg)
                #rpkm.append(utr.RPKM)

        ## Titles for the above variables
        #titles = ['Nr Supporting Reads', 'Downstream Coverage',
                  #'Upstream Coverage', 'RPKM']

        #array = [nr_supp_reads, covrg_down, covrg_up, rpkm]
        #p.triangleplot_scatter(array, titles)

    # For the ones with 1, 2, 3, 4, ... polyA sites:
        #1) box-plots of the ratio of downstream/upstream 
        #2) box-plots of the # of covering reads

    # Upstream/downstream ratios of polyA sites.

    for dset in dsets:
        # Get the maximum nr of polyA sites in one utr for max nr of plots
        # Give a roof to the number of plots; 5?
        max_cluster = max (utr.cluster_nr for utr in dset.utrs.itervalues())

        # Git yourself some arrays 
        # first second thrid and fourth
        # nr of clusters
        clrs_nr = 4
        clus_list = [{'ud_ratio':[], 'support':[]} for val in range(clrs_nr)]

        for (utr_id, utr) in dset.utrs.iteritems():

            if utr.cluster_nr < clrs_nr:
                continue

            if utr.strand == '+':
                clu = sorted(utr.clusters, key=attrgetter('polyA_coordinate'))
                clu = clu[-clrs_nr:] # get only the last clrs_nr
                # for + strand, reverse clu
                clu = clu[::-1]

            if utr.strand == '-':
                clu = sorted(utr.clusters, key=attrgetter('polyA_coordinate'))
                clu = clu[:clrs_nr] # get only the first clrs_nr

            # The order of clsters in clu is '1st, 2nd, 3rd...'

            eps_end = utr.beg + utr.eps_coord

            # only get UTRs that have the final polyA cluster close to the
            # coverage end!
            if not (eps_end - 50 < clu[0].polyA_coordinate < eps_end + 50):
                continue

            # Normalize the ratios by the largest absolute deviation from 1
            ud_ratios = []

            for cl in clu:
                # Make 0 into an almost-zero ...
                if cl.dstream_covrg == 0:
                    cl.dstream_covrg = 0.01

                # do log2 ratios (if ==1, twice as big)
                udratio = math.log(cl.ustream_covrg/cl.dstream_covrg,2)

                ud_ratios.append(udratio)

            #maxratio = max(max(ud_ratios), abs(min(ud_ratios)))
            #norm_ud_ratios = [rat/maxratio for rat in ud_ratios]
            norm_ud_ratios = ud_ratios

            # Append the normailzed ratios to the arrays
            for (indx, norm_rat) in enumerate(norm_ud_ratios):
                clus_list[indx]['ud_ratio'].append(norm_rat)

            # Normalize the read support
            read_supp = [cl.nr_support_reads for cl in clu]
            #maxsupp = max(read_supp)
            #norm_read_supp = [supp/maxsupp for supp in read_supp]
            norm_read_supp = read_supp

            # Append normalized support ratios to arrays
            for (indx, norm_supp) in enumerate(norm_read_supp):
                clus_list[indx]['support'].append(norm_supp)

        # Do teh plots
        p.last_three_clustersites(clus_list)

        # RESULT you see the trend you imagined.

def polyadenylation_comparison(settings, utrs):
    """
    * compare 3UTR polyadenylation in general
    * compare 3UTR polyadenylation UTR-to-UTR
    """
    # How to compare polyAdenylation between compartments?
    p = Plotter()

    # What is a good framework to call correlations from?
    # Is every 
    #
    # It depends on PAS-type, PAS-distance, polyA-number (1, 2, .., last)

    # Do we observe loss in coverage after sucessive polyadenylation sites?
    # Indeed, this will depend on the expression of each isomer, but still we
    # should be able to see something.

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

    # Get the dsets with utrs from the length and polyA files
    dsets = get_utrs(settings)

    output_control(settings, dsets)

    #utr_length_comparison(settings, utrs)

    polyadenylation_comparison(settings, dsets)

    reads(settings)

    data_annotation_correspondence(settings)

    UTR_processing(settings)

    other_methods_comparison(settings)

    rouge_polyA_sites(settings)

if __name__ == '__main__':
    main()


# DISCUSSION #

# You have to convince the reviewers that the difference you see in
# polyadenylation usage is not just down to low coverage or the biased effects
# of sequencing.

# Regardless, the datatype is poly(A)+, which should mean that we are only
# dealing with polyadenylated RNA sequences.

# XXX it seems clear that it's in your own interest to provide just one 3UTR
# class. You keep coming back to comparing UTRs -- not comparing polyA sites.

# How reproducable are the polyA reads?
# If for the same gene, for a similar before/after coverage, what is the
# variation within the number of polyA reads?

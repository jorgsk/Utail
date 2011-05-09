"""
Script for displaying and summarizing the results from utr_finder.py.
"""

from __future__ import division
import os
import ConfigParser
import sys
import itertools

import matplotlib.pyplot as plt
plt.ion() # turn on the interactive mode so you can play with plots

from operator import attrgetter
from operator import itemgetter
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

    def __init__(self, input_line, full_info = True):
        #
        (chrm, beg, end, ID, polyA_number, strand, polyA_coordinate,
         number_supporting_reads, dstream_covrg, ustream_covrg,
         annotated_polyA_distance, nearby_PAS, PAS_distance,
         rpkm) = input_line.split('\t')

        if full_info:
            self.chrm = chrm
            self.beg = str_to_intfloat(beg)
            self.end = str_to_intfloat(end)
            self.strand = strand
            self.ID = ID

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

    def scatterplot(self, dset1, dset2, label1, label2, title):#, xlim, ylim):
        """
        A basic scatter plot
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(dset1, dset2)
        ax.set_xlabel = label1
        ax.set_ylabel = label2
        ax.set_title(title)
        #ax.set_xlim = xlim
        #ax.set_ylim = ylim
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

        xlab = ["Poly(A) cluster {0} from 3' end".format(val) for
                val in range(1, cluster_nr+1)]

        for plotnr, plotarray in enumerate([ratios, supports]):

            # adjust to get subplot-index correct
            plotnr = plotnr+1

            # mean and std
            medians = [format(np.median(ar), '.2f') for ar in plotarray]
            stds = [format(np.std(ar), '.2f') for ar in plotarray]
            means = [format(np.mean(ar), '.2f') for ar in plotarray]

            ax = fig.add_subplot(2, 1, plotnr)

            labels = []
            for (med, std, mean) in zip(medians, stds, means):
                labels.append('median: '+med+'\nmean: '+mean+'\nstd: '+std)

            ax.boxplot(plotarray)

            n = str(len(plotarray[0]))


            # Set y limits depending on if log(ratio) or read count
            if plotnr == 1:
                ax.set_title("The 3'-most poly(A) cluster has most poly(A)"\
                             "reads and highest drop in coverage", size=25)
                ax.set_ylim(-3.2, 13.2)
                ax.set_xticks([])

                ## Plot text right onto the image
                for indx, lbl in enumerate(labels):
                    ax.text(0.55+float(indx), 10, lbl, size=13)

                ax.text(0.55, 8, 'n: '+n)

            if plotnr == 2:
                ax.set_ylim(-1,60)
                ax.set_xticklabels(xlab, size=15)

                ## Plot text right onto the image
                for indx, lbl in enumerate(labels):
                    ax.text(0.55+float(indx), 50, lbl, size=13)

                ax.text(0.55, 40, 'n: '+n)

            if plotnr == 1:
                ax.set_ylabel('Log2-ratio of upstream/downstream coverage', size=15)

            if plotnr == 2:
                ax.set_ylabel('Poly(A)-read count', size=15)

        plt.show()

    def cluster_count(self, cl_count):
        """
        Input is a matrix of the number of UTRs that have x nr of clusters,
        given a minimum y number of supporting reads at that cluster. Plot each
        row -- ignoring or keeping 0-counts as you wish.

        If some UTRs have very many clusters, the plot becomes long. Provide a
        '12 clusters or more'-feature. Take max-clusters, divide by two, and add
        20% of maxclusters.

        After that... look at the bottom of the script for things to do.
        """
        # include 0 or not?
        start_pos = 1

        # Slice matrix to remove 0s if set
        cl_count = cl_count[:, start_pos:]

        max_cluster = len(cl_count[0,:])
        read_limits = len(cl_count[:,0])

        # restrict to a certain maxcluster
        up_lim = True
        if up_lim:
            lim = int(math.floor(max_cluster/float(2)))
            lim = 5
            # Sum columns after lim to the lim-colum
            cl_count[:, lim] = cl_count[:, lim:].sum(axis=1)
            # Remove columns after the lim-column
            cl_count = cl_count[:, :lim+1]

            # Update max cluster
            max_cluster = len(cl_count[0,:])

        max_height = max(cl_count[:,0])

        fig = plt.figure()

        for lim in range(read_limits):
            row_nr = lim+1
            ax = fig.add_subplot(read_limits+1, 1, row_nr)

            ax.bar(range(start_pos, max_cluster+start_pos), cl_count[lim,:],
                   align = 'center', facecolor='#777777', width=0.5)

            if row_nr == 1:
                ax.set_title('The number of poly(A) clusters per 3UTR is stable')
                ax.set_ylabel('Min 1 read', rotation='horizontal',
                              horizontalalignment = 'right')
            else:
                ax.set_ylabel('Min {0} reads'.format(row_nr), rotation='horizontal')

            ax.set_xlim((start_pos-1, max_cluster+1))
            ax.set_ylim((0, max_height + 0.2*max_height))
            ax.set_yticks(range(0, int(math.ceil(max_height+0.2*max_height)), 2000))
            ax.yaxis.grid(True)

            if row_nr == read_limits:
                ax.set_xticks(range(start_pos,max_cluster+start_pos))
                ax.set_xlabel('Number of poly(A) cluster per 3UTR')

                # If you have limited the plot, say so in the last xtick
                if up_lim:
                    xticks = range(start_pos, max_cluster+start_pos)
                    xticks[-1] = ' > {0}'.format(max_cluster)
                    ax.set_xticklabels([str(tick) for tick in xticks])

            else:
                ax.set_xticks([])

        plt.show()

    def join_clusters(self, cluster_dicts, titles):
        """
        Show side-by-side the number of clusters with a given read count from
        the dicts where the number of clusters have been obtained by different
        means (all of them, just annotated ones, etc)
        """
        #cluster_dicts = cluster_dicts[:2]
        #titles = titles[:2]

        cutoff = 20

        dset_nr = len(cluster_dicts)
        counter = np.zeros([dset_nr, cutoff]) # 0-based

        percenter = np.zeros([dset_nr, cutoff]) # 0-based

        # Get the number of clusters with read count 1, 2, etc
        for (dset_indx, cl_dict) in enumerate(cluster_dicts):
            for (read_nr, clusters) in cl_dict.iteritems():
                if read_nr > cutoff-1:
                    counter[dset_indx, cutoff-1] += len(clusters) # add if > cutoff
                else:
                    counter[dset_indx, read_nr-1] = len(clusters)

        # Get the intersection of all clusters with the rest of the sets. That
        # intersection should be divided by the total clusters after.
        def intersect_percent(read_nr):
            """
            Take the union of all datasets (annotated, support in other
            datasets, svm ...), and finally do an intersection with ALL
            clusters. Return the percentage that intersect in the end.
            """
            cl_sets = [set(cl_dict[read_nr]) for cl_dict in cluster_dicts[1:]]
            uni = set.union(*cl_sets)
            uni = set(cluster_dicts[2][read_nr])
            intrsction = set.intersection(uni, set(cluster_dicts[0][read_nr]))
            return len(intrsction)/counter[0, read_nr-1]


        # TODO
        # Make a function that takes the pairwise intersection between all sets
        # and find a way to display the information.
        # Then do the svm thing.
        # Then figure out what this means for finding novel poly(A) sites.

        # What you did today: merged all clusters into a mega cluster
        # plotted the read count of two compartments against each other
        # 

        percentages = [intersect_percent(read_nr) for read_nr in range(1, cutoff+1)]

        form_perc = [format(val, '.2f') for val in percentages]

        debug()

        # Calculate the percentage of the last rows of the first one
        # The first row is always the all_clusters
        percentages = counter[1,:]/counter[0,:]
        form_perc = [format(val, '.2f') for val in percentages]

        # Get x-axis range
        ind = range(1,cutoff+1)
        max_height = counter.max()

        # get the plot
        (fig, ax) = plt.subplots()

        colors = ['#777777', '#FF00FF']
        axes = []

        for dset_ind in range(dset_nr):
            array = counter[dset_ind,:]
            title = titles[dset_ind]
            color = colors[dset_ind]
            ll = ax.bar(ind, array, align = 'center', facecolor=color, width=0.5)
            axes.append(ll)

        # put numbers on top of the bars
        for (rect_nr, rect) in enumerate(axes[0]):
            height = rect.get_height()
            #ax.text(xpos, ypos, your_text)
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height,
                     form_perc[rect_nr], ha='center', va='bottom')

        ax.set_xlim((0, cutoff+1))
        ax.set_ylim((0, max_height + 0.2*max_height))
        ax.set_yticks(range(0, int(math.floor(max_height+0.05*max_height)), 1000))
        ax.yaxis.grid(True)

        ax.set_title('Clusters with high coverage are more often in the annotation',
                     size=20)

        fig.legend((axes[0][0], axes[1][0]), ('All', 'Found in annotation'),
                   loc=10)

        # OK this stuff is working. Now do the same but with the number from
        # poly(A) db as well as the number from other datasets etc. Then
        # organize your functions better... this is a mess!

        plt.draw()

    def all_clusters(self, cl_read_counter, my_title):
        """
        Simply plot the number of clusters with X nr of reads. Make a cut-off at
        about 15
        """
        cutoff = 20
        counter = np.zeros(cutoff) # 0-based (# of 1-size clusters are in [0])

        for (read_nr, cluster_nr) in cl_read_counter.iteritems():
            if read_nr > cutoff-1:
                counter[cutoff-1] += cluster_nr # add up those bigger than cutoff
            else:
                counter[read_nr-1] = cluster_nr

        max_height = max(counter)
        (fig, ax) = plt.subplots()

        # the arbitrary x axis range
        ind = range(1,cutoff+1)

        ax.bar(ind, counter, align = 'center', facecolor='#777777', width=0.5)

        ax.set_title('Distribution of read counts of poly(A) cluster for {0}'\
                     .format(my_title))

        ax.set_xlim((0, cutoff+1))
        ax.set_ylim((0, max_height + 0.2*max_height))
        ax.set_yticks(range(0, int(math.ceil(max_height+0.2*max_height)), 1000))
        ax.yaxis.grid(True)

        # update the last value of the xtick
        ax.set_xticks(ind)
        ind[-1] = ' > {0}'.format(cutoff)
        ax.set_xticklabels([str(tick) for tick in ind])

        plt.draw()

def get_clusters(settings):

    cluster_files = settings.polyA_files()

    # Check if all length files exist or that you have access
    [verify_access(f) for f in cluster_files.values()]

    clusters = []
    for dset_name in settings.datasets:

        clusterfile = open(cluster_files[dset_name], 'rb')
        # Get headers
        c_header = clusterfile.next()

        dset_clusters = []

        for line in clusterfile:
            ln = line.split()
            hsh = '_'.join([ln[0], ln[6], ln[5]])
            #Index should be the location of each cluster: chrm_pos_strand
            dset_clusters.append(Cluster(line))

        clusters.append(PolyaCluster(dset_clusters, dset_name))

    return clusters

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
            utr_dict[line.split()[3]].clusters.append(Cluster(line,
                                                              full_info=False))

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


def polyadenylation_comparison(settings, dsets, clusters, super_clusters, dset_2super):
    """
    * compare 3UTR polyadenylation in general
    * compare 3UTR polyadenylation UTR-to-UTR
    """
    p = Plotter()

    ## Simply: the distribution of read count of the poly(A) clusters
    ##
    # Method: 
    #  for each dset, go through all 3utrs and count the number of clusters
    # i)    with a given read count
    # ii)   CLOSE TO ANNOTATION
    # 
    # 2) 

    # Create a mega-cluster of all the clusters in clusters...
    # On the form [cluster_coord][dsetx][read_count]

    for dset in dsets:

        all_read_counter = {} # cluster_size : read_count for all clusters
        annot_read_counter = {} # cluster_size : read_count for clusters w/annot
        other_dsets = {} # cluster_size : read_count for clusters in other dsets

        for (utr_id, utr) in dset.utrs.iteritems():

            if utr.clusters == []:
                continue

            for cls in utr.clusters:

                # key that uniquely defines this polyA_cluster
                # this key will be used to do unions and intersections
                keyi = dset.name+utr.chrm+utr.strand+str(cls.polyA_coordinate)

                # All clusters
                if cls.nr_support_reads in all_read_counter:
                    all_read_counter[cls.nr_support_reads].append(keyi)
                else:
                    all_read_counter[cls.nr_support_reads] = [keyi]

                # Annotated clusters
                if cls.nr_support_reads in annot_read_counter:
                    if cls.annotated_polyA_distance != 'NA':
                        annot_read_counter[cls.nr_support_reads].append(keyi)
                else:
                    if cls.annotated_polyA_distance != 'NA':
                        annot_read_counter[cls.nr_support_reads] = [keyi]

                # Clusters in other datasets
                all_key = dset_2super[keyi]

                if cls.nr_support_reads in other_dsets:
                    for (dn, sup_reads) in zip(*super_clusters[all_key]):
                        if dn != dset.name: # don't count yourself!!
                            if sup_reads > 1: # maybe set treshold?
                                other_dsets[cls.nr_support_reads].append(keyi)
                else:
                    other_dsets[cls.nr_support_reads] = [keyi]

        clusters = (all_read_counter, annot_read_counter, other_dsets)
        titles = ('All clusters', 'Annotated_clusters', 'Clusters_in_other_dsets')

        p.join_clusters(clusters, titles)

        # Now you have which ones, in each dataset, which are found in the
        # annotation.
        # I would also like to have a counter of how many ones of size 1,2, etc,
        # are found in either of the other two datasets with read count higher
        # than 2/3/4/5 something.
        #
        # make a scatter plot of the coverage of site 1 to site 2

    ###
    ### What is the correlation of read count for the same location in different
    ### compartments? NOTE must be extended to 3 datasets or more
    ###

    ## For all the clusters with more than one dataset supporting it, get a list
    ## of how many poly(A) reads are at that support.

    #count_dset = dict((dset.name, []) for dset in dsets)
    #for (chrpos, sup_cl) in super_clusters.iteritems():

        #if (len(set(sup_cl[0])) == 2) and (len(sup_cl[0]) == 2):
            #for (name, covrg) in zip(*sup_cl):
                #if len(zip(*sup_cl)) != 2:
                    #debug()
                #count_dset[name].append(covrg)

        ##if len(sup_cl[0]) > 2:
            ### TODO a bug. a polyA cluster has double representation in both
            ### datasets. whence this error?
            ##debug()


    ## get all pairwise combinations of the dsets
    #pairs = [pa for pa in itertools.combinations([ds.name for ds in dsets], 2)]
    #for pa in pairs:
        #(p1, p2) = (pa[0], pa[1])
    #title = 'PolyA-site read count variation'
    ##xlim = 
    #p.scatterplot(count_dset[p1], count_dset[p2], p1, p2, title)
    #debug()

    ##
    ## Comparing upstream/downstream coverage and read coverage for the 3-4 last
    ## polyA clusters in a 3UTR
    ##

    #for dset in dsets:
        ## Get the maximum nr of polyA sites in one utr for max nr of plots
        ## Give a roof to the number of plots; 5?
        #max_cluster = max (utr.cluster_nr for utr in dset.utrs.itervalues())

        ## Git yourself some arrays 
        ## first second thrid and fourth
        ## nr of clusters
        #clrs_nr = 3
        #clus_list = [{'ud_ratio':[], 'support':[]} for val in range(clrs_nr)]

        #for (utr_id, utr) in dset.utrs.iteritems():

            #if utr.cluster_nr < clrs_nr:
                #continue

            #if utr.strand == '+':
                #clu = sorted(utr.clusters, key=attrgetter('polyA_coordinate'))
                #clu = clu[-clrs_nr:] # get only the last clrs_nr
                ## for + strand, reverse clu
                #clu = clu[::-1]

            #if utr.strand == '-':
                #clu = sorted(utr.clusters, key=attrgetter('polyA_coordinate'))
                #clu = clu[:clrs_nr] # get only the first clrs_nr

            ## The order of clsters in clu is '1st, 2nd, 3rd...'

            #eps_end = utr.beg + utr.eps_coord

            ## only get UTRs that have the final polyA cluster close to the
            ## coverage end!
            #if not (eps_end - 50 < clu[0].polyA_coordinate < eps_end + 50):
                #continue

            ## Normalize the ratios by the largest absolute deviation from 1
            #ud_ratios = []

            #for cls in clu:
                ## Make 0 into an almost-zero ...
                #if cls.dstream_covrg == 0:
                    #cls.dstream_covrg = 0.01

                #if cls.ustream_covrg == 0:
                    #cls.ustream_covrg = 0.01

                ## do log2 ratios (if ==1, twice as big)
                #udratio = math.log(cls.ustream_covrg/cls.dstream_covrg,2)

                #ud_ratios.append(udratio)

            ##maxratio = max(max(ud_ratios), abs(min(ud_ratios)))
            ##norm_ud_ratios = [rat/maxratio for rat in ud_ratios]
            #norm_ud_ratios = ud_ratios

            ## Append the normailzed ratios to the arrays
            #for (indx, norm_rat) in enumerate(norm_ud_ratios):
                #clus_list[indx]['ud_ratio'].append(norm_rat)

            ## Normalize the read support
            #read_supp = [cl.nr_support_reads for cl in clu]
            ##maxsupp = max(read_supp)
            ##norm_read_supp = [supp/maxsupp for supp in read_supp]
            #norm_read_supp = read_supp

            ## Append normalized support ratios to arrays
            #for (indx, norm_supp) in enumerate(norm_read_supp):
                #clus_list[indx]['support'].append(norm_supp)

        ## Do teh plots
        #p.last_three_clustersites(clus_list)
         ##RESULT you see the trend you imagined.

    #debug()

    # 1) Investiage the sensitivity to the number of mapping poly(A) reads
    # Make plots like Pedro suggested AND check out how many of the 1-reads and
    # 2-reads clusters are found in the other compartments. This gives a measure
    # of how sensitive the poly(A) clusters are.
    # For the ones that don't get a re-confirmation, and those that do, compare
    # PAS presence and SVM presence?

    # Get number of utrs that have X clusters with minimum poly(A) read coverage
    # of Y
    #for dset in dsets:

        #min_covrg = [1, 2, 3, 4, 5, 6 ,7] # minimum coverage for each cluster
        #max_cluster = max(utr.cluster_nr for utr in dset.utrs.itervalues())

        ## number of clusters in each dset. matrix, use numpy.
        #cl_counter = np.zeros([len(min_covrg), max_cluster+1], dtype=int)

        #for (utr_id, utr) in dset.utrs.iteritems():
            ## If no clusters, add 1 to all the zeros
            #if utr.clusters == []:
                #cl_counter[:,0] += 1
                #continue

            ## you are here, so there are some 
            #for minc in min_covrg:
                ## For each cluster that has more or equal support to the minc
                ## value, give a +1 to the cluster count (x) at this mic level (y)
                #cl_count = 0
                #for cls in utr.clusters:
                    #if cls.nr_support_reads >= minc:
                        #cl_count += 1

                #if cl_count > 0:
                    ## minc must be adjusted with 1 to have 0-index of array
                    #cl_counter[minc-1, cl_count] += 1
                #else:
                    #cl_counter[minc-1, 0] += 1

        #p.cluster_count(cl_counter)

    #####################################################################
    # It depends on PAS-type, PAS-distance, polyA-number (1, 2, .., last)

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

def itatr(tup):
    """
    Yield the attribute polyA_coordinate from the first element in a tuple
    """
    f = itemgetter(0)
    g = attrgetter('polyA_coordinate')

    return g(f(tup))

def merge_clusters(clusters):
    """
    Merge all clusters from all datasets into one huge cluster.
    """

    all_clusters = {} # keys = cluster_centers
    each_cluster_2all = {}

    chrms = ['chr'+str(nr) for nr in range(1,23) + ['X','Y','M']]
    tsdict = dict((chrm, {'+':[], '-':[]}) for chrm in chrms)

    # Put all the polyAclusters in a dict with cluster pos and dset.name
    for dset in clusters:
        for cls in dset.pAclusters:
            tsdict[cls.chrm][cls.strand].append((cls, dset.name))

    for (chrm, strand_dict) in tsdict.items():
        if strand_dict == {'+': [], '-': []}:
            continue
        for (strand, cls) in strand_dict.iteritems():
            if cls == []:
                # or put a placeholder?
                continue

            # initialize the first mega_cluster
            clustsum = 0
            clustcount = 0
            this_mega_cluster = []
            mega_clusters = []
            mean = 0

            for (val, dset_name) in sorted(cls, key = itatr):
                ival = val.polyA_coordinate

                # If dist between new entry and cluster mean is < 40, keep in cluster
                if abs(ival - mean) < 20:
                    this_mega_cluster.append((val, dset_name))

                    # update cluster values
                    clustsum = clustsum + ival
                    clustcount += 1
                    mean = clustsum/clustcount

                else: # If not, start a new cluster, and save the old one
                    mega_clusters.append(this_mega_cluster)
                    clustsum = ival
                    clustcount = 1
                    this_mega_cluster = [(val, dset_name)]
                    mean = ival

            # Append the last cluster
            mega_clusters.append(this_mega_cluster)

            # if the first cluster is empty, remove it
            if mega_clusters[0] == []:
                mega_clusters.pop(0)

            # Get the mean of the clusters and save to the all_clusters dict
            for mega_cluster in mega_clusters:
                coordinates = [clu.polyA_coordinate for (clu, dn) in mega_cluster]

                mean = int(math.floor(sum(coordinates)/float(len(coordinates))))
                key = chrm + strand + str(mean)

                # The information you're looking for. You can add more later if
                # needed. If you need everything, you can add the whole object.
                coverages = [clu.nr_support_reads for (clu, dn) in mega_cluster]
                dsets = [dn for (clu, dn) in mega_cluster]

                # add entry to the all clusters
                all_clusters[key] = ((dsets, coverages))

                # store the key with the original polyA objects
                # for now, make a each_cluster -> all_cluster mapping
                # the idea is that you get the key for all_clusters by looking
                # in each_keys for your dataset-polyAcluster combination
                each_keys = [dset_name+cl.chrm+cl.strand+str(cl.polyA_coordinate) for
                             (cl, dset_name) in mega_cluster]

                for e_key in each_keys:
                    each_cluster_2all[e_key] = key

    return (all_clusters, each_cluster_2all)

def main():
    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    # Directory paths for figures and where the output lies
    (savedir, outputdir) = [os.path.join(here, d) for d in ('figures', 'output')]

    # Read UTR_SETTINGS
    settings = Settings(os.path.join(here, 'UTR_SETTINGS'), savedir, outputdir, here)

    # Get the SVM coordinates
    #svm = get_svm()

    # Get the dsets with utrs and their clusters from the length and polyA files
    dsets = get_utrs(settings)

    # Get just the clusters from the polyA files
    clusters = get_clusters(settings)

    # Merge all clusters in datset. Return interlocutor-dict for going from each
    # poly(A) cluster in each dataset to the
    (super_cluster, dset_2super) = merge_clusters(clusters)

    #output_control(settings, dsets)

    #utr_length_comparison(settings, utrs)

    polyadenylation_comparison(settings, dsets, clusters, super_cluster, dset_2super)
    debug()

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

# IDEA: measure of coverage-limitation: how many polyA clusters (>5 reads) are
# found with only 1 read in other compartments/datasets?

# Q: is there variation between datasets in polyA-cluster usage in utrs that are
# expressed in all datsets?

# Q: what should we call a reliable cluster?

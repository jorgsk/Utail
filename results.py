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

        # SVM info might be added
        self.svm_coordinates = []

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

        self.settings_file = settings_file

        # A bit ad-hoc: dicrectory of polyA read files
        self.polyAread_dir = os.path.join(here, 'polyA_files')

    # Return the paths of the length files
    def length_files(self):
        return dict((d, os.path.join(self.here, self.outputdir, 'length_'+d))
                    for d in self.datasets)

    # Return the paths of the polyA files
    def polyA_files(self):
        return dict((d, os.path.join(self.here, self.outputdir, 'polyA_' + d))
                    for d in self.datasets)

    # Return the paths of the epsilon files
    def epsilon_files(self):
        return dict((d, os.path.join(self.here, self.outputdir,'cumul_'+d+'.stat'))
                    for d in self.datasets)

    # Return the path of the svm-file. This file should be a bed-file where name
    # is the utr and the coordinate is the svm-coordinate. If the file does not
    # exist, create it by bedIntersecting the 3UTR file from
    # annotation_parser.py
    def svm_file(self, chr1=False):

        # Make the svm dir if it doesn't exist
        svmdir = os.path.join(self.here, 'svm')
        if not os.path.exists(svmdir):
            os.makedirs(svmdir)
        if chr1:
            svm_f = os.path.join(svmdir, 'svm_utr_intesection_chr1.bed')
        else:
            svm_f = os.path.join(svmdir, 'svm_utr_intesection.bed')

        # Get the settings file of the utr_analyzer; use this file to get the
        # correct path of the 3utr-bedfile.
        import utr_finder as finder

        finder_settings = finder.Settings\
                (*finder.read_settings(self.settings_file))

        if chr1:
            finder_settings.chr1 = True

        # Check if 3UTRfile has been made or provided; if not, get it from annotation
        beddir = os.path.join(self.here, 'source_bedfiles')
        utr_bed = finder.get_utr_path(finder_settings, beddir)

        # Check that it's accessible
        verify_access(utr_bed)

        # Run intersect bed on the utrfile with the svmfile
        svm_bed = os.path.join(svmdir, 'svm_final.bed')
        # If this file doesn't exist, you have to make it

        #if not os.path.isfile(svm_bed):
            #svm_from_source(svm_bed, svmdir, utr_bed, finder_settings, chr1)

        from subprocess import Popen, PIPE

        # Consider the strand; you have discarded the polyA reads with the
        # 'wrong' strand
        cmd = ['intersectBed', '-wb', '-s', '-a', svm_bed, '-b', utr_bed]
        f = Popen(cmd, stdout=PIPE)

        svm_dict = {}
        for line in f.stdout:
            (chrm, beg, end, d,d, strand, d,d,d, utr_exon_id, d,d) = line.split()
            # The bedfile has utr_exons in the index. You need the whole utr.
            # OBS 
            utr_id = '_'.join(utr_exon_id.split('_')[:-1])
            beg = int(beg)
            end = int(end)

            if utr_id in svm_dict:
                svm_dict[utr_id].append((chrm, beg, end, '0', '0', strand))
            else:
                svm_dict[utr_id] = [(chrm, beg, end, '0', '0', strand)]

        return svm_dict

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

    def join_clusters(self, cluster_dicts, titles, in_terms_of):
        """
        Show side-by-side the number of clusters with a given read count from
        the dicts where the number of clusters have been obtained by different
        means (all of them, just annotated ones, etc)
        """
        #cluster_dicts = cluster_dicts[:2]
        #titles = titles[:2]

        cutoff = 8

        dset_nr = len(cluster_dicts)
        counter = np.zeros([dset_nr, cutoff]) # 0-based

        # Get the number of clusters with read count 1, 2, etc
        for (dset_indx, cl_dict) in enumerate(cluster_dicts):
            for (read_nr, clusters) in cl_dict.iteritems():
                if read_nr > cutoff-1:
                    counter[dset_indx, cutoff-1] += len(clusters) # add if > cutoff
                else:
                    counter[dset_indx, read_nr-1] = len(clusters)

        dset_dict = dict(zip(titles, cluster_dicts))

        # pairw_all_annot contains is two tuples, each with two elements: a
        # matrix and a row-identifer for the matrix
        pairw_all_annot = pairwise_intersect(in_terms_of, dset_dict, cutoff)


        # Print two figures; one where all_clusters is main and one where
        # annotated_clusters are main
        for (dset_ind, (pairw_matrix, pairw_mrows)) in enumerate(pairw_all_annot):

            cols = ['#0000FF','#3333FF','#4C3380','#8A5CE6','#AD85FF','#AD39FF']

            # 1) Plot the first bars: the all_clusters ones.

            # Total number of bars in each complex
            bar_nr = len(pairw_matrix[:,0,0]) # actually =  + all_cl and - union
            # Set width of bars
            bar_width = 0.6 # total width of bar-cluster = wid*bar_nr
            # Get the width of the whole bar-compled
            complex_width = bar_width*bar_nr
            # Set how much space should be between the bar-complexes
            complex_interspace = complex_width/2
            # Total number of complexes is cutoff. Get the last x-coordinate.
            final_x = math.ceil((complex_width + complex_interspace)*cutoff)

            # Set your x-axis so that it will be wide enough for all complexes
            # this will be the leftmost x-position of the first bar
            ind = np.arange(1, final_x+1, complex_width+complex_interspace)
            # Shorten to make sure that this is as long as the data-points
            ind = ind[:cutoff]

            # Get max height of bars
            max_height = counter[dset_ind].max() # original clusters always highest

            # get the plot
            (fig, ax) = plt.subplots()

            # Plot the cluster counts (keep axis objects for later)
            ax_list = [ax.bar(ind, counter[dset_ind], facecolor=cols[0],
                          width=bar_width)]
            # Plot the union-counts on top of the cluster counts
            ax_list.append(ax.bar(ind, pairw_matrix[0, :, 0], facecolor=cols[1],
                                  width=bar_width))

            # Plot the rest of the bars.
            for int_ind in range(1, bar_nr):

                array = pairw_matrix[int_ind,:,0] # absolute numbers has dim 0
                clr = cols[int_ind+2]
                # ind+bar_width*(int_ind+1) adjusts the bars one 'bar_width' on the
                # x-axis
                ll = ax.bar(ind+bar_width*(int_ind), array, facecolor=clr,
                            width=bar_width)
                ax_list.append(ll)

            # format the union percentages nicely
            form_perc = [[format(el*100, '.0f')+'%' for el in pairw_matrix[ind,:,1]]
                         for ind in range(bar_nr)]

            # 4) put numbers on top of the bars (and skip the 'union' axis)
            myaxes = [ax_list[0]] + ax_list[2:]
            for (bars_nr, bars_axes) in enumerate(myaxes):
                for (rect_nr, rect) in enumerate(bars_axes):
                    height = rect.get_height()
                    #ax.text(xpos, ypos, your_text) and ax is the CURRENT axes
                    ax.text(rect.get_x()+rect.get_width()/2., 1.03*height,
                             form_perc[bars_nr][rect_nr], ha='center',
                            va='bottom', size='small')

            # Set the x-axis and y-axis
            ax.set_xlim((0, final_x+1))
            ax.set_ylim((0, max_height + 0.2*max_height))
            # Set the labels on the y axis
            ax.set_yticks(range(0, int(math.floor(max_height+0.05*max_height)), 1000))

            # set the positions of the ticks on the x axis; they should be the
            # centered on the bar clusters
            comp_center = complex_width/2
            complex_centers = np.arange(1+comp_center, final_x+comp_center,
                                         complex_width+complex_interspace)
            ax.set_xticks(complex_centers)

            # Set labels for those ticks (last one is > or more)
            ticklabels = [str(val) for val in range(1, cutoff+1)]
            ticklabels[-1] = ' > {0}'.format(cutoff)

            ax.set_xticklabels(ticklabels)

            # Put grids on the y-axis
            ax.yaxis.grid(True)

            # Set labels on the axes
            ax.set_xlabel('Number of reads covering poly(A) cluster')
            ax.set_ylabel('Number of poly(A) clusters')

            ax.set_title('Clusters with high coverage are often found in'
                         'supporting data ', size=20)

            legend_titles = [pairw_mrows[0][0]] + [tup[1] for tup in pairw_mrows]
            legend_axes = (ax[0] for ax in ax_list)
            fig.legend(legend_axes, legend_titles, loc=10)

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

def pairwise_intersect(in_terms_of, dset_dict, cutoff):
    """
    For 'all clusters' and 'annotated clusters', get the intersection with all
    other datasets; and the intersection with the union of all other datasets.
    """
    # Get the pairs: all with the rest
    # Get the union: all except the rest
    (all_cl, annot) = in_terms_of # all clusters and annotated TTS sites

    # Get the pairs of 'All_clusters' with the other datasets
    all_pairs = [(all_cl, dset_name) for dset_name in
                  dset_dict.keys() if dset_name != all_cl]

    # Get the names of all the datasets except 'All clusters'
    not_all = [dset_name for dset_name in dset_dict.keys() if dset_name !=
               all_cl]

    # Get the pairs of 'Annotated clusters' with the other datasets, except all
    annot_pairs = [(annot, dset_name) for dset_name in
                  dset_dict.keys() if (dset_name != annot) and (dset_name != all_cl)]

    # Get the names of all the datasets except 'Annotatedusters'
    not_annot = [dset_name for dset_name in dset_dict.keys() if (dset_name !=
               annot) and (dset_name != all_cl)]

    matr = []

    # Get the intersection matrices themselves. Rows are (main_dset, sub_dset)
    # intersection pars. Columns are the read count of the polyA sites.

    for (pairs, unions) in [(all_pairs, not_all), (annot_pairs,  not_annot)]:
        isec_matrix = get_intersection_matrix(pairs, unions, cutoff, dset_dict)

        # add "union"
        pairs = [(pairs[0][0], 'Union')] + pairs
        matr.append((isec_matrix, pairs))

    return matr

def get_intersection_matrix(pair_names, unions_names, cutoff, dset_dict):
    """
    Get the intersection of all the 'pairs' datasets.
    Create the union of all the 'union' datasets, and do the intersection of
    this unition with the leading pair in the pair datasets.

    Do all this for every read count up to cutoff.
    """

    dset_nr = len(pair_names)+1 #pairs and union

    # Counter is 3-dimensional for keeping both abs number of intersection AND
    # percentages. 

    counter = np.zeros([dset_nr, cutoff, 2]) # 0-based

    # Get the pairs 
    for (indx1, (main_name, sub_name)) in enumerate(pair_names):
        # Get the pair-dsets
        main_dset = dset_dict[main_name]
        sub_dset = dset_dict[sub_name]

        # Iterate through all (polyA-cluster, read_count) points in the
        # datasets, and add the polyA-clusters to two temporary lists, indexed
        # by the read count from 0 to cutoff-1.
        main_cls = [[] for val in range(cutoff)]
        sub_cls = [[] for val in range(cutoff)]

        for (dset, dset_l) in [(main_dset, main_cls), (sub_dset, sub_cls)]:

            for (read_nr, clusters) in dset.iteritems():
                if read_nr <= 0:
                    debug()
                if read_nr > cutoff-1:
                    dset_l[cutoff-1].append(clusters) # add if > cutoff
                else:
                    dset_l[read_nr-1] = clusters

                if dset_l[-1] != []:
                    # THE BUG a cluster [X, XX] are noted with 0 covering reads.
                    # This is a mistake.
                    debug()

        # Flatten the last arrays
        main_cls[-1] = sum(main_cls[-1], [])
        sub_cls[-1] = sum(sub_cls[-1], [])

        # Get number of intersections 
        isect_nrs = [len(set.intersection(set(main_cls[count]), set(sub_cls[count])))
                                      for count in range(0, cutoff)]

        # Get percent of intersection relative to 'main' dataset (will be all or annot)
        isect_pcnt = []
        for (indx, isect_nr) in enumerate(isect_nrs):

            # Only calculate percentage if more than 1 cluster with this read count
            if main_cls[indx] != 0:
                isect_pcnt.append(isect_nrs[indx]/len(main_cls[indx]))
            else:
                isect_pcnt.append(0)

        # Add the number and intersection to the array
        counter[indx1,:,0] = isect_nrs
        counter[indx1,:,1] = isect_pcnt

    # Now all the pairs have been added. Add the unions
    # Take the union of all dsets except
    all_cls = [[] for val in range(cutoff)]

    # add all the clusters from the union datasets to all_cls
    for u_name in unions_names:
        for (read_nr, clusters) in dset_dict[u_name].iteritems():

            if read_nr > cutoff-1:
                all_cls[cutoff-1].append(clusters) # add if > cutoff
            else:
                all_cls[read_nr-1].append(clusters)

    # flatten all_cls (which has all the clusters in the union dsets)
    # and take union at the same tim
    all_cls = [sum(el, []) for el in all_cls]

    # Get number of intersections 
    # (using main_cls from the previous for-loop -- dirty :S)
    all_isect_nrs = [len(set.intersection(set(main_cls[count]), set(all_cls[count])))
                                  for count in range(0, cutoff)]

    # Get percent of intersection relative to 'main' dataset (will be all or annot)
    all_isect_pcnt = []
    for (indx, isect_nr) in enumerate(isect_nrs):

        # Only calculate percentage if more than 1 cluster with this read count
        if main_cls[indx] != 0:
            all_isect_pcnt.append(all_isect_nrs[indx]/len(main_cls[indx]))
        else:
            all_isect_pcnt.append(0)

    # Add the number and intersection to the array
    counter[-1,:,0] = all_isect_nrs
    counter[-1,:,1] = all_isect_pcnt

    ### flip things around; put union row first. This is for better compliance
    # with downstream code

    newcount = np.zeros([dset_nr, cutoff, 2])
    newcount[0] = counter[-1]
    newcount[1:] = counter[0:-1]

    return newcount

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

def get_utrs(settings, svm=False):
    """
    Return a list of UTR instances. Each UTR instance is
    instansiated with a list of UTR objects and the name of the datset.

    1) Read through the length file, create UTR objects for each line
    2) If present, read through the polyA file, updating the UTR object created
    in step 1

    """
    # Do everything for the length files
    length_files = settings.length_files()
    cluster_files = settings.polyA_files()

    if svm:
        # Get a dictionary indexed by utr_id with the svm locations
        svm_dict = settings.svm_file()

    # Check if all length files exist or that you have access
    [verify_access(f) for f in length_files.values()]

    all_utrs = {}
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

        # Update the UTRs
        for (utr_id, utr_obj) in utr_dict.iteritems():
            # Update with poly(A) cluster info
            utr_dict[utr_id].cluster_nr = len(utr_obj.clusters)

            # If svm is not set, continue
            if not svm:
                continue

            # Don't add svm info if there is no svm
            if utr_id not in svm_dict:
                continue

            # Get SVM info
            SVMs = svm_dict[utr_id]

            for (chrm, svm_beg, svm_end, d, d, strand) in SVMs:
                # Check that strand and chrm of SVM are the same as UTR object ...
                assert (utr_dict[utr_id].chrm, utr_dict[utr_id].strand)\
                        == (chrm, strand), 'Mismatch'

                #Update UTR object with SVM info
                utr_dict[utr_id].svm_coordinates.append(svm_beg)

        all_utrs[dset_name] = UTRDataset(utr_dict, dset_name)

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


### PLOTTING METHODS ###

def cluster_size_sensitivity():
    """
    1) Investiage the sensitivity to the number of mapping poly(A) reads
    Make plots like Pedro suggested AND check out how many of the 1-reads and
    2-reads clusters are found in the other compartments. This gives a measure
    of how sensitive the poly(A) clusters are.
    """

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


def udstream_coverage_last_clusters(dsets):
    """
    Comparing upstream/downstream coverage and read coverage for the 3-4 last
    polyA clusters in a 3UTR
    """

    p = Plotter()

    for dset in dsets.values():

        # first second thrid and fourth
        # nr of clusters
        clrs_nr = 3
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

            for cls in clu:
                # Make 0 into an almost-zero ...
                if cls.dstream_covrg == 0:
                    cls.dstream_covrg = 0.01

                if cls.ustream_covrg == 0:
                    cls.ustream_covrg = 0.01

                # do log2 ratios (if ==1, twice as big)
                udratio = math.log(cls.ustream_covrg/cls.dstream_covrg,2)

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
         #RESULT you see the trend you imagined.

def correlate_polyA_coverage_counts(dsets, super_clusters):
    """
    What is the correlation of read count for the same location in different
    compartments? NOTE must be extended to 3 datasets or more
    """
    p = Plotter()

    # For all the clusters with more than one dataset supporting it, get a list
    # of how many poly(A) reads are at that support.
    count_dset = dict((dset_name, []) for dset_name in dsets)
    for (chrpos, sup_cl) in super_clusters.iteritems():

        if (len(set(sup_cl[0])) == 2) and (len(sup_cl[0]) == 2):
            for (name, covrg) in zip(*sup_cl):
                if len(zip(*sup_cl)) != 2:
                    debug()
                count_dset[name].append(covrg)

        #if len(sup_cl[0]) > 2:
            ## TODO a bug. a polyA cluster has double representation in both
            ## datasets. whence this error?
            #debug()

    # get all pairwise combinations of the dsets
    pairs = [pa for pa in itertools.combinations([ds_name for ds_name in dsets], 2)]
    for pa in pairs:
        (p1, p2) = (pa[0], pa[1])
    title = 'PolyA-site read count variation'
    p.scatterplot(count_dset[p1], count_dset[p2], p1, p2, title)

def compare_cluster_evidence(dsets, clusters, super_clusters, dset_2super):
    """
    Find the co-occurence for all the evidence for polyA clusters you have (in
    annotation, svm support, etc).

    Especially, what support do the annotated, alone, have?
    """

    p = Plotter()

    #for dset in dsets:
    for dset in dsets.values()[1:]:

        all_read_counter = {} # cluster_size : read_count for all clusters
        svm_support = {} # for all clusters: do they have SVM support?
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

                # SVM support
                if utr.svm_coordinates != []:
                    # Look for svm support for this cluster
                    found = False
                    for crd in utr.svm_coordinates:
                        if crd-30 < cls.polyA_coordinate < crd+30:
                            found = True
                            break

                    # If you found it, add to the svm_support cluster
                    if found:
                        if cls.nr_support_reads in svm_support:
                            svm_support[cls.nr_support_reads].append(keyi)
                        else:
                            svm_support[cls.nr_support_reads] = [keyi]

                # Annotated clusters
                if cls.annotated_polyA_distance != 'NA':
                    if cls.nr_support_reads in annot_read_counter:
                        annot_read_counter[cls.nr_support_reads].append(keyi)
                    else:
                        annot_read_counter[cls.nr_support_reads] = [keyi]

                # Clusters in other datasets
                if cls.nr_support_reads in other_dsets:
                    all_key = dset_2super[keyi] # the in-between key
                    for (dn, sup_reads) in zip(*super_clusters[all_key]):
                        if dn != dset.name: # don't count yourself!!
                            if sup_reads > 1: # maybe set treshold?
                                other_dsets[cls.nr_support_reads].append(keyi)
                else:
                    other_dsets[cls.nr_support_reads] = [keyi]

        cluster_dicts = (all_read_counter, svm_support, annot_read_counter)
        titles = ('All clusters', 'SVM_support', 'Annotated_TTS')
        #clusters = (all_read_counter, svm_support, annot_read_counter, other_dsets)
        #titles = ('All clusters', 'SVM_support', 'Annotated_TTS',
                  #'In_other_compartments')

        # make two figurses: one in terms of all clusters and on of annotated
        # TTS sites
        in_terms_of = (titles[0], titles[2])
        p.join_clusters(cluster_dicts, titles, in_terms_of)


def output_control(settings, dsets):
    """
    Control: what is the correlation between # supporting reads; change in
    coverage; PAS-type; PAS-distance; and rpkm?
    """
    p = Plotter()

    # TRIANGLE PLOT ##

    for dset in dsets.itervalues():

        nr_supp_reads = []
        covrg_down = []
        covrg_up = []
        rpkm = []

        for (utr_id, utr) in dset.utrs.iteritems():

            # For all polyA clusters get
            # 1) NR of supporting reads
            # 2) Coverage downstream
            # 3) Coverage upstream
            # 4) RPKM

            for cls in utr.clusters:
                nr_supp_reads.append(cls.nr_support_reads)
                covrg_down.append(cls.dstream_covrg)
                covrg_up.append(cls.ustream_covrg)
                rpkm.append(utr.RPKM)

        # Titles for the above variables
        titles = ['Nr Supporting Reads', 'Downstream Coverage',
                  'Upstream Coverage', 'RPKM']

        array = [nr_supp_reads, covrg_down, covrg_up, rpkm]
        p.triangleplot_scatter(array, titles)

     #For the ones with 1, 2, 3, 4, ... polyA sites:
        #1) box-plots of the ratio of downstream/upstream 
        #2) box-plots of the # of covering reads

     #Upstream/downstream ratios of polyA sites.


def polyadenylation_comparison(settings, dsets, clusters, super_clusters, dset_2super):
    """
    * compare 3UTR polyadenylation in general
    * compare 3UTR polyadenylation UTR-to-UTR
    """

    # Compare the compartments
    compare_cluster_evidence(dsets, clusters, super_clusters, dset_2super)

    debug()

    ## Correlate the coverage counts of common polyadenylation sites between
    ## clusters
    #correlate_polyA_coverage_counts(dsets, super_clusters)

    ## Get the change in coverage ratio of the 3 (or so) last polyadenylation
    ## sites from the 3 UTR
    #udstream_coverage_last_clusters(dsets)

    ## See the effects of excluding pA sites with 1, 2, 3, etc, coverage
    #cluster_size_sensitivity()

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

def get_reads_from_file(ds, dset, finder, pAread_file, utr_bed, utr_exons):
    """
    Go through
    """

    # Prepare a list of 'this_strand' clusters. They will serve as controls.
    this_strands = []

    # Get a dictionary for each utr_id with its overlapping polyA reads
    utr_polyAs = finder.get_polyA_utr(pAread_file, utr_bed)

    # Cluster the poly(A) reads for each utr_id.
    polyA_cls = finder.cluster_polyAs(utr_polyAs, utr_exons, polyA=True)

    # Add the read_count information to the UTR objects
    for (utr_exon_id, utr_exon_pAreads) in polyA_cls.iteritems():
        utr_id = '_'.join(utr_exon_id.split('_')[:-1])

        # Append the 'this_strand' clusters to the list
        if utr_exon_pAreads['this_strand'][1] != []:
            this_strands.append(utr_exon_pAreads['this_strand'][1])

        # Append the 'other_strand' clusters to the polyA clusters of the
        # UTR objects.
        if utr_exon_pAreads['other_strand'][0] != []:
            other_strand = utr_exon_pAreads['other_strand']
            # Go though all the clusters of this utr
            for cluster in dset.utrs[utr_id].clusters:
                # Go though all the polyA reads 
                for (cl_nr, cl_mean) in enumerate(other_strand[0]):
                    # If the match, update the cluster with the coverage
                    if cl_mean == cluster.polyA_coordinate:
                        cluster.all_pA_reads = other_strand[1][cl_nr]

    return this_strands

def classic_polyA_stats(settings, dsets):
    """
    * distances from polyA cluster to PAS site
    * PAS variant distribution
    * variance within polyA pAclusters
    * degree of usage of early and late polyA sites
    """

    # 1) Try to call rna_analyser methods to re-intersecting and
    # re-clustering the polyA reads

    import utr_finder as finder

    finder_settings = finder.Settings\
            (*finder.read_settings(settings.settings_file))

    chr1 = True
    if chr1:
        finder_settings.chr1 = True

    # Get the bedfile with 3UTR exons
    beddir = os.path.join(settings.here, 'source_bedfiles')
    utr_bed = finder.get_utr_path(finder_settings, beddir)
    polyA_dir = settings.polyAread_dir

    # Get the annotation class from utr_finder.py
    finder_annot = finder.Annotation(finder_settings.annotation_path)

    # Set the utrfile_path on finder's annotation class
    finder_annot.utrfile_path = utr_bed

    # Getting a dictionary of the 3UTRs from finder_annotation
    utr_exons = finder_annot.get_utrdict()

    # Do the analysis for all the polyA reads in the directory
    for (file_ds, pAread_file) in get_pA_files(polyA_dir):

        # Get the 'this_strands' clusters; and also update the dset class UTRs
        # with 'other_strand' info

        # only proceed if the file_dataset has already been imported in dsets
        if file_ds in dsets:
            dset = dsets[file_ds] # get the dset class instance
        else:
            continue

        this_strands = get_reads_from_file(file_ds, dset, finder, pAread_file,
                                           utr_bed, utr_exons)

        debug()
        # Do the actual analysis on the 'other_strand'

    # Basically: give me the UTR dict
    # TODO:

        # 2) Calculate the distribution of polyA variation about a site
        # 3) Calculate the distance from polyA-site to 1) closest 2) best PAS
        # 4) Intersect with the UTR objects to get the PAS type
        #   -- Get the PAS distribution
        # 4) Get the variation about a SINGLE polyA-site. It tells us something
        # about the abundance of different types of transcript? Or just some
        # artifact of the PCR process? It possibly points to the minimum number
        # of unique transcripts in the sample. 

        # 5) The same stuff but for 'other_strand'.

    pass

def get_pA_files(polyA_dir):
    """
    Return a (datset, read_file) tuple for the files in polyA_dir
    """
    return [(f.lstrip('polyA_reads_').rstrip('.bed'), os.path.join(polyA_dir, f))
                      for f in os.listdir(polyA_dir)]

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

def svm_from_source(svm_path, svm_dir, utr_path, hgfasta_path, chr1):
    print('SVM file not found. Generating from SVM program ...')
    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    sys.path.append(os.path.join(here, 'modules'))
    sys.path.append(os.path.join(here, 'modules/pyfasta'))

    if chr1:
        fasta_file = os.path.join(svm_dir, 'utr_exon_fastas_chr1.fa' )
    else:
        fasta_file = os.path.join(svm_dir, 'utr_exon_fastas.fa' )

    fasta_handle = open(fasta_file, 'wb')

    from fasta import Fasta

    # 1) For all the utr-exons in utr_path, get the genomic sequence, and save
    # in a fasta file

    # The fasta object
    f = Fasta(hgfasta_path)
    for line in open(utr_path, 'rb'):

        (chrm, beg, end, utr_exon_id, ex_nr, strand) = line.split()
        beg = int(beg)
        end = int(end)

        if strand == '+':
            seq = f.sequence({'chr':chrm,'start':beg+1, 'stop':end-1,
                                          'strand':strand}).upper()
        if strand == '-':
            seq = f.sequence({'chr':chrm,'start':beg+2, 'stop':end,
                                          'strand':strand}).upper()

        fasta_handle.write('>{0}\n{1}\n'.format(utr_exon_id, seq))

    fasta_handle.close()

    # 2) Run the SVM machine on the fasta file

    # 3) Parse the output into a nicely formatted bedfile in svm_patth

def main():
    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    # Directory paths for figures and where the output lies
    (savedir, outputdir) = [os.path.join(here, d) for d in ('figures', 'output')]

    # Read UTR_SETTINGS
    settings = Settings(os.path.join(here, 'UTR_SETTINGS'), savedir, outputdir, here)

    # Get the dsets with utrs and their clusters from the length and polyA files
    # Optionally get SVM information as well
    dsets = get_utrs(settings, svm=True)

    # Get just the clusters from the polyA files
    #clusters = get_clusters(settings)

    # Merge all clusters in datset. Return interlocutor-dict for going from each
    # poly(A) cluster in each dataset to the
    #(super_cluster, dset_2super) = merge_clusters(clusters)

    #output_control(settings, dsets)

    #utr_length_comparison(settings, utrs)

    #polyadenylation_comparison(settings, dsets, clusters, super_cluster, dset_2super)

    ## The classic polyA variation distributions
    classic_polyA_stats(settings, dsets)

    debug()

    #reads(settings)

    #data_annotation_correspondence(settings)

    #UTR_processing(settings)

    #other_methods_comparison(settings)

    #rouge_polyA_sites(settings)


if __name__ == '__main__':
    main()

# DISCUSSION #

# You have to convince the reviewers that the difference you see in
# polyadenylation usage is not just down to low coverage or the biased effects
# of sequencing.

# Regardless, the datatype is poly(A)+, which should mean that we are only
# dealing with polyadenylated RNA sequences.

# How reproducable are the polyA reads?
# If for the same gene, for a similar before/after coverage, what is the
# variation within the number of polyA reads?

# IDEA: measure of coverage-limitation: how many polyA clusters (>5 reads) are
# found with only 1 read in other compartments/datasets?

# Q: is there variation between datasets in polyA-cluster usage in utrs that are
# expressed in all datsets?

# Q: what should we call a reliable cluster?

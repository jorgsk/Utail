"""
Script for displaying and summarizing the results from utr_finder.py.
"""

from __future__ import division
import os
import ConfigParser
import sys
from itertools import combinations as combins

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
    Class that contains helper functions to work on datasets. A dataset is
    typically an individual cell compartment.
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
            for (utr_id, utr) in self.utrs.iteritems():
                if (utr.eps_rel_size != 'NA') and (utr.RPKM > minRPKM):
                    if utr.ID in IDs:
                        eps_lengths.append(utr.eps_rel_size)

        # Include all utr objects
        else:
            for (utr_id, utr) in self.utrs.iteritems():
                if utr.eps_rel_size != 'NA' and utr.RPKM > minRPKM:
                    eps_lengths.append(utr.eps_rel_size)

        return eps_lengths

    def expressed_IDs(self):
        """
        Return list (or set) of IDs of utrs that are expressed in this dataset. The
        criterion for 'expressed' is that eps_rel_size is not 'NA'.
        """
        return [utr.ID for (utr_id, utr) in self.utrs.iteritems()
                if utr.eps_rel_size != 'NA']


class UTR(object):
    """
    For UTR objects from the 'length' output file in the 'output' directory.
    """

    def __init__(self, input_line):

        # Read all the parameters from line
        (chrm, beg, end, utr_extended_by, strand, ID, epsilon_coord,
         epsilon_rel_size, epsilon_downstream_covrg, epsilon_upstream_covrg,
         annotTTS_dist, annot_downstream_covrg, annot_upstream_covrg,
         epsilon_PAS_type, epsilon_PAS_distance, epsilon_beyond, RPKM,
         avrg_covrg ) = input_line.split('\t')

        self.chrm = chrm
        self.beg = str_to_intfloat(beg)
        self.end = str_to_intfloat(end)
        self.length = self.end-self.beg
        self.extended_by = str_to_intfloat(utr_extended_by)
        self.strand = strand
        self.ID = ID
        self.eps_coord = str_to_intfloat(epsilon_coord)
        self.eps_rel_size = str_to_intfloat(epsilon_rel_size)
        self.epsilon_beyond_aTTS = str_to_intfloat(epsilon_beyond)
        # Get the _absoulte_ epsilong length and the distance from annotanted
        # end
        if self.eps_rel_size != 'NA':
            self.eps_abs_size = math.ceil(self.eps_rel_size*self.length)
            self.eps_remainder = self.length - self.eps_abs_size
        else:
            self.eps_abs_size = 'NA'
            self.eps_remainder = 'NA'
        #
        self.eps_downstream_covrg = str_to_intfloat(epsilon_downstream_covrg)
        self.eps_upstream_covrg = str_to_intfloat(epsilon_upstream_covrg)
        self.annot_downstream_covrg = str_to_intfloat(annot_downstream_covrg)
        self.annot_upstream_covrg = str_to_intfloat(annot_upstream_covrg)
        self.annotTTS_dist = annotTTS_dist

        # PAS type and PAS distance are space-delimited
        PAS_type = epsilon_PAS_type.split(' ')
        self.eps_PAS_type = [str_to_intfloat(pas) for pas in PAS_type]

        PAS_distance = epsilon_PAS_distance.split(' ')
        self.eps_PAS_distance = [str_to_intfloat(dist) for dist in PAS_distance]

        self.RPKM = str_to_intfloat(RPKM)
        self.avrg_covrg = str_to_intfloat(avrg_covrg.rstrip())

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

        # Later, you might add a list of genomic coordinates of the polyA-reads
        # for this cluster
        self.all_pA_coords = 'NA'

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
        fig.draw()

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
        fig.draw()

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

        fig.draw()

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

        plt.draw()

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

        plt.draw()

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

            ax.set_title('Clusters with high coverage are often found in'\
                         ' supporting data ', size=20)

            legend_titles = [pairw_mrows[0][0]] + [tup[1] for tup in pairw_mrows]
            legend_axes = (ax[0] for ax in ax_list)
            fig.legend(legend_axes, legend_titles, loc=10)

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

    def distance_histogram(self, distance_dict):
        """
        Plot the normalized histogram of the three types of distances
        superimposed and look for distances
        """

        cols = ['#0000FF', '#3333FF', '#4C3380']
        col_dict = dict(zip(distance_dict.keys(), cols))

        for (title, distances) in distance_dict.items():
            (fig, ax) = plt.subplots()
            # TODO I think the normalization is a bit off. The 0-bar stretches
            # all the way to 1.
            ax.hist(distances, bins=200, normed=True, label = title, color =
                    col_dict[title])

            # Put a fitted normal distribution on top
            mean = np.mean(distances)
            std = np.std(distances)

            ax.set_ylim((0,1))
            ax.set_xlim((-30, 30))
            ax.legend()

        plt.draw()

        # NOTE the plots work well. You need to find a way to combine them into
        # one plot (which is first etc?)

    def cluster_size_distribution(self, cluster_sizes, cutoff):
        """
        Make bar plots of the (normalized) distribution of cluster sizes.
        """
        cols = ['#4C3380', '#0000FF', '#3333FF']
        col_dict = dict(zip(cluster_sizes.keys(), cols))

        (fig, axes) = plt.subplots(len(cluster_sizes), sharex=True, sharey=True)

        x_coords = range(1, cutoff+1)

        for indx, (name, size_dist) in enumerate(cluster_sizes.items()):
            # Normalize size dist
            size_sum = size_dist.sum()
            ax = axes[indx]
            size_dist = [val/size_sum for val in size_dist]
            ax.bar(x_coords, size_dist, color = col_dict[name],
                           align='center', label = name)

            ax.set_ylim((0, 1.1))
            ax.legend()
            ax.set_xlabel('PolyA cluster sizes', size=20)
            ax.set_ylabel('Frequency of occurrence', size=20)

            ax.text(5, 0.8, 'Total clusters: {0}'.format(size_sum))

        # Hide xticks for all but the last plot
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        axes[-1].set_xticks(x_coords)
        x_coords[-1] = ' > {0}'.format(cutoff)
        axes[-1].set_xticklabels([str(tick) for tick in x_coords])
        axes[0].set_title('PolyA clusters in opposite strand are fewer, have'\
                          ' mostly 1-size clusters, and have few large clusters',
                         size=23)
        # Make subplots close
        fig.subplots_adjust(hspace=0)

        plt.draw()

        # OK but I think it's best to print them on top of each other to really
        # show how many more there are of the cis-strand one.


    def rec_sensitivity(self, sensitivities, intervals, attributes):
        """
        Plot how the false negative poly(A) cluster discovery rate varies with
        the minimum value of the different attributes of the 3UTR.
        """
        # One plot for each attribute; but in that plot all datasets should be
        # present. Potentially make different ranges for the different
        # attributes.

        colors = ['c', 'r', 'b', 'g']
        for atr in attributes:
            int_ranges = intervals[atr]

            (fig, ax) = plt.subplots()

            utr_count = []

            col_dict = dict(zip(sensitivities.keys(), colors))


            for (dset, atr_dict) in sensitivities.items():

                # unzip the false positive rate and nr of utrs
                (fals_pos, utr_nrs) = zip(*atr_dict[atr])

                # save utr_nrs for making a bar plot later
                utr_count.append(utr_nrs)

                # Get the x_coordinates
                x_coords = range(1,len(fals_pos)+1)

                # Make a line-plot of the fals_pos
                ax.plot(x_coords, fals_pos, label=dset, c=col_dict[dset], lw=2)

            # Set y-ticks
            ax.set_ylim((0,1))
            yticks = np.arange(0,1.1,0.1)
            ax.set_yticks(yticks)
            ax.set_yticklabels([format(val*100,'.0f')+'%' for val in yticks])
            ax.yaxis.grid(True)

            # Make a bar plot of utr_nrs on yaxis nr 2
            # Get the counts of 3UTRs
            mean_counts = np.mean(utr_count, axis=0)
            std_counts = np.std(utr_count, axis=0)

            # Create a 'x-twinned' y axis.
            ax2 = ax.twinx()
            x_coords = range(1,len(mean_counts)+1)
            ax2.bar(x_coords, mean_counts, color='#4C3380', yerr=std_counts,
                   width=0.6, align='center', label='# of 3UTRs')
            ax2.set_ylabel('Number of 3UTRs', size=15)

            # Set the colors and fontsizes of the ticks
            for tl in ax2.get_yticklabels():
                tl.set_color('#4C3380')
                tl.set_fontsize(12)
                tl.set_fontweight('bold')

            # Some hack to get the line-plot in front
            ax.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
            ax.patch.set_visible(False) # hide the 'canvas'

            # Set x-ticks
            ax.set_xticks(x_coords)
            xlabels = ['('+str(v[0])+', '+str(v[1])+']' for v in int_ranges]
            ax.set_xticklabels(xlabels)
            ax.legend(loc='center right')
            ax.set_ylabel('Sensitivity of poly(A) recovery', size=20)
            ax.set_xlabel('RPKM ranges for 3UTRs', size=20)

            # Set xlim so that 
            ax.set_xlim((0.5, max(x_coords)+0.5))

            title = 'We detect most poly(A) clusters for high-RPKM 3UTRs'
            ax.set_title(title, size=20)


    def rpkm_dependent_epsilon(self, distances, rpkm_intervals, titles, order):
        """
        The distances dict is like this
        distances[dset][outp_var] = [[interv1], [intrv2], ...,[intrvN]]
        where each of the intrv are the values of those output variables for all
        UTRs in the intervals given in the input rpkm_intervals.

        2 results:
            1) UTR-length is independent of RPKM (short and long UTRs are not
            differently expressed)
            2) The epsilon-length parameter should only be trusted for UTRs with
            RPKM more than 1.
        """
        for (dset_name, outp_dict) in distances.items():

            nr_outp_var = len(order)
            (fig, axes) = plt.subplots(nr_outp_var, sharex=True)
            for (indx, outp_variable) in enumerate(order):
                ax = axes[indx]
                ax.boxplot(outp_dict[outp_variable])
                ax.set_title(titles[outp_variable])

            # Reduce the space between subplots
            fig.subplots_adjust(hspace=0.1)
            # Hide ticks for all but lowest axis
            plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

            # Set x-tick labels
            x_range = range(1, nr_outp_var+1)
            xlabels = ['('+str(v[0])+', '+str(v[1])+']' for v in rpkm_intervals]
            ax.set_xticklabels(xlabels)
            ax.set_xlabel('RPKM ranges for 3UTRs', size=20)

            # Set 'super title'
            fig.suptitle('3UTR length is independent of RPKM, but relative '\
                         'length is highly variable for low RPKM for {0}'\
                         .format(dset_name), size=23)

    def wc_compartment(self, co_occurence):
        """
        Bar-plot of the co-occurence of polyA clusters in different compartments
        """
        # title translation dict
        title_transf = {'Cytoplasm Whole_Cell': 'C + WC',
                        'Cytoplasm Nucleus Whole_Cell': 'C + N + WC',
                        'Cytoplasm': 'C',
                        'Whole_Cell': 'WC',
                        'Nucleus': 'N',
                        'Nucleus Whole_Cell': 'N + WC',
                        'Cytoplasm Nucleus': 'C + N'
                       }
        # make one plot per cell line
        for (cell_line, comp_sets) in co_occurence.items():
            # Skip non-included cell lines
            if comp_sets == {}:
                continue

            # Modify the thing to remove the single-compartment ones
            for name, count in comp_sets.items():
                if len(name.split()) == 1:
                    comp_sets.pop(name)

            for_plot = sorted([(v,k) for (k,v) in comp_sets.items()], reverse=True)
            (heights, labels) = zip(*for_plot)
            transf_labels = [title_transf[lab] for lab in labels]
            x_coords = range(1, len(for_plot)+1)

            fig, ax = plt.subplots()
            ax.bar(x_coords, heights, align='center', width=0.5)
            ax.set_xticks(x_coords)
            ax.set_xticklabels(transf_labels)

            ax.set_ylabel('Poly(A) cluster count', size=20)
            ax.set_xlabel('Poly(A) clusters found in the same position in more than'\
                          ' one compartment', size=20)
            fig.suptitle('Few poly(A) clusters are common to only cytosol and '\
                         'nucleus -- as expected', size=20)

    def utr_length(self, for_pie, cell_line):
        """
        For one cell-line, plot the difference between the compartments, and
        whole cell + MAX(compartment) for estimation of data integrity.
        # REMEMBER: WHOLE CELL IS THE SUM OF NUCLEUS AND CYTOPLASM
        # WHOLE CELL = NUCLEUS + CYTOPLASM
        # THUS YOU ALWAYS EXPECT len(WHOLE_CELL) = MAX(len(NUCLEUS),len(CYTOPLASM)
        # FOR HIGH RPKM -- CAN YOU TEST FOR THIS? IT WOULD BE A BOON TO THE
        # ACERTATION OF THE VERACITY OF YOUR DATA

        Note with plots: are you comparing relevant enteties? Maybe two plots?
        1st is Expressed vs Not expressed and lowly expressed (explain what they
        are).

        Remember that this is what you select from the annotation -- after
        removing for overlapping stuffz! That number should be included here in
        the first pie-plot. There should be a 'removed due to overlap' slice of
        pie.
        """

        # Two subplots :)
        # TODO aspect ratio is not good. The figures don't come out good. Color
        # must be changed. Maybe not colors but stripes? "Screened out" should
        # appear. 
        # TODO INCLUDE absoulte numbers
        (fig, axxes) = plt.subplots(1,2)

        # Make count and frequency dictionaries 
        count_dict = dict((name, len(entry)) for (name, entry) in
                          for_pie.items())
        total_count = sum(count_dict.values())
        freq_dict = dict((n, l/total_count) for (n, l) in count_dict.items())

        # First plot will only have the expressed vs nonexpressed parts
        expressed = ['Undetermined', 'Same length', 'Different length']
        non_expressed = ['Lowly expressed', 'Not expressed']

        # Get the expressed and non_expressed parts
        # What I need are the fractions and labels of expressed and nonexpressd
        expr_namefreq = [('Expressed', sum([freq_dict[exp] for exp in expressed]))]
        non_expr_namefreq = [(name, freq_dict[name]) for name in non_expressed]
        namefreq = expr_namefreq + non_expr_namefreq
        labels, freqs = zip(*namefreq)

        # Sort labels and freqs
        freqs, labels = zip(*sorted(zip(freqs, labels)))

        # Frequencies of the above
        ax = axxes[0]
        ax.pie(freqs, labels=labels, autopct='%1.1f%%', shadow=True)
        ax.set_title('All 3UTRs', size=20)

        # Second plot is for the details of the expressed 3UTRs
        ax = axxes[1]

        exprsum = sum([count_dict[name] for name in expressed])
        expr_namefreqs = [(name, count_dict[name]/exprsum) for name in expressed]
        labels, piefracs = zip(*expr_namefreqs)

        # Sort piefracs, explode, and labels according to piefracs
        piefracs, labels = zip(*sorted(zip(piefracs, labels), reverse=True))
        ax.pie(piefracs, labels=labels, autopct='%1.1f%%', shadow=True)
        ax.set_title('Expressed 3UTRs', size=20)

        # Set figure super-title
        fig.suptitle('Fraction of 3UTRs with different length in {0}'\
                     .format(cell_line), size=23)

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

                #if dset_l[-1] != []:
                    #debug()

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

    I would like a structure like this:

        dsets['cell lines'] = ([HeLa][comp1, comp2, comp3], [K562][com[1, ...]])

        and

        dsets['compartments'] = ([comp1][HeLa, K562, ...], [comp2][HeLa, K563, ...] )

    I will use the following terminology
    dsets[dset][dset_inverse] = utr_dict
    because the second level is the inverse of the first level.

    """
    # Define your 'cell lines' and 'compartments' variable names: the datasets
    # will be scanned for these
    cell_lines = ['K562', 'K562_Minus', 'HeLa', 'GM12878']
    compartments = ['Whole_Cell', 'Nucleus', 'Cytoplasm', 'Chromatin']

    dsets = {'cell lines': dict((cln, {}) for cln in cell_lines),
             'compartments': dict((comp, {}) for comp in compartments)}

    # Do everything for the length files
    length_files = settings.length_files()
    cluster_files = settings.polyA_files()

    if svm:
        # Get a dictionary indexed by utr_id with the svm locations
        svm_dict = settings.svm_file()

    # Check if all length files exist or that you have access
    [verify_access(f) for f in length_files.values()]

    for dset_name in settings.datasets:

        lengthfile = open(length_files[dset_name], 'rb')
        clusterfile = open(cluster_files[dset_name], 'rb')

        # Get headers
        l_header = lengthfile.next()
        c_header = clusterfile.next()

        # Create the utr objects from the length file
        utr_dict = dict((line.split()[5], UTR(line)) for line in lengthfile)

        # Add cluster objects to the utr objects (line[3] is UTR_ID
        for line in clusterfile:
            utr_dict[line.split()[3]].clusters.append(Cluster(line))

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

        # Add the utr_dict for this cellLine-compartment dataset to BOTH the
        # cell-line and compartment dictionaries.
        for cln in cell_lines:
            if cln in dset_name:
                for comp in compartments:
                    if comp in dset_name:
                        dsets['cell lines'][cln][comp] = utr_dict

        for comp in compartments:
            if comp in dset_name:
                for cln in cell_lines:
                    if cln in dset_name:
                        dsets['compartments'][comp][cln] = utr_dict

    return dsets

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

def utr_length_diff(dsets):
    """
    Compare the epsilon length of 3UTRs in the different datasets.
    """
    # And what about poly(A) support? You should call a function to check if
    # each epsilon site has poly(A) support within +/- 100 nt or so.

    # Estimation of sensitivity: you expect results to be similar between whole
    # cell and eah of the compartments; you get get sensitivity for length
    # difference and for presence/nonpresence by comparing these for high RPKMs.

    # RESULT: most lenghts are the same. You need to find a way that highlights
    # those that DO differ, both in terms of presence-nonpresence or
    # significantly different expression.

    # New idea: classify cell compartment differences in 3 categories:
        # 1) Same length
        # 2) Indeterminable (whole-cell is too different)
        # 3) Different length (by whole-cell criteria)
        # 4) Mutually exclusive. Found in WC and one compartment, but not in the
        # other.

    # In the long term, do this for each cell line -- and make a second pie of
    # which have different length in only 1 cell line, 2 cell lines, or in 3 or
    # more cell lines. Then investigate those that are found in more than 3 cell
    # lines case by case (should not be too many).

    # If one cell-line is particularily weak, you may leave it out of the
    # comparison (3utr coverage in the nucleus compartment will be my bane)
    # Why is coverage of 3UTR so bad in nucleus?????????????
    # Is the loss of polyadenylation in nucleus solely due to low rpkm?

    p = Plotter()

    # Set for checking of you have enough compartments
    WC_N_C = set(['Whole_Cell', 'Nucleus', 'Cytoplasm'])
    length_dict = {}
    for cell_line, comp_dict in dsets['cell lines'].items():

        # Skip empty ones
        if comp_dict == {}:
            continue

        # check if this cell line has at least (WC, Nuc, and Cyt)
        if not WC_N_C.issubset(comp_dict.keys()):
            print('Not enough compartments for comparison. Skipping.')
            continue

        # Get all combinations of the compartments and whole-cell compartments
        all_combs = list(combins(comp_dict.keys(), 2))
        for (compartment, utrs) in comp_dict.iteritems():

            # Set the name ... should be a more standardized way of doing that.
            dset_name = cell_line +' '+ compartment

            for (utr_id, utr) in utrs.iteritems():

                # Get the absoulte 3UTR length determined by epsilon
                if utr.eps_rel_size != 'NA':
                    eps_length = utr.eps_rel_size*utr.length
                else:
                    eps_length = 'NA'

                # Add compartment name, length, and rpkm information
                if utr_id not in length_dict:
                    length_dict[utr_id] = {compartment: (eps_length, utr.RPKM,
                                                       utr.eps_rel_size)}
                else:
                    length_dict[utr_id][compartment] = (eps_length, utr.RPKM,
                                                      utr.eps_rel_size)

        # IDEA: the difference is going to scale with length to a certain
        # degree, simply because the epsilon end gets less accurate for long
        # 3UTRs (I assume).
        # SOLUTION: make a function that shows how the length of WC, N, and C
        # co-vary for increasing RPKM.

        # Question? How many have 'NA' in whole cell but presence in a
        # compartment?  How does this scale for RPKM? Another measure of
        # correctness.

        # Set the minimum rpkm
        rpkm_min = 1
        v_low_rpkm = 0.1 # if below this, consider as 'not expressed'

        # Make another dict with the length-differences themselves 
        for_pie = {'Same length': [],
                   'Undetermined': [],
                   'Different length': [],
                   'Not expressed': [],
                   'Lowly expressed': []}


        # Get the index of the whole cell datapoint in the tuple
        for (utr_id, lpm) in length_dict.iteritems():

            # Get lenghts and rpkms for checking min rpkm
            (lengths, rpkms, eps_rel) = zip(*lpm.values())

            # Define some shortcut parameters to save space
            WClen = lpm['Whole_Cell'][0]
            Nlen = lpm['Nucleus'][0]
            Clen = lpm['Cytoplasm'][0]

            WCrpkm = lpm['Whole_Cell'][1]
            Nrpkm = lpm['Nucleus'][1]
            Crpkm = lpm['Cytoplasm'][1]

            WCepsrel = lpm['Whole_Cell'][2]
            Nepsrel = lpm['Nucleus'][2]
            Cepsrel = lpm['Cytoplasm'][2]

            # Skip those with too small an rpkm
            # ignore rpkms if the length is 'NA'
            noNArpkms = [rp for (i, rp) in enumerate(rpkms) if lengths[i] != 'NA']

            # If non are expressed anywhere, or have VERY low expression 
            # -> 'Not expressed'
            if noNArpkms == []:
                for_pie['Not expressed'].append((utr_id, lpm))
                continue

            # If the maximum rpkm is lower than the set minimum
            # -> 'Low RPKM'
            # the checks so far
            if max(noNArpkms) < rpkm_min:
                for_pie['Lowly expressed'].append((utr_id, lpm))
                continue

            # If one is NA and the rest are 6 and 9?
            # => 'Different length'.
            if (Clen == 'NA') and (Nrpkm > rpkm_min and WCrpkm > rpkm_min):
                for_pie['Different length'].append((utr_id, lpm))
                continue

            if (Nlen == 'NA') and (Crpkm > rpkm_min and WCrpkm > rpkm_min):
                for_pie['Different length'].append((utr_id, lpm))
                continue

            # If any of the comps have 'NA' and have not been picked up by now, or
            # whole cell rpkm is so low that it cannot be used for comparison
            # -> 'Undetermined'
            # add the rpkms if you want to check it out later
            if 'NA' in eps_rel or WCrpkm < v_low_rpkm:
                for_pie['Undetermined'].append((utr_id, lpm))
                continue

            # All differences in relative and absolute length
            absdiffs = [abs(lpm[comb[0]][0] - lpm[comb[1]][0]) for comb in all_combs]
            reldiffs = [abs(lpm[comb[0]][2] - lpm[comb[1]][2]) for comb in all_combs]

            # If they have a relative length within 10% or 100nt variation
            # -> 'same length'
            if max(reldiffs) < 0.1 or max(absdiffs) < 100:
                for_pie['Same length'].append((utr_id, lpm))
                continue
            # If not, check if whole cell supports the longest.
            else:
                # Finally check if some of them have too low RPKM to be caught by
                # the checks so far
                if min(noNArpkms) < rpkm_min:
                    for_pie['Lowly expressed'].append((utr_id, lpm))
                    continue

                # The distance from WC to the smallest length must be 80% or more
                # of the distance from the smallest to the largest length
                # -> 'Different length'
                wc_compmin_reldist = abs((WCepsrel - min(Cepsrel, Nepsrel))) # eg 0.05
                comp_reldist = abs(Cepsrel - Nepsrel) # eg 0.4
                # hey, did I lose my favorite case???
                if wc_compmin_reldist > 0.8*comp_reldist:
                    for_pie['Different length'].append((utr_id, lpm))
                    continue

                # WC could be low by accident of PCR. If the other two are 200 nt or
                # more apart and have high enough RPKM, keep them
                comp_absdist = abs(Clen - Nlen) #absolute distance
                comp_minrpkms = min(Crpkm, Nrpkm) # minimum rpkm
                if (comp_absdist > 200 or comp_reldist > 0.2) and comp_minrpkms < 3:
                    for_pie['Different length'].append((utr_id, lpm))
                    continue

                # If not close enough, you cannot determine
                # -> 'Undetermined'
                else:
                    for_pie['Undetermined'].append((utr_id, lpm))
                    continue
                    # eg 0.08, the upper limit for wc

            # If you make it here, I'd be very interested in how to file you! :)
            debug()

        # only send the lp part of for_pie to the plotter
        for_pie_plot = {}
        for (category, id_lpmlist) in for_pie.items():
            for_pie_plot[category] = [entr[1] for entr in id_lpmlist]

        p.utr_length(for_pie_plot, cell_line)

        # TODO compare the 'for_pie' dict for the different cell lines
        # Save the UTR ID as part of the set and compare them; are they in
        # common? Do they go the same way?

def utr_length_comparison(settings, dsets):
    """
    * compare 3UTR length in general (with boxplot)
    * compare 3UTR length UTR-to-UTR (with boxplot)
    """
    # What we expect with length:
    # 1) WholeCell is an average of cytoplasm and nucleus:
        # i) Difference between cytoplasm and nucleus is expected
        # ii) Difference between WC and each of C and N is not expected
        # iii) Should see more difference between N and N and C and C within
        # compartmens but between cell lines, than when comparing WC with WC,
        # because WC can even out some of those differences.
    # 2) Differences are only significant for high-RPKM and 'long' 3UTRS:
        # i) Short 3UTRs are more susceptible to coverage variation

    # Show the fraction of 3UTRs with different length across the cell
    # compartments Cytosol and Nucleus
    utr_length_diff(dsets)

    # From the datasets you need to differentiate between the different cell
    # lines, the compartments, and 'whole cell'.

    # NOTE from the presentattion Sarah sent you, it looks as if RPKM values for
    # 'all biotypes' is more or less similar between nucl and cytol.

    # Pedro wants to see this for the kidney datasets. Running now.
    # IDEA: juxtapositioned columns in plot: one with different utrs, another
    # with different utrs, supported by poly(A) reads!!!11

    # Make a simple boxplot of all the relative lengths above a certain rpkm
    #simple_length_boxplot(dsets, minRPKM=2)

    # PLOT 1 Box plot of distances from annotated end as function of RPKM
    #rpkm_dependent_distance_from_end(dsets)

    # PLOT 2 similar thing but for changes in upstream/downsteam coverage

def rpkm_dependent_distance_from_end(dsets):
    """
    Make a series of box plots showing the distribution of distances from the
    annoated end in clusters of RPKM, just as for the poly(A)s. Here, I expect
    to see little difference.

    Show the distributions both of absolute distance from end, and relative
    distance from end. Low RPKM 3UTRs might be longer/shorter. Also show the
    distribution of lenghts? (3 plots- length, abs, and norm dist from end)
    """
    p = Plotter()

    def dist_fromend(utrs, rpkm, intvals, outp_variables):
        """
        input:
         1) The string "RPKM" or another delimiting variable
         2) the attribute ranges (min, max) in [min, max) notation
         3) the dataset
         4) the output variables of the datset
        output:
         1) The variables for that datset in that RPKM range
        """
        rpkm_attr = attrgetter(rpkm)
        outp_dict = dict((outp_var, []) for outp_var in outp_variables)
        attribute = dict((o_var, attrgetter(o_var)) for o_var in outp_variables)

        # min_val is on the form [(0,1), (1,5), (5,10), (10,20), (20,40) ... 
        for (intrvl_min, intrvl_max) in intvals:
            intrv_dict = dict((outp_var, []) for outp_var in outp_variables)

            # If we are on the last interval, set interval_max to be 1 billion
            # For RPKM and rna-seq based values, this is as good as max
            if intrvl_max == 'inf':
                intrvl_max = 10**9

            for (utr_name, utr) in utrs.iteritems():
                # Skip those with no epsilon length
                if utr.eps_coord == 'NA':
                    continue

                # Skip those that fall outside this interval
                if (rpkm_attr(utr) < intrvl_min) or (rpkm_attr(utr) > intrvl_max):
                    continue

                # For each of the outp_variables, add to outp_dict
                for o_var in outp_variables:
                    intrv_dict[o_var].append(attribute[o_var](utr))

            for (op_var, intrv_list) in intrv_dict.items():
                outp_dict[op_var].append(intrv_list)

        return outp_dict

    distances = {}
    # The utrs that fall in these bins will be reported
    rpkm_intervals = [(0,1), (1,3), (3,6), (6,10), (10,20), (20,40), (40,80),
                               (80,120), (120,'inf')]
    #rpkm_intervals = [(0,1), (1,'inf')]
    # The output variables
    outp_vars = set(['eps_rel_size', 'eps_abs_size', 'length', 'eps_remainder'])
    titles = {'eps_rel_size': 'Relative epsilon length', 'eps_abs_size':
              'Absolute epsilon length', 'length': 'Length of 3UTR',
              'eps_remainder': 'Distance from epsilon_end to annotated end'}
    # set the order in which you want the plots
    #order = ['length', 'eps_abs_size', 'eps_remainder', 'eps_rel_size']
    order = ['length',  'eps_rel_size']
    # The output variables

    # Go through all dsets and get the respective output variables for the
    # different RPKM clusters 
    for (cell_line, compartment_dict) in dsets['cell lines'].items():
        for (compartment, utrs) in compartment_dict.items():
            dset_name = cell_line +' '+ compartment
            distances[dset_name] = dist_fromend(utrs, 'RPKM', rpkm_intervals,
                                                outp_vars)

    p.rpkm_dependent_epsilon(distances, rpkm_intervals, titles, order)

def simple_length_boxplot(dsets, minRPKM):
    """
    Simply box plots of length-epsilon values.
    """
    # Get the plotter object
    p = Plotter()
    # 1) Box plot of 3UTR lengths in all utrs
    # get dset names and an array of arrays [[], [],...,[]], of lengths
    names = dsets.keys()
    eps_lengths = [dset.get_eps_length(minRPKM=2) for dset in dsets.values()]

    ylim = (0, 1.2)
    #p.boxplot(eps_lengths, names, ylim)

    # 2) Box plot of 3UTR lengths for utrs expressed in all datasets.
    # Get the IDs of utrs that are expressed in each datset
    utrID_sets = [set(dset.expressed_IDs()) for dset in dsets.values()] # list of sets

    # Do the intersection of these sets to get only the ones that are expressed
    # in all datsets
    common = set.intersection(*utrID_sets)

    # Get the lengths of only the commonly expressed utrs
    eps_lengths_c = [dset.get_eps_length(minRPKM=50, IDs=common)
                     for dset in dsets.values()]

    p.boxplot(eps_lengths_c, names, ylim)


def cluster_size_sensitivity(dsets):
    """
    1) Investiage the sensitivity to the number of mapping poly(A) reads
    Make plots like Pedro suggested AND check out how many of the 1-reads and
    2-reads clusters are found in the other compartments. This gives a measure
    of how sensitive the poly(A) clusters are.
    """

    # Get number of utrs that have X clusters with minimum poly(A) read coverage
    # of Y
    p = Plotter()
    for dset in dsets:

        min_covrg = [1, 2, 3, 4, 5, 6 ,7] # minimum coverage for each cluster
        max_cluster = max(utr.cluster_nr for utr in dset.utrs.itervalues())

        # number of clusters in each dset. matrix, use numpy.
        cl_counter = np.zeros([len(min_covrg), max_cluster+1], dtype=int)

        for (utr_id, utr) in dset.utrs.iteritems():
            # If no clusters, add 1 to all the zeros
            if utr.clusters == []:
                cl_counter[:,0] += 1
                continue

            # you are here, so there are some 
            for minc in min_covrg:
                # For each cluster that has more or equal support to the minc
                # value, give a +1 to the cluster count (x) at this mic level (y)
                cl_count = 0
                for cls in utr.clusters:
                    if cls.nr_support_reads >= minc:
                        cl_count += 1

                if cl_count > 0:
                    # minc must be adjusted with 1 to have 0-index of array
                    cl_counter[minc-1, cl_count] += 1
                else:
                    cl_counter[minc-1, 0] += 1

        p.cluster_count(cl_counter)


def udstream_coverage_last_clusters(dsets):
    """
    Comparing upstream/downstream coverage and read coverage for the 3-4 last
    polyA clusters in a 3UTR
    """

    p = Plotter()

    #for dset in dsets:
    for (cell_line, compartment_dict) in dsets['cell lines'].items():

        for (compartment, utrs) in compartment_dict.items():

            # first second thrid and fourth
            # nr of clusters
            clrs_nr = 3
            clus_list = [{'ud_ratio':[], 'support':[]} for val in range(clrs_nr)]

            for (utr_id, utr) in utrs.iteritems():

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

                eps_end = utr.eps_coord

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
    pairs = [pa for pa in combins([ds_name for ds_name in dsets], 2)]
    for pa in pairs:
        (p1, p2) = (pa[0], pa[1])
    title = 'PolyA-site read count variation'
    p.scatterplot(count_dset[p1], count_dset[p2], p1, p2, title)

def compare_cluster_evidence(dsets, super_clusters, dset_2super):
    """
    Find the co-occurence for all the evidence for polyA clusters you have (in
    annotation, svm support, etc).

    Especially, what support do the annotated, alone, have?
    """

    p = Plotter()

    #for dset in dsets:
    for (cell_line, compartment_dict) in dsets['cell lines'].items():

        for (compartment, utrs) in compartment_dict.items():

            all_read_counter = {} # cluster_size : read_count for all clusters
            svm_support = {} # for all clusters: do they have SVM support?
            annot_read_counter = {} # cluster_size : read_count for clusters w/annot
            other_dsets = {} # cluster_size : read_count for clusters in other dsets

            dset_name = cell_line +' '+ compartment

            for (utr_id, utr) in utrs.iteritems():

                if utr.clusters == []:
                    continue

                for cls in utr.clusters:

                    # key that uniquely defines this polyA_cluster
                    # this key will be used to do unions and intersections
                    keyi = dset_name+utr.chrm+utr.strand+str(cls.polyA_coordinate)

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
                            if dn != dset_name: # don't count yourself!!
                                if sup_reads > 1: # maybe set treshold?
                                    other_dsets[cls.nr_support_reads].append(keyi)
                    else:
                        other_dsets[cls.nr_support_reads] = [keyi]

            cluster_dicts = (all_read_counter, svm_support, annot_read_counter)
            titles = ('All clusters', 'SVM_support', 'Annotated_TTS')

            #clusters = (all_read_counter, svm_support, annot_read_counter,
                        #other_dsets)
            #titles = ('All clusters', 'SVM_support', 'Annotated_TTS',
                      #'In_other_compartments')

            # make two figurses: one in terms of all clusters and on of annotated
            # TTS sites

            in_terms_of = (titles[0], titles[2])
            p.join_clusters(cluster_dicts, titles, in_terms_of)

    plt.draw()


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

def recovery_sensitivity(dsets):
    """
    Gauge the false negative rate of discovery of poly(A) clusters. The false
    negative ratio is the fraction of aTTS that have an epsilon-coord close but
    do not have a polyA cluster close.

    Show this function for different cell lines for different compartments.

    A false positive version should also be done: do this on regions of the
    genome that are unlikely to contain 3UTR sequences, such as 5'UTR sequences
    that do not overlap any other genomic region.
    """

    p = Plotter()

    def rec_sens(utrs, utr_attribute, intvals):
        """
        input:
         1) The attribute (rpkm, u/dstream coverage)
         2) the attribute ranges (min, max) in [min, max) notation
         3) the dataset
        output:
         1) The false-negative rate for each of the minimum values
        """
        get_attr = attrgetter(utr_attribute)
        output = []

        # min_val is on the fo [(0,1), (1,5), (5,10), (10,20), (20,40) ... 
        # A neat approach is to iterate through this list and check if the
        # attribute value is SMALLER than the max! this works ONLY because the
        # list has contiguous values; otherwise something would be lost.
        for (interval_min, interval_max) in intvals:
            has_polya = []

            # If we are on the last interval, set interval_max to be 1 billion
            # For RPKM and rna-seq based values, this is as good as max
            if interval_max == 'inf':
                interval_max = 10**9

            for (utr_name, utr) in utrs.iteritems():
                # Skip those that don't have nearby annotated TTS
                if utr.annotTTS_dist == 'NA':
                    continue

                # Skip those that fall outside this interval
                if get_attr(utr) > interval_max:
                    continue

                if get_attr(utr) < interval_min:
                    continue

                # If the utr has no clusters, add a 1 and skip the rest
                if utr.clusters == []:
                    has_polya.append(0)
                    continue

                found = False
                # Check if any of the polyA clusters are around the epsilon_end site
                for cls in utr.clusters:
                    if utr.eps_coord-100 < cls.polyA_coordinate < utr.eps_coord+100:
                        found = True

                if found:
                    has_polya.append(1)
                else:
                    has_polya.append(0)

            # Add the ratio of present clusters to nonpresent clusters
            output.append((np.mean(has_polya), len(has_polya)))

        return output

    sensitivity = {}
    attributes = ['RPKM']
    # The utrs that fall in these bins will be reported
    intervals = {'RPKM': [(0,1), (1,3), (3,6), (6,10), (10,20), (20,40), (40,80),
                               (80,120), (120,'inf')]}

    # Go through all dsets and give 
    for (cell_line, compartment_dict) in dsets['cell lines'].items():
        for (compartment, utrs) in compartment_dict.items():
            dset_name = cell_line +' '+ compartment

            # Ceate empty dict for each d_set combination
            sensitivity[dset_name] = {}

            for atr in attributes:
                intrv = intervals[atr]
                sensitivity[dset_name][atr] = rec_sens(utrs, atr, intrv)

    # Plot the false negative ratio as a function of the attributes
    # NOTE for comparing between cell lines or between compartments, it could be
    # wise to move this plotting function into the previous for-loop and re-make
    # the sensitbity dict for each either 'cell line' or 'compartments'.
    p.rec_sensitivity(sensitivity, intervals, attributes)

def wc_compartment_reproducability(dsets, super_clusters, dset_2super):
    """
    Compare the overlap of polyA sites between compartment-compartment and
    between whole_cell-compartment. You expect higher reproducibility between
    whole_cell-compartment than for compartment-compartment.

    Should you add a min-RPKM or a min-coverage?
    """

    # single cluster minimum read coverage
    min_covr = 2

    mis_aligned_clusters = 0

    # make a co_occurence dict for all cell lines
    co_occurence = {}
    # Get all combinations of compartments of this cell line
    for (cell_line, comp_dict) in dsets['cell lines'].items():
        combs = []
        for cn in range(1, len(comp_dict)+1):
            all_combs = [p for p in combins([co for co in comp_dict.keys()], cn)]
            combs.append(all_combs)

        # Sort the combinations and make a hash-entry by combining them
        # Later you can then sort these things and look them up with certainty.
        co_occurence[cell_line] = {}
        for co in combs:
            for ki in co:
                co_occurence[cell_line][' '.join(sorted(ki))] = 0

    # Go through the super_clusters
    for ikey, clst in super_clusters.iteritems():
        (cl_comp, covrg) = clst

        # if the clusters have been mis-aligned on the second run, skip them,
        # but keep count!
        if len(set(cl_comp)) < len(cl_comp):
            mis_aligned_clusters += 1
            continue

        # Skip the single-hit ones if they have length 1
        if len(covrg) == 1 and covrg[0] < min_covr:
            continue

        cluster_stuff = {}
        # Split up the cell_line--compartment information
        for cc in cl_comp:
            (cl, comp) = cc.split()
            # attribute to each cell line all the compartment in this cluster
            if cl in cluster_stuff:
                cluster_stuff[cl].append(comp)
            else:
                cluster_stuff[cl] = [comp]

        # add a count to the co_occurence list depending on compartmentcomp
        for cline, compments in cluster_stuff.items():
            co_occurence[cline][' '.join(sorted(compments))] += 1

    p = Plotter()
    p.wc_compartment(co_occurence)


def polyadenylation_comparison(settings, dsets, super_clusters, dset_2super):
    """
    * compare 3UTR polyadenylation in general
    * compare 3UTR polyadenylation UTR-to-UTR
    """

    #### Compare the compartments; how many annotated do we find? etc.
    compare_cluster_evidence(dsets, super_clusters, dset_2super)

    #### Get the sensitivity rates for of recovery of poly(A) clusters
    recovery_sensitivity(dsets)

    #### Get the reproducibility rate of poly(A) clusters for each compartment
    #### compared to whole cell (if it is in on of the compartment, it should
    #### presumably be in whole cell as well, albeit dilluted, that is, coverage
    #### should be lower).
    #wc_compartment_reproducability(dsets, super_clusters, dset_2super)

    #### Correlate the coverage counts of common polyadenylation sites between
    #### clusters
    #correlate_polyA_coverage_counts(dsets, super_clusters) NOT DSETS READY

    #### Get the change in coverage ratio of the 3 (or so) last polyadenylation
    #### sites from the 3 UTR
    #udstream_coverage_last_clusters(dsets)

    #### See the effects of excluding pA sites with 1, 2, 3, etc, coverage
    #cluster_size_sensitivity(dsets) NOT DSETS READY

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

def get_reads_from_file(ds, dsets, finder, pAread_file, utr_bed, utr_exons):
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
            #### TEMPORARY CODE WHILE WORKING WITH OLD FILES
            this_strands.append(utr_exon_pAreads['this_strand'][1])

        # Append the 'other_strand' clusters to the polyA clusters of the
        # UTR objects.
        if utr_exon_pAreads['other_strand'][0] != []:
            other_strand = utr_exon_pAreads['other_strand']
            # Go though all the clusters of this utr

            ####
            # TEMPORARY FIX UNTIL YOU RUN YOUR DATASETS AGAIN
            if utr_id in dsets[ds].utrs:
            # REMOVE THIS WHEN FINISHED XXX
            ####

                for cluster in dsets[ds].utrs[utr_id].clusters:
                    # Go though all the polyA reads 
                    for (cl_nr, cl_mean) in enumerate(other_strand[0]):
                        # If they match, update the cluster with all the genomic
                        # positions of the reads at this cluster
                        if cl_mean == cluster.polyA_coordinate:
                            cluster.all_pA_coords = other_strand[1][cl_nr]

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

    chr1 = False
    #chr1 = True
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
            pass
        else:
            continue

        # Update the datasets in dsets with 'other_strand' polyA information
        # As well, get a list 'this_strands' which contains the polyA clusters
        # in the other strands.
        this_strands = get_reads_from_file(file_ds, dsets, finder, pAread_file,
                                           utr_bed, utr_exons)

        size_cutoff = 30
        # Do the 'this_strands' analysis. Return distribution of distances from
        # mean and sizes of clusters (sizes are up to a cutoff value)
        # TODO get distribution of how many are close to annotated sites
        t_strnd_dists = siz_dis_this_strands(this_strands)

        # Get distribution of cluster sizes and distances from cluster means, as
        # well as the number of clusters close to annotated ones, and, further,
        # the size and distance distribution from these annotated ones.
        (o_strnd_dists, o_strnd_dists_anot) = siz_dis_other_strands(file_ds,
                                                                    dsets)

        # Compare the distributions to each other
        compare_this_other(t_strnd_dists, o_strnd_dists, o_strnd_dists_anot,
                           size_cutoff)

        # 2) Calculate the distribution of polyA variation about a site
        # 3) Calculate the distance from polyA-site to 1) closest 2) best PAS
        # 4) Intersect with the UTR objects to get the PAS type
        #   -- Get the PAS distribution
        # 4) Get the variation about a SINGLE polyA-site. It tells us something
        # about the abundance of different types of transcript? Or just some
        # artifact of the PCR process? It possibly points to the minimum number
        # of unique transcripts in the sample. 

        # 5) The same stuff but for 'other_strand'.

def compare_this_other(this_strnd, oth_strnd, oth_strnd_anot, cutoff):
    """
    Plot the distributions of distances of poly_A reads as well and the
    distribution of the sizes of the polyA_clusters. Compare this_strand to
    other_strand and for other_strand, compare those that are annotated.
    """

    p = Plotter()

    (this_sizes, this_dists) = this_strnd
    (other_sizes, other_dists) = oth_strnd
    (annot_other_sizes, annot_other_dists) = oth_strnd_anot

    # These are all dictionaries. Compare the (normalized) distribution of
    # lenghts from all of them
    #sizes = {'this': this_sizes, 'other': other_sizes,
             #'annot_other':annot_other_sizes}
    sizes = {'Opposite strand': this_sizes, 'Annotated strand': other_sizes}

    distances = {'this': this_dists, 'other': other_dists,
                 'annot_other':annot_other_dists}

    ## Get all dists, irrespective of cluster size
    merged_dists = {}
    for (dist_name, dist_dict) in distances.items():
        merged_dists[dist_name] = sum(dist_dict.itervalues(), [])

    #p.distance_histogram(merged_dists)

    ## Create a zero-array for the max sizes
    all_sizes = {}
    for (size_name, size_dict) in sizes.items():
        this_size = np.zeros(cutoff)
        for (size, size_count) in size_dict.iteritems():
            if size < cutoff:
                this_size[size-1] = size_count
            if size >= cutoff:
                this_size[-1] += size_count

        all_sizes[size_name] = this_size

    p.cluster_size_distribution(all_sizes, cutoff)

    debug()

def siz_dis_other_strands(file_ds, dsets):
    """
    Return distributions of the the sizes of clusters and the distances from the
    cluster mean. Also return the same distributions but only for the clusters
    that are near annotated TTS.
    """

    dist_from_mean = {}
    sizes_count = {}

    dist_from_mean_annot = {}
    sizes_count_annot = {}

    utrs = dsets[file_ds].utrs

    for (utr_id, utr) in utrs.iteritems():
        for cls in utr.clusters:
            cluster = cls.all_pA_coords

            # Skip NA clusters that are here from some reason
            if cluster == 'NA':
                continue

            cl_len = len(cluster)
            cl_mean = np.mean(cluster)

            # Unique list of distances from the mean
            cl_dists_from_mean = list(set([cl_mean - pos for pos in cluster]))

            # Save the cluster distancs from mean
            if cl_len not in dist_from_mean:
                dist_from_mean[cl_len] = cl_dists_from_mean
            else:
                for dist in cl_dists_from_mean:
                    dist_from_mean[cl_len].append(dist)

            # Save the cluster read-count
            if cl_len not in sizes_count:
                sizes_count[cl_len] = 1
            else:
                sizes_count[cl_len] += 1


            # Do the same if this is an annotated cluster
            if cls.annotated_polyA_distance != 'NA':

                # Cluster distancs from mean
                if cl_len not in dist_from_mean_annot:
                    dist_from_mean_annot[cl_len] = cl_dists_from_mean
                else:
                    for dist in cl_dists_from_mean:
                        dist_from_mean_annot[cl_len].append(dist)

                # Cluster read-count
                if cl_len not in sizes_count_annot:
                    sizes_count_annot[cl_len] = 1
                else:
                    sizes_count_annot[cl_len] += 1

    return ((sizes_count, dist_from_mean), (sizes_count_annot,
                                            dist_from_mean_annot))

def siz_dis_this_strands(this_strands):
    """
    Input is a list of lists of poly(A) events coming from what I expect is the
    other strand. Oh damn. An obvious check: which of the other strand UTRs are
    close to annotated ones???? :S:S Shit. Maybe I can do an ad-hoc test in this
    script? I can print out bed_clusters of the regions +/- 40 nt around these
    clusters, and I can run them against the annotated poly(A) sites-bedfile I
    have lying around somewhere.

    1) Get the distribution about mean and the distribution of sizes (how many
    with 1, how many with 2, etc) -- this is for comparing with the
    'other_strands' distributions

    2) Get how many of these are close to annotated TTS sites. This is important
    for the control.
    """

    # 1) Go through all clusters and return a list with distances from the mean
    # of each cluster
    dist_from_mean = {}
    sizes_count = {}

    for clusters in this_strands:
        for cluster in clusters:
            cl_len = len(cluster)
            cl_mean = np.mean(cluster)
            # Unique list of distances from the mean
            cl_dists_from_mean = list(set([cl_mean - pos for pos in cluster]))

            # Save the cluster distancs
            if cl_len not in dist_from_mean:
                dist_from_mean[cl_len] = cl_dists_from_mean
            else:
                for dist in cl_dists_from_mean:
                    dist_from_mean[cl_len].append(dist)

            # Save the cluster count
            if cl_len not in sizes_count:
                sizes_count[cl_len] = 1
            else:
                sizes_count[cl_len] += 1


    return (sizes_count, dist_from_mean)


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

def merge_clusters(dsets):
    """
    Merge all clusters from all datasets into one huge cluster.
    """

    all_clusters = {} # keys = cluster_centers
    each_cluster_2all = {} # mapping from the individual clusters to super-cluster

    chrms = ['chr'+str(nr) for nr in range(1,23) + ['X','Y','M']]
    tsdict = dict((chrm, {'+':[], '-':[]}) for chrm in chrms)

    # Put all the polyAclusters in a dict with cluster pos and dset.name
    for (cell_line, compartment_dict) in dsets['cell lines'].items():
        for (compartment, utrs) in compartment_dict.items():
            dset_name = cell_line +' '+ compartment
            for (utr_id, utr) in utrs.iteritems():
                if utr.clusters != []:
                    for cls in utr.clusters:
                        tsdict[cls.chrm][cls.strand].append((cls, dset_name))

    for (chrm, strand_dict) in tsdict.items():
        if strand_dict == {'+': [], '-': []}:
            continue
        for (strand, cls) in strand_dict.iteritems():

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

def write_verified_TTS(clusters):
    """
    Output in bed-format the 3UTR clusters that are verified by both annotation
    and polyA reads.
    """

    # Doing it for He-La
    sure_clusters = {}
    for dset in clusters:
        for k in dset.pAclusters:
            if k.annotated_polyA_distance != 'NA':

                k_beg = k.polyA_coordinate
                bed_entry = (k.chrm, k_beg , k_beg+1, k.ID, '0', k.strand)

                if k.ID in sure_clusters:
                    # If another cluster is nearby, don't add it.
                    for entry in sure_clusters[k.ID]:
                        if entry[1]-30 < k_beg < entry[1] + 30:
                            break
                        else:
                            sure_clusters[k.ID].append(bed_entry)
                else:
                    sure_clusters[k.ID] = [bed_entry]

    # Write to file
    with open('polyA-verified_TTS_sites_HeLa.bed', 'wb') as f:
        for (cls, ent) in sure_clusters.iteritems():
            for entry in ent:
                str_entry = [str(en) for en in entry]
                f.write('\t'.join(str_entry) + '\n')

def before_after_ratio(dsets):
    """
    Give emipirical credence to the before/after ratios
    """

    for (cell_line, compartment_dict) in dsets['cell lines'].items():
        for (compartment, utrs) in compartment_dict.items():
            for rpkm_lim in [0, 1, 5, 20]:
                ratios = []
                dset_name = cell_line + ' ' + compartment
                for (utr_id, utr) in utrs.iteritems():
                    #save before/after ratio if epsilon end is close to annotated end
                    if utr.RPKM < rpkm_lim:
                        continue
                    # screen out those that extend well beyond the annotated end
                    if utr.epsilon_beyond_aTTS > 60:
                        continue
                    if utr.annotTTS_dist != 'NA':
                        us_covrg = utr.eps_upstream_covrg
                        if us_covrg == 0:
                            us_covrg = 0.04
                        ds_covrg = utr.eps_downstream_covrg
                        if ds_covrg == 0:
                            ds_covrg = 0.04
                        ratios.append(math.log(us_covrg/ds_covrg, 2))

                #plt.ioff()
                fig, ax = plt.subplots()
                #(N, bins, patches) = ax.hist(ratios, bins=200) # for log
                (N, ains, patches) = ax.hist(ratios, bins=50) # for +1
                ax.set_xlabel('Log2 ratios of upstream/downstream coverage'\
                              ' for RPKM > {0}'\
                              .format(rpkm_lim), size=20)
                ax.set_title('$Log2(upstream/downstream)$ {0}'.format(dset_name),
                             size=20)

                # Set some vertical lines
                vert_lines = [-1, 0, 1] # for log ratios
                line_nr = len(vert_lines)
                for x_pos in vert_lines:
                    ax.axvline(x=x_pos, c='r')

                # Calculate percentaes the before 0, between 0 and 2, and after 2
                partition = [0 for v in range(line_nr+1)]
                tot_ratios = len(ratios)
                for val in ratios:
                    found = False
                    for (indx, x_pos) in enumerate(vert_lines):
                        if val < x_pos:
                            partition[indx] += 1
                            found = True
                            break
                    if not found:
                        partition[-1] +=1

                percentages = [part/tot_ratios for part in partition]

                # Place text with percentages
                max_yval = ax.get_ylim()[1]
                y_coord = max_yval - math.floor(max_yval*0.05)

                # Get x-coordinates depending on the vertical lines
                (xmin, xmax) = ax.get_xlim()
                x_coords = []
                for (indx, x_pos) in enumerate(vert_lines):
                    if indx == 0:
                        x_coords.append((xmin-x_pos)/2)
                    else:
                        x_coords.append((vert_lines[indx-1]+x_pos)/2 -0.3)
                # Ad hoc at the end add the last x_coord
                x_coords.append((xmax-x_pos)/2)

                for (indx, x_coord) in enumerate(x_coords):
                    percent = format(percentages[indx]*100, '.0f')+'%'
                    ax.text(x_coord, y_coord, percent, size=15)

def rpkm_dist(dsets):
    """
    Simply plot the rpkms of the different datasets
    """
    thresh = 50
    axis_list = []
    for (cell_line, compartment_dict) in dsets['cell lines'].items():
        for (compartment, utrs) in compartment_dict.items():
            rpkms_norm = []
            rpkms_log = []
            dset_name = cell_line +' '+ compartment
            for (utr_id, utr) in utrs.iteritems():
                if utr.eps_coord != 'NA':
                    # SKip those with rpkm = 0
                    if utr.RPKM == 0:
                        continue
                    if utr.RPKM > thresh:
                        rpkms_norm.append(thresh)
                    else:
                        rpkms_norm.append(utr.RPKM)
                    rpkms_log.append(math.log(utr.RPKM, 2))

            fig, ax = plt.subplots()
            axis_list.append(ax)
            (N, bins, blah) = ax.hist(rpkms_log, bins=200)
            ax.set_xlabel('RPKM, log2 transformed', size=20)
            ax.set_ylabel('3UTR count', size=20)
            ax.set_title('Distribution of RPKM values for {0}'.format(dset_name))

            # Percentages above/below 0
            less_zero = len([1 for v in rpkms_norm if v>1])/float(len(rpkms_norm))
            print('\n{0}: Percentage of UTRs with RPKM larger than 1: {1}%'\
                  .format(dset_name, format(less_zero*100, '.0f')))

    # Set the x-axis according to the largest x-interval found
    # Strangely, these things are only correct after a final plt.draw()
    min_x = min([ax.get_xlim()[0] for ax in axis_list])
    max_x = max([ax.get_xlim()[1] for ax in axis_list])
    for ax in axis_list:
        ax.set_xlim((min_x, max_x))
        ax.set_xticks(np.arange(min_x, max_x, 2))
        ax.grid()
        plt.draw()

def main():
    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    # Directory paths for figures and where the output lies
    (savedir, outputdir) = [os.path.join(here, d) for d in ('figures', 'output')]

    # Read UTR_SETTINGS (there are two -- for two different annotations!)
    settings = Settings(os.path.join(here, 'UTR_SETTINGS'), savedir, outputdir, here)
    #settings = Settings(os.path.join(here, 'OLD_ENCODE_SETTINGS'),
                        #savedir, outputdir, here)

    # Get the dsets with utrs and their clusters from the length and polyA files
    # Optionally get SVM information as well
    dsets = get_utrs(settings, svm=True)

    #### RPKM distribution
    #### What is the distribution of 3UTR RPKMS? run this function and find out!
    #rpkm_dist(dsets)

    #### Extend beyond distribution and number
    #### How many 3UTRs extend significantly? beyond the aTTS
    #### NOTE that, with current code, they must be checked by hand, as the
    #### extended area can overlap a genomic feature
    # Give feedback to Simone, then start documenting.
    #beyond_aTTS(dsets) TODO

    ##### Write a bedfile with all the polyA sites confirmed by both annotation and
    ##### 3UT poly(A) reads (For Sarah)
    #write_verified_TTS(clusters) broken, clusters replaced with dsets

    ##### get the before/after coverage ratios for trusted epsilon ends
    #before_after_ratio(dsets)

    # RESULT the before/after ratio is not as informative as we had hoped.
    # I should print out a good few of them that have ratios in the gray zone
    # and check them in the browser. However, this could just mean that the
    # ratio is not important for classifying ends.
    # TODO check it in the browser to see if our intutition is correct ... !:)
    # For the non-beyond_aTTS ones.

    ##### Merge all clusters in datset in the dictionary super_cluster. Also
    ##### return the dset_2super dict, which acts as a translator dict from the
    ##### key of each individual cluster to its corresponding key in super_cluster
    #(super_cluster, dset_2super) = merge_clusters(dsets)

    #output_control(settings, dsets)

    utr_length_comparison(settings, dsets)

    #polyadenylation_comparison(settings, dsets, super_cluster, dset_2super)

    ## The classic polyA variation distributions
    #classic_polyA_stats(settings, dsets)

    debug()

    #reads(settings)

    #data_annotation_correspondence(settings)

    #UTR_processing(settings)

    #other_methods_comparison(settings)

    #rouge_polyA_sites(settings)


if __name__ == '__main__':
    main()

# DISCUSSION #

# On internal representation:
    # Your main dataset will be the cell compartments; therefore, it makes sense
    # to organize your internal datasets to facilitate comparison within
    # cell-lines as well as between compartments of different cell lines.
    # However, forcing this structure upon the internal representation would
    # make it impossible to analyze datasets that do not come from the cell
    # lines (human genome project, drosophila, wtf).
    # An alternative could be to define the two dimensions of comparison in a
    # pre-written list. One dimension could be ('whole_cell' ,'cy..', 'nu..),
    # and another dimension could be ('HeLa', 'K562', ...) and make your plots
    # take these dimensions into account. Then these lists can be modified for
    # each type of datset (gencode, 1000 genomes, drosophila...)
    # I want to do this so that I can quickly compare across compartments or
    # across cell lines without too much fuss.

# You have to convince the reviewers that the difference you see in
# polyadenylation usage is not just down to low coverage or the biased effects
# of sequencing.
# UPDATE: you have shown that you readily detect poly(A) clusters for high RPKM

# How reproducable are the polyA reads?
# If for the same gene, for a similar before/after coverage, what is the
# variation within the number of polyA reads?

# IDEA: measure of coverage-limitation: how many polyA clusters (>5 reads) are
# found with only 1 read in other compartments/datasets?

# Q: is there variation between datasets in polyA-cluster usage in utrs that are
# expressed in all datsets?

# Q: what should we call a reliable cluster?

# Q: poly(A) cluster coverage in each 3UTR -- does it show that 1 3UTR isomer is
# preferentially used?

# Q: In the Whole Cell -- what percentage of sample is cytoplasm and what
# percentage is nucleus? I suspect that there is more cytoplasm than whole cell.

# TODO IDEA TO PROVE THE VERACITY OF THE EPSILON PARAMETER:
    # show the distribution FOR THOSE that fall within 100 nt of annotated
    # end: observe as this distribution narrows as you approach the optimal
    # value. This is nothing else than optimization. However, you will not
    # implement anything to make it rigorous.
    # The proper way to do it would be on 2 cell lines, all compartments,
    # and epsilon values of 0.99, 0.995, 0.998, 0.999, 0.9995,
    # 0.99999 (you need one that is too high!)
    # count the number of UTRs that fall in that rage and also show the
    # distribution about the annotated end value.
    # This is not a repeatable experiment so make sure you do it right. You
    # should run all of them at each epsilon at the same time. Then copy the
    # output dir and label the dir with that epsilon value.
    # ALSO!!!!! print the distance from the poly(A) ends, if found!!! :)
    # THEN!!! The method is backed up as soundly as possible.
    # however: increasing epsilon DEcreases the before/after coverage ratio.
    # There is a trade-off: become more epsilon-accurate, or have better
    # before/after ratios that make you comfortable that a change has really
    # happened.
    # IDEA XXX TODO XXX what happens if you take random before/after ratios
    # in the coverage vector? And what happens if you increase the
    # average-area to 100?

# code snippet for generating colors; this way you don't have to know in advance
# how many colors you need.
#def gen_color():
    #"""generator for getting n of evenly distributed colors using
    #hsv color and golden ratio. It always return same order of colors
 #
    #:returns: RGB tuple
    #"""
    #import colorsys
    #golden_ratio = 0.618033988749895
    #h = 0.22717784590367374
 #
    #while 1:
        #h += golden_ratio
        #h %= 1
        #HSV_tuple = [h, 0.95, 0.95]  # this defines how "deep" are the colors
        #RGB_tuple = colorsys.hsv_to_rgb(*HSV_tuple)
        #yield map(lambda x:str(int(x * 256)), RGB_tuple)

"""
Script for displaying and summarizing the results from utail.py.
"""

from __future__ import division
import os
import ConfigParser
import sys
from itertools import combinations as combins
from copy import deepcopy

from subprocess import Popen, PIPE, STDOUT

import matplotlib.pyplot as plt
#import matplotlib.cm as cm
from matplotlib import lines

#plt.ion() # turn on the interactive mode so you can play with plots
plt.ioff() # turn off interactive mode for working undisturbed

from operator import attrgetter
from operator import itemgetter
import math

import numpy as np
from scipy import stats

import time
import glob
import cPickle as pickle

# For making nice tables
#from TableFactory import *

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
        print ('Warning: debugging mark present.')
        pass
########################################

# Horrible, global variables

first_pas = 'AATAAA'
second_pas = 'ATTAAA'
top_pas = set([first_pas, second_pas])
lower_pas = set(['AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA',
             'AATACA', 'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA',
             'AATAGA'])

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

class BasicUtr(object):
    """ The most simple UTR object: genomic coordinates plus a clusters list
    """
    def __init__(self, chrm, beg, end, strand, ID, line):

        self.chrm = chrm
        self.beg = beg
        self.end = end
        self.strand = strand
        self.ID = ID

        # Create the first clusters object. the rest will be added later
        self.clusters = [Only(line)]

class UTR(object):
    """
    For UTR objects from the 'length' output file in the 'output' directory.
    """

    def __init__(self, input_line):

        # Read all the parameters from line
        (chrm, beg, end, utr_extended_by, strand, ID, epsilon_coord,
         epsilon_rel_size, epsilon_downstream_covrg, epsilon_upstream_covrg,
         annotTTS_dist, epsilon_PAS_type, epsilon_PAS_distance, RPKM,
         avrg_covrg) = input_line.split('\t')

        self.chrm = chrm
        self.beg = str_to_intfloat(beg)
        self.end = str_to_intfloat(end)
        self.length = self.end-self.beg
        self.extended_by = str_to_intfloat(utr_extended_by)
        self.strand = strand
        self.ID = ID
        self.eps_coord = str_to_intfloat(epsilon_coord)
        self.eps_rel_size = str_to_intfloat(epsilon_rel_size)
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

    def __repr__(self):
        return self.ID[-8:]

    def __str__(self):
        return "\nChrm\t{0}\nBeg\t{1}\nEnd\t{2}\nStrand\t{3}\n"\
                .format(self.chrm, self.beg, self.end, self.strand)

class Super(object):
    """
    A simplified class for the super objects
    """

    def __init__(self, dsets, center, covrg, tail_types, tail_infos, strand,
                 nearby_PASes, PAS_distances, annot_distances):

        self.dsets = dsets
        self.polyA_coordinate = center
        self.nr_support_reads = covrg
        self.tail_types = tail_types
        self.tail_infos = tail_infos
        self.strand = strand

        annotdist = [int(d) for d in annot_distances if d != 'NA']

        if annotdist == []:
            self.annotated_polyA_distance = 'NA'
        else:
            self.annotated_polyA_distance = int(np.mean(annotdist))

        # get the average A/T
        # this is a study in it's own!
        As = tail_types.count('A')
        Ts = tail_types.count('T')

        if As > Ts:
            self.tail_type = 'A'
        else:
            self.tail_type = 'T'

        # add the first tail
        self.tail_info = tail_infos[0]

        # if there are more, sum them
        if len(tail_infos) > 1:
            nuc_dict = dict([b.split('=') for b in self.tail_info.split(':')])
            nuc_dict = dict((k, float(v)) for k,v in nuc_dict.items())
            for tail_i in tail_infos[1:]:
                for comsplit in tail_i.split(':'):
                    nuc, count = comsplit.split('=')
                    nuc_dict[nuc] += float(count)

            # avererage with number of tails
            nuc_dict = dict((k, str(v/len(tail_infos))) for k,v in nuc_dict.items())

            self.tail_info = ':'.join(['='.join(item) for item in nuc_dict.items()])

        # add the PAS type and distances
        # associate to each PAS type a distance and take the average
        pases = {}
        for pa_nr, pa in enumerate(nearby_PASes):
            pasplit = pa.split('--')
            displit = PAS_distances[pa_nr].split('--')

            for split_nr, pas in enumerate(pasplit):
                if pas not in pases:
                    pases[pas] = [displit[split_nr]]
                else:
                    pases[pas].append(displit[split_nr])

        self.nearby_PAS = []
        self.PAS_distance = []

        # add average distance to pas list
        for pas, pasdistlist in pases.items():
            # skip the NA things
            if pas == 'NA':
                continue
            else:
                self.nearby_PAS.append(pas)
                distance = np.mean([int(d) for d in pasdistlist if d!='NA'])
                self.PAS_distance.append(int(math.floor(distance)))

        # if it's empty, means all PAS were 'NA'
        if self.nearby_PAS == []:
            self.nearby_PAS.append('NA')
            self.PAS_distance.append('NA')


class Only(object):
    """ For the onlypolyA files. They don't have coverage and have fewer
    parameters than the othes.
    """

    def __init__(self, input_line):

        (chrm, beg, end, utr_ID, strand, polyA_coordinate,
         polyA_coordinate_strand, tail_info, annotated_polyA_distance, nearby_PAS,
         PAS_distance, number_supporting_reads, nr_unique_supporting_reads,
         unique_reads_spread) = input_line.split('\t')

        self.chrm = chrm
        self.beg = str_to_intfloat(beg)
        self.end = str_to_intfloat(end)
        self.ID = utr_ID
        self.strand = polyA_coordinate_strand
        self.tail_info = tail_info

        # Get if it is an 'A' or a 'T' tail
        ga = [g.split('=') for g in self.tail_info.split(':')]
        self.tail_type = sorted([(float(g[1]), g[0]) for g in ga])[-1][-1]

        self.polyA_coordinate = str_to_intfloat(polyA_coordinate)
        self.annotated_polyA_distance = annotated_polyA_distance

        # PAS type and PAS distance are space-delimited
        self.nearby_PAS = [str_to_intfloat(pas) for pas in nearby_PAS.split('#')]

        self.PAS_distance = [str_to_intfloat(dist) for dist in
                             PAS_distance.split('#')]

        self.nr_support_reads = str_to_intfloat(number_supporting_reads)
        self.nr_unique_support_reads = str_to_intfloat(nr_unique_supporting_reads)
        self.uniquet_reads_spread = str_to_intfloat(unique_reads_spread)


class Cluster(object):
    """
    For polyA cluster objects from the 'polyA' output file in the 'output' directory
    """

    def __init__(self, input_line, full_info = True):
        #
        (chrm, beg, end, ID, polyA_number, strand, polyA_coordinate,
         polyA_strand, number_supporting_reads, dstream_covrg, ustream_covrg,
         annotated_polyA_distance, nearby_PAS, PAS_distance,
         rpkm) = input_line.split('\t')

        if full_info:
            self.chrm = chrm
            self.beg = str_to_intfloat(beg)
            self.end = str_to_intfloat(end)
            self.strand = strand
            self.ID = ID

        # to make it compatible
        self.tail_info = 'NA'
        self.tail_type = 'NA'

        self.cluster_nr = str_to_intfloat(polyA_number)
        self.polyA_coordinate = str_to_intfloat(polyA_coordinate)
        self.nr_support_reads = str_to_intfloat(number_supporting_reads)

        self.dstream_covrg = str_to_intfloat(dstream_covrg)
        self.ustream_covrg = str_to_intfloat(ustream_covrg)
        self.annotated_polyA_distance = str_to_intfloat(annotated_polyA_distance)

        # PAS type and distance are space-delimited
        PAS_type = nearby_PAS.split('#')
        self.nearby_PAS = [str_to_intfloat(pas) for pas in PAS_type]

        PAS_distance = PAS_distance.split('#')
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
    Convenience class. One instance is created from this class: it the pertinent
    settings parameters obtained from the UTR_SETTINGS file.
    """
    def __init__(self, settings_file, savedir, outputdir, here, chr1):


        conf = ConfigParser.ConfigParser()
        conf.optionxform = str
        conf.read(settings_file)

        self.datasets = conf.get('PLOTTING', 'datasets').split(':')

        self.savedir = savedir
        self.outputdir = outputdir
        self.here = here

        # which genomic regions are to be investigated
        # is this how you should get them out? I think you should get them from
        # the files you have loaded
        # Return the paths of the onlypolyA files

        self.chr1 = chr1

        self.settings_file = settings_file

        # A bit ad-hoc: dicrectory of polyA read files
        self.polyAread_dir = os.path.join(here, 'polyA_files')

        # The hg19 genome
        self.hg19_path = os.path.join(self.here, 'ext_files', 'hg19')

    def get_annot_polyas(self):
        """
        Simply return the annotated poly(A) sites
        """
        beddir = os.path.join(self.here, 'source_bedfiles')

        import utail as utail

        utail_settings = utail.Settings\
                (*utail.read_settings(self.settings_file))

        return utail.get_a_polyA_sites_path(utail_settings, beddir)

    # Return the paths of the length files
    def length_files(self, region):
        return dict((d, os.path.join(self.here, self.outputdir,
                                     'length_'+d+'_'+region))
                    for d in self.datasets)

    # Return the paths of the polyA files
    def polyA_files(self, region):
        return dict((d, os.path.join(self.here, self.outputdir,
                                     'polyA_' + d+'_'+region))
                    for d in self.datasets)

    # Return the paths of the onlypolyA files
    def only_files(self, region):
        return dict((d, os.path.join(self.here, self.outputdir,
                                     'onlypolyA_'+d+'_'+region))
                    for d in self.datasets)

    # Return the paths of the polyAstats files
    def polyAstats_files(self, region):
        return dict((d, os.path.join(self.here, self.outputdir,
                                     d+'_'+region+'_polyA_statistics'))
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

    def last_three_clustersites(self, clus_list, dset_name):
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
                             "reads and highest drop in coverage\n{0}"\
                             .format(dset_name), size=25)
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

    def join_clusters(self, cluster_dicts, titles, in_terms_of, dset_name):
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
            # Skip the ones where annotated clusters are main
            if dset_ind == 1:
                continue

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
            #myaxes = [ax_list[0]] + ax_list[2:]
            myaxes = ax_list[1:]
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
                         ' supporting data\n{0}'.format(dset_name), size=22)

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


    def rec_sensitivity(self, nrs, mean_clus_fracs, std_clus_fracs, intervals,
                        here):
        """
        Plot how relative # of found intervals correlates with RPKM
        """

        colors = ['AntiqueWhite', 'Aquamarine', 'BlueViolet', 'Brown', 'Coral',
                  'CornflowerBlue', 'Cornsilk', 'Crimson', 'Cyan', 'DarkBlue',
                  'DarkCyan', 'DarkGoldenRod', 'Red']

        int_ranges = intervals

        (fig, ax) = plt.subplots()

        utr_count = nrs

        # Get the x_coordinates
        x_coords = range(1,len(intervals)+1)

        #fraclabel = 'Average polyadenylation sites found relative to annotation'
        # Make a line-plot of the fals_pos
        #ax.plot(x_coords, mean_clus_fracs, label=fraclabel, c='Green', lw=2)
        ax.plot(x_coords, mean_clus_fracs, c='Green', lw=2)
        ax.errorbar(x_coords, mean_clus_fracs, yerr=std_clus_fracs, c='Green',
                    lw=1, fmt=None)

        # Set y-ticks
        ax.set_ylim((0,3))
        yticks = np.arange(0,3.5,0.5)
        ax.set_yticks(yticks)
        ax.set_yticklabels([val for val in yticks])
        ax.yaxis.grid(True)
        # Set the colors and fontsizes of the ticks
        for t in ax.get_yticklabels():
            t.set_color('Green')
            t.set_fontsize(10)

        # Create a 'x-twinned' y axis.
        ax2 = ax.twinx()
        x_coords = range(1,len(utr_count)+1)
        ax2.bar(x_coords, utr_count, color='Blue', width=0.6,
                align='center', label='# of 3UTRs in interval')
        ax2.set_ylabel('Number of 3UTRs', size=13)

        # Set the colors and fontsizes of the ticks
        for tl in ax2.get_yticklabels():
            tl.set_color('Blue')
            tl.set_fontsize(10)
            #tl.set_fontweight('bold')

        # Some hack to get the line-plot in front
        ax.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
        ax.patch.set_visible(False) # hide the 'canvas'

        # Set x-ticks
        ax.set_xticks(x_coords)
        xlabels = ['('+str(v[0])+', '+str(v[1])+']' for v in int_ranges]
        xlabels[-1] = '('+str(int_ranges[-1][0])+', inf)'
        ax.set_xticklabels(xlabels)
        ax.legend(loc='upper right')
        ax.set_ylabel('Discovered/annotated poly(A) sites in 3UTRs', size=13)
        ax.set_xlabel('RPKM ranges for 3UTRs', size=13)

        # Set xlim so that 
        ax.set_xlim((0.5, max(x_coords)+0.5))

        title = 'More poly(A) clusters are found for high-RPKM 3UTRs'
        ax.set_title(title, size=15)

        output_dir = os.path.join(here, 'Results_and_figures', 'GENCODE_report',
                                  'Figures')
        filename = 'More_polyA_clusters_for_high_RPKM_3UTRS'
        filepath = os.path.join(output_dir, filename+'.pdf')
        fig.savefig(filepath, format='pdf')
        filepath = os.path.join(output_dir, filename+'.eps')
        fig.savefig(filepath, format='eps', papertype='A4')

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
                ax.set_ylabel(titles[outp_variable])

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

    def utr_length(self, for_pie, cell_line, index_dict):
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

        (fig, axxes) = plt.subplots(1,2)

        # First plot will only have the expressed vs nonexpressed parts
        expressed = ['Undetermined', 'Same length', 'Different length']
        non_expressed = ['Lowly expressed', 'Not expressed']

        # Make count and frequency dictionaries 
        count_dict = dict((name, len(entry)) for (name, entry) in
                          for_pie.items())
        # Coutn dict for updating labels
        labelcount = dict((name, len(entry)) for (name, entry) in
                          for_pie.items())
        labelcount['Expressed'] = sum([count_dict[exp] for exp in expressed])
        #
        total_count = sum(count_dict.values())
        freq_dict = dict((n, l/total_count) for (n, l) in count_dict.items())

        # Get the expressed and non_expressed parts
        # What I need are the fractions and labels of expressed and nonexpressd
        expr_namefreq = [('Expressed', sum([freq_dict[exp] for exp in expressed]))]
        non_expr_namefreq = [(name, freq_dict[name]) for name in non_expressed]
        namefreq = expr_namefreq + non_expr_namefreq
        labels, freqs = zip(*namefreq)

        # Sort labels and freqs
        freqs, labels = zip(*sorted(zip(freqs, labels)))
        # Update labels to hold the absolute numbers as well
        labels = [lab +' ('+str(labelcount[lab])+')' for lab in labels]

        # Frequencies of the above
        ax = axxes[0]
        cols = ['r', 'g', 'b']
        ax.pie(freqs, labels=labels, autopct='%1.1f%%', colors=cols, shadow=True)
        ax.set_title('All 3UTRs', size=20)

        # Second plot is for the details of the expressed 3UTRs
        ax = axxes[1]

        exprsum = sum([count_dict[name] for name in expressed])
        expr_namefreqs = [(name, count_dict[name]/exprsum) for name in expressed]
        labels, piefracs = zip(*expr_namefreqs)

        # Sort piefracs, explode, and labels according to piefracs
        piefracs, labels = zip(*sorted(zip(piefracs, labels), reverse=True))
        col = {'Same length': 'b', 'Different length':'#FFFF00',
                'Undetermined': '#000000'}
        cols = [col[lab] for lab in labels]
        # Update labels to hold the absolute numbers as well
        labels = [lab +' ('+str(labelcount[lab])+')' for lab in labels]
        ax.pie(piefracs, labels=labels, colors=cols, autopct='%1.1f%%', shadow=True)
        ax.set_title('Expressed 3UTRs', size=20)

        # Set figure super-title
        fig.suptitle(cell_line, size=22)
        # can you add the total number after the labels? WholeCell (2001)

    def compartment_congruence_table(self, isect_utrs, random_isects, super_pie, bias):
        """
        For utrs that have different length in the different compartments,
        compare the lengths by plotting.
        """

        # TABLE 1 JUST THE CELL LINES -- NO TALK OF INTERSECTION
        # First make a table of the different cell line lengths
        alones = RowSpec(
            ColumnSpec('cell_line', 'Cell line', width=1),
            ColumnSpec('diflen_nr', 'Total 3UTRs w/different length', width=1),
            ColumnSpec('long_nuc', 'Percent longest in nucleus', width=1),
            ColumnSpec('long_cyt', 'Percent longest in cytosol', width=1),
            ColumnSpec('rpkm_nuc', 'Percent highest RPKM in nucleus', width=1),
            ColumnSpec('rpkm_cyt', 'Percent highest RPKM in cytosol', width=1))

        rows = []
        # Getting the table row-values
        for cell_line, stat_dict in bias.items():
            diflen_nr = len(super_pie[cell_line]['Different length'])

            cyt_frac = format(stat_dict['Cytosol-longest']*100, '.0f')+ ' %'
            nuc_frac = format(stat_dict['Nucleus-longest']*100, '.0f')+ ' %'

            cyt_rpkm = format(stat_dict['Cytosol RPKM largest']*100, '.0f')+ ' %'
            nuc_rpkm = format(stat_dict['Nucleus RPKM largest']*100, '.0f') + ' %'

            rows.append(alones({'cell_line': cell_line,
                                'diflen_nr': diflen_nr,
                                'long_nuc': nuc_frac,
                                'long_cyt': cyt_frac,
                                'rpkm_nuc': nuc_rpkm,
                                'rpkm_cyt': cyt_rpkm}))

        # Writing the table itself
        outfile = open('table1.pdf', 'wb')
        mytable = PDFTable('Tab1', 'test', headers = [alones]).render(rows)
        outfile.write(mytable)
        outfile.close()

        # TABLE 2: The UTRS that intersect
        # Make a the columns for the table
        tog_rows = []
        together = RowSpec(
            ColumnSpec('cell_combo', 'Intersection', width=2),
            ColumnSpec('common_nr', 'Total 3UTRs in common', width=1),
            ColumnSpec('samplr_nr', 'Common 3UTRs by sampling', width=1),
            ColumnSpec('same_way', 'Same direction of lengthening', width=1))

        al_rows = []
        # TABLE 3: The UTRS that intersect -- but each alone
        alone = RowSpec(
            ColumnSpec('cell_line', 'Cell line', width=2),
            ColumnSpec('long_nuc', 'Longest in nucleus', width=1),
            ColumnSpec('long_cyt', 'Longest in cytosol', width=1),
            ColumnSpec('rpkm_nuc', 'Highest RPKM in nucleus', width=1),
            ColumnSpec('rpkm_cyt', 'Highest RPKM in cytosol', width=1))

        for (cl_combo, common_utrs) in isect_utrs.items():
            c_lines = cl_combo.split('+')

            cell_combo = cl_combo
            common_nr = len(common_utrs)
            randd = random_isects[cl_combo]
            samplr_nr = format(randd['mean'], '.1f') + ' +/- '\
            + format(randd['std'], '.1f')

            # Get how many go the same way, and how many go different ways
            same_way = 0
            for utr in common_utrs:
                direction = []
                for cline in c_lines:
                    difflen = super_pie[cline]['Different length']
                    if difflen[utr]['Cytoplasm'][0] > difflen[utr]['Nucleus'][0]:
                        direction.append(1)
                    else:
                        direction.append(-1)

                if len(set(direction)) == 1:
                    same_way +=1

            tog_rows.append(together({'cell_combo': cell_combo,
                                'common_nr': common_nr,
                                'samplr_nr': samplr_nr,
                                'same_way': same_way}))

            # Go through the cell lines in this combo, adding their values
            for cline in c_lines:
                cyto_longer = 0
                nucl_longer = 0
                cyto_rpkm_bigger = 0
                nucl_rpkm_bigger = 0

                difflen = super_pie[cline]['Different length']

                for utr in common_utrs:

                    if difflen[utr]['Cytoplasm'][0] >= difflen[utr]['Nucleus'][0]:
                        cyto_longer +=1
                    else:
                        nucl_longer +=1

                    if difflen[utr]['Cytoplasm'][1] >= difflen[utr]['Nucleus'][1]:
                        cyto_rpkm_bigger +=1
                    else:
                        nucl_rpkm_bigger +=1

                cyfrac = format(100*cyto_longer/common_nr, '.0f') +'%'
                nufrac = format(100*nucl_longer/common_nr, '.0f') +'%'

                cyrpkm = format(100*cyto_rpkm_bigger/common_nr, '.0f') +'%'
                nurpkm = format(100*nucl_rpkm_bigger/common_nr, '.0f') +'%'

                al_rows.append(alone({'cell_line': cline,
                                    'long_nuc': cyfrac,
                                    'long_cyt': nufrac,
                                    'rpkm_nuc': cyrpkm,
                                    'rpkm_cyt': nurpkm}))

        debug()
        # Writing the table itself
        outfile = open('table2.pdf', 'wb')
        mytable = PDFTable('Tab2', 'tust', headers=[together]).render(tog_rows)
        outfile.write(mytable)
        outfile.close()

        outfile = open('table3.pdf', 'wb')
        mytable = PDFTable('Tab3', 'tast', headers=[alone]).render(al_rows)
        outfile.write(mytable)
        outfile.close()
            # RESULT: they certainly don't follow the same pattern: they almost
            # follow the opposite pattern. Is there any dataset bias that cause
            # HeLa to have cyto longer than nucl, or K562 to have nucl longer
            # than cyto? Is it part of the over-all trend in the two
            # differentially expressed ones? You need to make the same table for
            # each individual cell line (is there a shif)

            # Finally, make plots of the lengths of those in common? Plot the
            # length of whole cell, nucleus, and cytoplasm for all

    def isect_utrs_lengths(self, isect_utrs, super_pie):
        """
        Plot the lenghts of those 3UTRs that intersect. Show 10 random 3UTRs if
        the number of intersecting 3UTRs is above 10.
        """

        # make one plot per intersection
        for (cl_combo, utrs) in iseect_utrs.items():

            debug()
            (fig, ax) = plt.subplots()

            ## the arbitrary x axis range
            ind = range(1, len(isect_utrs)+1)

            #ax.bar(ind, counter, align = 'center', facecolor='#777777', width=0.5)

            #ax.set_title('Distribution of read counts of poly(A) cluster for {0}'\
                         #.format(my_title))

            #ax.set_xlim((0, cutoff+1))
            #ax.set_ylim((0, max_height + 0.2*max_height))
            #ax.set_yticks(range(0, int(math.ceil(max_height+0.2*max_height)), 1000))
            #ax.yaxis.grid(True)
        debug()

    def compartment_congruence_barplot(self, isect_utrs, super_pie):
        """
        For each intersection of diff-len 3UTRS, plot the fraction of
        cytosol-longest 3UTRs for each compartment + those only in the fraction.
        """
        pass

    def polydist_plot(self, polydist, annotdist, count, rel_usage):
        """
        Plot the distribution of polyA site usage from 3' to 5'.
        polydist are the count of polyA sites read coverage
        annotdist coutns only those sites that are annotated
        count just counts how many of each 3utr you have

        make 3 plots, one for each compartment (including whole cell). each plot
        has 3 subplots. the subplots are the length 2, length 3, and length 4
        polyA-count 3UTRs. simply plot the count of each with plot(). Should you
        normalize though?
        """

        polycoverage = AutoVivification()
        annotcount = AutoVivification()
        utrcount = AutoVivification()

        # make helper-dict for reorganizing dicts in one loop.
        helper = {'covr': (polydist, polycoverage),
                  'annot': (annotdist, annotcount),
                  'utrcout': (count, utrcount)}

        keepers = [2,3,4]

        # re-organize the dicts to a [comp][count][cell_line] = [1,55,77]
        for (name, changeme) in helper.items():
            for (c_line, comp_dict) in changeme[0].items():
                for comp, count_dict in comp_dict.items():
                    for nr, nrlist in count_dict.items():
                        if nr in keepers:
                            changeme[1][comp][nr][c_line] = nrlist

        # For each cell line, create an averaged poly(A) count for all the
        # '1,2,3,...'. Keep the std.
        avrg_annot = AutoVivification()

        for (compartment, count_dict) in annotcount.items():
            avrg_annot[compartment] = dict((ke, []) for ke in keepers)
            for (count, cl_dict) in count_dict.items():
                for (cl, countlist) in cl_dict.items():
                    avrg_annot[compartment][count].append(countlist)

                # averge depending on the number of cell lines
                mean = np.mean(avrg_annot[compartment][count], axis=0)
                stds = np.std(avrg_annot[compartment][count], axis=0)
                avrg_annot[compartment][count] = (mean, stds)

        # re-structure the rel_usage to a [cell_line][count][compartment]
        # structure
        rel_usage_restruct = AutoVivification()
        for (cell_line, comp_dict) in rel_usage.items():
            for (comp, count_dict) in comp_dict.items():
                for (nr, nrlist) in count_dict.items():
                    rel_usage_restruct[cell_line][nr][comp] = nrlist

        # Get average and std of the relative usage of polyA sites to the
        # 3'-one. If the variation is too big, make box-plots?
        # AND do log-2 plots
        # AND make one plot per cell line ... ... ... it's too much information
        # alrady. Probably you'll only keep one cell line anyway.

        # old version of this plot.. 
        #old_plot()

        for (cell_line, nr_dict) in rel_usage_restruct.items():
            if cell_line != 'GM12878':
                continue

            # shorten the coutn_dict
            newcount = dict((nr, nr_dict[nr]) for nr in nr_dict if nr in
                            keepers)

            #fig, axes = plt.subplots(1, len(newcount))
            fig, axes = plt.subplots(len(newcount))

            # set some scaling variables
            tweenboxes = 0.5
            tweengroups = tweenboxes*3

            boxwidth = 0.25

            # x_nr is the subplot number (the polyA cluster number)
            for x_n, (nr, comp_dict) in enumerate(newcount.items()):

                ax = axes[x_n]

                ax.set_title('{0} poly(A) clusters'.format(nr))

                # set coordinates according to scaling variables
                xcoords = np.arange(1, nr, tweengroups)
                if len(xcoords) < nr-1:
                    xcoords = np.append(xcoords, xcoords[-1]+tweengroups)

                # ticklist for makign the ticklabels
                ticklist = []

                comp_dict.pop('Whole_Cell') # Remove whole cell from comparison

                # y_n is the number of boxplot-clusters in the subplot
                for y_n, (comp, nrlist) in enumerate(comp_dict.items()):

                    # make nr boxplots of the log2 transformed values.
                    log2list = np.log2(nrlist)[:,:-1]

                    # Scale for space between the boxes
                    xpos = xcoords + tweenboxes*y_n

                    wtf = ax.boxplot(log2list, widths=boxwidth, positions=xpos)

                    # add the xpos to the ticklist
                    for pos in xpos:
                        ticklist.append((pos, comp))

                    # add the annotated ratio
                    annots = annotdist[cell_line][comp][nr]

                    # skip those with 0 annotated in one of them
                    if 0 not in annots:
                        log2an = np.log2([an/annots[-1] for an in annots])[:-1]

                        for (indx, ypos) in enumerate(log2an):
                            x,y = ([xpos[indx]-0.15, xpos[indx]+0.15], [ypos, ypos])
                            line = lines.Line2D(x, y, lw=3, color='g')
                            ax.add_line(line)

                ax.set_xlim(xcoords[0]-tweenboxes,
                            xcoords[-1]+len(comp_dict)*tweenboxes)

                # Get the ticks and their labels in order
                tickpos, labs = zip(*sorted(ticklist))

                ax.set_xticks(tickpos)
                ax.set_xticklabels(labs)

                ax.set_yticks(range(-3, 4))
                ax.set_yticklabels(range(-3,4))
                ax.set_ylim(-5, 5)
                ax.yaxis.grid(True)

            fig.suptitle('Relative usage of proximal and distal '\
                         'polyA-sites\n{0}'.format(cell_line))

            #fig.subplots_adjust(hspace=0.1)
            #fig.subplots_adjust(wspace=0.2)

        debug()

        # TODO include the annotation coutn somehow. ticklabels 5', 3', and
        # middle? better than 1, 2, 3.
        # Then get the log2 of the average annotation and insert on the plot as
        # a bar of some width.

        # Then re-do the novel PA sites thing for each compartment

        # Then go back and re-examine the 3UTR lengths that have
        # different length in the compartments. Can you support your claims with
        # poly(A) reads?

        # Finally, redo all t
        # This must lead you to reconsidering how you calculate he epsilon
        # value. You can no longer base yourself on the annotated value -- you
        # have to base yourself on the extension. I see a mountain of bugs
        # before me, but it needs must be done.

    def oldplot(self, polycoverage, avrg_annot):

        # create 1 subplot for each compartment, and 1 subplot for each
        # '2,3,4' polyA cluster counts. With 3 compartments, you get 9 subplots.
        comp_nr = len(polycoverage)
        clscount_nr = len(polycoverage.values()[0].keys())
        fig, axes = plt.subplots(comp_nr, clscount_nr)

        #normalize = False
        normalize = True

        # IDEA TODO IDEA
        # From Pedro: make everything normalized to the 3'-proximal one. Thus
        # the 3' would be 1 and the other would be relative to that one. You
        # would do the normaization on each 3UTR.

        for x_n, (compartment, clscount_dict) in enumerate(polycoverage.items()):
            for y_n, (cls_count, cell_linedict) in enumerate(clscount_dict.items()):
                # The base-level x-coords for this cls_count. place the bars
                # incrimentally after this
                xcoords = np.arange(1, cls_count+1)
                # create the axes
                ax = axes[x_n, y_n]
                for z_n, (cell_line, actual_list) in enumerate(cell_linedict.items()):

                    # Normalize the actual_list
                    if normalize:
                        lsum = sum(actual_list)
                        actual_list = [ac/lsum for ac in actual_list]

                    ## Plot the polyA-read counts
                    # update the x-coords for this plot
                    # the space of the bars depend on the size of the x-axis.
                    plotcoords = xcoords + 0.2*z_n
                    #cols = [cm.summer(val) for val in range(cls_count)]
                    cols = ['g', 'r', 'k']
                    ax.bar(plotcoords, actual_list, label=cell_line, width=0.2,
                          align='center', color=cols)

                ## Plot the frequency with which they are annotated
                # get the average and std of the annotated coutns
                anavrg = avrg_annot[compartment][cls_count][0]
                anstd = avrg_annot[compartment][cls_count][1]

                # normalize the annotated counts
                ansum = sum(anavrg)
                anavrg = [av/ansum for av in anavrg]
                anstd = [ast/ansum for ast in anstd]

                # Create a 'x-twinned' y axis.
                ax2 = ax.twinx()
                a_coords = plotcoords + 0.2
                ax2.bar(a_coords, anavrg, color='#4C3380', width=0.2,
                        yerr=anstd, label='Annotation frequency',
                        align ='center')

                ax2.set_ylim(0,1)

                # Set the colors and fontsizes of the ticks
                for tl in ax2.get_yticklabels():
                    tl.set_color('#4C3380')
                    tl.set_fontsize(12)

                # Some hack to get the line-plot in front
                ax.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
                ax.patch.set_visible(False) # hide the 'canvas'

                ## Set axis labels and that stuff
                ax.set_xlim(1-cls_count*0.4, cls_count+cls_count*0.4)
                if x_n == comp_nr-1:
                    ax.set_xlabel('{0} poly(A) clusters'.format(cls_count))
                    ax.set_xticks(xcoords)
                else:
                    ax.set_xticklabels([], visible=False)

                if normalize:
                    ax.set_ylim(0,1)
                if y_n == 0:
                    ax.set_ylabel(compartment)
                else:
                    if normalize:
                        ax.set_yticklabels([], visible=False)

        fig.suptitle('Usage of proximal and distal polyadenylation sites',
                     size=20)

        # Fine-tune: remove space between subplots
        #if normalize:
            #adj = 0.07
        #else:
            #adj = 0.15
        #fig.subplots_adjust(hspace=adj)
        #fig.subplots_adjust(wspace=adj)
        debug()

    def region_figure(self, plot_clusters, normlized, dsets, everything,
                      settings, clr_type):
        """ Create the poly(A) reads/clusters-inregions figure!
        """

        reg_nr = len(everything)

        if reg_nr < 3:
            #(fig, axes) = plt.subplots(int(math.ceil(reg_nr/2))+1, 2, sharey=True)
            (fig, axes) = plt.subplots(int(math.ceil(reg_nr/2))+1, 2)
        else:
            #(fig, axes) = plt.subplots(int(math.ceil(reg_nr/2)), 2, sharey=True)
            (fig, axes) = plt.subplots(int(math.ceil(reg_nr/2)), 2)

        if plot_clusters:
            if clr_type == 'all':
                clreads = 'all clusters'
            if clr_type == '2+':
                clreads = 'clusters with 2+ reads or annotated'
        else:
            clreads = 'reads'

        if normlized:
            reltv = ', normalized to region size. '
        else:
            reltv = ''

        fig.suptitle('Poly(A) {0} in PolyA+ and polyA- samples for '\
                     ' different '\
                     ' genomic regions{1}{2}'.format(clreads, reltv, ' K562'),
                     size=23)

        # The order in which you wish to plot the data. This list will be reduced if
        # not all datasets are present
        plot_order = ['Chromatin', 'Nucleoplasm', 'Nucleus PolyA-', 'Nucleus',
                      'PolyA+', 'Cytoplasm PolyA-', 'Cytoplasm PolyA+']

        item_sum = dict()

        # the max y, so that the plots are comparable
        ymax = 0

        for reg_ind, (region, reg_dict) in enumerate(sorted(everything.items())):

            for cell_line, cl_dict in reg_dict.items():

                # where you keep the height of the bars
                height_dict = AutoVivification()

                for compartment, comp_dict in cl_dict.items():

                    # Get values for the bars for the different plots
                    for replicate, rep_dict in comp_dict.items():
                        for polypl, polyA_statistics in rep_dict.items():

                            if plot_clusters:
                                if normlized:
                                    treads = polyA_statistics['total_sites_normalized']
                                else:
                                    if clr_type == 'all':
                                        treads = polyA_statistics['total_sites']
                                    if clr_type == '2+':
                                        treads = polyA_statistics['total_sites2']
                            else:
                                if normlized:
                                    treads = polyA_statistics['total_reads_region_normalized']
                                else:
                                    treads = polyA_statistics['total_reads']

                            key = ':'.join([compartment, polypl])

                            # 1) Absolute nr of reads
                            height_dict[key][replicate] = treads
                            # 2) sum the total number of treads
                            if key in item_sum:
                                item_sum[key] += treads
                            else:
                                item_sum[key] = treads

                # crop plot order according to the datasets that you actually have
                for item in plot_order[:]:
                    if len(item.split()) == 1:
                        key = ':'.join([item, 'PolyA+'])
                    else:
                        key = ':'.join(item.split())

                    if key not in height_dict.keys():
                        plot_order.remove(item)

                # Note: it would be best if each compartment had the same color; you
                # separate poly(A) + and - with alpha = 0.5
                colors = ['b', 'g', 'c', 'm', 'y', 'r']

                color_map = dict(zip(plot_order, colors))
                cols = colors[:len(color_map)]

                bar_cl_nr = len(plot_order)
                #bar_nr = bar_cl_nr*2
                x_arr = np.arange(1, bar_cl_nr+1) # bar at 1, 2, 3, etc

                # Get the axis
                column = reg_ind % 2
                row = int(math.ceil(reg_ind//2))
                ax = axes[row, column]

                ax.set_title(region, size=16)
                bar_width = 0.25

                means = []
                stds = []

                # loop through the data in the order you like
                for position, data_name in enumerate(plot_order):
                    # chromatin and nucleoplasm don't have polyA+/-. use + as dummy.
                    if len(data_name.split()) == 1:
                        key = ':'.join([data_name, 'PolyA+'])
                    else:
                        key = ':'.join(data_name.split())

                    repl_dict = height_dict[key]

                    # Add the mean and std of the replicates
                    means.append(np.mean([repl_dict['replicate'],
                                          repl_dict['not_replicate']]))
                    stds.append(np.std([repl_dict['replicate'],
                                          repl_dict['not_replicate']]))

                ax.bar(x_arr, means, width=bar_width, yerr=stds, color =
                       cols, alpha=0.9, edgecolor = 'k')

                ax.yaxis.grid(True)
                ax.set_xticks(x_arr)
                ax.set_xticklabels(plot_order, rotation=15, size=10)
                ax.set_xlim((min(x_arr)-0.5, max(x_arr)+0.5))

                ymax = max(ymax, ax.get_ylim()[1])

        # finally, set ymax again for all windows
        for axpair in axes:
            for ax in axpair:
                ax.set_ylim((0, ymax))

        # Finally, plot some summary statistics for each compartment: how many
        # poly(A) reads do we find and how many clusters do we find.
        if reg_nr >= 3:
            lastax = axes[-1][-1]
            sums = []
            maxy = 0
            for position, data_name in enumerate(plot_order):
                # chromatin and nucleoplasm don't have polyA+/-. use + as dummy.
                if len(data_name.split()) == 1:
                    key = ':'.join([data_name, 'PolyA+'])
                else:
                    key = ':'.join(data_name.split())

                sums.append(item_sum[key])
                maxy = max(item_sum[key], maxy)

            lastax.bar(x_arr, sums, width=bar_width, color = cols, edgecolor = 'k')
            lastax.set_ylim((0, maxy+maxy*0.2))
            lastax.yaxis.grid(True)
            lastax.set_xticks(x_arr)
            lastax.set_xticklabels(plot_order, rotation=15, size=10)
            lastax.set_xlim((min(x_arr)-0.5, max(x_arr)+0.5))

            lastax.set_title('Sum of {0} for whole genome'.format(clreads), size=13)

        # remove y-axis except for the last one
        plt.setp([a.get_yticklabels() for a in axes[:-1,1]], visible=False)
        fig.subplots_adjust(wspace=0.1)
        fig.subplots_adjust(hspace=0.7)

        return (fig, axes)

    def all_lying_bar(self, data_dict, regions, title, here):
        """ Plot the poly(A) sites from the different regions
        """

        fractions = ['plus', 'minus']

        titls = {'plus': 'P(A)+',
                 'minus':'P(A)-'}

        # The nr and names of bars in the plot
        #plot_keys = ['all', 'T', 'PAS']
        plot_keys = ['PAS', 'all']
        colors = {'all': 'm', 'PAS': 'b'}

        labels = {'all': 'All', 'PAS': 'With downstream PAS'}

        (fig, axes) = plt.subplots(1,2, sharex=True)
        #plt.ion()
        #plt.ioff()

        for frac_nr, frac in enumerate(fractions):

            plotme = {'all': [], 'PAS': []}

            # get the height of the bars from the input
            for region in regions:

                for k in plot_keys:
                    plotme[k].append(data_dict[region][frac][k])

            # you want to place the left foot of all at 1,2,3, etc
            # you want to place the left foot of T at 1.25
            # you want to place the left foot of PAS at 1.375
            heights = {'all': 0.25, 'T': 0.125, 'PAS': 0.125}

            # number of data-points
            dpoints = len(plotme.values()[0])

            # where you want the plotting to start
            start = 1

            # automated calculation of bar positions given the
            # height/width. this one's a keeper!
            pos = dict()
            for knr, k in enumerate(plot_keys):
                if knr == 0:
                    pos[k] = np.arange(start, dpoints+start)
                else:
                    adjust = sum([heights[plot_keys[x]] for x in
                               range(knr)])
                    pos[k] = np.arange(start+adjust, dpoints+start)

            ax = axes[frac_nr]
            rects = dict() # save the rect information

            # make the actual plots
            for pkey in plot_keys:
                rects[pkey] = ax.barh(bottom=pos[pkey],
                                      width=plotme[pkey],
                                      height=heights[pkey],
                                      color=colors[pkey],
                                      label=labels[pkey])

            # print either the number or percentage
            for pkey, rs in rects.items():
                for r_nr, rect in enumerate(rs):
                    width = int(rect.get_width())
                    xloc = width + 100
                    yloc = rect.get_y()+rect.get_height()/2.0
                    clr = 'black'
                    align = 'left'

                    if pkey == 'all':
                        txt = width
                        fsize=10
                    else:
                        divby = plotme['all'][r_nr]
                        try:
                            txt = format(width/divby, '.2f')
                        except ZeroDivisionError:
                            txt = '0'

                        fsize = 8.5
                        yloc = yloc - 0.03

                    # ylocation, centered at bar

                    ax.text(xloc, yloc, txt,
                             horizontalalignment=align,
                             verticalalignment='center', color=clr,
                             weight='bold', fontsize=fsize)

            # print the total number for 'all', and the percentage of
            # 'all' for the other two
            # specify xticks if needeed

            # put some titles here and there
            # get the y-ticks. they should centered
            center = sum(heights.values())/2.0
            yticks = np.arange(start+center, dpoints+start)

            if frac_nr == 0:
                ax.set_ylabel('Genomic regions', size=20)
                ax.set_yticks(yticks) # set the 3utr-exonic etc
                ax.set_yticklabels(regions) # set the 3utr-exonic etc
            else:
                ax.set_yticklabels([])

            ax.set_ylim(start-0.5, dpoints+1) # extend the view

            ax.set_title(titls[frac], size=21)

            # put the legend only in the top-left corner plot
            ax.legend(loc='upper right')

        # Set xlim (it's shared)
        xlm = ax.get_xlim()
        stepsize = 5000

        ax.set_xlim((0, xlm[1]+stepsize))
        #xticks = range(0, xlm[1]+stepsize, stepsize)
        #ax.set_xticks(xticks)
        #f = lambda x: '' if x%(stepsize*2) else x
        #ticklabels = [f(tick) for tick in xticks]
        #ax.set_xticklabels(ticklabels)
        ax.set_xticklabels([]) # remove xticks

        fig.subplots_adjust(wspace=0.1)
        fig.subplots_adjust(hspace=0.2)
        #fig.suptitle(title+ 'for {0}'.format(titles[key1]), size=20)
        fig.set_size_inches(16,19)

        fig.suptitle('Poly(A) sites for WC, N, and C merged', size=20)

        output_dir = os.path.join(here, 'Results_and_figures', 'GENCODE_report',
                                  'Figures')

        filename = 'All_polyA_different_compartmentsd'
        filepath = os.path.join(output_dir, filename+'.pdf')
        fig.savefig(filepath, format='pdf')
        filepath = os.path.join(output_dir, filename+'.eps')
        fig.savefig(filepath, format='eps', papertype='A4')

    def intersect_lying_bar(self, data_dict, regions, title, here):
        """ Plot the poly(A) sites from the different regions
        """

        compartments = ['Whole_Cell', 'Cytoplasm', 'Nucleus']
        fractions = ['plus_sliced', 'intersection', 'minus_sliced']

        titls = {'plus_sliced': 'P(A)+ unique',
                 'intersection': 'P(A)+/P(A)- intersection',
                 'minus_sliced':'P(A)- unique'}

        # The nr and names of bars in the plot
        #plot_keys = ['all', 'T', 'PAS']
        plot_keys = ['PAS', 'all']
        colors = {'all': 'm', 'T': 'g', 'PAS': 'b'}

        labels = {'all': 'All', 'T': 'Mapped with poly(T)',
                  'PAS': 'With downstream PAS'}

        sidelabels = {'Whole_Cell': 'Whole cell', 'Cytoplasm': 'Cytoplasm',
                      'Nucleus': 'Nucleus'}

        (fig, axes) = plt.subplots(3,3, sharex=True)
        #plt.ion()
        #plt.ioff()

        for comp_nr, comp in enumerate(compartments):
            for frac_nr, frac in enumerate(fractions):

                plotme = {'all': [], 'PAS': []}

                # get the height of the bars from the input
                for region in regions:

                    for k in plot_keys:
                        plotme[k].append(data_dict[comp][region][frac][k])

                # you want to place the left foot of all at 1,2,3, etc
                # you want to place the left foot of T at 1.25
                # you want to place the left foot of PAS at 1.375
                heights = {'all': 0.25, 'T': 0.125, 'PAS': 0.125}

                # number of data-points
                dpoints = len(plotme.values()[0])

                # where you want the plotting to start
                start = 1

                # automated calculation of bar positions given the
                # height/width. this one's a keeper!
                pos = dict()
                for knr, k in enumerate(plot_keys):
                    if knr == 0:
                        pos[k] = np.arange(start, dpoints+start)
                    else:
                        adjust = sum([heights[plot_keys[x]] for x in
                                   range(knr)])
                        pos[k] = np.arange(start+adjust, dpoints+start)

                ax = axes[comp_nr, frac_nr]
                rects = dict() # save the rect information

                # make the actual plots
                for pkey in plot_keys:
                    rects[pkey] = ax.barh(bottom=pos[pkey],
                                          width=plotme[pkey],
                                          height=heights[pkey],
                                          color=colors[pkey],
                                          label=labels[pkey])

                # print either the number or percentage
                for pkey, rs in rects.items():
                    for r_nr, rect in enumerate(rs):
                        width = int(rect.get_width())
                        xloc = width + 100
                        yloc = rect.get_y()+rect.get_height()/2.0
                        clr = 'black'
                        align = 'left'

                        if pkey == 'all':
                            txt = width
                            fsize=10
                        else:
                            divby = plotme['all'][r_nr]
                            try:
                                txt = format(width/divby, '.2f')
                            except ZeroDivisionError:
                                txt = '0'

                            fsize=8.5
                            yloc = yloc - 0.03

                        # ylocation, centered at bar

                        ax.text(xloc, yloc, txt,
                                 horizontalalignment=align,
                                 verticalalignment='center', color=clr,
                                 weight='bold', fontsize=fsize)

                # print the total number for 'all', and the percentage of
                # 'all' for the other two
                # specify xticks if needeed

                # put some titles here and there
                # get the y-ticks. they should centered
                center = sum(heights.values())/2.0
                yticks = np.arange(start+center, dpoints+start)

                if frac_nr == 0:
                    ax.set_ylabel(sidelabels[comp], size=20)
                    ax.set_yticks(yticks) # set the 3utr-exonic etc
                    ax.set_yticklabels(regions) # set the 3utr-exonic etc
                else:
                    ax.set_yticklabels([])

                ax.set_ylim(start-0.5, dpoints+1) # extend the view

                if comp_nr == 0:
                    ax.set_title(titls[frac], size=21)

                # put the legend only in the top-left corner plot
                if frac_nr == 1 and comp_nr == 0:
                    ax.legend(loc='upper right')

        # Set xlim (it's shared)
        xlm = ax.get_xlim()
        stepsize = 5000

        ax.set_xlim((0, xlm[1]+stepsize))
        #xticks = range(0, xlm[1]+stepsize, stepsize)
        #ax.set_xticks(xticks)
        #f = lambda x: '' if x%(stepsize*2) else x
        #ticklabels = [f(tick) for tick in xticks]
        #ax.set_xticklabels(ticklabels)
        ax.set_xticklabels([]) # remove xticks

        fig.subplots_adjust(wspace=0.1)
        fig.subplots_adjust(hspace=0.2)
        #fig.suptitle(title+ 'for {0}'.format(titles[key1]), size=20)
        fig.set_size_inches(14,19)

        output_dir = os.path.join(here, 'Results_and_figures', 'GENCODE_report',
                                  'Figures')

        filename = 'Intersected_nr_of_polyA_different_compartments_non_stranded'
        filepath = os.path.join(output_dir, filename+'.pdf')
        fig.savefig(filepath, format='pdf')
        filepath = os.path.join(output_dir, filename+'.eps')
        fig.savefig(filepath, format='eps', papertype='A4')

    def lying_bar_regions(self, data_dict, regions, title, ID, here):
        """ Plot the poly(A) sites from the different regions
        """

        compartments = ['Whole_Cell', 'Cytoplasm', 'Nucleus']
        fractions = ['+', '-']

        # The nr and names of bars in the plot
        #plot_keys = ['all', 'T', 'PAS']
        plot_keys = ['PAS', 'T', 'all']
        colors = {'all': 'm', 'T': 'g', 'PAS': 'b'}

        labels = {'all': 'All', 'T': 'Mapped with poly(T)',
                  'PAS': 'With downstream PAS'}

        titles = {'2+': 'two or more reads or annotated',
                  '1': 'maximum one read'}

        sidelabels = {'Whole_Cell': 'Whole cell', 'Cytoplasm': 'Cytoplasm',
                      'Nucleus': 'Nucleus'}

        # Make plots both for 2+/annot reads or for 1/reads
        # more honest: 2/1 and show % of annot ...
        for key1 in ['3+', '2+', '1']:

            (fig, axes) = plt.subplots(3,2, sharex=True)
            #plt.ion()
            #plt.ioff()

            for comp_nr, comp in enumerate(compartments):
                for frac_nr, frac in enumerate(fractions):

                    plotme = {'all': [], 'PAS': [], 'T': []}

                    # get the height of the bars from the input
                    for region in regions:
                        thiskey = ':'.join([comp, frac, region])
                        for k in plot_keys:
                            plotme[k].append(data_dict[thiskey][key1][k])

                    # you want to place the left foot of all at 1,2,3, etc
                    # you want to place the left foot of T at 1.25
                    # you want to place the left foot of PAS at 1.375
                    heights = {'all': 0.25, 'T': 0.125, 'PAS': 0.125}

                    # number of data-points
                    dpoints = len(plotme.values()[0])

                    # where you want the plotting to start
                    start = 1

                    # automated calculation of bar positions given the
                    # height/width. this one's a keeper!
                    pos = dict()
                    for knr, k in enumerate(plot_keys):
                        if knr == 0:
                            pos[k] = np.arange(start, dpoints+start)
                        else:
                            adjust = sum([heights[plot_keys[x]] for x in
                                       range(knr)])
                            pos[k] = np.arange(start+adjust, dpoints+start)

                    ax = axes[comp_nr, frac_nr]
                    rects = dict() # save the rect information

                    # make the actual plots
                    for pkey in plot_keys:
                        rects[pkey] = ax.barh(bottom=pos[pkey],
                                              width=plotme[pkey],
                                              height=heights[pkey],
                                              color=colors[pkey],
                                              label=labels[pkey])

                    # print either the number or percentage
                    for pkey, rs in rects.items():
                        for r_nr, rect in enumerate(rs):
                            width = int(rect.get_width())
                            if key1 =='2+':
                                xloc = width + 100
                            else:
                                xloc = width + 500
                            yloc = rect.get_y()+rect.get_height()/2.0
                            clr = 'black'
                            align = 'left'
                            if pkey == 'all':
                                txt = width
                                fsize=10
                            else:
                                divby = plotme['all'][r_nr]
                                txt = format(width/divby, '.2f')
                                fsize=8.5
                                yloc = yloc - 0.03

                            # ylocation, centered at bar

                            ax.text(xloc, yloc, txt,
                                     horizontalalignment=align,
                                     verticalalignment='center', color=clr,
                                     weight='bold', fontsize=fsize)

                    # print the total number for 'all', and the percentage of
                    # 'all' for the other two
                    # specify xticks if needeed

                    # put some titles here and there
                    # get the y-ticks. they should centered
                    center = sum(heights.values())/2.0
                    yticks = np.arange(start+center, dpoints+start)

                    if frac_nr == 0:
                        ax.set_ylabel(sidelabels[comp], size=20)
                        ax.set_yticks(yticks) # set the 3utr-exonic etc
                        ax.set_yticklabels(regions) # set the 3utr-exonic etc

                    ax.set_ylim(start-0.5, dpoints+1) # extend the view

                    if frac_nr == 1:
                        ax.set_yticklabels([])

                    if comp_nr == 0:
                        ax.set_title('P(A){0}'.format(frac), size=21)

                    # put the legend only in the top-left corner plot
                    if frac_nr == 1 and comp_nr == 0:
                        ax.legend(loc='upper right')

            # Set xlim (it's shared)
            if ID == 'sense':
                ax.set_xticks(np.arange(0,1.1,0.1))
                ax.set_xlim((0,1.0))
            else:
                xlm = ax.get_xlim()
                if key1 == '2+':
                    stepsize = 5000
                elif key1 == '3+':
                    stepsize = 4000
                else:
                    stepsize = 10000
                ax.set_xlim((0, xlm[1]+stepsize))
                xticks = range(0, xlm[1]+stepsize, stepsize)
                ax.set_xticks(xticks)
                f = lambda x: '' if x%(stepsize*2) else x
                ticklabels = [f(tick) for tick in xticks]
                ax.set_xticklabels(ticklabels)

            fig.subplots_adjust(wspace=0.1)
            #fig.suptitle(title+ 'for {0}'.format(titles[key1]), size=20)
            fig.set_size_inches(14,19)

            output_dir = os.path.join(here, 'Results_and_figures', 'GENCODE_report',
                                      'Figures')

            filename = 'nr_of_polyA_different_compartments_non_stranded_ABRIDGED'
            filename += '_{0}'.format(key1)
            filepath = os.path.join(output_dir, filename+'.pdf')
            fig.savefig(filepath, format='pdf')
            filepath = os.path.join(output_dir, filename+'.eps')
            fig.savefig(filepath, format='eps', papertype='A4')

    def non_PAS_difference(self, ratio_dict, regions, title, here):
        compartments = ['Whole_Cell', 'Cytoplasm', 'Nucleus']
        fractions = ['+', '-']

        # The nr and names of bars in the plot
        plot_keys = ['T percentage', 'non_gPAS_Ts']
        plot_keys = ['T percentage', 'non_gPAS_Ts', 'non_PAS_Ts']
        colors = {'non_gPAS_Ts': 'm', 'non_PAS_Ts': 'g', 'T percentage': 'b'}

        labels = {'non_gPAS_Ts': 'Non-canonical-PAS clusters',
                  'non_PAS_Ts': 'Non-PAS clusters',
                  'T percentage': 'T percentage'}

        sidelabels = {'Whole_Cell': 'Whole cell', 'Cytoplasm': 'Cytoplasm',
                      'Nucleus': 'Nucleus'}

        # Make plots both for 2+/annot reads or for 1/reads
        # more honest: 2/1 and show % of annot ...
        for key1 in ['morethan1OA', 'only1']:

            (fig, axes) = plt.subplots(3,2, sharex=True)
            #plt.ion()
            #plt.ioff()

            for comp_nr, comp in enumerate(compartments):
                for frac_nr, frac in enumerate(fractions):

                    plotme = {'non_gPAS_Ts': [],
                              'non_PAS_Ts': [],
                              'T percentage': []}

                    # get the height of the bars from the input
                    for region in regions:
                        thiskey = ':'.join([comp, frac, region])
                        for k in plot_keys:
                            plotme[k].append(ratio_dict[thiskey][key1][k])

                    # you want to place the left foot of all at 1,2,3, etc
                    # you want to place the left foot of T at 1.25
                    # you want to place the left foot of PAS at 1.375
                    heights = {'T percentage': 0.25,
                               'non_gPAS_Ts': 0.125,
                               'non_PAS_Ts': 0.125}

                    # number of data-points
                    dpoints = len(plotme.values()[0])

                    # where you want the plotting to start
                    start = 1

                    # automated calculation of bar positions given the
                    # height/width. this one's a keeper!
                    pos = dict()
                    for knr, k in enumerate(plot_keys):
                        if knr == 0:
                            pos[k] = np.arange(start, dpoints+start)
                        else:
                            adjust = sum([heights[plot_keys[x]] for x in
                                       range(knr)])
                            pos[k] = np.arange(start+adjust, dpoints+start)

                    ax = axes[comp_nr, frac_nr]
                    rects = dict() # save the rect information

                    # make the actual plots
                    #for pkey in plot_keys:
                    for pkey in plot_keys[:-1]:
                        rects[pkey] = ax.barh(bottom=pos[pkey],
                                              width=plotme[pkey],
                                              height=heights[pkey],
                                              color=colors[pkey],
                                              label=labels[pkey])

                    # print either the number or percentage
                    for pkey, rs in rects.items():
                        for r_nr, rect in enumerate(rs):
                            width = rect.get_width()
                            xloc = width + 0.01
                            # ylocation, centered at bar
                            yloc = rect.get_y()+rect.get_height()/2.0
                            clr = 'black'
                            align = 'left'
                            #if pkey == 'all':
                            txt = format(width, '.2f')
                            fsize = 9

                            ax.text(xloc, yloc, txt,
                                     horizontalalignment=align,
                                     verticalalignment='center', color=clr,
                                     weight='bold', fontsize=fsize)

                    # print the total number for 'all', and the percentage of
                    # 'all' for the other two
                    # specify xticks if needeed

                    # put some titles here and there
                    # get the y-ticks. they should centered
                    center = sum(heights.values())/2.0
                    yticks = np.arange(start+center, dpoints+start)

                    if frac_nr == 0:
                        ax.set_ylabel(sidelabels[comp], size=20)
                        ax.set_yticks(yticks) # set the 3utr-exonic etc
                        ax.set_yticklabels(regions) # set the 3utr-exonic etc

                    ax.set_ylim(start-0.5, dpoints+1) # extend the view

                    if frac_nr == 1:
                        ax.set_yticklabels([])

                    if comp_nr == 0:
                        ax.set_title('P(A){0}'.format(frac), size=21)

                    # put the legend only in the top-left corner plot
                    if frac_nr == 1 and comp_nr == 0:
                        ax.legend(loc='upper right')

            # Set xlim (it's shared)
            ax.set_xticks(np.arange(0,1.1,0.1))
            ax.set_xlim((0,1.1))

            fig.subplots_adjust(wspace=0.1)
            fig.set_size_inches(14.4,19)

            output_dir = os.path.join(here, 'Results_and_figures', 'GENCODE_report',
                                      'Figures')

            filename = 'non-PAS ratios'
            filename += '_{0}'.format(key1)
            filepath = os.path.join(output_dir, filename+'.pdf')
            fig.savefig(filepath, format='pdf')
            filepath = os.path.join(output_dir, filename+'.eps')
            fig.savefig(filepath, format='eps', papertype='A4')

    def non_PAS_difference_separate(self, ratio_dict, regions, title, here):
        """
        Plot the difference in non-PAS usage
        """

        # figure out how many subplots you will need
        plot_order = ['Nucleus', 'Cytoplasm', 'Whole_Cell']
        x_order = ['non_gPAS_Ts']

        figs = {}
        axes = {}

        for region in regions:
            fig, axs = plt.subplots(2)

            # save for future reference
            figs[region] = fig
            axes[region] = axs

            for str_nr, strand in enumerate(['-', '+']):
                # get out only the dsets for this strand for this region
                these_dsets = {}
                for key, val in ratio_dict.items():
                    if region in key.split(':') and strand in key.split(':'):
                        these_dsets[key.split(':')[0]] = val

                ax = axs[str_nr]

                #heights = [these_dsets[k]['non_gPAS_Ts'] for k in plot_order]
                heights = [these_dsets[k]['non_gPAS_Ts'] for k in plot_order]
                x_pos = range(1, len(heights)+1)

                ax.bar(x_pos, heights, align='center', width=0.6)

                ax.set_xticks(x_pos)
                ax.set_xticklabels(plot_order)

                ax.set_title('Poly(A){0}'.format(strand), size=20)
                ax.set_ylabel('Non-PAS-mediated-polyadenylation-index (NMP-index)',
                              size=10)

                ax.set_ylim(0,1)

            fig.suptitle(region)

            # save figure
            output_dir = os.path.join(here, 'Results_and_figures', 'GENCODE_report',
                                      'Figures')

            filename = 'non-PAS-index'
            filename += '_{0}'.format(region)
            filepath = os.path.join(output_dir, filename+'.pdf')
            fig.savefig(filepath, format='pdf')
            filepath = os.path.join(output_dir, filename+'.eps')
            fig.savefig(filepath, format='eps', papertype='A4')


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
        isect_nrs = [len(set.intersection(set(main_cls[count]),
                                          set(sub_cls[count]))) for count in
                     range(0, cutoff)]

        # Get percent of intersection relative to 'main' dataset (will be all or
        # annot)
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
    # Take the union of all dsetsxcept
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
    all_I_nrs = [len(set.intersection(set(main_cls[count]),
                                          set(all_cls[count]))) for count in
                     range(0, cutoff)]

    # Get percent of intersection relative to 'main' dataset (will be all or annot)
    all_I_pcnt = []
    for (indx, isect_nr) in enumerate(isect_nrs):

        # Only calculate percentage if more than 1 cluster with this read count
        if main_cls[indx] != 0:
            all_I_pcnt.append(all_I_nrs[indx]/len(main_cls[indx]))
        else:
            all_I_pcnt.append(0)

    # Add the number and intersection to the array
    counter[-1,:,0] = all_I_nrs
    counter[-1,:,1] = all_I_pcnt

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
        # Skip headers
        clusterfile.next()

        dset_clusters = []

        for line in clusterfile:
            ln = line.split()
            hsh = '_'.join([ln[0], ln[6], ln[5]])
            #Index should be the location of each cluster: chrm_pos_strand
            dset_clusters.append(Cluster(line))

        clusters.append(PolyaCluster(dset_clusters, dset_name))

    return clusters

def super_falselength(settings, region, batch_key, subset=[], speedrun=False):
    """ Return the super clusters with information only from the poly(A) files
    without using length (coverage) information.

    Region is a list of genomic regions to be treated. Datasets are obtained
    from 'settings'. You can leave out some datasets by putting their names in
    the subset-list.
    """

    dsets = AutoVivification()
    super_3utr = AutoVivification()

    # if sending in just one region, make it into a string
    if type(region) == str:
        regions = [region]

    for region in regions:

        # get only files from this region
        regionfiles = settings.only_files(region)

        # Check if all length files exist or that you have access
        for dset, fpath in regionfiles.items():
            if subset == []:
                verify_access(fpath)
            else:
                if dset not in subset:
                    continue
                else:
                    verify_access(fpath)

        # limit the nr of reads from each dset if you want to be fast
        if speedrun:
            maxlines = 100

        linenr = 0

        for dset_name in settings.datasets:

            # If you have specified a subset of the datasets, skip those that
            # are in that subset
            if subset != [] and dset_name not in subset:
                continue

            # work with files from this region only
            onlyfile = open(regionfiles[dset_name], 'rb')

            # Skip headers
            onlyfile.next()

            # Create the utr objects from the length file
            utr_dict = {}

            for (linenr, line) in enumerate(onlyfile):

                # If a speedrun, limit the nr of 3UTRs to read.
                if speedrun:
                    if linenr > maxlines:
                        break

                # key = utr_ID, value = UTR instance
                (chrm, beg, end, utr_id, strand) = line.split()[:5]

                if utr_id in utr_dict:
                    utr_dict[utr_id].clusters.append(Only(line))
                else:
                    utr_dict[utr_id] = BasicUtr(chrm, beg, end, strand, utr_id,
                                               line)

            # Add the utr_dict for this cellLine-compartment 
            this_cl = dset_name.split('_')[0]
            this_comp = '_'.join(dset_name.split('_')[1:])

            dsets[region][this_cl][this_comp] = utr_dict

        # You send stuff to get super 3utr. the dsets themselves and the merged
        # clusters of the dsets
        #super_3utr[region] = get_super_3utr(dsets[region],
                                            #*merge_clusters(dsets[region]))

        # TODO implement the same thing with bed-tools and obeserve the same
        # results but much faster and with mich simpler code. more memory
        # efficient too. RESULT you save more memory than speed.
        super_3utr[region] = get_super_regions(dsets[region], settings,
                                               batch_key)

    return dsets, super_3utr

def get_super_regions(dsets, settings, batch_key):
    """
    Use bed-tools to cluster the 3utrs

    Let each 3UTR know about all the cell lines and compartments where its
    poly(A) clusters are found.

    Update all the 3UTR clusters in dsetswith their super-cluster key
    At the same time create a 3UTR dict. Each 3UTR has information about the
    compartments and cell lines where it resides.

    dsets has a dsets[cell_line][compartment][feature_dict] structure

    super_3utr[region]

    1) for all cell_lines and compartments in dsets, go through all utrs and
    save all clusters to a bed file. in the val field put number of covering
    reads, and in the name field put the feature_id so we can go back to it
    later. Also put the transcript type and transcript info in name.

    2) slopBed with 15 nt in both directions
    """
    super_features = {}

    # 1) super cluster key: 'chr1-10518929', for example
    # 2) the other datasets where this cluster appears (hela nucl etc)

    cluster_file = os.path.join(settings.here, 'merge_dir', 'merged_pA_' + batch_key)
    temp_cluster_file = cluster_file + '_temp'

    out_handle = open(cluster_file, 'wb')

    for cell_line, comp_dict in dsets.items():
        for comp, utr_dict in comp_dict.items():
            for feature_id, feature in utr_dict.iteritems():

                # store every foudn feature in super_features for lalter
                if feature_id not in super_features:
                    feature.super_clusters = [] # superclusters later
                    super_features[feature_id] = feature

                # Create the super_key and store the super_cluster information
                for cls in feature.clusters:
                    fromwhere = '_'.join([cell_line, comp])

                    PAS = '--'.join(cls.nearby_PAS)
                    PAS_distance = '--'.join([str(d) for d in cls.PAS_distance])

                    name = '%'.join([fromwhere,
                                     feature_id,
                                     cls.tail_type,
                                     cls.tail_info,
                                     PAS,
                                     PAS_distance,
                                     str(cls.annotated_polyA_distance)])
                    # also add PAS_distance and PAS_type
                    out_handle.write('\t'.join([cls.chrm,
                                                str(cls.polyA_coordinate),
                                                str(cls.polyA_coordinate+1),
                                                name,
                                                str(cls.nr_support_reads),
                                                cls.strand+'\n'
                                               ]))
    out_handle.close()

    junctions = os.path.join(settings.here, 'junctions',
                             'splice_junctions_merged.bed')

    # i.5) cut away the exon-exon noise
    cmd = ['subtractBed', '-a', cluster_file, '-b', junctions]
    #p1 = Popen(cmd, stdout=open(temp_cluster_file, 'wb'), stderr=PIPE)
    p1 = Popen(cmd, stdout=open(temp_cluster_file, 'wb'))
    p1.wait()
    #err1 = p1.stderr.read()
    #if err1:
        #debug()

    # ii) expand the entries
    cmd = ['slopBed', '-i', temp_cluster_file, '-g', settings.hg19_path, '-b', '15']
    p2 = Popen(cmd, stdout=open(cluster_file, 'wb'))
    p2.wait()
    #err2 = p2.stderr.read()
    #if err2:
        #debug()

    # iii) merge and max scores
    cmd = ['mergeBed', '-nms', '-s', '-scores', 'sum', '-i', cluster_file]
    #cmd = ['mergeBed', '-nms', '-s', '-scores', 'max', '-i', cluster_file]
    #p3 = Popen(cmd, stdout=PIPE, stderr=PIPE)
    p3 = Popen(cmd, stdout=PIPE)

    # iv) read therough the centered output and update the super_cluster
    for line in p3.stdout:
        (chrm, beg, end, names, maxcovrg, strand) = line.split()

        center = int((int(beg)+int(end))/2)
        polyA_coverage = int(float(maxcovrg))

        dsets = []
        f_ids = []
        tail_types = []
        tail_infos = []
        nearby_PASes = []
        PAS_distances = []
        annot_distances = []

        name_list = names.split(';')
        for names2 in name_list:
            (dset, f_id, t_type, t_info, nearby_PAS, PAS_distance, annot_dist)\
                    = names2.split('%')
            dsets.append(dset)
            f_ids.append(f_id)
            tail_types.append(t_type)
            tail_infos.append(t_info)
            nearby_PASes.append(nearby_PAS)
            PAS_distances.append(PAS_distance)
            annot_distances.append(annot_dist)

        feature_id = set(f_ids).pop()

        # create cluster object
        cls = Super(dsets, center, polyA_coverage, tail_types, tail_infos,
                    strand, nearby_PASes, PAS_distances, annot_distances)

        # add the cluster to super_featuress
        super_features[feature_id].super_clusters.append(cls)

    # any errors?
    #err3 = p3.stderr.read()
    #if err3:
        #debug()

    return super_features

def get_utrs(settings, region, dset_subset=[], speedrun=False):
    """
    Return a list of UTR instances. Each UTR instance is
    instansiated with a list of UTR objects and the name of the datset.

    1) Read through the length file, create UTR objects for each line
    2) If present, read through the polyA file, updating the UTR object created
    in step 1

        dsets[cell_line][comp1][utr1])

    Finally merge the dsets in dsets into a super_3utr structure.

    If dset_subset is != [], restrict the analysis to those dsets in the subset

    """
    # parse the datasets to get the cell lines and compartments in this dataset
    cell_lines = list(set([ds.split('_')[0] for ds in settings.datasets]))

    # Initialize the dsets
    dsets = dict((cln, {}) for cln in cell_lines)

    # Do everything for the length files
    # XXX you must update these things to take the region (UTR etc) into
    # account, since it's now part of the filename
    length_files = settings.length_files(region)
    cluster_files = settings.polyA_files(region)

    # Check if all length files exist or that you have access
    # if you specify a certain subset, don't check those that are not in the
    # subset, since they will not be processed anyway
    for dset, fpath in length_files.items():
        if dset_subset != []:
            if dset not in dset_subset:
                continue
        else:
            verify_access(fpath)

    # limit the nr of reads from each dset if you want to be fast
    if speedrun:
        maxlines = 300

    linenr = 0

    for dset_name in settings.datasets:

        # Filter against dsets not in subset
        if dset_subset != []:
            if dset_name not in dset_subset:
                continue

        lengthfile = open(length_files[dset_name], 'rb')
        clusterfile = open(cluster_files[dset_name], 'rb')

        # Skip headers
        lengthfile.next()
        clusterfile.next()

        # Create the utr objects from the length file
        utr_dict = {}

        for (linenr, line) in enumerate(lengthfile):
            # If a speedrun, limit the nr of 3UTRs to read.
            if speedrun:
                if linenr > maxlines:
                    break
            # key = utr_ID, value = UTR instance
            utr_dict[line.split()[5]] = UTR(line)

        # Add cluster objects to the utr objects (line[3] is UTR_ID
        for line in clusterfile:
            if line.split()[3] in utr_dict:
                utr_dict[line.split()[3]].clusters.append(Cluster(line))

        # Update the UTRs
        for (utr_id, utr_obj) in utr_dict.iteritems():

            # Update with poly(A) cluster info
            utr_dict[utr_id].cluster_nr = len(utr_obj.clusters)

        # Add the utr_dict for this cellLine-compartment 
        this_cl = dset_name.split('_')[0]
        this_comp = '_'.join(dset_name.split('_')[1:])

        dsets[this_cl][this_comp] = utr_dict

    ## Merge all clusters in datset in the dictionary super_cluster. Also
    ## return the dset_2super dict, which acts as a translator dict from the
    ## key of each individual cluster to its corresponding key in super_cluster
    # send them to get_super_2utr to make super_3utrs


    super_3utr = get_super_3utr(dsets, *merge_clusters(dsets))

    return dsets, super_3utr

def get_super_3utr(dsets, super_cluster, dset_2super):
    """
    Let each 3UTR know about all the cell lines and compartments where its
    poly(A) clusters are found.

    Update all the 3UTR clusters in dsetswith their super-cluster key
    At the same time create a 3UTR dict. Each 3UTR has information about the
    compartments and cell lines where it resides.

    3utr[cell_line][compartment].clusters
    3utr.super_cluster_coords = super_key: (chr1, 223432, '-')
    3utr.super_in[cell_line][compartment] = True/False?
    """
    super_3utr = {}

    # 1) super cluster key: 'chr1-10518929', for example
    # 2) the other datasets where this cluster appears (hela nucl etc)

    for cell_line, comp_dict in dsets.items():
        for comp, utr_dict in comp_dict.items():
            for utr_id, utr in utr_dict.iteritems():

                # a [super_coord][cell_line][comp] = covr dict
                utr.super_cover = {}

                # Create the super_key and store the super_cluster information
                for cl in utr.clusters:

                    # superkey for dset_2super
                    superkey = cell_line+' '+comp + utr.chrm + cl.strand +\
                            str(cl.polyA_coordinate)

                    supr_covr = AutoVivification()

                    # Add a in[cell_line][compartment] = covr
                    for dset_name, cov in zip(*super_cluster[dset_2super[superkey]]):
                        (cl_line, compart) = dset_name.split(' ')
                        supr_covr[cl_line][compart] = cov

                    utr.super_cover[dset_2super[superkey]] = supr_covr

                # Add this 3UTR to the super-3utr dict
                if utr_id not in super_3utr:
                    super_3utr[utr_id] = utr
                    super_3utr[utr_id].rpkms = [utr.RPKM]

                # if it's there, add any super-cluster: [cl][comp][covr] that aren't
                # there already
                else:
                    for (s_key, covr_dict) in utr.super_cover.items():
                        if s_key not in super_3utr[utr_id].super_cover:
                            super_3utr[utr_id].super_cover[s_key] = covr_dict

                    ## Add the RPKMs for later
                    super_3utr[utr_id].rpkms.append(utr.RPKM)

    # go through all the super_clusters and add them to super_clusters
    for (utrid, utr) in super_3utr.items():
        utr.super_clusters = []
        for superkey, cl_dict in utr.super_cover.items():
            this_cls = 0
            for cl, comp_dict in cl_dict.items():
                for comp, clstr in comp_dict.items():
                    # add the first with no adding
                    if this_cls == 0:
                        this_cls = clstr
                        this_cls.from_dsets = ['_'.join([cl, comp])]
                        this_cls.tail_types = [clstr.tail_type]
                        this_cls.tail_infos = [clstr.tail_info]
                        this_cls.strands = [clstr.strand]
                    else:
                        this_cls.from_dsets.append('_'.join([cl, comp]))
                        this_cls.nr_support_reads += clstr.nr_support_reads
                        this_cls.tail_types.append(clstr.tail_type)
                        this_cls.tail_infos.append(clstr.tail_info)
                        this_cls.strands.append(clstr.strand)

            if '+' in superkey:
                this_cls.polyA_coordinate = superkey.split('+')[-1]
            else:
                this_cls.polyA_coordinate = superkey.split('-')[-1]

            utr.super_clusters.append(this_cls)

    return super_3utr

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

    # Later, you will add (eps_length, utr.RPKM, utr.eps_rel_size) as a tuple.
    # Make a dictionary for indexing this

    index_dict = {'eps_len':0, 'rpkm':1, 'rel_size':2}

    p = Plotter()

    # Set for checking of you have enough compartments
    WC_N_C = set(['Whole_Cell', 'Nucleus', 'Cytoplasm'])

    # a super-pie-dict for each cell_line
    super_pie = {}

    for cell_line, comp_dict in dsets.items():

        # Skip empty ones
        if comp_dict == {}:
            continue

        # check if this cell line has at least (WC, Nuc, and Cyt)
        if not WC_N_C.issubset(comp_dict.keys()):
            print('Not enough compartments for comparison for {0}. Skipping.'\
                 .format(cell_line))
            continue

        length_dict = {}

        # Get all combinations of the compartments and whole-cell compartments
        all_combs = list(combins(comp_dict.keys(), 2))

        for (compartment, utrs) in comp_dict.iteritems():

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
        for_pie = {'Same length': {},
                   'Undetermined': {},
                   'Different length': {},
                   'Not expressed': {},
                   'Lowly expressed': {}}

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
                for_pie['Not expressed'][utr_id] = lpm
                continue

            # If the maximum rpkm is lower than the set minimum
            # -> 'Low RPKM'
            # the checks so far
            if max(noNArpkms) < rpkm_min:
                for_pie['Lowly expressed'][utr_id] = lpm
                continue

            # If one is NA and the rest are 6 and 9?
            # => 'Different length'.
            if (Clen == 'NA') and (Nrpkm > rpkm_min and WCrpkm > rpkm_min):
                for_pie['Different length'][utr_id] = lpm
                continue

            if (Nlen == 'NA') and (Crpkm > rpkm_min and WCrpkm > rpkm_min):
                for_pie['Different length'][utr_id] = lpm
                continue

            # If any of the comps have 'NA' and have not been picked up by now, or
            # whole cell rpkm is so low that it cannot be used for comparison
            # -> 'Undetermined'
            # add the rpkms if you want to check it out later
            if 'NA' in eps_rel or WCrpkm < v_low_rpkm:
                for_pie['Undetermined'][utr_id] = lpm
                continue

            # All differences in relative and absolute length
            absdiffs = [abs(lpm[comb[0]][0] - lpm[comb[1]][0]) for comb in all_combs]
            reldiffs = [abs(lpm[comb[0]][2] - lpm[comb[1]][2]) for comb in all_combs]

            # If they have a relative length within 10% or 100nt variation
            # -> 'same length'
            if max(reldiffs) < 0.1 or max(absdiffs) < 100:
                for_pie['Same length'][utr_id] = lpm
                continue
            # If not, check if whole cell supports the longest.
            else:
                # Finally check if some of them have too low RPKM to be caught by
                # the checks so far
                if min(noNArpkms) < rpkm_min:
                    for_pie['Lowly expressed'][utr_id] = lpm
                    continue

                # The distance from WC to the smallest length must be 80% or more
                # of the distance from the smallest to the largest length
                # -> 'Different length'
                wc_compmin_reldist = abs((WCepsrel - min(Cepsrel, Nepsrel))) # eg 0.05
                comp_reldist = abs(Cepsrel - Nepsrel) # eg 0.4
                # hey, did I lose my favorite case???
                if wc_compmin_reldist > 0.8*comp_reldist:
                    for_pie['Different length'][utr_id] = lpm
                    continue

                # WC could be low by accident of PCR. If the other two are 200 nt or
                # more apart and have high enough RPKM, keep them. Choose rpkm=2
                # for more hits ...
                comp_absdist = abs(Clen - Nlen) #absolute distance
                comp_minrpkms = min(Crpkm, Nrpkm) # minimum rpkm
                if (comp_absdist > 200 or comp_reldist > 0.2) and comp_minrpkms < 3:
                    for_pie['Different length'][utr_id] = lpm
                    continue

                # If not close enough, you cannot determine
                # -> 'Undetermined'
                else:
                    for_pie['Undetermined'][utr_id] = lpm
                    continue
                    # eg 0.08, the upper limit for wc

            # If you make it here, I'd be very interested in how to classify you! :)
            debug()

        # Add to super-pie for across-cell_line comparison
        super_pie[cell_line] = for_pie

        # only send the lp part of for_pie to the plotter
        for_pie_plot = {}
        for (category, lpm_dict) in for_pie.items():
            for_pie_plot[category] = lpm_dict.values()

        # plot each pie
        #p.utr_length(for_pie_plot, cell_line, index_dict)

    # Quick -- what is the percent of overlap between the sections, including
    # 'expressed'?
    #classification_overlap(super_pie) TODO

    # If you have more than 1 cell line, start comparing them
    #compare_pies(super_pie) TODO

    debug()

def wc_is_n_plus_c(dsets):
    """
    Your assumption is that whole cell is n + c. This means that the length of
    whole cell should often be between n or c. Is this true?
    """
    lengths = {}
    for (cell_line, comp_dict) in dsets.items():
        lengths[cell_line] = {}
        for (compartment, utrs) in comp_dict.items():
            for (utr_id, utr) in utrs.iteritems():

                # Skip the small-rpkm utrs
                if utr.RPKM < 1:
                    continue

                # Get the absolute 3UTR length determined by epsilon
                if utr.eps_rel_size != 'NA':
                    eps_length = utr.eps_rel_size*utr.length
                else:
                    eps_length = 'NA'

                # Add compartment name, length, and rpkm information
                if utr_id not in lengths[cell_line]:
                    lengths[cell_line][utr_id] = {compartment: eps_length}
                else:
                    lengths[cell_line][utr_id][compartment] = eps_length

    # Go through all the compartments. Count if whole cell length is lower,
    # bigger, or in the middle of the two compartments.

    #countdict = {}
    #abs_dist_dict = {} # also get the absolute differences for the different cases

    #for (cell_line, utrs) in lengths.items():

        #countdict[cell_line] = {'low':0, 'mid':0, 'big':0}
        #abs_dist_dict[cell_line] = {'low':[], 'mid':[], 'big':[]}

        #for (utr_id, c_d) in utrs.iteritems():
            ##Skip those that don't have all compartments
            #if len(c_d) < 3:
                #continue

            #if c_d['Whole_Cell'] < min(c_d['Nucleus'], c_d['Cytoplasm']):
                #countdict[cell_line]['low'] += 1

                #abs_dist = c_d['Whole_Cell'] - min(c_d['Nucleus'], c_d['Cytoplasm'])
                #abs_dist_dict[cell_line]['low'].append(abs_dist)

            #elif c_d['Whole_Cell'] > max(c_d['Nucleus'], c_d['Cytoplasm']):
                #countdict[cell_line]['big'] += 1

                #abs_dist = c_d['Whole_Cell'] - max(c_d['Nucleus'], c_d['Cytoplasm'])
                #abs_dist_dict[cell_line]['big'].append(abs_dist)

            #else:
                #countdict[cell_line]['mid'] += 1

                #mid = (2*c_d['Whole_Cell']-c_d['Nucleus']-c_d['Cytoplasm'])/2
                #abs_dist_dict[cell_line]['mid'].append(mid)

    #avrg_dist_dict = {}
    #for (cl, ddict) in abs_dist_dict.items():
        #avrg_dist_dict[cl] = dict((n, sum(o)/len(o)) for (n, o) in
                                  #ddict.iteritems())

    # NOTE you stopped to investigate extensions. What purpose is there to
    # measure this when almost all of them have abs_length 1 or 0.999?

    #ALSO! what are the distances for the min, max, and medium lengths??
    minmax = {} # ALSO! get min max distances
    for (cell_line, utrs) in lengths.items():
        minmax[cell_line] = []
        for (utr_id, c_d) in utrs.iteritems():
            #Skip those that don't have all compartments
            if len(c_d) < 3:
                continue

            #minmax[cell_line].append(max(c_d.values())-min(c_d.values()))
            minmax[cell_line].append(abs((c_d['Nucleus']-c_d['Cytoplasm'])))

    for (cl, mm) in minmax.items():
        print(cl)
        print('Mean min-max distance: {0} +/- {1}'.format(np.mean(mm),
                                                         np.std(mm)))
        fig, ax = plt.subplots()
        ax.boxplot(mm)
        ax.set_title('Max-min distance for cytosol and nucleus\n{0}'.format(cl))

    # RESULT Mean min-max distance: 29.9979385274 +/- 66.2790973153
    # For RPKM > 1. This lends credance to not accepting 3UTRs as of different
    # length when the distance is less than 100 nt :)
    # RESULT: min-max for CYTOPLASM AND NUCLEUS only are 22. +/- 58
    debug()

def classification_overlap(super_pie):
    """
    3UTRs are classified into 'Expressed', 'Different length', etc. What is the
    degree of overlap between each of these classifications between the cell
    lines?
    """
    debug()

def lenbias(diff_lens, super_pie):
    """
    Check if the cell lines with different lengths between cytosol and nucleus
    have any bias in which is longer.
    """
    bias = {}
    for (cell_line, utrs) in diff_lens.items():
        cyto_longer = 0
        nucl_longer = 0

        cyto_rpkm_bigger = 0
        nucl_rpkm_bigger = 0

        difflen = super_pie[cell_line]['Different length']
        total_diff = len(difflen)
        for utr in utrs:
            if difflen[utr]['Cytoplasm'][0] > difflen[utr]['Nucleus'][0]:
                cyto_longer +=1
            if difflen[utr]['Cytoplasm'][0] < difflen[utr]['Nucleus'][0]:
                nucl_longer +=1

            if difflen[utr]['Cytoplasm'][1] > difflen[utr]['Nucleus'][1]:
                cyto_rpkm_bigger +=1
            if difflen[utr]['Cytoplasm'][1] < difflen[utr]['Nucleus'][1]:
                nucl_rpkm_bigger +=1

        cyfrac = cyto_longer/total_diff
        nufrac = nucl_longer/total_diff

        cyrpkmfrac = cyto_rpkm_bigger/total_diff
        nurpkmfrac = nucl_rpkm_bigger/total_diff

        bias[cell_line] = {'Cytosol-longest': cyfrac,
                           'Nucleus-longest': nufrac,
                           'Cytosol RPKM largest': cyrpkmfrac,
                           'Nucleus RPKM largest': nurpkmfrac,
                          }

        print('{0}\tcytosol rpkm larger: {1}\n\tnucleus rpkm larger: {2}\n'\
              .format(cell_line, cyrpkmfrac, nurpkmfrac, total_diff))
        print('{0}\tcytosol longer: {1}\n\tnucleus longer: {2}\ntotal:\t{3}\n'.\
              format(cell_line, cyfrac, nufrac, total_diff))

    return bias

def compare_pies(super_pie):
    """
    Look for 3UTRs that have different length in all the cell compartments
    """

    # Dictionary that holds the 3UTRs of different length for each cell line
    diff_lens = {}
    for (cell_line, cl_pie) in super_pie.items():
        diff_lens[cell_line] = set(cl_pie['Different length'].keys())

    # Check the diff_len utrs for bias in the length-difference (is nucleus
    # always shorter than cytosol, for example).
    bias = lenbias(diff_lens, super_pie)

    # RESULT there is a clear bias in which compartment is longer between the
    # two cell lines. What is good, however, is that the RPKM distribution is
    # the same for both datasets. Conclusion, difference in RPKM distribution
    # between nucl and cyt don't explain the bias in which compartment is
    # longest.

    # Method: make sets from the UTRs in each pie_dicts, and simply intersect
    # the dicts. See what comes out.
    # Get all intersection numbers of 2 or more couples. Also get the number of
    # intersections that would be expected from random.

    # Do it the monte carlo way. Sample randomly from the whole pool of 3UTRs
    # the number that is in each cell line, and do the same intersection 1000
    # times, and take a mean and std.

    # Get all combinations of the compartments and whole-cell compartments
    all_combs = []
    for i in range(2, len(super_pie.keys()) + 1):
        all_combs += list(combins(super_pie.keys(), i))

    # Go through all combinations and get the intersection of UTRs
    isect_nrs = {} # the number of intersecting utsr
    isect_utrs = {} # the intersecting 3UTRs themselves
    for comb in all_combs:
        # list of sets of diff-length utrs from the combination
        comb_sampl = [set(diff_lens[co]) for co in comb]
        # for each combination, add the number of intersecting utrs
        isect_nrs['+'.join(comb)] = len(set.intersection(*comb_sampl))
        isect_utrs['+'.join(comb)] = set.intersection(*comb_sampl)

    # See how many intersection you would expect between the cell lines by
    # random sampling of the commonly expressed 3UTRs
    random_isects = random_intersections(super_pie, diff_lens, all_combs)

    # RESULT YES! They have 10-20 times more UTRS in common than expected, and
    # further, all 3 cell lines have some of these in common!
    # The final test is if they vary in the same direction...!

    p = Plotter()

    # TODO make a table of the distribution of same length, different length,
    # together with the overall bias for the compartments. Dashing on the
    # RPKM-bias between the compartments.

    #p.compartment_congruence_table(isect_utrs, random_isects, super_pie, bias)
    #p.compartment_congruence_barplot(isect_utrs, super_pie)

    # If you intersect 3 or more cell lines, specifically plot those that
    # intersect for all cell lines. Plot WC, N, and C, including RPKM somewhere,
    # for all of them!
    p.isect_utrs_lengths(isect_utrs, super_pie)

    debug()

def random_intersections(super_pie, diff_lens, all_combs):
    """
    Return the number of intersections expected by random. Chance from all 3UTRs
    or chance from all expressed 3UTRs?

    Result: by comparing with all 3UTRs the average is very low, just 1.8. Maybe
    you have to compare with all expressed ones instead to get a realistic
    picture?
    """
    import random

    counts = dict((cl, len(utr_ids)) for (cl, utr_ids) in diff_lens.items())

    # Those that can be considered expressed
    expressed = ['Undetermined', 'Same length', 'Different length']

    # Now you're getting all the expressed from one of them. You should only get
    # those that are expressed in common.

    # Get the expressed utrs in all cell lines. Then intersect them so you only
    # sample from the utrs that are commonly expressed
    expressed_utrs = {}
    for (cell_line, class_dict) in super_pie.items():
        cl_expressed_utrs = set()

        for (gog, lendict) in class_dict.items():
            if gog in expressed:
                cl_expressed_utrs.update(set(lendict.keys()))

        expressed_utrs[cell_line] = cl_expressed_utrs

    all_utrs = set.intersection(*expressed_utrs.values())

    random_nr = 1000

    r_isects = dict(('+'.join(comb), []) for comb in all_combs)

    for i in range(random_nr):

        #draw as many as 'counts' utrs for each cell line and intersect them
        r_sampl = dict((cl, random.sample(all_utrs, nr)) for (cl, nr) in
                    counts.items())

        for comb in all_combs:
            comb_sampl = [set(r_sampl[co]) for co in comb]
            # for each combination, add the number of intersecting utrs
            r_isects['+'.join(comb)].append(len(set.intersection(*comb_sampl)))

    # Return a dictionary of the mean and std of the number of intersections
    # between the datasets resulting from the random sampling
    ret_dict = {}
    for combo, isects in r_isects.items():
        ret_dict[combo] = {'mean': np.mean(isects), 'std': np.std(isects)}

    return ret_dict

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
    #simple_length_boxplot(dsets, minRPKM=2) # BUSTED

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
    titles = {'eps_rel_size': 'Relative length by read coverage', 'eps_abs_size':
              'Absolute epsilon length', 'length': 'Length of 3UTR in annotation',
              'eps_remainder': 'Distance from epsilon_end to annotated end'}
    # set the order in which you want the plots
    #order = ['length', 'eps_abs_size', 'eps_remainder', 'eps_rel_size']
    order = ['length',  'eps_rel_size']
    # The output variables

    # Go through all dsetsand get the respective output variables for the
    # different RPKM clusters 
    for (cell_line, compartment_dict) in dsets.items():
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
    for (cell_line, compartment_dict) in dsets.items():

        for (compartment, utrs) in compartment_dict.items():

            # first second thrid and fourth
            # nr of clusters
            clrs_nr = 3
            clus_list = [{'ud_ratio':[], 'support':[]} for val in range(clrs_nr)]

            dset_name = cell_line +' '+ compartment

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
            p.last_three_clustersites(clus_list, dset_name)
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
    annotation, etc).

    Especially, what support do the annotated, alone, have?
    """

    p = Plotter()

    #for dset in dsets:
    for (cell_line, compartment_dict) in dsets.items():

        for (compartment, utrs) in compartment_dict.items():

            all_read_counter = {} # cluster_size : read_count for all clusters
            annot_read_counter = {} # cluster_size : read_count for clusters w/annot
            other_dsets= {} # cluster_size : read_count for clusters in other dsets

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

                    # Annotated clusters
                    if cls.annotated_polyA_distance != 'NA':
                        if cls.nr_support_reads in annot_read_counter:
                            annot_read_counter[cls.nr_support_reads].append(keyi)
                        else:
                            annot_read_counter[cls.nr_support_reads] = [keyi]

                    # Clusters in other datasets TODO broken somehow.
                    if cls.nr_support_reads in other_dsets:
                        all_key = dset_2super[keyi] # the in-between key
                        for (dn, sup_reads) in zip(*super_clusters[all_key]):
                            if dn != dset_name: # don't count yourself!!
                                if sup_reads > 1: # maybe set treshold?
                                    other_dsets[cls.nr_support_reads].append(keyi)
                    else:
                        other_dsets[cls.nr_support_reads] = [keyi]

            cluster_dicts = (all_read_counter, annot_read_counter)
            titles = ('All clusters', 'Annotated_TTS')

            # make two figurses: one in terms of all clusters and on of annotated
            # TTS sites

            in_terms_of = (titles[0], titles[2])
            p.join_clusters(cluster_dicts, titles, in_terms_of, dset_name)

    plt.draw()


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

    # Go through all dsetsand give 
    for (cell_line, compartment_dict) in dsets.items():
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
    p.rec_senitivity(sensitivity, intervals, attributes)

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
    for (cell_line, comp_dict) in dsets.items():
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

def novel_polyAs(dsets, super_clusters, dset_2super):
    """
    For each cell line, what are the novel poly(A) sites?

    Accept 2> reads OR 1 read if supported by PAS or SVM
    # 1) For each cell line, make a utr_novel_PA dict:
        # utr_novelPA[utr] = {super_cluster} = {'max_coverage', 'support in # of
        # compartments', 'in_annotation', 'has_PAS', 'has_SVM'}
    # 2) Go through each compartment
    # 3) Go through all utrs
    # 4) Store the utr's dset_2super-linked cluster, with read_coverage
    # 5) If the dset_2super is already there, update read coverage if higher,
    # and also update that it's found in 2 or more compartments
    """

    p = Plotter()

    # prepare the novel_PA dict. It looks like this
    super_cluster_nr = 0
    novel_PA = {}
    for (cell_line, comp_dict) in dsets.items():
        novel_PA[cell_line] = {}
        for (compartment, utrs) in comp_dict.items():
            novel_PA[cell_line][compartment] = {}

            # make a shortcut to the previous dict
            compcls = novel_PA[cell_line][compartment]

            dset_name = cell_line+' '+compartment
            for (utr_id, utr) in utrs.iteritems():

                # SKip those without clusters
                if utr.clusters == []:
                    continue

                # Make dict if utr not found before
                if utr_id not in compcls:
                    compcls[utr_id] = {}

                for cls in utr.clusters:
                    # Get the various things from the cluster you want to
                    # save
                    if cls.nearby_PAS[0] != 'NA':
                        has_PAS = True
                    else:
                        has_PAS = False

                    max_covrg = cls.nr_support_reads

                    if cls.annotated_polyA_distance != 'NA':
                        has_annotation = True
                    else:
                        has_annotation = False

                    # Get the super ID that links you with the super-cluster
                    super_key = dset_name+cls.chrm+cls.strand+\
                            str(cls.polyA_coordinate)
                    super_id = dset_2super[super_key]

                    # Check if this super_id has already been found; if so,
                    # update it. If not, create it anew.

                    if super_id not in compcls[utr_id]:

                        super_cluster_nr += 1
                        compcls[utr_id][super_id] =\
                                {'max_covrg': max_covrg,
                                 'compartment': 1,
                                 'has_annotation': has_annotation,
                                 'has_PAS': has_PAS,
                                 'RPKM': utr.RPKM}
                    else:
                        # Update the entry (through g as shortcut)
                        # PS normally you shuoldn't arrive here, but some spook
                        # of the re-clustering causes two originally distinct
                        # polyA clusters to cluster as one.
                        g = compcls[utr_id][super_id]
                        if g['max_covrg'] < max_covrg:
                            g['max_covrg'] = max_covrg

                        # True or False = True
                        g['has_annotation'] = g['has_annotation'] or has_annotation
                        g['has_PAS'] = g['has_PAS'] or has_PAS
                        g['RPKM'] = max(utr.RPKM, g['RPKM'])

                        g['compartment'] += 1


    # Include in the count if the poly(A) site is annotated. From this you can
    # also get the distribution of annotated vs novel coverage distributions.


    #minRPKM = 30 # skip UTRs that have lower RPKM than this
    # each UTR entry will
    #polyAcounter = count_polyAs(novel_PA, dsets, minRPKM)
    debug()

    # plot the table you get out
    #p.polyAcounterplot(polyAcounter, minRPKM)

    # RESULT for minRPKM = 50
    #
    #{'GM12878': {'Cytoplasm': 372, 'Nucleus': 321},
     #'HeLa': {'Cytoplasm': 500, 'Nucleus': 322},
     #'K562': {'Cytoplasm': 312, 'Nucleus': 205}}
    # To me this correlates with the sensitivity graphs.

    # for each 3UTR, count the coverage of proximal to distal poly(A) sites.
    # NOTE! This is independent of the coutn-analysis. Even if that doesn't
    # work, this one should be good.
    polydist, annotdist, count, rel_usage = polyAdistribution(novel_PA)

    # Shit. You know, this information would have been better if you could
    # compare with all other dsets. Then you could add 1 if found in other
    # dsetsfor example.

    #p.polydist_plot(polydist, annotdist, count, rel_usage)
    # RESULT the fact is that the 3'-most site is expressed the most. There is
    # no difference between cell compartments.
    # I recomend that you ditch the cases with 4 clustesr. Then make it clearer
    # from the context that they are on average less expressed than the 3'-most
    # poly(A) site, but they are more expressed than the ratio of
    # annotated/non-annotated would hint.

    # What you you need to start doing:
    # 1) Using novel_PA, Print a small table of the novel poly(A) sites you
    # find. First print for all compartments; then merge compartments and print
    # for cell lines; then merge for cell lines and print for whole dataset.
    # 2) Look again at the polydist_plot. You find no differences in the
    # distributions of usage of poly(A) sites. You need to present this. Perhaps
    # you can present just the cases with 3 poly(A) sites and show there is no
    # difference. Also, make it more clear that the 3'-most one is mostly higher
    # expressed.
    # maybe it would be clear with distributions of the boxplots of the counts
    # themselves? Don't add the coutns, but keep the distribution. Plot the
    # distributions using box plots. I think that's the simplest, clearest way
    # of showing the differencel.
    # 3) Go back to the diff-len 3UTRs that Roderic wasn't so impresed by.
    # Characterize better which compartment they lengthen in, and check if any
    # of the cases are backed up by poly(A) evidence.
    # 4) Start working on detecting 3UTRs that are longer than the annotation.
    # Reshape the entire system for rel-length then ... or, keep the rel-length
    # according to annotation, but introduce another variable (why not!?) that
    # can be 1.3 1.5 whatever, and include the stop site for this one. Then you
    # are in a better position to determine where things stop.

    debug()

def polyAdistribution(novel_PA):
    """
    For each compartment, for each gene, get how many polyA clusters the gene
    has. Make a dict for each compartment '1, 2, 3,4' as keys. The key is the
    max nr of "goood" clusters defined by your criteria. Each key contains the
    number of read counts at those internal positions in the 3UTR. 1 is 5' and
    the last value is 3' (If only 1 they are the same). Thus, when you find
    clusters from the - strand, you need to take that into account.
    """
    polydist = AutoVivification()
    annotdist = AutoVivification()
    count = AutoVivification()

    # pedro's idea: relative useage according to 3'-most
    rel_usage = AutoVivification()

    getreverse = {'-': True, '+': False}

    for (cell_line, comp_dict) in novel_PA.items():
        for (compartment, utr_dict) in comp_dict.iteritems():
            for (utr_id, cls_dict) in utr_dict.iteritems():

                trusted = {}
                strand = 0
                # Get the trusted polyA sites from cls_dict
                for chmStrandCoord, cls in cls_dict.items():

                    # if this UTR's RPKM is low, skip the whole thing.
                    if cls['RPKM'] < 10:
                        break

                    # Demapd, pas, annot, or SVM, and covr > 1.
                    if ((cls['has_PAS'] or\
                       cls['has_annotation'] ) or\
                       cls['max_covrg']) > 1:

                        # Extract the strand and coordinate
                        if len(chmStrandCoord.split('-')) == 2:
                            strand = '-'
                            coord = chmStrandCoord.split('-')[1]
                        elif len(chmStrandCoord.split('+')) == 2:
                            strand = '+'
                            coord = chmStrandCoord.split('+')[1]

                        # key by coordinate
                        trusted[int(coord)] = cls
                    else:
                        # go to next utr if one polyA site fails
                        break

                # Skip if no clusters are trusted
                if len(trusted) == 0:
                    continue

                # nr of clusters
                polyA_nr = len(trusted)

                # Create [0, 0, .. ,0] if doesn't exist for this len
                if polyA_nr not in polydist[cell_line][compartment]:
                    polydist[cell_line][compartment][polyA_nr] =\
                            [0 for i in range(polyA_nr)]

                # same for annotdist
                if polyA_nr not in annotdist[cell_line][compartment]:
                    annotdist[cell_line][compartment][polyA_nr] =\
                            [0 for i in range(polyA_nr)]

                # same for count
                if polyA_nr not in count[cell_line][compartment]:
                    count[cell_line][compartment][polyA_nr] = 0

                # Create a list for appending for rel_usage
                if polyA_nr not in rel_usage[cell_line][compartment]:
                    rel_usage[cell_line][compartment][polyA_nr] = []

                # add a count of how many you have
                count[cell_line][compartment][polyA_nr] +=1

                # Go through trusted in the 5->3 order by sort. If strand is -,
                # use reverse = True
                rev = getreverse[strand]

                usage = []
                # Iterate through dict in sorted (5->3) manner
                for (indx, pos) in enumerate(sorted(trusted, reverse=rev)):
                    clstr = trusted[pos]
                    polycov = clstr['max_covrg']
                    polydist[cell_line][compartment][polyA_nr][indx] += polycov

                    # append for 'usage'
                    usage.append(polycov)

                    # add for annotdist if this cluster is in annotation
                    if clstr['has_annotation']:
                        annotdist[cell_line][compartment][polyA_nr][indx] += 1

                # normalize according to the last index in usage
                rel_us = [us/usage[-1] for us in usage]
                rel_usage[cell_line][compartment][polyA_nr].append(rel_us)

    return polydist, annotdist, count, rel_usage

def count_polyAs(novel_PA, dsets, minRPKM):
    """
    For each compartment, count the number of sure polyA sites above a certain
    minimum RPKM.
    """

    polyAcompare = {}
    for (cell_line, comp_dict) in novel_PA.items():
        polyAcompare[cell_line] = {}
        for (compartment, utr_dict) in comp_dict.iteritems():
            for (utr_id, cls_dict) in utr_dict.iteritems():

                #  if utr not present for this cell line, add it
                if utr_id not in polyAcompare[cell_line]:
                    polyAcompare[cell_line][utr_id] = {compartment: cls_dict}

                else:
                    polyAcompare[cell_line][utr_id][compartment] = cls_dict


    polyAcounter = {}
    for (cell_line, utrs) in polyAcompare.items():
        all_comps = dsets[cell_line].keys()
        all_comps.remove('Whole_Cell')

        # initialize dict for all compartments
        polyAcounter[cell_line] = dict((comp, 0) for comp in all_comps)

        for (utr_id, comp_dict) in utrs.iteritems():
            # Skip those that are not present in all compartments
            if not set(comp_dict.keys()).issuperset(set(all_comps)):
                continue

            # Skip those that have low RPKM
            rpkms = []
            for (comp, cls_dict) in comp_dict.items():
                if comp in all_comps:
                    rpkms.append(cls_dict.items()[0][1]['RPKM'])

            if min(rpkms) < minRPKM:
                continue

            for (comp, cls_dict) in comp_dict.items():

                # Skip 'whole cell'
                if comp not in all_comps:
                    continue

                # add the polyA sites in the cls_dict supported by at least one
                # line of evidence:
                for (super_coord, cls) in cls_dict.items():
                    if cls['has_PAS'] or\
                       cls['has_annotation'] or\
                       cls['max_covrg'] > 1:

                        polyAcounter[cell_line][comp] += 1

    return polyAcounter

def novel_polyAs_2(dsets, super_3utr, settings):
    """
    1) How many total annotated poly(A) sites
    wc -l on the annot poly(A) file

    2) How many annotated poly(A) sites in "my" 3UTRs
    wc -l on the same file after intersection with "my" 3UTRs

    3) How many novel poly(A) sites in "my" 3UTRs
    loop through super_3utr, and get those with 2> reads

    4) How many of these have been detected with ESTs/cDNA?
    expand the polyA_db.bed-file and unique-intersect with a bedfile you write
    of all your poly(A) sites
    """

    my_novel_polyAs = []

    print('\n'+'-'*70+'\n')
    ## 1) How many total annotated poly(A) sites (and tian tb sites)
    #cmd1 = ['slopBed', '-b', '200', '-i', settings.annot_polyA_path, '-g',
           #settings.hg19_path]

    #f1 = Popen(cmd1, stdout=PIPE, bufsize=-1)

    #cmd2 = ['mergeBed', '-s', '-i', 'stdin'] # read from stdin

    #f2 = Popen(cmd2, stdin=f1.stdout, stdout=PIPE, bufsize=-1)

    all_annot_polyA = sum(1 for i in open(settings.annot_polyA_path))
    tian_db = sum(1 for i in open(settings.polyA_DB_path))

    print('Number of GENCODE annotated poly(A) sites: '\
          '{0}\n'.format(all_annot_polyA))
    print('Number of poly(A) sites in polyA_db: '\
          '{0}\n'.format(tian_db))

    # 2) How many annotated poly(A) sites in "my" 3UTRs. How many of these do I
    # recover? How many new do I find?
    cmd = ['intersectBed', '-u', '-a', settings.annot_polyA_path, '-b',
           settings.utr_exons_path]

    f = Popen(cmd, stdout=PIPE)
    # The number of annotated poly(A) sites in YOUR 3utrs
    my_annot = sum(1 for i in f.stdout)
    print('Number of GENCODE annotated poly(A) sites in "our" 3UTRS: '\
          '{0}\n'.format(my_annot))

    # 3)
    # The number of annotated pol(A) sites that I recover
    my_recover = 0
    # The number of novel poly(A) sites in 'my' 3UTRs
    my_novel = 0
    for utr_id, utr in super_3utr.iteritems():
        for super_polyA, cell_dict in utr.super_cover.items():
            has_annot = False
            enough_covrg = False
            # Trust if you find pA in more than 1 cell line
            if len(cell_dict) > 1:
                enough_covrg = True
            for c_line, comp_dict in cell_dict.items():
                # Trust if you find pA in more than 1 compartment
                if len(comp_dict) > 1:
                    enough_covrg = True
                for comp, polyA in comp_dict.items():
                    if polyA.annotated_polyA_distance != 'NA':
                        has_annot = True
                    #  Trust if you find pA with more than 1 read
                    if polyA.nr_support_reads > 1:
                        enough_covrg = True

            if has_annot:
                my_recover += 1
            if enough_covrg:
                my_novel += 1

                # Save the novel polyA sites to a list
                if len(super_polyA.split('+')) == 2:
                    my_novel_polyAs.append(tuple(super_polyA.split('+') + ['+']))
                else:
                    my_novel_polyAs.append(tuple(super_polyA.split('-') + ['-']))

    print('Number of "our" poly(A) sites that are also in the GENCODe annotation: '\
          '{0}\n'.format(my_recover))
    print('Number of "our" poly(A) sites that are new to GENCODE: '\
          '{0}\n'.format(my_novel))

    #4) How many of these have been detected with ESTs/cDNA?
    #expand the polyA_db.bed-file and unique-intersect with a bedfile you write
    #of all your poly(A) sites
    # i) write these novel polyAs to bedfile
    out_dir = os.path.dirname(settings.annot_polyA_path)
    temp_path = os.path.join(out_dir, 'novel_polyA.bed')
    temp_handle = open(temp_path, 'wb')

    for nov in my_novel_polyAs:
        out = '\t'.join([nov[0], nov[1], str(int(nov[1])+1), nov[2]]) + '\n'
        temp_handle.write(out)

    temp_handle.close()

    # ii) Expand all entries in annotated polyAs with 2 nt, merge these entries,
    # and feed this into an overlapper with novel_polyA.bed; count the number of
    # overlaps
    cmd1 = ['slopBed', '-b', '20', '-i', settings.polyA_DB_path, '-g',
           settings.hg19_path]

    f1 = Popen(cmd1, stdout=PIPE, bufsize=-1)

    cmd2 = ['mergeBed', '-s', '-i', 'stdin'] # read from stdin

    f2 = Popen(cmd2, stdin=f1.stdout, stdout=PIPE, bufsize=-1)

    cmd3 = ['intersectBed', '-a', temp_path, '-b', 'stdin']

    f3 = Popen(cmd3, stdin=f2.stdout, stdout=PIPE, bufsize=-1)

    # the number of my novel poly(A) sites that are in Tian's database
    in_tian = sum(1 for i in f3.stdout)

    print('Number of "our" novel poly(A) sites that are also in polyA_db: '\
          '{0}\n'.format(in_tian))

    print('-'*70+'\n')


def polyadenylation_comparison(dsets, super_3utr, settings):
    """
    * compare 3UTR polyadenylation in general
    * compare 3UTR polyadenylation UTR-to-UTR
    """

    #### Novel polyA sites in annotated 3UTRS
    #novel_polyAs(dsets, super_clusters, dset_2super)
    # The new version that makes stuff easy
    novel_polyAs_2(dsets, super_3utr, settings)

    #### Compare the compartments; how many annotated do we find? etc.
    #compare_cluster_evidence(dsets, super_clusters, dset_2super)

    #### Get the sensitivity rates for of recovery of poly(A) clusters
    #recovery_sensitivity(dsets)

    #### Get the reproducibility rate of poly(A) clusters for each compartment
    #### compared to whole cell (if it is in on of the compartment, it should
    #### presumably be in whole cell as well, albeit dilluted, that is, coverage
    #### should be lower).
    #wc_compartment_reproducability(dsets, super_clusters, dset_2super)

    #### Correlate the coverage counts of common polyadenylation sites between
    #### clusters
    #correlate_polyA_coverage_counts(dsets, super_clusters) NOT dsetsREADY

    #### Get the change in coverage ratio of the 3 (or so) last polyadenylation
    #### sites from the 3 UTR
    #udstream_coverage_last_clusters(dsets)

    #### See the effects of excluding pA sites with 1, 2, 3, etc, coverage
    #cluster_size_sensitivity(dsets) NOT dsetsREADY


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

    import utail as finder

    finder_settings = finder.Settings\
            (*finder.read_settings(settings.settings_file))

    chr1 = False
    #chr1 = True
    if chr1:
        finder_settings.chr1 = True

    # Get the bedfile with 3UTR exons
    beddir = os.path.join(settings.here, 'source_bedfiles')
    utr_bed = finder.get_utr_path(finder_settings, beddir, False)
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

        # Update the datasets in dsetswith 'other_strand' polyA information
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

    takes into account the strand of the cluster
    """

    all_clusters = {} # keys = cluster_centers
    each_cluster_2all = {} # mapping from the individual clusters to super-cluster

    chrms = ['chr'+str(nr) for nr in range(1,23) + ['X','Y','M']]
    tsdict = dict((chrm, {'+':[], '-':[]}) for chrm in chrms)

    # Put all the polyAclusters in a dict with cluster pos and dset.name
    for (cell_line, compartment_dict) in dsets.items():
        for (compartment, utrs) in compartment_dict.items():
            dset_name = cell_line +' '+ compartment
            for (utr_id, utr) in utrs.iteritems():
                if utr.clusters != []:
                    for cls in utr.clusters:
                        tsdict[cls.chrm][cls.strand].append((cls, dset_name))

    # iterate through each strand for each chromosome and cluster pA-sites
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
                #coverages = [clu.nr_support_reads for (clu, dn) in mega_cluster]
                # NOTE: adding whole object!
                coverages = [clu for (clu, dn) in mega_cluster]
                dsets= [dn for (clu, dn) in mega_cluster]

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

    for (cell_line, compartment_dict) in dsets.items():
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
    Simply plot the rpkms of the different datasets. Each row is a cell line and
    each column is a compartment.
    """
    #thresh = 0
    cl_max = len(dsets)
    comp_max = max([len(comps) for cl, comps in dsets.items()])

    (fig, axes) = plt.subplots(comp_max, cl_max, sharex=True, sharey=True)

    for cl_nr, (cell_line, comp_dict) in enumerate(dsets.items()):
        for comp_nr, (compartment, utrs) in enumerate(comp_dict.items()):

            rpkms_norm = []
            rpkms_log = []
            dset_name = cell_line +' '+ compartment

            for (utr_id, utr) in utrs.iteritems():
                if utr.eps_coord != 'NA':
                    # SKip those with rpkm = 0
                    if utr.RPKM == 0:
                        continue
                    #if utr.RPKM > thresh:
                        #rpkms_norm.append(thresh)
                    #else:
                    rpkms_norm.append(utr.RPKM)
                    rpkms_log.append(math.log(utr.RPKM, 2))

            ax = axes[comp_nr, cl_nr]
            (N, bins, blah) = ax.hist(rpkms_log, bins=200)
            ax.grid()

            if cl_nr == 0:
                ax.set_ylabel('{0}'.format(compartment), size=20)
            if comp_nr == comp_max-1:
                ax.set_xlabel('{0}'.format(cell_line), size=20)

            # Percentages above/below 0
            less_zero = len([1 for v in rpkms_norm if v>1])/float(len(rpkms_norm))
            print('\n{0}: Percentage of UTRs with RPKM larger than 1: {1}%'\
                  .format(dset_name, format(less_zero*100, '.0f')))

    fig.suptitle('log2-transformed distribution of RPKM values', size=20)

    # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
    plt.setp([a.get_xticklabels() for a in axes[0:comp_max-1,:].flatten()],
             visible=False)
    plt.setp([a.get_yticklabels() for a in axes[:,1:].flatten()],
             visible=False)

    # Fine-tune: remove space between subplots
    fig.subplots_adjust(hspace=0.07)
    fig.subplots_adjust(wspace=0.07)

def gen_color():
    """generator for getting n of evenly distributed colors using
    hsv color and golden ratio. It always return same order of colors

    :returns: RGB tuple
    """
    import colorsys
    golden_ratio = 0.618033988749895
    h = 0.22717784590367374

    while 1:
        h += golden_ratio
        h %= 1
        HSV_tuple = [h, 0.95, 0.95]  # this defines how "deep" are the colors
        RGB_tuple = colorsys.hsv_to_rgb(*HSV_tuple)
        yield map(lambda x:str(int(x * 256)), RGB_tuple)

def beyond_aTTS(dsets):
    """
    How much do 3UTRs extend beyong the aTTS?
    """
    # RESULT it seems that 95% are not extended beyond the annotation. For those
    # that remain, a few might overlap some other annotatated area. A few might
    # be genuinly longer. How to get them? They would slightly aid the diff-len
    # compartment data.

    for (cell_line, comp_dict) in dsets.items():
        for (compartment, utrs) in comp_dict.items():
            biggerthan = []
            not_bigger = 0
            for (utr_id, utr) in utrs.iteritems():
                if utr.epsilon_beyond_aTTS == 0:
                    not_bigger += 1
                else:
                    biggerthan.append(utr.epsilon_beyond_aTTS)

        fig, ax = plt.subplots()
        ax.boxplot(biggerthan)
        print('{0}\nNot longer than annotation: {1}\n'.format(cell_line,
                                                              not_bigger/len(utrs)))

    debug()

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def visualizedict(dictionary):
    """
    Print out the structire of a dictionary
    """
    maxdepth = 5
    print dictionary.keys()[:maxdepth]
    for nr1, (keys, subs) in enumerate(dictionary.iteritems()):
        if nr1 > maxdepth:
            continue
        if type(subs) is not dict:
            continue
        print subs.keys()[:maxdepth]
        for nr2, (k, su) in enumerate(subs.iteritems()):
            if nr2 > maxdepth:
                continue
            if type(su) is not dict:
                continue
            print su.keys()[:maxdepth]
            for nr3, (k3, su3) in enumerate(su.iteritems()):
                if nr3 > maxdepth:
                    continue
                if type(su3) is not dict:
                    continue
                print su3.keys()[:maxdepth]
                for nr4, (k4, su4) in enumerate(su3.iteritems()):
                    if nr4 > maxdepth:
                        continue
                    if type(su4) is not dict:
                        continue
                    print su4.keys()[:maxdepth]

def pipeline_utrpath(settings, chr1=False):
    """
    Return the path of the 3UTR bed file that was used for simulating the
pipeline
    """

    import utail as finder

    finder_settings = finder.Settings\
            (*finder.read_settings(settings.settings_file))

    if chr1:
        finder_settings.chr1 = True

    beddir = os.path.join(settings.here, 'source_bedfiles')

    return finder.get_utr_path(finder_settings, beddir, False)

def nucsec(dsets, super3UTR, settings, here):
    """ Get the sequences of all the 3UTR objects (a few hundred megabytes of
    memory?). Then intersect the 3UTR regions with the nucleosome bedGraph
    files. TODO your chromatin runs are done. Read them (fast) and start reading
    sequences and other stuffs.

    IDEA: for each extra dataset, how many new (high quality) poly(A) sites are
    found? Can we draw a nice curve to saturation? :) To do this properly, you
    should add each 'lane' file for both biological replictates of one of the
    cell lines. Maybe even add cytosolic samples too. This can be done easily!
    Just define datasets: 1lane, 2lanes, 3lanes, 4lanes etc! Then also output
    the number of reads in each (alternatively you get get this later with wc -l
    on the raw reads file in the temp-dirs.)

    Another idea: a measure for sensitivity: how many of the datasets have a
    1-read poly(A) site that is supported by >2 in other datasets?

    """
    # 1) Get the path of the 3UTR-bedfile used in the pipeline
    # 2) Get the path of the nucleosomal reads
    nuclpaths = [os.path.join(here, 'nucleosomes', dset,
                 'nucl_'+dset+'_sample.bed') for dset in dsets]

    utr_path = pipeline_utrpath(settings, chr1=True)

    # 3) run bed-intersect on these two files

    # HEY! you are not getting the coverages 'out' of the file. The coverages
    # are in the 'val' lane -- this is lost. I think you want intersectBed!
    # Either that, or you need  to run coverageBed on the ORIGINAL raw reads.
    for nupath in nuclpaths:
        cmd = ['intersectBed', '-wa', '-wb', '-a', nupath, '-b', utr_path]

        p = Popen(cmd, stdout=PIPE)
        # 4) Loop through the output and create a nucleosome coverage 
        for line in p.stdout:
            (chrmP, beP, enP, covr, chrmU, beU, enU, utr_id, val,
             strand) = line.split()
            (tot_exons, extendby) = val.split('+')

            # For now, ignore the 3UTRs with more than 1 exon
            if int(tot_exons) != 1:
                continue

            # OK, there's something you need to know about UTR_ID
            # UTR_ID = emsembID + _UTRnr + _EXONnr in this utr
            local_utrID = '_'.join(utr_id.split('_')[:-1])

            # DEBUGGING XXX
            if local_utrID not in super3UTR:
                continue
            # XXX DEBUGGING
            # XXX what does it mean that an utr_ID is not in super3UTR? maybe
            # because you only get the first 1000 lines from each file!

            # Get the slice of intersection and fill this slice in the coverage
            # vector
            vec_beg = max(int(beP), int(beU)) - int(beU)
            vec_end = min(int(enP), int(enU)) - int(beU)

            # Check if 
            super3UTR[local_utrID].has_nucl = True
            super3UTR[local_utrID].nucl_covr[vec_beg:vec_end] = int(covr)

    # OK, now you are covering the 3UTRs with their nucleosomes. Now you just
    # need to do some correlation.

    # How to do that exactly? First: plot. Compute the average of all nucl_covr
    # vector +/- 400 nt around each poly(A) site

    # 1) Go through each poly(A) site. If there is 400 on 'both sides' of it,
    # and if has_nucl-coverage, and if there is at least 500 nt to the next
    # poly(A) site, add this coverage vector to a coverage matrix
    coverages = []
    for utr_id, utr in super3UTR.iteritems():
        # Skip those without nucleosome coverage
        if not utr.has_nucl:
            continue
        # Skip those that don't have any poly(A) clusters
        if utr.super_cover == {}:
            continue
        debug()

    # XXX Do the below with 1000-extended 3UTRs using biological replicates.

    # TODO the below! Now you are turning to making "Roderic's table ... "
    # NOTE this is just basics. What you really should do is go through
    # super_3UTR and classify poly(A) sites in different ways. Then you simply
    # write out as many bedfiles as you have classifications with +/- 1000 nt
    # extensions. Then you intersect those extensions with the genome-covered
    # file. Then you make graphs based on those extensions. That's more or less
    # what Hagen did. When you do the above, also get controls: 1) non-expressed
    # transcripts' poly(A) sites, 2) random AATAAA sites in 5'utr and cds
    # (intergenic regions in molecular Cell paper).

    debug()

def get_region_sizes():
    region_dir =\
    '/users/rg/jskancke/phdproject/3UTR/the_project/genomic_regions_hg19'

    # don't do the same job again ... only change if your ananlyssi changes!
    sizes = {'3UTR-exonic': 29901395,
             '3UTR-intronic': 11911101,
             '5UTR-exonic': 8176282,
             '5UTR-intronic': 146019292,
             'CDS-exonic': 35392076,
             'CDS-intronic': 976171790,
             'Intergenic': 1607232798,
             'Nocoding-exonic': 32682004,
             'Noncoding-intronic': 247883855}

    return sizes


    cmd = 'ls '+ os.path.join(region_dir)
    p = Popen(cmd, stdout=PIPE, shell=True)
    for filename in p.stdout:
        filename = filename.rstrip()

        # skip the 'chr1' directory
        if filename.startswith('chr'):
            continue

        reg = filename.split('_')[0]

        sizes[reg] = 0

        for line in open(os.path.join(region_dir, filename), 'rb'):
            (chrm, beg, end, d,d, strand) = line.split()
            sizes[reg] += int(end) - int(beg)

    return sizes

def polyA_summary(dsets, super_3utr, polyAstats, settings, side):
    """ Show that the polyA minus are fishy

    1) Number of total polyA sites
    2) Number and % of polyA sites in annotation
    3) Number and % of polyA sites with PAS
    4) Total number of poly(A) reads

    Now that you're getting lots of extra information, what about splitting this
    up into

    i) just the poly(A) reads between compartments (an anlysis of
    fitness of the experiment) and

    ii) a comparason between the strands for each compartment, cl, etc

    NEW: you now have sites on both strand of the annotated one. I think you
    should give a report for both.

    """

    # Get the sizes of the different genomic regions
    region_sizes = get_region_sizes()
    tot_size = sum(region_sizes.values())


    # make a [region][cl][compartment][polyA+/-][replica] structure; each region
    # will have its own plots. cl will be K562 only. Then compare compartments
    # first in terms of replicas, then in terms of polyA+/-
    everything = AutoVivification()
    unique_sites = AutoVivification()

    for region, reg_dict in dsets.items():

        print('-'*60+'\n\nResults for region: {0} {1}\n\n'.format(region, side)\
              +'-'*60+'\n')

        for cl, cl_dict in reg_dict.items():

            minus = AutoVivification()

            for comp, comp_dict in cl_dict.items():

                minus[comp] = {'Cls nr':0, 'Annot pA':0, 'PAS pA':0, 'read_nr':0,
                               'both_an_and_PAS':0, 'Cls min2':0}

                for utr_id, utr in comp_dict.iteritems():
                    if utr.clusters != []:
                        for cls in utr.clusters:
                            # skip according to which side you want to asess
                            if side == 'annotated' and cls.strand != utr.strand:
                                continue

                            if side == 'opposite' and cls.strand == utr.strand:
                                continue

                            minus[comp]['Cls nr'] += 1
                            minus[comp]['read_nr'] += cls.nr_support_reads

                            # if more than 1 read supporting
                            if cls.nr_support_reads > 1 or\
                               cls.annotated_polyA_distance != 'NA':

                                minus[comp]['Cls min2'] += 1

                            if cls.annotated_polyA_distance != 'NA':
                                minus[comp]['Annot pA'] += 1

                            if cls.nearby_PAS != ['NA']:
                                minus[comp]['PAS pA'] += 1

                            if (cls.annotated_polyA_distance != 'NA') \
                               and (cls.nearby_PAS != ['NA']):
                                minus[comp]['both_an_and_PAS'] += 1

            for (comp, count_dict) in minus.items():
                total = count_dict['Cls nr']
                total2 = count_dict['Cls min2']

                total2_frac = format(total2/total, '.2f')

                annot = count_dict['Annot pA']
                annot_frac = format(annot/total, '.2f')

                pas = count_dict['PAS pA']
                pas_frac = format(pas/total, '.2f')

                read_nr = count_dict['read_nr']
                read_per_site = format(read_nr/total, '.2f')

                anYpas = count_dict['both_an_and_PAS']

                if annot == 0:
                    anYpas_rate = 'NA'
                else:
                    anYpas_rate = format(anYpas/annot, '.2f')

                # PLOT THIS STATISTIC WHEN RUNNING ONLY FOR ANNOTATED 3UTR REGIONS
                # YOU KNOW 3UTRs!
                oppos_readnr = polyAstats[region][cl][comp]['Other strand count']
                this_readnr = polyAstats[region][cl][comp]['Annotated strand count']
                both_readnr = polyAstats[region][cl][comp]['Both strands count']

                print cl
                print comp
                print("Total poly(A) sites: {0}".format(total))
                print("Total poly(A) sites with 2+ reads or annotated: {0} ({1})"\
                      .format(total2, total2_frac))
                print("Total poly(A) reads (per site): {0} ({1})"\
                      .format(read_nr, read_per_site))
                print("Annotated poly(A) sites: {0} ({1})"\
                      .format(annot, annot_frac))
                print("poly(A) sites with PAS: {0} ({1})".format(pas, pas_frac))
                print("Annotated poly(A) sites with PAS: {0} ({1})"\
                      .format(anYpas, anYpas_rate))
                print("")

                # get poly(A) + or - and get if it is a replicate or not it's in
                # comp
                realcomp = comp # successively remove 'replicate' or 'minus'

                if 'Replicate' in comp:
                    replicate = 'replicate'
                    realcomp = comp.partition('Replicate')[0]
                else:
                    replicate = 'not_replicate'

                if 'Minus' in comp:
                    polyMinus = 'PolyA-'
                    realcomp = realcomp.partition('Minus')[0]
                else:
                    polyMinus = 'PolyA+'

                mdict = {}

                mdict['total_sites'] = total
                mdict['total_sites2'] = total2 # sits with 2+ reads or annotated
                mdict['total_reads'] = read_nr

                # ad hoc solution: when you run 3UTR only (from the careful stuff)
                # How many reads normalized by region size?
                if region == '3UTR':
                    mdict['total_reads_region_normalized'] =\
                    read_nr*(1-(region_sizes['3UTR-exonic']/tot_size))
                    mdict['total_sites_normalized'] =\
                    total*(1-(region_sizes['3UTR-exonic']/tot_size))
                    mdict['sites_with_pas_fraction'] = pas_frac
                else:
                    mdict['total_reads_region_normalized'] =\
                    read_nr*(1-(region_sizes[region]/tot_size))
                    mdict['total_sites_normalized'] =\
                    total*(1-(region_sizes[region]/tot_size))
                    mdict['sites_with_pas_fraction'] = pas_frac

                # finally, add all to this:
                everything[region][cl][realcomp][replicate][polyMinus] = mdict

        total_unique = 0
        for utr_id, utr in super_3utr[region].iteritems():
            total_unique += len(utr.super_cover)

        unique_sites[region]['Total_poly(A)_sites'] = total_unique
        if region == '3UTR':
            unique_sites[region]['Total_normalized_poly(A)_sites'] =\
                    total_unique/region_sizes['3UTR-exonic']
        else:
            unique_sites[region]['Total_normalized_poly(A)_sites'] =\
                    total_unique/region_sizes[region]

        print("Total unique poly(A) sites for {1}: {0}".format(total_unique,
                                                               region))

    return
    p = Plotter()
    #For each region, for each cell type, compute the read counts etc as in
    #'everything' for poly(A) minus and for replicates. Make one multi-plot per
    #region.
    #p.region_figure(True, False, dsets, everything, settings, 'all')

    #debug()
    for plot_clusters in [True, False]:
        #for normlized in [True, False]:
        normlized = False

        if plot_clusters:
            for clr_type in ['all', '2+']:
                p.region_figure(plot_clusters, normlized, dsets, everything,
                                settings, clr_type)
        else:
            p.region_figure(plot_clusters, normlized, dsets, everything,
                            settings, 'all')

def get_polyA_stats(settings):
    """ Return a dictionary that for each [region][cell_line][compartment]
    contains some output statistics for the poly(A) reads for that dataset.
    """

    polyA_stats = AutoVivification()
    regions = settings.regions

    for region in regions:

        regionfiles = settings.polyAstats_files(region)

        # Check if all length files exist or that you have access
        [verify_access(f) for f in regionfiles.values()]

        for dset_name in settings.datasets:

            statsfile = open(regionfiles[dset_name], 'rb')

            # Create the utr objects from the length file
            stat_dict = {}

            for (linenr, line) in enumerate(statsfile):
                (key, value) = line.split('\t')

                stat_dict[key] = value.rstrip()


            # Add the utr_dict for this cellLine-compartment 
            this_cl = dset_name.split('_')[0]
            this_comp = '_'.join(dset_name.split('_')[1:])

            polyA_stats[region][this_cl][this_comp] = stat_dict

    return polyA_stats

def get_dsetreads(settings, region):
    """ For each dataset, get the number of total reads. The region doesn't
    matter, because the number of reads are dataset-specific.
    """

    dsetreads = {}
    polyA_files = settings.polyAstats_files(region)
    for dset, dsetpath in polyA_files.items():

        filedict = dict((line.split('\t')[0], line.split('\t')[1])
                        for line in open(dsetpath, 'rb'))

        dsetreads[dset] = int(filedict['Total number of reads'].rstrip())

    return dsetreads

def avrg_tail(new_tail, sum_tail):
    """ Add new tail to sum tail. Return sum tail.
    """
    # Add all nucleotides to the sum
    # is this an a-tail or a t-tail?
    nuc_dict = dict([b.split('=') for b in new_tail.split(':')])

    def sortbyme(val):
        return float(val[1])
    tail_type = sorted(nuc_dict.items(), key=sortbyme, reverse=True)[0][0]

    sum_tail[tail_type][0] += 1 # count each one
    sum_tail[tail_type][1] += float(nuc_dict['G'])
    sum_tail[tail_type][2] += float(nuc_dict['A'])
    sum_tail[tail_type][3] += float(nuc_dict['T'])
    sum_tail[tail_type][4] += float(nuc_dict['C'])

    return sum_tail

def data_scooper(cls, keyw, this_dict):
    """ Get data from the merged cls clusters
    """

    # Count all clusters
    this_dict['All']['info_dict'][keyw] += 1

    # Count tails
    taildict = this_dict['All']['tail_lens'][keyw]
    taildict = avrg_tail(cls.tail_info, taildict)

    if cls.PAS_distance[0] != 'NA':
        this_dict['wPAS']['info_dict'][keyw] += 1

        taildict = this_dict['wPAS']['tail_lens'][keyw]
        taildict = avrg_tail(cls.tail_info, taildict)

        if 'AATAAA' in cls.nearby_PAS or 'ATTAAA' in cls.nearby_PAS:
            this_dict['goodPAS']['info_dict'][keyw] += 1

            taildict = this_dict['goodPAS']['tail_lens'][keyw]
            taildict = avrg_tail(cls.tail_info, taildict)

        if 'AATAAA' in cls.nearby_PAS:
            this_dict['bestPAS']['info_dict'][keyw] += 1

            taildict = this_dict['bestPAS']['tail_lens'][keyw]
            taildict = avrg_tail(cls.tail_info, taildict)

    if cls.annotated_polyA_distance != 'NA':
        this_dict['annotated']['info_dict'][keyw] += 1

        taildict = this_dict['annotated']['tail_lens'][keyw]
        taildict = avrg_tail(cls.tail_info, taildict)

        if cls.PAS_distance[0] != 'NA':
            this_dict['annotated_wPAS']['info_dict'][keyw] += 1

            taildict = this_dict['annotated_wPAS']\
                    ['tail_lens'][keyw]
            taildict = avrg_tail(cls.tail_info, taildict)

    return this_dict


def get_dsetclusters(subset, region, settings, speedrun, batch_key):
    """ Get counts for all clusters and 'good' clusters (2 or more reads or
    annotated).
    """

    # count if the below variables are in same or in opposite strand: in the end
    # sum them. This is only valid for those genomic regions where you know the
    # strand.

    # just the total reads. this is something separate.
    total_reads = {'same': 0, 'opposite': 0}

    # info on the number of 
    info_dict = {'same': 0, 'opposite': 0}

    tail_lens = {'same': {'A': [0,0,0,0,0], 'T': [0,0,0,0,0]},
                 'opposite':{'A': [0,0,0,0,0], 'T': [0,0,0,0,0]}} # total,g,a,t,c

    # method: each category in categories1 have each of the subcategories in
    # subcategories. these subcategories will have one of two dicsts: info_dict
    # and tail_lens.

    categories1 = ['Total clusters', 'morethan1', 'morethan1OA', 'only1']
    subcategories = ['All', 'annotated', 'wPAS', 'annotated_wPAS', 'goodPAS',
                     'bestPAS']

    bigcl = {}
    for cat1 in categories1:
        bigcl[cat1] = {}
        bigcl['total_reads'] = total_reads
        for cat2 in subcategories:
            bigcl[cat1][cat2] = {}
            bigcl[cat1][cat2]['info_dict'] = deepcopy(info_dict)
            bigcl[cat1][cat2]['tail_lens'] = deepcopy(tail_lens)

    dsets, super_3utr = super_falselength(settings, region, batch_key, subset,
                                          speedrun)

    for utr_name, utr in super_3utr[region].iteritems():

        for cls in utr.super_clusters:

            if cls.strand == utr.strand:
                keyw = 'same'
            else:
                keyw = 'opposite'

            total_reads[keyw] += cls.nr_support_reads

            bigcl['Total clusters'] = data_scooper(cls, keyw, bigcl['Total clusters'])

            # Count clusters with 2 or more reads
            if cls.nr_support_reads > 1:

                bigcl['morethan1'] = data_scooper(cls, keyw, bigcl['morethan1'])

            # Count clusters with 2 or more reads or annotated
            if cls.nr_support_reads > 1 or\
               cls.annotated_polyA_distance != 'NA':

                bigcl['morethan1OA'] = data_scooper(cls, keyw, bigcl['morethan1OA'])

            # Count clusters with only 1 read
            if cls.nr_support_reads == 1:

                bigcl['only1'] = data_scooper(cls, keyw, bigcl['only1'])


    return bigcl

def super_cluster_statprinter(dsetclusters, region, thiskey, settings, filename):

    statdict = dsetclusters[thiskey]

    keys = ['Total clusters', 'morethan1OA', 'morethan1', 'only1']

    subkeys =  ['All', 'wPAS', 'goodPAS', 'bestPAS', 'annotated',
                'annotated_wPAS']

    datakeys = ['info_dict', 'tail_lens']

    headers = {'Total clusters': '### All clustes ###',
               'morethan1': '### Clusters with 2 or more coverage ###',
               'only1': '### Clusters with only 1 coverage ###',
               'morethan1OA': '### Clusters with 2 or more or annotated ###'}

    subheaders = {'wPAS': 'With PAS',
                  'All': 'All',
                  'goodPAS': 'AATAAA or ATTAAA',
                  'bestPAS': 'AATAAA',
                  'annotated': 'Annotated',
                  'annotated_wPAS': 'Annotated with PAS'}

    reads_same = statdict['total_reads']['same']
    reads_opposite = statdict['total_reads']['opposite']
    reads_sum = reads_same + reads_opposite
    reads = (reads_sum, reads_same, reads_opposite)

    # you need it for percentages
    total_sum = statdict['Total clusters']['All']['info_dict']['opposite']\
            +\
            statdict['Total clusters']['All']['info_dict']['same']

    # just initating this value
    local_sum = total_sum

    output_path = os.path.join(settings.here,
                          'Results_and_figures/GENCODE_report/region_stats')

    output_file = os.path.join(output_path, region+'_'+filename+'.stats')
    handle = open(output_file, 'wb')

    handle.write('########################################################\n')
    handle.write(region+'\n')

    handle.write('Reads:{0} (same: {1}, opposite: {2})\n'.format(*reads))

    def divformat(numer, denom):
        if numer == 0:
            return 0
        else:
            return format(numer/float(denom), '.2f')

    for key in keys:

        handle.write('\n'+headers[key]+'\n')

        for dkey in datakeys:

            for subkey in subkeys:

                if dkey == 'info_dict':

                    same = statdict[key][subkey][dkey]['same']
                    opposite = statdict[key][subkey][dkey]['opposite']
                    so_sum = same + opposite

                    # All clusters
                    if subkey == 'All':
                        so_pcnt = divformat(so_sum, local_sum)
                        # must store the local sum for when not total clusters
                        local_sum = so_sum
                    else:
                        so_pcnt = divformat(so_sum, local_sum)

                    same_pcnt = divformat(same, so_sum)
                    oppo_pcnt = divformat(opposite, so_sum)

                    so = (so_sum, so_pcnt, same, same_pcnt, opposite, oppo_pcnt)
                    handle.write(subheaders[subkey]+':\t{0} ({1})\tsame {2} ({3})'\
                          '\topposite {4} ({5})\n'.format(*so))

                if dkey == 'tail_lens':
                    same = statdict[key][subkey][dkey]['same']
                    osite = statdict[key][subkey][dkey]['opposite']
                    keys = ['A', 'T']
                    so_sum = {}
                    for k in keys:
                        so_sum[k] = [same[k][i]+osite[k][i] for i in range(5)]

                    handle.write(subheaders[subkey]+'\n')

                    def divme(a,b):
                        try:
                            return format(a/b, '.2f')
                        except ZeroDivisionError:
                            return '0'

                    indx = {'T':3, 'A':2}
                    for k in keys:
                        snr = same[k][0]
                        ssnr = str(snr)

                        #sprint = ssnr+' '+str([divme(v, snr) for v in same[k][1:]])
                        sprint = ssnr+' '+divme(same[k][indx[k]], snr)
                        br = osite[k][0]
                        sor = str(br)
                        #oprint = sor+' '+str([divme(v, br) for v in osite[k][1:]])
                        oprint = sor+' '+divme(osite[k][indx[k]], br)

                        sr = so_sum[k][0]
                        ssr = str(sr)

                        #smprint = ssr+' '+str([divme(v, sr) for v in so_sum[k][1:]])
                        smprint = ssr+' '+divme(so_sum[k][indx[k]], sr)

                        #print('\tsum:\t\t' +smprint+'\n')
                        handle.write(k+'\tsame: '+sprint+'\topposite: '\
                                     +oprint+'\tsum: '+smprint+'\n')

    handle.write('########################################################\n')
    handle.close()

def clusterladder(settings, speedrun):
    """
    The more reads, the more polyAs, up to a point.
    """

    #1) Make a dictionary: dataset-> nr of total reads
    dsetreads = get_dsetreads(settings, region='3UTR')

    #2) Make super-clusters for your datasets of choice

    wc_c = [ds for ds in settings.datasets if (('Cytoplasm' in ds) or
                   ('Whole_Cell' in ds) or ('Nucleus' in ds)) and (not 'Minus' in ds)]

    wc_c_minus = [ds for ds in settings.datasets if (('Cytoplasm' in ds) or
                   ('Whole_Cell' in ds) or ('Nucleus' in ds)) and 'Minus' in ds]

    #data_grouping = {'Poly(A) plus': c,
                    #'Poly(A) minus': c_minus}
    data_grouping = {'Poly(A) plus': wc_c,
                     'Poly(A) minus': wc_c_minus}

    # small color dictionary
    colors = ['m', 'r', 'b', 'g', 'k']
    cols = {}
    for indx, title in enumerate(data_grouping.keys()):
        cols[title] = colors[indx]

    # keep a dictionary with reference to all the plots
    plotdict = {}

    #speedrun = True
    speedrun = False
    if speedrun:
        data_grouping['Poly(A) plus'] = data_grouping['Poly(A) plus'][:2]
        data_grouping['Poly(A) minus'] = data_grouping['Poly(A) minus'][:2]

    region = 'whole'
    for title, dsets in data_grouping.items():

        # sort the dsets in cell_lines by # of reads
        def mysorter(dset):
            return get_dsetreads(settings, region='3UTR')[dset]
        all_dsets = sorted(dsets, key=mysorter, reverse=True)
        #all_dsets = sorted(dsets, key=mysorter)

        # add more and more datasets
        subsets = [all_dsets[:end] for end in range(1, len(all_dsets)+1)]

        subsetcounts = {}

        for subset in subsets:

            # Get the number of 'good' and 'all' clusters
            key = ':'.join(subset)
            batch_key = 'first_ladder'
            dsetclusters = get_dsetclusters(subset, region, settings,
                                            speedrun, batch_key)

            subsetcounts[key] = count_clusters(dsetclusters, dsetreads)

        plotdict[title] = subsetcounts

    # loop through and do the plotting
    (fig, axes) = plt.subplots(1,2)

    for plot_nr, (title, dictsum) in enumerate(plotdict.items()):
        ax = axes[plot_nr]

        read_counts = []
        cluster_counts = []
        PAScluster_counts = []

        # you must loop throug dependnig on lenghr
        def sorthelp(tup):
            return len(tup[0])
        for dsets, countdict in sorted(dictsum.items(), key=sorthelp):
            # the the sum of rpeads from these datasets
            x = [get_dsetreads(settings, '3UTR')[ds] for ds in dsets.split(':')]
            read_counts.append(sum(x))

            # get all cluster counts
            cluster_counts.append(countdict['All'])
            PAScluster_counts.append(countdict['PAS'])

        ax.plot(read_counts, cluster_counts, color=cols[title],
                      linewidth=2, label='All sites')[0]
        ax.plot(read_counts, PAScluster_counts, ls='--', color=cols[title],
                      linewidth=2, label='Sites with PAS')[0]

        ax.set_xlabel('Billons of reads', size=18)
        ax.set_ylabel('Polyadenylation sites', size=18)
        ax.set_title('Polyadenylation site discovery saturates fast', size=20)

        # Sort the legends to your preference
        ax.legend(loc=0)

        # Set a grid on the y-axis
        ax.yaxis.grid(True)
        ax.xaxis.grid(True)

        ax.set_title(title, size=15)

    output_dir = os.path.join(settings.here, 'Results_and_figures',
                              'GENCODE_report', 'Figures')

    fig.set_size_inches(26,12)
    filename = 'Saturation_plot'
    filepath = os.path.join(output_dir, filename+'.pdf')
    fig.savefig(filepath, format='pdf')
    filepath = os.path.join(output_dir, filename+'.eps')
    fig.savefig(filepath, format='eps', papertype='A4')

def count_clusters(dsetclusters, dsetreads):
    """
    Parse through and count the number of clusters with +1 and number with PAS
    """

    countdict = {
        'All': sum(dsetclusters['morethan1']['All']['info_dict'].values()),
        'PAS': sum(dsetclusters['morethan1']['wPAS']['info_dict'].values())}

    return countdict


def merge_regions(clusterstorage):
    """ You get in all the individual clusterstorage things for each region.
    [[1,2,3], [1,2,3]] the 1,2,3 are compatible, so merge them with each other.
    they should have the same keys, so you need to make a new external one, and
    iterate over the keys for each one and simply add the internal keys.

    you should output the a dict that can be read in this form. actually the
    form of those individual 1 2 and 3.

        for names, numbers in sorted(dsetclusters.items(), key=something):

    """
    outlist = []

    for dset_indx, dsetcluster in enumerate(clusterstorage[0]):

        output = dict((key, val) for key, val in dsetcluster.items())

        # go over the other region(s)
        for nam, num in clusterstorage[1][dset_indx].items():

            for k in output[nam].keys():
                output[nam][k] += num[k]

        outlist.append(output)

    return outlist

def super_bed(handle, super_3utr, region):
    """ Go through the super beds and save them to a bedfile if >1 or annotated,
    and create summary statistics
    """

    one_stats = {'one_with_PAS': 0, 'one_with_good_PAS': 0,
                 'one_with_annot': 0, 'one_with_annot_and_PAS': 0,
                 'one_with_annot_and_good_PAS': 0, 'total': 0,
                 'one_with_PAS_sans_annot': 0,
                 'one_with_good_PAS_sans_annot': 0}

    two_stats = {'two_with_PAS': 0, 'two_with_good_PAS': 0,
                 'two_with_annot': 0, 'two_with_annot_and_PAS': 0,
                 'two_with_annot_and_good_PAS': 0, 'total': 0,
                 'two_with_PAS_sans_annot': 0,
                 'two_with_good_PAS_sans_annot': 0}

    goodpas = set(['ATTAAA', 'AATAAA'])

    total = 0

    for utr_name, utr in super_3utr[region].iteritems():

        for cls in utr.super_clusters:
            total +=1

            # Write to file
            if cls.nr_support_reads>1 or cls.annotated_polyA_distance!='NA':

                beg = cls.polyA_coordinate

                entry = '\t'.join([utr.chrm, str(beg), str(beg+1), utr.ID,
                                   str(cls.nr_support_reads), utr.strand])

                # write only those with PAS!
                if cls.nearby_PAS[0] != 'NA':
                    handle.write(entry + '\n')

            # Make statistics

            # for PAS with 1 read
            if cls.nr_support_reads==1:
                one_stats['total'] += 1

                # with annot
                if cls.annotated_polyA_distance != 'NA':
                    one_stats['one_with_annot'] += 1

                    if cls.nearby_PAS[0] != 'NA':
                        one_stats['one_with_annot_and_PAS'] +=1

                # without annot
                else:
                    if cls.nearby_PAS[0] != 'NA':
                        one_stats['one_with_PAS_sans_annot'] += 1

                        if set.intersection(set(cls.nearby_PAS),goodpas)!=set([]):
                            one_stats['one_with_good_PAS_sans_annot'] += 1

                if cls.nearby_PAS[0] != 'NA':
                    one_stats['one_with_PAS'] += 1

                    if set.intersection(set(cls.nearby_PAS), goodpas) != set([]):
                        one_stats['one_with_good_PAS'] += 1

                        if cls.annotated_polyA_distance != 'NA':
                            one_stats['one_with_annot_and_good_PAS'] +=1

            # for PAS with 2 reads or more
            if cls.nr_support_reads > 1:
                two_stats['total'] += 1

                # with annot
                if cls.annotated_polyA_distance != 'NA':
                    two_stats['two_with_annot'] += 1

                    if cls.nearby_PAS[0] != 'NA':
                        two_stats['two_with_annot_and_PAS'] +=1

                # without annot
                else:
                    if cls.nearby_PAS[0] != 'NA':
                        two_stats['two_with_PAS_sans_annot'] += 1

                        if set.intersection(set(cls.nearby_PAS),goodpas)!=set([]):
                            two_stats['two_with_good_PAS_sans_annot'] += 1

                if  cls.nearby_PAS[0] != 'NA':
                    two_stats['two_with_PAS'] += 1

                    if set.intersection(set(cls.nearby_PAS), goodpas) != set([]):
                        two_stats['two_with_good_PAS'] += 1

                        if cls.annotated_polyA_distance != 'NA':
                            two_stats['two_with_annot_and_good_PAS'] +=1

    return one_stats, two_stats, total


def venn_polysites(settings, speedrun):
    """ Output all clustered poly(A) sites for 3UTR exons and for the rest of
    the genome.

    Can you somehow restrict poly(A)DB and GENCODEv7 sites to those 3UTRs that
    are actually expressed in your datasets? That would simplify matters.
    """

    all_ds = [ds for ds in settings.datasets if ((('Cytoplasm' in ds) or
                                                ('Whole_Cell' in ds) or
                                                ('Nucleus' in ds)) and
                                               (not 'Minus' in ds))]
    #all_ds = [ds for ds in settings.datasets if (('Cytoplasm' in ds) and
                                               #(not 'Minus' in ds))]
    #speedrun = True
    speedrun = False

    outdir = os.path.join(settings.here,
                          'Results_and_figures/GENCODE_report/venn_diagram')

    # will hold the paths of all the files that result from merging sites with
    # one another. begin by adding the two sources of poly(A) sites
    paths = {'gencode':\
             os.path.join(settings.here, 'annotated_polyAsites/gencode_polyA.bed'),
             'polyAdb':\
             os.path.join(settings.here, 'annotated_polyAsites/polyA_db_proper.bed')}

    region = 'whole'

    subset = all_ds
    if speedrun:
        subset = subset[:2]

    batch_key = 'venn'
    dsets, super_3utr = super_falselength(settings, region, batch_key,
                                          subset, speedrun)

    # create 'outfile' and fill up the stats dicts
    paths['whole'] = whole_tobed(super_3utr, outdir, region)

    # get all paths to merged polyA files and their linenumbers
    (paths, merged_stats) = intersect_polyAs(paths, outdir, region)

    # make a venn-diagram of the intersection of poly(A) sites between each of
    # the below regions and the polyAdb and the GENCODE poly(A) annotation.

    make_venn(paths, merged_stats, outdir, region, settings)


def make_venn(paths, wcdict, outdir, region, settings):
    """ Call upon R to summon forth the elusive Venn diagrams
    """

    outdir = os.path.join(settings.here, 'analysis', 'R_code')
    outpath = os.path.join(outdir, 'venn.r')
    outfile = open(outpath, 'wb')

    not_weighted = os.path.join(outdir,'gencode_polyAdb_discovered_not_weighted.pdf')
    weighted = os.path.join(outdir, 'gencode_polyAdb_discovered_weighted.pdf')

    # write the commands to an R script
    line00 = 'library(Vennerable, lib.loc='\
            '"/users/rg/jskancke/R/x86_64-unknown-linux-gnu-library/2.13/")'
    line0 = "pdf('{0}')".format(not_weighted)
    line1 = 'Vcombo <- Venn(SetNames = c("Discovered", "GENCODE V7", "polyAdb")'\
            ', Weight= c(0, {0}, {1}, {2}, {3}, {4}, {5}, {6}))'\
            .format(wcdict[region], #only polyAdb
                    wcdict['gencode'], # only GENC
                    wcdict['_I_'.join(sorted([region,'gencode']))],
                    wcdict['polyAdb'], # only in region
                    wcdict['_I_'.join(sorted([region,'polyAdb']))],
                    wcdict['_I_'.join(sorted(['gencode','polyAdb']))], # pAdb + GENC
                    wcdict['_I_'.join(sorted(['gencode','polyAdb', region]))])

    line2 = 'plot(Vcombo, doWeights=FALSE)'
    line3 = 'dev.off()'  # close the plot to access the file
    line4 = "pdf('{0}')".format(weighted)
    line5 = 'plot(Vcombo, doWeights=TRUE)'
    line6 = 'dev.off()'  # close the plot to access the file

    for line in [line00, line0, line1, line2, line3, line4, line5, line6]:
        outfile.write(line+'\n')

    outfile.close()

    # call r on the outfile
    cmd = ['Rscript', outpath]
    p = Popen(cmd)
    p.wait()

    fig_dir = os.path.join(settings.here, 'Results_and_figures',
                              'GENCODE_report', 'Figures')

    # move the pdfs to output
    import shutil
    for pdf in glob.glob(outdir+'/*.pdf'):
        dirname, filename = os.path.split(pdf)
        shutil.copy(pdf, os.path.join(fig_dir, filename))


def intersect_polyAs(paths, outdir, region):
    """ You have regions 3UTR, a-3UTR, genome (A1, A2, A3) where polyA reads
    land and PAdb (B) and GENCODE (C)

    For each of these 3 regions, compute all intersections with B and C:
     U = union
     I = intersection
     \ = complement

     Compute with bedTools:
     AIBIC
     AIB
     BIC
     AIC
     From these 4 intersections you can compute all the numbers you need:

     A\(BUC) = A - AIB - AIC - AIBIC (Just A)
     B\(AUC) = B - BIC - BIA - AIBIC (Just B)
     C\(BUC) = C - CIB - CIA - AIBIC (Just C)

     (AIB)\C = AIB - AIBIC (Just A and B)
     (AIC)\B = AIC - AIBIC (Just A and C)
     (CIB)\A = CIB - AIBIC (Just C and B)

     AIBIC (A, B and C) """

    # Compute with bedTools:
    # AIBIC
    # AIB
    # BIC
    # From these 3 intersections you can compute all the numbers you need

    # 4) intersect the two annotations
    intersecters = [region] + ['gencode', 'polyAdb']

    # Send in A and [B,C]
    paths = intersect_wrap(paths, intersecters, outdir, extendby=10)

    # Do a wc of all the paths in paths
    wc = {}
    for pathname, path in paths.items():
        wc[pathname] = sum((1 for line in open(path, 'rb')))

    # make a printout of wc for comparing with the final plots
    for name, count in wc.items():
        print name + ':\t' + str(count)
    print('\n')

    A = region
    B = 'gencode'
    C = 'polyAdb'

    AIBIC = '_I_'.join(sorted([A,B,C]))

    venncount = {}
    venncount[AIBIC] = wc[AIBIC] # we don't need to compute this one

    # A\(BUC) = A - AIB - AIC - AIBIC (Just A)
    for let in [A, B, C]:
        compl = [A, B, C]
        compl.remove(let)
        compl.sort()

        # initialize with all poly(A) sites
        venncount[let] = wc[let]
        # remove the two "BIC" or "AIC" types
        for k in compl:
            venncount[let] -= wc['_I_'.join(sorted([let, k]))]
        # add the union of the three types, since it was subtracted twice
        venncount[let] += wc[AIBIC]

    # (AUB)\C = AIB - AIBIC (Just A and B)
    for let in [A, B, C]:
        compl = [A, B, C]
        compl.remove(let)
        (c1, c2) = sorted(compl)

        key = '{0}_I_{1}'.format(c1, c2) # one of (BUC)\A
        venncount[key] = wc['_I_'.join(sorted([c1, c2]))] - wc[AIBIC]

    return paths, venncount

def intersect_wrap(paths, intersecters, outdir, extendby=20):
    """
    1) extend each of the joiner files with some value
    2) intersect them
    3) return paths dict with link to the intersected path
    """
    #hg19 = '/home/jorgsk/work/3UTR/ext_files/hg19'
    hg19 = '/users/rg/jskancke/phdproject/3UTR/the_project/ext_files/hg19'

    # function used only here
    def extender(extendby, path, hg19, ext_path):
        cmd1 = ['slopBed', '-b', str(extendby), '-i', path, '-g', hg19]
        p1 = Popen(cmd1, stdout = open(ext_path, 'wb'))
        p1.wait()

    # 1) extend all files
    for pathname, path in paths.items():
        # skip those you don't want
        if pathname not in intersecters:
            continue

        ext_name = pathname+'_extended_{0}'.format(extendby)
        ext_path = os.path.join(outdir, ext_name)

        # if you have extended already, don't do it again
        if ext_name in paths:
            continue

        extender(extendby, path, hg19, ext_path)
        paths[ext_name] = ext_path

    # there are 4 intersections to be performed:
        # AIB, AIC, BIC, and AIBIC

    for r in [2, 3]:
        combinations = list(combins(intersecters, r))
        for combo in combinations:
            combo = sorted(combo)

            isect_name = '_I_'.join(combo) # sort to access easier
            isect_path = os.path.join(outdir, isect_name)

            if len(combo) == 2:
                # refer to the extended files
                (filea, fileb) = [paths[combo[0]+'_extended_{0}'.format(extendby)],
                                  paths[combo[1]+'_extended_{0}'.format(extendby)]]

            if len(combo) == 3:
                # 1) extend the intersection that exists of the first two
                orig_name = '_I_'.join(combo[:-1])
                orig_path = os.path.join(outdir, orig_name)
                ext_name = '_I_'.join(combo[:-1])+'_extended_{0}'.format(extendby)
                ext_path = os.path.join(outdir, ext_name)
                extender(extendby, orig_path, hg19, ext_path)
                # 2) Use the extended intersection as file a and the normal as b
                # Why is this one so sensitive to extendby?
                (filea, fileb) = (ext_path,
                                  paths[combo[2]+'_extended_{0}'.format(extendby)])

            cmd2 = ['intersectBed', '-wa', '-s', '-a', filea, '-b', fileb]
            p2 = Popen(cmd2, stdout=PIPE )

            out_handle = open(isect_path, 'wb')
            for line in p2.stdout:
                # re-center the intersection
                if len(line.split()) == 6:
                    (chrm, beg, end, name, val, strnd) = line.split()

                if len(line.split()) == 4:
                    (chrm, beg, end, strnd) = line.split()

                center_dist = math.ceil((int(end)-int(beg))/2)
                center = int(int(end) - center_dist)

                # write the center of the new mergers
                out_handle.write('\t'.join([chrm, str(center), str(center+1), '0',
                                           '0', strnd]) + '\n')
            out_handle.close()
            paths[isect_name] = isect_path

    return paths

def whole_tobed(super_3utr, outdir, region):
    """ Merge the paths in 'only_these' into one file, and then run mergeBed on
    them, and finally report the centers of the mergings obtained.  Return
    'paths' with the path to the merged version.
    """

    superbed_path = os.path.join(outdir, region+'.bed')
    handle = open(superbed_path, 'wb')


    ### lines in jouned_path
    #goodpas = set(['ATTAAA', 'AATAAA'])

    annotated = 0
    for utr_name, utr in super_3utr[region].iteritems():

        for cls in utr.super_clusters:

            # Write to file
            if cls.nr_support_reads>1 or cls.annotated_polyA_distance!='NA':

                beg = cls.polyA_coordinate

                entry = '\t'.join([utr.chrm, str(beg), str(beg+1), utr.ID,
                                   str(cls.nr_support_reads), utr.strand])

                handle.write(entry + '\n')
                # write only those with PAS!
                if cls.nearby_PAS[0] != 'NA':
                    annotated +=1

    handle.close()

    print 'annotated', annotated

    return superbed_path


def get_dsetnames(settings, compartments, ignorers, demanders):
    all_dsets = []

    for comp in compartments:
        for ds in settings.datasets:
            # if not discarding any
            if ignorers == [] and demanders == []:
                if comp in ds:
                    all_dsets.append(ds)

            # you must dicard some
            elif ignorers != [] and demanders == []:
                not_throw = True

                for ign in ignorers:
                    if ign in ds:
                        not_throw = False

                if comp in ds and not_throw:
                    all_dsets.append(ds)

            # some must be there
            elif ignorers == [] and demanders != []:
                not_throw = True

                for dem in demanders:
                    if dem not in ds:
                        not_throw = False

                if comp in ds and not_throw:
                    all_dsets.append(ds)

            elif ignorers != [] and demanders != []:
                not_throw = True

                for ign in ignorers:
                    if ign in ds:
                        not_throw = False

                for dem in demanders:
                    if dem not in ds:
                        not_throw = False

                if comp in ds and comp in demanders and not_throw:
                    all_dsets.append(ds)

    return all_dsets

def cumul_stats_printer_all(settings, speedrun):

    #speedrun = True
    speedrun = False

    region = 'whole'
    for fraction in ['Plus', 'Minus']:

        if fraction == 'Plus':
            all_dsets = [ds for ds in settings.datasets if (('Cytoplasm' in ds) or
                                                            ('Whole_Cell' in ds) or
                                                            ('Nucleus' in ds)) and
                         ('Minus' not in ds)]

        if fraction == 'Minus':
            all_dsets = [ds for ds in settings.datasets if (('Cytoplasm' in ds) or
                                                            ('Whole_Cell' in ds) or
                                                            ('Nucleus' in ds)) and
                         ('Minus' in ds)]

        if speedrun:
            all_dsets = all_dsets[:2]

        dsetclusters = {}

        batch_key = 'all_stats_3utr'
        dsetclusters[region] = get_dsetclusters(all_dsets, region, settings,
                                                speedrun, batch_key)

        filename = '_'.join([region, fraction])

        super_cluster_statprinter(dsetclusters, region, region, settings, filename)

def cumul_stats_printer_genome(settings, speedrun):

    # RARE LOOP JUST FOR GETTING THE CYTOPLASM/NUCLEUS INTEL
    region = 'whole'
    for compartments in [['Cytoplasm'], ['Nucleus'], ['Whole_Cell']]:
        for ignorers in [[], ['Minus']]:
            if ignorers == []:
                demanders = ['Minus']
            else:
                demanders = []

            if ignorers == []:
                filename = '+'.join(compartments+demanders)
            else:
                filename = '+'.join(compartments+demanders)+'-'+'-'.join(ignorers)

            # gets all dsets with the above settings
            all_dsets = get_dsetnames(settings, compartments, ignorers, demanders)

            dsetclusters = {}

            batch_key = 'cumul_stats_genome'
            dsetclusters['region'] = get_dsetclusters(all_dsets, region,
                                                    settings, speedrun, batch_key)

            # print the two regions
            super_cluster_statprinter(dsetclusters, region, region, settings,
                                      filename)

def cumul_stats_printer(settings, speedrun):

    regions = ['3UTR-exonic', 'CDS-exonic', 'CDS-intronic', 'Intergenic', 'whole']

    for region in regions:

        for compartments in [['Cytoplasm'], ['Nucleus'], ['Whole_Cell']]:
            for ignorers in [[], ['Minus']]:
                if ignorers == []:
                    demanders = ['Minus']
                else:
                    demanders = []

                if ignorers == []:
                    filename = '+'.join(compartments+demanders+region)
                else:
                    filename = '+'.join(compartments+demanders+region)\
                            +'-'+'-'.join(ignorers)

                subset = get_dsetnames(settings, compartments, ignorers, demanders)

                dsetclusters = {}


                # Get the number of 'good' and 'all' clusters
                key = ':'.join(subset)
                batch_key = 'cuml_st_printer'
                dsetclusters[key] = get_dsetclusters(subset, region, settings,
                                                     speedrun, batch_key)

                # NEW! Print out some statistics like this for each clustering.
                super_cluster_statprinter(dsetclusters, region, key, settings,
                                          filename)


def dsetcluster_join(dsetclusters):
    dsetclusters['genome'] = deepcopy(dsetclusters['3UTR'])
    for key1, dict1 in dsetclusters['anti-3UTR'].items():
        for key2, dict2 in dict1.items():
            if type(dict2) is int:
                dsetclusters['genome'][key1][key2] += dict2
            else:
                for key3, dict3 in dict2.items():
                    for key4, dict4 in dict3.items():
                        if type(dict4) is int:
                            dsetclusters['genome'][key1][key2][key3][key4] += dict4
                        else:
                            for key5, dict5 in dict4.items():
                                for indx, val in enumerate(dict5):
                                    dsetclusters['genome'][key1][key2][key3][key4]\
                                            [key5][indx] += val
    return dsetclusters['genome']

def apa_dict(settings):
    """
    Fetch dictionary of annotated poly(A) sites

    OBS! Will use whatever genomic region is in UTR_SETTINGS!
    The method will force 3UTRs from annotation, so this path must be provided
    in the UTR_SETTINGS file
    """
    import utail as utail

    utail_settings = utail.Settings(*utail.read_settings(settings.settings_file))
    utail_annotation = utail.Annotation(utail_settings.annotation_path,
                                        utail_settings.annotation_format,
                                        utail_settings.annotated_polyA_sites)

    chr1 = False
    beddir = os.path.join(settings.here, 'source_bedfiles')
    rerun_annotation_parser = False
    region_file = 'NA'
    # You need to set region_file_provided to false
    utail_settings.regionfile_provided = False
    utail_annotation.utrfile_path = utail.get_utr_path(utail_settings, beddir,
                                           rerun_annotation_parser, region_file)
    utail_annotation.feature_coords = utail_annotation.get_utrdict()
    return utail_annotation.get_polyA_dict(chr1)

def rpkm_polyA_correlation(settings, speedrun):
    """
    Gather the poly(A) sites in 3UTRs for all cytoplasmic compartments.
    Actually, don't bother. Just show the scatter for two compartments. You'll
    need to read both length and non-length input files. Make a simple
    dict[utr_id] = (#poly(A)_reads, #poly(A)_sites, RPKM)
    """

    # obs nr1. you don't see what you thought you'd see. if this scatter plot is
    # to believed, we don't see the tip of the ice-berg: rather, we see casual
    # glipses, as patches of land spotted through shifting coulds as a plane is
    # ascending into the skies.

    # Could it be that the number of poly(A) sites in the 3UTR is independent of
    # the RPKM? Maybe you should compare with the % of "recaptured" poly(A)
    # sites. What is the easiest way to find out about this? 
    # hey, you weren't getting all the p(A) sites dude ...

    a_polyA_sites_dict = apa_dict(settings)

    # sum the poly(A) sites for the murky few multi-exons out there
    apa_sites = {}
    for key, pAsites in a_polyA_sites_dict.iteritems():
        core_key = '_'.join(key.split('_')[:-1])
        if core_key not in apa_sites:
            apa_sites[core_key] = pAsites
        else:
            if pAsites == []:
                continue
            else:
                for pA in pAsites:
                    apa_sites[core_key].append(pA)

    # use the good old reader
    #speedrun = True
    speedrun = False
    region = '3UTR'

    dset_subsets = [ds for ds in settings.datasets if (('Cytoplasm' in ds) or
                   ('Whole_Cell' in ds) or ('Nucleus' in ds)) and (not 'Minus' in ds)]

    #dset_subsets = [ds for ds in settings.datasets if ('Cytoplasm' in ds) and
                   #(not 'Minus' in ds)]

    #dset_subsets = [ds for ds in settings.datasets if not 'Minus' in ds]

    # Now going through every single cytoplasm and whole cell alone
    go = [] # rpkm, total_nr_reads, nr_clusters
    go2 = [] # rpkm, % of clusters covered

    for dubset in dset_subsets:

        dsets, super_3utr = get_utrs(settings, region, [dubset], speedrun)

        for utr_id, utr in super_3utr.iteritems():
            if len(utr.rpkms) >1:
                debug()
            utr.RPKM = np.mean(utr.rpkms) # redefine the rpkm
            if len(utr.rpkms) >1:
                debug()
            if utr.RPKM > 0.5:
                clu_nr = len(utr.super_clusters)

                if clu_nr > 0:
                    covnr = sum([clu.nr_support_reads for clu in utr.super_clusters])
                else:
                    covnr = 0
                go.append((utr.RPKM, clu_nr, covnr))

                if utr_id in apa_sites:
                    a_clu_nr = len(apa_sites[utr_id])
                    if a_clu_nr > 0:
                        clu_frac = clu_nr/a_clu_nr
                    go2.append((utr.RPKM, clu_frac))

    go = np.array(go)
    rpkms = go[:,0]
    cluster_nrs = go[:,1]
    read_nrs = go[:,2]

    go2 = np.array(go2)
    rpkms2 = go2[:,0]
    clus_frac = go2[:,1]

    print('\nMean discovered/annotated ratio: {0:.2f}'.format(np.mean(clus_frac)))
    print('Mean discovered: {0:.2f}\n'.format(np.mean(cluster_nrs)))

    print(stats.spearmanr(rpkms2, clus_frac))

    # Divide the RPKMs into intervals
    #intervals = [(0,1), (1,3), (3,6), (6,10), (10,20), (20,40), (40,80),
                               #(80, max(rpkms2)+100)]
    intervals = [(0,5), (5,10), (10,15), (15,20), (20,25), (25,30), (30,35),
                               (35, max(rpkms2)+100)]

    clus_fracs = [[] for o in range(len(intervals))]

    for (rpkm, cl_frac) in go2:
        for ival_index, ival in enumerate(intervals):

            if ival[0] < rpkm < ival[1]:
                clus_fracs[ival_index].append(cl_frac)

    clus_fracs = [np.array(cls_fr) for cls_fr in clus_fracs]

    mean_clus_fracs = [np.mean(cl_frac) for cl_frac in clus_fracs]
    std_clus_fracs = [np.std(cl_frac) for cl_frac in clus_fracs]
    nrs = [len(cl_frac) for cl_frac in clus_fracs]

    # USe this good old plot :) code-reuse :)
    p = Plotter()
    p.rec_sensitivity(nrs, mean_clus_fracs, std_clus_fracs, intervals,
                      settings.here)

def barsense_counter(super_3utr, region):
    """
    Count for making bar plots and strand plots!

    UPDATE the count needs to be like this:
        [plus1][all, PAS, with T]
        [only1][all, PAS, with T]
    """
    subdict = {'all': 0, 'PAS': 0, 'T': 0, 'goodPAS': 0}
    bardict = {'3+': deepcopy(subdict),
               '2+': deepcopy(subdict),
               '1': deepcopy(subdict)}

    stranddict = {'3+': deepcopy(subdict),
                  '2+': deepcopy(subdict),
                  '1': deepcopy(subdict)}

    for utr_id, utr in super_3utr[region].iteritems():

        for cls in utr.super_clusters:

            # 1- sites
            if cls.nr_support_reads == 1:
                key = '1'

                bardict[key]['all'] += 1
                if cls.strand == utr.strand:
                    stranddict[key]['all'] += 1

                # any PAS
                if cls.nearby_PAS[0] != 'NA':
                    bardict[key]['PAS'] += 1
                    if cls.strand == utr.strand:
                        stranddict[key]['PAS'] += 1

                # for one of the two canonical PAS
                if 'AATAAA' in cls.nearby_PAS or 'ATTAAA' in cls.nearby_PAS:
                    bardict[key]['goodPAS'] += 1
                    if cls.strand == utr.strand:
                        stranddict[key]['goodPAS'] += 1

                # Get if this was an A or a T cluster
                if cls.tail_type == 'T':
                    bardict[key]['T'] += 1
                    if cls.strand == utr.strand:
                        stranddict[key]['T'] += 1

            # trusted sites
            if cls.nr_support_reads > 1:
                key = '2+'

                bardict[key]['all'] += 1
                if cls.strand == utr.strand:
                    stranddict[key]['all'] += 1

                # any PAS
                if cls.nearby_PAS[0] != 'NA':
                    bardict[key]['PAS'] += 1
                    if cls.strand == utr.strand:
                        stranddict[key]['PAS'] += 1

                # for one of the two canonical PAS
                if 'AATAAA' in cls.nearby_PAS or 'ATTAAA' in cls.nearby_PAS:
                    bardict[key]['goodPAS'] += 1
                    if cls.strand == utr.strand:
                        stranddict[key]['goodPAS'] += 1

                # Get if this was an A or a T cluster
                if cls.tail_type == 'T':
                    bardict[key]['T'] += 1
                    if cls.strand == utr.strand:
                        stranddict[key]['T'] += 1

            # more trusted sites
            if cls.nr_support_reads > 2:
                key = '3+'

                bardict[key]['all'] += 1
                if cls.strand == utr.strand:
                    stranddict[key]['all'] += 1

                # any PAS
                if cls.nearby_PAS[0] != 'NA':
                    bardict[key]['PAS'] += 1
                    if cls.strand == utr.strand:
                        stranddict[key]['PAS'] += 1

                # for one of the two canonical PAS
                if 'AATAAA' in cls.nearby_PAS or 'ATTAAA' in cls.nearby_PAS:
                    bardict[key]['goodPAS'] += 1
                    if cls.strand == utr.strand:
                        stranddict[key]['goodPAS'] += 1

                # Get if this was an A or a T cluster
                if cls.tail_type == 'T':
                    bardict[key]['T'] += 1
                    if cls.strand == utr.strand:
                        stranddict[key]['T'] += 1

    normstranddict = {}
    # normalize the numbers in stranddict with those in bardict
    # XXX SKIPPING FOR NOW, not sure you want to use this again
    #for key, val in bardict.items():
        #normstranddict[key] = stranddict[key]/val

    return bardict, normstranddict

def side_sense_plot(settings, speedrun):
    """
    Side plot of poly(A) reads in different regions, as well as the
    sense/antisense debacle. You should probably run the region-stuff with
    non_stranded and the sense/antisense with stranded.

    Initial results show some contradicting results for the stranded regions.
    Maybe it's the overlap that's killing you? Maybe you need an overlap free
    region.

    IDEA: to simplify, should you show just the Exonic areas for 3UTR and 5UTR?
    Maybe show CDS intronic as well.
    """

    #1 Gather the data in a dictlike this [wc/n/c|/+/-]['region'] = [#cl, #cl w/pas]

    compartments = ['Whole_Cell', 'Cytoplasm', 'Nucleus']

    fractions = ['+', '-']

    #regions = ['5UTR-exonic', '5UTR-intronic', '3UTR-exonic', '3UTR-intronic',
               #'CDS-exonic', 'CDS-intronic', 'Nocoding-exonic',
               #'Noncoding-intronic', 'Intergenic']
    #regions = ['3UTR-exonic', 'CDS-exonic', 'CDS-intronic',
               #'Nocoding-exonic', 'Noncoding-intronic', 'Intergenic']
    regions = ['3UTR-exonic', 'CDS-exonic', 'CDS-intronic', 'Intergenic']

    # Get one dict for the bar plot and one dict for the sense-plot
    bar_dict = {}
    sense_dict = {}

    for comp in compartments:
        for frac in fractions:

            if frac == '+':
                subset = [ds for ds in settings.datasets if (comp in ds) and
                          (not 'Minus' in ds)]
            if frac == '-':
                subset = [ds for ds in settings.datasets if (comp in ds) and
                          ('Minus' in ds)]

            for region in regions:

                batch_key = 'side_sense'
                dsets, super_3utr = super_falselength(settings, region,
                                                      batch_key, subset,
                                                      speedrun=speedrun)

                key = ':'.join([comp, frac, region])

                # count the number clusters with +1, of those with PAS/good_PAS
                bar_dict[key], sense_dict[key] = barsense_counter(super_3utr,
                                                                  region)

    #pickfile = 'TEMP_PICKLE'

    #if not os.path.isfile(pickfile):
        #pickle.dump(bar_dict, open(pickfile, 'wb'))
    #else:
        #bar_dict = pickle.load(open(pickfile, 'rb'))

    p = Plotter()

    title = 'Polyadenlyation in different regions for different'\
            ' cellular compartments'
    ID = 'side'
    p.lying_bar_regions(bar_dict, regions, title, ID, settings.here)


def ratio_counter(dsetclusters):
    """ Count the numbers important for plotting difference between non-PAS
    ratios in datasets
    """
    thisdict = {'only1': {}, 'morethan1OA': {}}

    for dkey in ['morethan1OA', 'only1']:

        total_As = 0
        total_Ts = 0
        PAS_Ts = 0
        gPAS_Ts = 0

        for kw in ['opposite', 'same']:
            total_Ts += dsetclusters[dkey]['All']['tail_lens'][kw]['T'][0]
            total_As += dsetclusters[dkey]['All']['tail_lens'][kw]['A'][0]
            PAS_Ts += dsetclusters[dkey]['wPAS']['tail_lens'][kw]['T'][0]
            gPAS_Ts += dsetclusters[dkey]['goodPAS']['tail_lens'][kw]['T'][0]

        thisdict[dkey]['non_PAS_Ts'] = (total_Ts - PAS_Ts)/total_Ts
        thisdict[dkey]['non_gPAS_Ts'] = (total_Ts - gPAS_Ts)/total_Ts
        thisdict[dkey]['T percentage'] =  total_Ts/(total_Ts+total_As)
        thisdict[dkey]['TA ratio'] =  total_Ts/total_As

    return thisdict

def non_PAS_polyA(settings, speedrun):
    """ Make a measure from the relative amounts of non-PAS T-reads to T-reads.
    """
    compartments = ['Whole_Cell', 'Cytoplasm', 'Nucleus']
    fractions = ['+', '-']
    #regions = ['5UTR-exonic', '5UTR-intronic', '3UTR-exonic', '3UTR-intronic',
               #'CDS-exonic', 'CDS-intronic', 'Nocoding-exonic',
               #'Noncoding-intronic', 'Intergenic']
    #regions = ['3UTR-exonic', 'CDS-exonic', 'CDS-intronic',
               #'Nocoding-exonic', 'Noncoding-intronic', 'Intergenic']
    regions = ['3UTR-exonic', 'CDS-exonic', 'CDS-intronic', 'Intergenic',
               'anti-3UTR-exonic']

    #regions = ['anti-3UTR-exonic']

    #speedrun = True
    speedrun = False

    # Get one dict for the bar plot and one dict for the sense-plot
    ratio_dict = {}
    for comp in compartments:
        for frac in fractions:

            if frac == '+':
                subset = [ds for ds in settings.datasets if (comp in ds) and
                          (not 'Minus' in ds)]
            if frac == '-':
                subset = [ds for ds in settings.datasets if (comp in ds) and
                          ('Minus' in ds)]

            for region in regions:

                batch_key = 'nonPAS_polyAz'
                dsetclusters = get_dsetclusters(subset, region, settings,
                                                speedrun, batch_key)

                key = ':'.join([comp, frac, region])

                # count the number clusters with +1, of those with PAS/good_PAS
                ratio_dict[key] = ratio_counter(dsetclusters)

    p = Plotter()
    title = 'Non-PAS polyadenylation'
    p.non_PAS_difference(ratio_dict, regions, title, settings.here)
    # below prints each region separate, not updated. to new level in ratio_dict
    #p.non_PAS_difference_separate(ratio_dict, regions, title, settings.here)


def isect(settings, merged, extendby, temp_dir, comp, reg):
    """
    Expand file A and B file by 15 nt, intersect, and return
    A - intersection
    intersection
    B - intersection
    """

    # function used only here
    def extender(extendby, path, hg19, ext_path):
        cmd1 = ['slopBed', '-b', str(extendby), '-i', path, '-g', hg19]
        p1 = Popen(cmd1, stdout = open(ext_path, 'wb'))
        p1.wait()

    # 1) extend all files
    expanded = {}
    for pathname, path in merged.items():
        # skip those you don't want

        ext_name = pathname+'{0}_{1}_REextended_{2}'.format(comp, reg, extendby)
        ext_path = os.path.join(temp_dir, ext_name)

        extender(extendby, path, settings.hg19_path, ext_path)
        expanded[pathname] = ext_path

    filea, fileb = expanded.values()
    cmd2 = ['intersectBed', '-wa', '-wb', '-s', '-a', filea, '-b', fileb]
    p2 = Popen(cmd2, stdout=PIPE )

    isect_name = str(comp)+'_'+str(reg)+'_'+'_I_'.join(expanded.keys())
    isect_path = os.path.join(temp_dir, isect_name)

    out_handle = open(isect_path, 'wb')

    # save the intersected path
    for line in p2.stdout:
        (achrm, abeg, aend, aname, aval, astrnd,
         bchrm, bbeg, bend, bname, bval, bstrnd) = line.split('\t')

        name = '#'.join(aname.split(' '))

        # re-center the intersection
        if int(abeg) < int(bbeg):
            center_dist = math.ceil((int(bend)-int(abeg))/2)
            center = int(int(bend) - center_dist)
        else:
            center_dist = math.ceil((int(aend)-int(bbeg))/2)
            center = int(int(aend) - center_dist)

        # write the center of the new mergers
        out_handle.write('\t'.join([achrm, str(center), str(center+1), name,
                                   aval, astrnd]) + '\n')

    out_handle.close()

    # save the output in an output dict
    output_dict = {'intersection': isect_path}
    # remove the intersection from + and -

    # 1) expand the new intersection path!
    isect_exp_name = isect_name + '_ExPanDeD{0}'.format(extendby)
    isect_exp_path = os.path.join(temp_dir, isect_exp_name)
    extender(extendby, isect_path, settings.hg19_path, isect_exp_path)

    # 2) subtract from the two origial files that which overlaps the extended
    # intersection
    for pathname, path in merged.items():
        cmd = ['subtractBed', '-s', '-a', path, '-b', isect_exp_path]
        outpath = path + '_Sliced'
        p3 = Popen(cmd, stdout = open(outpath, 'wb'))
        p3.wait()

        output_dict[pathname+'_sliced'] = outpath

    return output_dict


def merge_paths(settings, paths, temp_dir, key, min_covr, comp, reg):
    # the file you will keep expanding and subtracting

    extendpath = os.path.join(temp_dir,
                              '{0}_{1}_extendify_{2}'.format(comp, reg, key))
    juncfreepath = os.path.join(temp_dir,\
                                '{0}_{1}_extendify_juncfree{2}'.format(comp, reg,
                                                                      key))
    tempextendpath = os.path.join(temp_dir,
                                  '{0}_{1}_tempextendify_juncfree{2}'.format(comp,
                                                                            reg,
                                                                            key))
    # append to the extend-path
    with open(extendpath, 'wb') as ext_handle:
        # i) concatenate the files
        for path in paths:
            handle = open(path, 'rb')
            handle.next() # skip header

            for line in handle:
                (chrm, beg, end, utr_ID, strand, polyA_coordinate,
                 polyA_coordinate_strand, polyA_average_composition,
                 annotated_polyA_distance, nearby_PAS, PAS_distance,
                 number_supporting_reads, number_unique_supporting_reads,
                 unique_reads_spread) = line.split('\t')

                # write a bed format to disc
                ext_handle.write('\t'.join([chrm,
                                          polyA_coordinate,
                                          str(int(polyA_coordinate)+1),
                                          nearby_PAS,
                                          number_supporting_reads,
                                          polyA_coordinate_strand +'\n']))

    # i.5) cut away the exon-exon noise
    junctions = os.path.join(settings.here, 'junctions',
                             'splice_junctions_merged.bed')
    cmd = ['subtractBed', '-a', extendpath, '-b', junctions]
    p = Popen(cmd, stdout=open(juncfreepath, 'wb'))
    p.wait()

    # ii) expand the entries
    cmd = ['slopBed', '-i', juncfreepath, '-g', settings.hg19_path, '-b', '15']
    p = Popen(cmd, stdout=open(tempextendpath, 'wb'))
    p.wait()

    # iii) merge and sum scores
    cmd = ['mergeBed','-s', '-nms', '-scores', 'sum', '-i', tempextendpath]
    #cmd = ['mergeBed','-s', '-nms', '-scores', 'max', '-i', tempextendpath]
    p = Popen(cmd, stdout=PIPE)

    # iv) center the new transcript
    finalpath = os.path.join(temp_dir,
                             '{0}_{1}_extendify_juncfree_merged{2}'.format(comp,
                                                                          reg,
                                                                          key))
    out_handle = open(finalpath, 'wb')
    for line in p.stdout:
        (chrm, beg, end, pases, maxcovrg, strand) = line.split('\t')

        meanbeg = int((int(beg)+int(end))/2)

        ps0 = pases.split(';')
        ps1 = pases.split(' ')
        ps2 = pases.split('#')

        go = ps0 + ps1 + ps2
        supergo = [g for g in go if (len(g)==6 or len(g)==2)]

        ps = list(set(supergo))
        if len(ps) > 1 and 'NA' in ps:
            ps.remove('NA')

        pas = '#'.join(ps)

        # filter on coverage
        if int(float(maxcovrg)) > min_covr:

            out_handle.write('\t'.join([chrm,
                                       str(meanbeg),
                                       str(meanbeg+1),
                                       pas,
                                       maxcovrg,
                                       strand.rstrip()+'\n']))
    out_handle.close()

    return finalpath

def all_sideplot(settings):
    """
    Intersect poly(A)+ and poly(A)- for the all regions in all compartments.
    Finally merge the poly(A) sites for each region and plot them sidebarwise.
    """
    co = ['Whole_Cell', 'Cytoplasm', 'Nucleus']
    #co = ['Whole_Cell', 'Whole_Cell', 'Whole_Cell']

    #regions = ['5UTR-exonic', '5UTR-intronic', '3UTR-exonic', '3UTR-intronic',
               #'CDS-exonic', 'CDS-intronic', 'Nocoding-exonic',
               #'Noncoding-intronic', 'Intergenic']
    regions = ['3UTR-exonic', 'CDS-exonic', 'CDS-intronic', 'Intergenic']
    #regions = ['3UTR-exonic', 'anti-3UTR-exonic']

    # Get one dict for the bar plot and one dict for the sense-plot

    # file paths for all regions
    dsetdict = dict((reg, settings.only_files(reg)) for reg in regions)

    temp_dir = os.path.join(settings.here, 'temp_files')
    min_covr = 1

    # Merge the plus and minus subsets for each region
    reg_dict = {}
    for reg in regions:

        plus_subset = [path for ds, path in dsetdict[reg].items() if
                       (not 'Minus' in ds) and ((co[0] in ds) or
                                                (co[1] in ds) or
                                                (co[2] in ds))]

        minus_subset = [path for ds, path in dsetdict[reg].items() if
                        ('Minus' in ds) and ((co[0] in ds) or
                                             (co[1] in ds) or
                                             (co[2] in ds))]
        # FOR SPEED
        #plus_subset = plus_subset[:1]
        #minus_subset = minus_subset[:1]

        path_dict = {'plus': plus_subset, 'minus': minus_subset}
        merged = {}

        comp ='All'
        for key, paths in path_dict.items():
            merged[key] = merge_paths(settings, paths, temp_dir, key,
                                      min_covr, comp, reg)

        reg_dict[reg] = merged

    ## save the bedfiles for the different regions
    #save_dir = os.path.join(settings.here, 'analysis', 'all_sideplot')
    #save_pure(comp_dict, save_dir)

    data_dict = count_these_all(reg_dict)

    p = Plotter()

    title = 'Polyadenylation in different regions'

    p.all_lying_bar(data_dict, regions, title, settings.here)

def intersection_sideplot(settings):
    """
    Intersect poly(A)+ and poly(A)- for the different regions and
    """
    compartments = ['Whole_Cell', 'Cytoplasm', 'Nucleus']
    #compartments = ['Whole_Cell']

    #regions = ['5UTR-exonic', '5UTR-intronic', '3UTR-exonic', '3UTR-intronic',
               #'CDS-exonic', 'CDS-intronic', 'Nocoding-exonic',
               #'Noncoding-intronic', 'Intergenic']
    #regions = ['3UTR-exonic', 'CDS-exonic', 'CDS-intronic',
               #'Nocoding-exonic', 'Noncoding-intronic', 'Intergenic']
    regions = ['3UTR-exonic', 'CDS-exonic', 'CDS-intronic', 'Intergenic']
    #regions = ['3UTR-exonic']

    # Get one dict for the bar plot and one dict for the sense-plot

    # file paths for all regions
    dsetdict = dict((reg, settings.only_files(reg)) for reg in regions)

    temp_dir = os.path.join(settings.here, 'temp_files')
    min_covr = 1
    extendby = 15

    comp_dict = {}
    for comp in compartments:

        # Merge the plus and minus subsets for each region
        reg_dict = {}
        for reg in regions:

            plus_subset = [path for ds, path in dsetdict[reg].items() if 
                           (not 'Minus' in ds) and comp in ds]
            minus_subset = [path for ds, path in dsetdict[reg].items() if 
                            ('Minus' in ds) and comp in ds]

            path_dict = {'plus': plus_subset, 'minus': minus_subset}
            merged = {}

            for key, paths in path_dict.items():
                merged[key] = merge_paths(settings, paths, temp_dir, key,
                                          min_covr, comp, reg)

            # return a dictionary: separated['plus'/'minus'/'separated]
            separated = isect(settings, merged, extendby, temp_dir, comp, reg)

            reg_dict[reg] = separated

        comp_dict[comp] = reg_dict

    # save the bedfiles for the different regions
    save_dir = os.path.join(settings.here, 'analysis', 'pure')
    save_pure(comp_dict, save_dir)

    data_dict = count_these(comp_dict)

    p = Plotter()

    title = 'Polyadenlyation in different regions for different'\
            ' cellular compartments'

    p.intersect_lying_bar(data_dict, regions, title, settings.here)


def save_pure(comp_dict, save_dir):
    """
    Save as bed-files the poly(A) sites from the pure regions
    """

    for compartment, reg_dict in comp_dict.items():
        for region, slice_dict in reg_dict.items():
            for sl, slice_path in slice_dict.items():

                outname = '__'.join([compartment, region, sl])
                outpath = os.path.join(save_dir, outname)
                outhandle = open(outpath, 'wb')

                for line in open(slice_path, 'rb'):
                    outhandle.write(line)

                outhandle.close()

def count_these(some_dict):
    """
    Count PAS and all. You know its +1 already
    """

    goodpas = set(['AATAAA', 'ATTAAA'])

    allpas = set(['AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA',
                 'AATACA', 'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA',
                 'AATAGA'])

    data_dict = AutoVivification()

    for comp, comp_dict in some_dict.items():
        for reg, reg_dict in comp_dict.items():
            for key, keypath in reg_dict.items():

                # get how many of the isect have PAS
                pas = 0
                goodPas = 0
                All = 0

                for line in open(keypath, 'rb'):
                    (chrm, beg, end, PAS, covr, strand) = line.split('\t')

                    PAS = '#'.join(PAS.split(' '))

                    All += 1

                    has_pas = False
                    has_good_pas = False

                    for pa in PAS.split('#'):
                        if pa in allpas:
                            has_pas = True
                        if pa in goodpas:
                            has_good_pas = True

                    if has_pas:
                        pas += 1
                    if has_good_pas:
                        goodPas +=1

                data_dict[comp][reg][key]['all'] = All
                data_dict[comp][reg][key]['PAS'] = pas

    return data_dict

def count_these_all(some_dict):
    """
    Count PAS and all. You know its +1 already
    """

    goodpas = set(['AATAAA', 'ATTAAA'])

    allpas = set(['AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA',
                 'AATACA', 'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA',
                 'AATAGA'])

    data_dict = AutoVivification()

    for reg, reg_dict in some_dict.items():
        for key, keypath in reg_dict.items():

            # get how many of the isect have PAS
            pas = 0
            goodPas = 0
            All = 0

            for line in open(keypath, 'rb'):
                (chrm, beg, end, PAS, covr, strand) = line.split('\t')

                PAS = '#'.join(PAS.split(' '))

                All += 1

                has_pas = False
                has_good_pas = False

                for pa in PAS.split('#'):
                    if pa in allpas:
                        has_pas = True
                    if pa in goodpas:
                        has_good_pas = True

                if has_pas:
                    pas += 1
                if has_good_pas:
                    goodPas +=1

            data_dict[reg][key]['all'] = All
            data_dict[reg][key]['PAS'] = pas

    return data_dict

def get_genc3_paths(carrie, my_cell_lines, regions):

    #0) make a directory -> cell line dictionary
    dir_2_cell_line = {}
    for line in open(os.path.join(carrie, 'Dataset.info.txt')):
        (dirinfo, cell_l, comp, longshort, rplica) = line.split()
        dirname = dirinfo[4:]
        dir_2_cell_line[dirname] = cell_l

    #1 collect file paths for nuclear and cytoplasmic exons
    regions = ['nuclear', 'cytoplasmic']
    paths = dict((reg, []) for reg in regions)

    for dirname, dirnames, filenames in os.walk(carrie):
        if dirname.endswith('exons'):

            main_dir = os.path.split(os.path.split(dirname)[0])[1]

            # include only those cell lines you have been using
            if main_dir not in dir_2_cell_line:
                print main_dir
                continue

            if dir_2_cell_line[main_dir] not in my_cell_lines:
                continue

            compartment = main_dir[3:]

            # we're only looking at nuclear and cytoplasm
            if compartment not in ['N','C']:
                continue

            for f in filenames:
                if f == 'All.exon.rpkm.pooled.txt.gz':

                    mypath = os.path.join(dirname, f)

                    if compartment == 'C':
                        paths['cytoplasmic'].append(mypath)
                    if compartment == 'N':
                        paths['nuclear'].append(mypath)

    return paths


def cytonuclear_rpkms_genc3(settings):
    """
    See other function
    """
    temp_dir = os.path.join(settings.here, 'temp_files')
    carrie = '/users/rg/projects/NGS/Projects/ENCODE/hg19main/'\
            'GingerasCarrie'

    my_cell_lines = ('K562', 'GM12878', 'HUVEC', 'HELAS3', 'HEPG2')

    #0) get paths to cytoplasmic and nuclear rpkms
    paths = get_genc3_paths(carrie, my_cell_lines)


    #1) Get dicts of the exons' rpkms
    #  Can you keep all exons in memory? I think so.

    ex_dict = {'nuclear': {}, 'cytoplasmic': {}}

    for region in ['nuclear', 'cytoplasmic']:

        for f in paths[region]:
            cmd = ['zcat', f]
            p = Popen(cmd, stdout=PIPE)
            for line in p.stdout:
                coord, rpkm, d = line.split()

                if coord in ex_dict[region]:
                    ex_dict[region][coord].append(rpkm)
                else:
                    ex_dict[region][coord] = [rpkm]

    # 2) Write a bed-file with the average rpkm of the exons
    strandmaker = {'1':'+', '-1':'-'}

    rpkm_paths = {}
    # save all rpkms for getting the average differene between nucleus and
    # cytoplasm
    check_rpkms = {'nuclear': [], 'cytoplasmic': []}

    for region in ['nuclear', 'cytoplasmic']:

        outfile = region+'_all_exons_rpkm.bed'
        outpath = os.path.join(temp_dir, outfile)
        outhandle = open(outpath, 'wb')

        for coord, rpkms in ex_dict[region].iteritems():
            (chrm, beg, end, strnd) = coord.split('_')

            mean_rpkm = sum([float(fl) for fl in rpkms])/len(rpkms)
            check_rpkms[region].append(mean_rpkm)

            outhandle.write('\t'.join([chrm, beg, end, '0', str(mean_rpkm),
                                       strandmaker[strnd]]) + '\n')

        rpkm_paths[region] = outpath

        outhandle.close()

    # 2.5) As a control, get the average RPKM difference between nucleus and
    # cytoplasm
    average_ratio = sum(check_rpkms['cytoplasmic'])/sum(check_rpkms['nuclear'])

    # 3) intersect the two big files with the poly(A) reads you have in the
    #nucleus (optionally, do this only for exonic poly(A) reads).
    sliced_reads = os.path.join(settings.here, 'analysis', 'pure', '*')

    rpkms = {'nuclear': {}, 'cytoplasmic': {}}

    for f in glob.glob(sliced_reads):
        fname = os.path.split(f)[1]
        (compartment, region, fraction) = fname.split('__')

        for mycompartment in ['nuclear', 'cytoplasmic']:

            #if compartment == 'Nucleus' and region.endswith('exonic') and \
               #fraction == 'minus_sliced':

            if compartment == 'Nucleus' and fraction == 'minus_sliced':

                cmd = ['intersectBed', '-wa', '-wb', '-a', f, '-b',
                       rpkm_paths[mycompartment]]

                cmd = ['intersectBed', '-u', '-a', rpkm_paths[mycompartment],
                       '-b', f]

                p = Popen(cmd, stdout=PIPE)

                for line in p.stdout:
                    #(chma, bega, enda, pasa, covrga, stranda, chrm, begb, endb,
                    #d, rpkmb, strandb) = line.split('\t')
                    (chma, bega, enda, d, rpkma,  stranda) = line.split('\t')

                    if region in rpkms[mycompartment]:
                        rpkms[mycompartment][region].append(float(rpkma))
                    else:
                        rpkms[mycompartment][region] = [float(rpkma)]

    results = {}
    for region in rpkms['nuclear'].keys():
        results[region] = sum(rpkms['cytoplasmic'][region])/\
                sum(rpkms['nuclear'][region])

    debug()

    # would the results be better with gencode 7? who knoes!
    #Again compare the RPKM values. Do we see a difference?
    # There is no such difference.
    #ipdb> pp average_ratio 
    # results:

        #1.9767955597874043
        #ipdb> pp results
        #{'3UTR-exonic': 2.3708455693225732,
         #'5UTR-exonic': 1.4557833987822002,
         #'5UTR-intronic': 0.18392529867006108,
         #'CDS-exonic': 2.1470952631719835,
         #'CDS-intronic': 1.4636101747123955,
         #'Intergenic': 3.7702052091564511,
         #'Nocoding-exonic': 1.7442351505710341,
         #'Noncoding-intronic': 1.0318510528636997}


    debug()


def cytonuclear_rpkms(settings):
    """
    Observation: there are more poly(A)- reads in exonic regions in the nucleus
    than in the cytoplasmic regions (this is even more obvious in the intergenic
    region actually. wonder why.)

    Needs must: 1) Get the exons. They can be found in
    /users/rg/projects/NGS/Projects/ENCODE/hg19main/gencode7/CSHL_longRNAseq/PolyA/parts/*/exons/K5621LNS.exon.rpkm.pooled.txt.gz
    etc..

    It seems that the pipeline has not been run for all cell lines.
    Nevertheless, go through the poly(A)+ libraries and crate 2 big bedfiles
    from the nuclear and cytoplasmic regions. Do mergeBed and use the mean
    value. Compare the exon's rpkms so you know the overall rate.
    # NOTE merge-bed would collapse the exons into each other if they overlap.
    # You'd better do it yourself: making a hash table for each exon and adding
    # the rpkms you find

    Then intersect the two big files with the poly(A) reads you have in the
    nucleus (optionally, do this only for exonic poly(A) reads).

    Again compare the RPKM values. Do we see a difference?
    """

    temp_dir = os.path.join(settings.here, 'temp_files')
    polyAplus = '/users/rg/projects/NGS/Projects/ENCODE/hg19main/'\
            'gencode7/CSHL_longRNAseq/PolyA/parts/'
    #1 collect file paths for nuclear and cytoplasmic exons
    regions = ['nuclear', 'cytoplasmic']
    paths = dict((reg, []) for reg in regions)

    for dirname, dirnames, filenames in os.walk(polyAplus):
        if dirname.endswith('exons'):

            main_dir = os.path.split(os.path.split(dirname)[0])[1]
            compartment = main_dir.lstrip('CSHL').rstrip('Plus')[3:]

            if compartment not in ['N', 'WC', 'C']:
                print('Wrong parsing')
                debug()

            # we're only looking at nuclear and cytoplasm
            if compartment == 'WC':
                continue

            for f in filenames:
                if f.endswith('exon.rpkm.pooled.txt.gz'):

                    if compartment == 'C':
                        paths['cytoplasmic'].append(f)
                    if compartment == 'N':
                        paths['nuclear'].append(f)

    # XXX you're stopping here because the gencode v.7 run is not finished
    debug()


            #debug()

        #if 'HELAS31LNS.exon.rpkm.pooled.txt.gz' in filenames:
            #debug()
        #if myfile in filenames:
            #readpath = os.path.join(dirname, myfile)

    #1 Get exons
    #ex_dict = {'nuclear': {}, 'cytoplasmic': {}}

    #for region in ['nuclear', 'cytoplasmic']:

        #outfile = region+'_all_exons_rpkm.bed'
        #outpath = os.path.join(temp_dir, outfile)

        #handle = open(outpath, 'wb')

def clusterladder_cell_lines(settings):
    """
    Same as custerladder but now for each cell line
    """
    #1) Make a dictionary: dataset-> nr of total reads
    dsetreads = get_dsetreads(settings, region='3UTR')

    cell_lines = ['HUVEC', 'GM12878', 'HeLa-S3', 'K562', 'HEPG2']
    # cell_line color dictionary
    colors = ['m', 'r', 'b', 'g', 'k']
    cell_cols = {}

    #speedrun = True
    speedrun = False

    for indx, cline in enumerate(cell_lines):
        cell_cols[cline] = colors[indx]

    wc_c = [ds for ds in settings.datasets if ((('Cytoplasm' in ds) or
                                                ('Whole_Cell' in ds) or
                                                ('Nucleus' in ds)) and
                                               (not 'Minus' in ds))]

    wc_c_minus = [ds for ds in settings.datasets if ((('Cytoplasm' in ds) or
                                                      ('Whole_Cell' in ds)
                                                      or ('Nucleus' in ds))
                                                     and ('Minus' in ds))]
    # set figure and axes
    (fig, axes) = plt.subplots(1,2)

    region = 'whole'
    for cell_line in cell_lines:

        wc_c = [ds for ds in wc_c if cell_line in ds]
        wc_c_minus = [ds for ds in wc_c_minus if cell_line in ds]

        if speedrun:
            wc_c = wc_c[:2]
            wc_c_minus = wc_c_minus[:2]

        data_grouping = {'Poly(A) plus': wc_c,
                         'Poly(A) minus': wc_c_minus}

        # keep a dictionary with reference to all the plots
        plotdict = {}

        for title, dsets in data_grouping.items():

            # sort the dsets in cell_lines by # of reads
            def mysorter(dset):
                return get_dsetreads(settings, region='3UTR')[dset]
            all_dsets = sorted(dsets, key=mysorter, reverse=True)
            #all_dsets = sorted(dsets, key=mysorter)

            # add more and more datasets
            subsets = [all_dsets[:end] for end in range(1, len(all_dsets)+1)]

            # A dictionary with all clusters and +2 or annot clusters
            subsetcounts = {}

            for subset in subsets:

                # Get the number of 'good' and 'all' clusters
                key = ':'.join(subset)
                batch_key = 'cell_line_ladder' # a unique key for super
                dsetclusters = get_dsetclusters(subset, region, settings,
                                                speedrun, batch_key)

                subsetcounts[key] = count_clusters(dsetclusters, dsetreads)

            plotdict[title] = subsetcounts

        # go through the sums and plot them
        for plot_nr, (title, dictsum) in enumerate(plotdict.items()):
            ax = axes[plot_nr]

            read_counts = []
            cluster_counts = []
            PAScluster_counts = []

            # you must loop throug dependnig on length
            def sorthelp(tup):
                return len(tup[0])
            for dsets, countdict in sorted(dictsum.items(), key=sorthelp):
                # the the sum of reads from these datasets
                x = [get_dsetreads(settings, '3UTR')[ds] for ds in dsets.split(':')]
                read_counts.append(sum(x))

                # get all cluster counts
                cluster_counts.append(countdict['All'])
                PAScluster_counts.append(countdict['PAS'])

            ax.plot(read_counts, cluster_counts, color=cell_cols[cell_line],
                          linewidth=2, label=cell_line)[0]
            ax.plot(read_counts, PAScluster_counts, ls='--',
                          color=cell_cols[cell_line], linewidth=2)[0]

        # set the axes
        for ax in axes:
            ax.set_xlabel('Billons of reads', size=16)
            ax.set_ylabel('Polyadenylation sites', size=16)
            ax.set_title('Polyadenylation site discovery saturates fast', size=18)

            # Sort the legends to your preference
            ax.legend(loc=0)

            # Set a grid on the y-axis
            ax.yaxis.grid(True)
            ax.xaxis.grid(True)

            ax.set_title(title, size=13)

    output_dir = os.path.join(settings.here, 'Results_and_figures',
                              'GENCODE_report', 'Figures')

    fig.set_size_inches(26,12)
    filename = 'Saturation_plot_cell_lines'
    filepath = os.path.join(output_dir, filename+'.pdf')
    fig.savefig(filepath, format='pdf')
    filepath = os.path.join(output_dir, filename+'.eps')
    fig.savefig(filepath, format='eps', papertype='A4')

def gencode_report(settings, speedrun):
    """ Make figures for the GENCODE report. The report is an overview of
    evidence for polyadenylation from the Gingeras RNA-seq experiments.

    NOTE: you should estimate a some kind of maximum value. get the total number
    of expressed genes and use the average number of poly(A) sites per gene.

    Two main figures:
        1) Total number of new poly(A) sites in 3UTRs and outside (+2 reads or
        annotated) both with and without normalization to region size
        2) of new poly(A) sites obtained with increasing number of reads
    """
    #speedrun = True
    speedrun = False

    # 0.1) Core stats selected regions + genome
    #cumul_stats_printer(settings, speedrun)
    # 0.2 Core stats All regions All cell lines
    #cumul_stats_printer_all(settings, speedrun)

    # 1) nr of polyA sites obtained with increasing readnr
    # TODO store the results in 'max'. The run everything for 'sum'. Then decide
    # which you want to keep.
    #clusterladder(settings, speedrun)
    #clusterladder_cell_lines(settings)
    # RESULT: make two plots: one from each of the 'data_groupign' below. One
    # shows best how discovery tapers out, the other shows that poly(A) minus
    # controls have few poly(A) reads. Also list how many of those poly(A) minus
    # are annotated ones. Use the remainder as a false positive.
    # maybe the gencode people are mostly interested in a kind of validation of
    # their poly(A) minus datasets: I can show that there are very few poly(A)
    # sites in poly(A)+ compared to poly(A)-, and further, that those that are
    # there are mostly noise from the poly(A)- dataset. It is a way a measure of
    # the robustness of the poly(A)- dataset.

    # 2) venn diagram of all the poly(A) sites from the whole genome for 3UTR and
    #venn_polysites(settings, speedrun)
    # TODO! the venn isn't correct. When you intersect with only annotated
    # sites, you still get 16k hits that are not annotated. need to look at
    # that.

    # 3) plot of correlation between number of poly(A) sites expressed
    # transcript 3UTRs and the RPKM of the 3UTRs.
    #rpkm_polyA_correlation(settings, speedrun)# wait: re-reun length files
    # NOTE: this is not right. you must compare with RPKM of the exon in which
    # the poly(A) reads are found. Your best chance is the gencode annotation.

    # 4) Sidewise bar plots of the number of poly(A) incidents in the different
    # genomic regions, for poly(A)+ and poly(A)-. Can also make the plot of the
    # 'sensedness' of the strands, but for this you must change the files in the
    # output directory to come from stranded_sequences.
    #side_sense_plot(settings, speedrun) # DONE!

    # 5) Poly(A)+ pure, intersection, poly(A)-
    #intersection_sideplot(settings)

    # 6) All side plot!
    #all_sideplot(settings) # just this one left, then copy to 'sum'. then
    #compare. then choose.

    # 6) Cytoplasmic vs. nuclear RPKMS in 3UTR exonic exons
    #cytonuclear_rpkms(settings)
    #cytonuclear_rpkms_genc3(settings) # negative result

    # Within: PAS, Good PAS, Annotated, T/A ratio for these?
    # It's not going to be easiy. Maybe just 3UTR exonic/non-3UTR exonic?
    # is it worth it to make this big table? Maybe you can represent the T %
    # from 0 to 100, just like with the strandedness. Then you could make a plot
    # of that! Would be nice if you could represent your key plots with PAS, and
    # T%! Better than the ratio.

    # 6) Comparison of A/T PAS annotated etc for poly(A)- C/N show that there is
    # enriched expression of short polyA(A) reads. Can potentially be
    # spontaneous premature annotation, or degradation related polyadenylation.
    #non_PAS_polyA(settings, speedrun)

    # 7) Demonstrate that you can predict the strand with 3UTR

    # From working at home: it's most convincing if you give a bar plot that
    # contains 5X2 bars, each one a comparison of the A/T ratio in each region.
    # This information will already be in the big figure, but re-make it in a
    # standing-up way for emphasis

    # split-mapped reads:     2.94e+07 (0.010)

    #5UTR exons:              8.18e+06 (0.003)
    #3UTR introns:            1.19e+07 (0.004)
    #3UTR exons:              2.99e+07 (0.010)
    #Nocoding exons:          3.27e+07 (0.011)
    #CDS exons:               3.54e+07 (0.011)
    #5UTR introns:            1.46e+08 (0.047)
    #Noncoding introns:       2.48e+08 (0.080)
    #CDS introns:             9.76e+08 (0.315)
    #Intergenic:              1.61e+09 (0.519)

    #Total:                   3.10e+09

    # 16000  in 3UTr and 13000 outside.
    # the 3UTR regions is 30 million basepairs, vs 3.1 billion for the whole genome.
    # this leaves. 3UTR is 1 percent of the genome

    # XXX 16/0.01 = 1600. 13/0.99 = 13.13. 1600/13.13 = 121 fold enrichment of
     #poly(A) sites in annotated 3UTRs.


class Cufflink_exon(object):

    def __init__(self, chrm, beg, end, strand, transcript_id, gene_id, exon_nr,
                Utail_ID):
        #self.chrm = chrm
        #self.beg = beg
        #self.end = end
        self.strand = strand
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.exon_nr = exon_nr

        self.Utail_ID = Utail_ID

        # will be added later
        #self.cluster_coords = []
        #self.annot_list = []

def cufflink_super(settings, speedrun, region):
    """ Return a super-3UTR for the cufflink model (poly(A)+ all cell lines
    cytoplasm, nucleus, and whole cell)
    """
    subset = [ds for ds in settings.datasets if (('Cytoplasm' in ds) or
                   ('Whole_Cell' in ds)) and (not 'Minus' in ds)]

    if speedrun:
        subset = [ds for ds in settings.datasets if ('Cytoplasm' in ds) and
                                                           (not 'Minus' in ds)]

    batch_key = 'cuff_super'
    dsets, super_3utr = super_falselength(settings, region, batch_key, subset,
                                          speedrun)

    return super_3utr

def join_antiexonic(exonic, anti, genome_dir):
    """ concatenate exonic with anti-exonic to get the whole genome
    """
    genome = {}
    for dset in exonic:
        outpath = os.path.join(genome_dir, dset+'genome_nonmerged')

        ## don't re-concatenate if you've done it already
        # XXX remove for debugging
        #if os.path.isfile(outpath):
            #genome[dset] = outpath
            #continue

        outhandle = open(outpath, 'wb')

        # write so that together it's for the whole genome
        for path in [exonic[dset], anti[dset]]:
            # save the things with val = nr of reads
            file_handle = open(path, 'rb')
            file_handle.next() # skip header
            for line in file_handle:
                (chrm, beg, end, utr_ID, strand, polyA_coordinate,
                 polyA_coordinate_strand, polyA_average_composition,
                 annotated_polyA_distance, nearby_PAS, PAS_distance,
                 number_supporting_reads, number_unique_supporting_reads,
                 unique_reads_spread) = line.split('\t')

                # write a bed format to disc
                outhandle.write('\t'.join([chrm,
                                          polyA_coordinate,
                                          str(int(polyA_coordinate)+1),
                                          nearby_PAS,
                                          number_supporting_reads,
                                          polyA_coordinate_strand +'\n']))
        outhandle.close()

        genome[dset] = outpath

    return genome

def merge_polyAs(settings, min_covr, minus, subset):

    # 1) go through all datasets, and concatenate the onlypolya files for the two
    # regions (they will be more or less independent, a few k overlap, fix later
    # if necessary)

    exonic = settings.only_files('3UTR-exonic')
    anti = settings.only_files('anti-3UTR-exonic')

    if subset == 'All_Cell_Lines':
        if minus:
            exonic = dict((k,v) for k,v in exonic.items() if 'Minus' in k)
            anti = dict((k,v) for k,v in anti.items() if 'Minus' in k)
        else:
            exonic = dict((k,v) for k,v in exonic.items() if 'Minus' not in k)
            anti = dict((k,v) for k,v in anti.items() if 'Minus' not in k)
    else:
        if minus:
            exonic = dict((k,v) for k,v in exonic.items() if ('Minus' in k) and
                          (subset in k))
            anti = dict((k,v) for k,v in anti.items() if ('Minus' in k) and
                        (subset in k))
        else:
            exonic = dict((k,v) for k,v in exonic.items() if ('Minus' not in k) and
                        (subset in k))
            anti = dict((k,v) for k,v in anti.items() if ('Minus' not in k) and
                        (subset in k))


    genome_dir = os.path.join(settings.here, 'genome_wide_dir')
    genome = join_antiexonic(exonic, anti, genome_dir)

    # merge and center each of the dsets and keep the path of the merged copy
    genomem = merge_center(genome, settings)

    # 2) Extend and merge all poly(A) files in genomem, summing the number of
    # annotated poly(A) files
    ext, firstpath = genomem.items()[0]

    # the file you will keep expanding and subtracting
    extendpath = os.path.join(os.path.dirname(firstpath), 'extendify')
    tempextendpath = os.path.join(os.path.dirname(firstpath), 'tempextendify')

    # XXX REMOVE IF YOU START CHANGIGN STUFF
    #if os.path.isfile(extendpath):
        #return extendpath

    # append to the extend-path
    ext_handle = open(extendpath, 'wb')

    # i) concatenate the files
    for dset, dpath in genomem.items():
        ext_handle.write(open(dpath, 'rb').read())

    ext_handle.close()

    # ii) expand the entries
    cmd = ['slopBed', '-i', extendpath, '-g', settings.hg19_path, '-b', '10']
    p = Popen(cmd, stdout=open(tempextendpath, 'wb'))
    p.wait()

    # iii) merge and sum scores
    cmd = ['mergeBed','-s', '-nms', '-scores', 'sum', '-i', tempextendpath]
    #cmd = ['mergeBed','-s', '-nms', '-scores', 'max', '-i', tempextendpath]
    p = Popen(cmd, stdout=PIPE)

    # iv) center the new transcript
    out_handle = open(extendpath, 'wb')
    for line in p.stdout:
        (chrm, beg, end, pases, maxcovrg, strand) = line.split()

        meanbeg = int((int(beg)+int(end))/2)

        ps = pases.split(';')
        ps2 = [pa.split('#') for pa in ps if pa != 'NA']
        if ps2 == []:
            pas = 'NA'
        else:
            uniq = list(set(sum(ps2, [])))
            pas = '#'.join(uniq)

        # filter on coverage
        if int(float(maxcovrg)) > min_covr:

            out_handle.write('\t'.join([chrm,
                                       str(meanbeg),
                                       str(meanbeg+1),
                                       pas,
                                       maxcovrg,
                                       strand+'\n']))
    out_handle.close()

    return extendpath

def get_gencode_models(model):
    """
    Get the gencode models. Fast.
    """
    transcripts = dict()
    genes = set([])
    exons = dict()

    # skip the first few lines
    gencode_handle = open(model, 'rb')
    for skipme in range(5):
        gencode_handle.next()

    t2 = time.time()
    for line_nr, line in enumerate(gencode_handle):
        (chrm, d, feature_type, beg, end, d, strand, d,d, gene_id,d, transcript_id,
         d, gene_type,d,d,d,d,d,ts_type) = line.split()[:20]

        # skip non-exon ones
        if feature_type != 'exon':
            continue

        transcript_id = transcript_id.rstrip(';').strip('"')
        gene_id = gene_id.rstrip(';').strip('"')
        ts_type = ts_type.rstrip(';').strip('"')

        # beg_transcript_ID_exonNr
        Utail_ID = '_'.join([beg, transcript_id])

        exons[Utail_ID] = (chrm, beg, end, transcript_id, gene_id)

        genes.add(transcript_id)

        if transcript_id in transcripts:
            transcripts[transcript_id][0] +=1 # nr of exons for this transcript
        else:
            transcripts[transcript_id] = [1, strand, gene_id, ts_type]

    print 'Made exon and gene dict: {0} seconds'.format(time.time()-t2)
    print 'Nr of exons: {0}'.format(len(exons))
    print 'Nr of genes: {0}'.format(len(genes))

    return genes, transcripts, exons

def get_gentrex(novel_ext, model):
    """
    transcripts[transcript_id] = (1/2, nr_exons, strand)
    1/2 is for novel or exteded transcript

    genes[gene_id] = set([transcript_id])

    exons[Utail_ID] = [beg, end, exon_nr, transcript_id, gene_id, []]
    the empty list is for the polyA sites with coverage

    """

    transcripts = dict()

    t1 = time.time()
    # this doesn't seem to cover all transcripts
    # (also it doesn't take up a lot of memory or time)
    for line_nr, line in enumerate(open(novel_ext, 'rb')):
        ts_id, novex = line.split()
        transcripts[ts_id] = [novex, 0, 0, 0]

    print 'Made transcript dict: {0} seconds'.format(time.time()-t1)
    print 'Nr of transcripts: {0}'.format(len(transcripts))
    # gene_id -> transcritpt ID
    genes = dict()
    exons = dict()

    t2 = time.time()
    # store the exons in an exon dict
    for line_nr, line in enumerate(open(model, 'rb')):
        (chrm, d,d, beg, end, d, strand, d,d,d,d, transcript_id, d, exon_nr)\
                = line.split()[:14]
        transcript_id = transcript_id.rstrip(';').strip('"')
        gene_id = '.'.join(transcript_id.split('.')[:-1])
        exon_nr = exon_nr.rstrip(';').strip('"')

        # beg_transcript_ID_exonNr
        Utail_ID = '_'.join([beg, transcript_id, exon_nr])
        #exons[Utail_ID] = Cufflink_exon(chrm, beg, end, strand, transcript_id,
                                        #gene_id, exon_nr, Utail_ID)
        exons[Utail_ID] = [chrm, beg, end, exon_nr, transcript_id, gene_id, []]

        if gene_id in genes:
            if not transcript_id in genes[gene_id]:
                genes[gene_id].add(transcript_id)
        else:
            genes[gene_id] = set([transcript_id])

        if transcript_id in transcripts:
            transcripts[transcript_id][1] +=1 # nr of exons for this transcript
            transcripts[transcript_id][2] = strand
            transcripts[transcript_id][3] = gene_id
        else:
            transcripts[transcript_id] = [0, 1, strand, gene_id] # 0 for not in the txt
            pass

    print 'Made exon and gene dict: {0} seconds'.format(time.time()-t2)
    print 'Nr of exons: {0}'.format(len(exons))
    print 'Nr of genes: {0}'.format(len(genes))

    return genes, transcripts, exons

def write_gencode_ends(settings, transcripts):

    gencends = os.path.join(settings.here, 'source_bedfiles', 'gencode_ends.bed')

    outhandle = open(gencends, 'wb')

    total_gencode = {'all': set([])}

    for ts_id, transcript in transcripts.iteritems():

        # coutn the total number
        total_gencode['all'].add(ts_id)

        if transcript.t_type in total_gencode:
            total_gencode[transcript.t_type].add(ts_id)
        else:
            total_gencode[transcript.t_type] = set([ts_id])

        # the 3' end depends on the strand
        if transcript.strand == '-':
            point = transcript.exons[0][1]
        if transcript.strand == '+':
            point = transcript.exons[-1][2]

        outhandle.write('\t'.join([transcript.chrm,
                                   str(point-40),
                                   str(point+40),
                                   transcript.t_type,
                                   transcript.ts_id,
                                   transcript.strand + '\n']))
    outhandle.close()

    return gencends, total_gencode

def gencode_cufflinks_report(settings, subset):
    """
    Same as for cufflinks, but for gencode
    """

    # 1) and 2) merge all polyA files (not poly(A) minus)
    min_covr = 1
    minus = True
    #minus = False
    merged_polyA_path = merge_polyAs(settings, min_covr, minus, subset)

    juncfree_pA_path = merged_polyA_path + '_juncfree'

    # remove areas around exon junctions, because it causes biases
    junctions = os.path.join(settings.here, 'junctions',
                             'splice_junctions_merged.bed')

    # i.5) cut away the exon-exon noise
    cmd = ['subtractBed', '-a', merged_polyA_path, '-b', junctions]
    p = Popen(cmd, stdout=open(juncfree_pA_path, 'wb'))
    p.wait()

    merged_polyA_path = juncfree_pA_path

    #chr1 = True
    chr1 = False

    model = '/users/rg/jskancke/phdproject/3UTR/gencode7/gencode7_annotation.gtf'
    if chr1:
        model = '/users/rg/jskancke/phdproject/3UTR/'\
                'gencode7/gencode7_annotation_chr1.gtf'

    an_frmt = 'GENCODE'
    # get dicts for the cufflinks model
    import annotation_parser as annparse

    (transcripts, genes) = annparse.make_transcripts(model, an_frmt)

    gencends, total_gencode = write_gencode_ends(settings, transcripts)

    # count how many of each type you found
    gencode_nrs = dict((key, len(v)) for key, v in total_gencode.items())

    # 3) intersect the bed_model with the merged polyA sites
    cmd = ['intersectBed', '-wo', '-a', gencends, '-b', merged_polyA_path]
    p = Popen(cmd, stdout=PIPE)

    goodpas = set(['AATAAA', 'ATTAAA'])

    allpas = set(['AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA',
                 'AATACA', 'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA',
                 'AATAGA'])

    # a dictionary that counts the number of transcripts and if there is a PAS
    # or not. When a new transcript type is identified, a new dictionary is made
    # for it
    subdict = {'PAS': set([]), 'Good PAS': set([]), 'No PAS': set([]), 'Cov':[]}

    event_counter = {'all': deepcopy(subdict) }

    for line in p.stdout:
        (hrm, beg, end, ts_type, ts_id, strand, d,d,d, pas, cov, d,d) = line.split()


        pases = pas.split('#')

        # add a dict for this transcript type
        if ts_type in event_counter:
            pass
        else:
            event_counter[ts_type] = deepcopy(subdict)

        event_counter['all']['Cov'].append(cov)
        event_counter[ts_type]['Cov'].append(cov)

        # check for PAS presence
        has_PAS = False
        has_good_PAS = False
        for pa in pases:
            if pa in allpas:
                has_PAS = True
            if pa in goodpas:
                has_good_PAS = True

        if has_PAS:
            event_counter['all']['PAS'].add(ts_id)
            event_counter[ts_type]['PAS'].add(ts_id)

        if has_good_PAS:
            event_counter['all']['Good PAS'].add(ts_id)
            event_counter[ts_type]['Good PAS'].add(ts_id)

        if not has_PAS:
            event_counter['all']['No PAS'].add(ts_id)
            event_counter[ts_type]['No PAS'].add(ts_id)

    #### summarize the above $
    event_sum = {}
    for key, subdict in event_counter.items():
        event_sum[key] = {}
        for subkey, ts_set in subdict.items():
            if subkey == 'Cov':
                event_sum[key][subkey] = np.mean([float(s) for s in ts_set])
            else:
                event_sum[key][subkey] = len(ts_set)

    # use gencode_nrs and event_sum to print out neat statistics about this now
    output_dir = os.path.join(settings.here, 'Results_and_figures', 'GENCODE_report',
                                  'csv_files')

    output_path = 'gencode_transcript_types_{0}_{1}_Minus{2}\
                          '.format(min_covr, subset, str(minus))

    gencsum_handle = open(os.path.join(output_dir, output_path), 'wb')

    gencsum_handle.write('Type\t#found\t% of total\twith PAS'\
                  '\twith good PAS\tmean_covr\n')

    # sorted by how many you find.
    def sorthelp(tupadict):
        adict = tupadict[1]
        return adict['PAS'] + adict['No PAS']

    key_tups = sorted(event_sum.items(), key=sorthelp, reverse=True)
    for key in (k[0] for k in key_tups):

        # skip those for which we don't find any
        if key not in event_sum:
            print('Not found: {0}'.format(key))
            continue

        nr_found = event_sum[key]['PAS'] + event_sum[key]['No PAS']
        pcnt_of_an = nr_found/gencode_nrs[key]

        pcnt_PAS = event_sum[key]['PAS']/nr_found
        pcnt_Good_PAS = event_sum[key]['Good PAS']/nr_found

        mean_cov =event_sum[key]['Cov']

        gencsum_handle.write('{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4:.2f}\t{5}\n'.format(key,
                                                           nr_found,
                                                           pcnt_of_an,
                                                           pcnt_PAS,
                                                           pcnt_Good_PAS,
                                                            mean_cov))
    gencsum_handle.close()


def new_cufflinks_report(settings, speedrun=False):
    """ Basic statistics on the cufflinks data.
    1) How many p(A)? How many overlap annotated?
    2) How many novel?
    3) How many fall very close to 3' ends?
    4) We verify that XXX of these genes are polyadenylated (min 1 transcript)
    5)

    1) Get poly(A)s for the whole genome for all the cell lines (cat those
    for anti-3UTR and 3UTR-exonic?)
    2) For all those files, merge using bed tools. Should be faster than your
    own. With a for-loop it will be quick. (Actually you could have implemented
    this in your own loops no? the thing was that clustering was never a cpu and
    memory hog until I got overlapping exons so damn many of them)
    3) From the cufflinks models, extract +/- 40 nt around the last exons
    4) Intersect what's around the last exons with the polyA file

    # OK, done. Now roderic wants to know about
    1) PAS
    2) How this compares to GENCODE
    3) How the whole analysis goes if you use poly(A) minus instead
    """

    # 1) and 2) merge all polyA files (not poly(A) minus)
    min_covr = 1
    minus = False
    subset = 'All_Cell_Lines'
    merged_polyA_path = merge_polyAs(settings, min_covr, minus, subset)

    juncfree_pA_path = merged_polyA_path + '_juncfree'

    # remove areas around exon junctions, because it causes biases
    junctions = os.path.join(settings.here, 'junctions',
                             'splice_junctions_merged.bed')

    # i.5) cut away the exon-exon noise
    cmd = ['subtractBed', '-a', merged_polyA_path, '-b', junctions]
    p = Popen(cmd, stdout=open(juncfree_pA_path, 'wb'))
    p.wait()

    merged_polyA_path = juncfree_pA_path
    debug()

    cuff_dict = '/users/rg/jskancke/phdproject/3UTR/CUFF_LINKS'
    model = cuff_dict + '/CuffCode_CSHL_exons_sans_GENCODE7.gtf'
    #bed_model = os.path.join(settings.here, 'source_bedfiles', 'cufflinks_exons.bed')

    novel_ext = cuff_dict + '/novel_or_extended.txt'
    # get dicts for the cufflinks model
    genes, transcripts, exons = get_gentrex(novel_ext, model)

    cuffends = os.path.join(settings.here, 'source_bedfiles', 'cufflinks_ends.bed')
    outhandle = open(cuffends, 'wb')

    for ex_id, (chrm, beg, end, exon_nr, ts_id, gene_id, lst) in exons.iteritems():

        # get information from the transcript
        nov_ext, ts_exon_nr, ts_strand, gene_id = transcripts[ts_id]

        if int(exon_nr) == ts_exon_nr:
            if ts_strand == '+':

                beg = str(int(end)-40)
                end = str(int(end)+40)

                outhandle.write('\t'.join([chrm,
                                           beg,
                                           end,
                                           ts_id,
                                           str(nov_ext),
                                           ts_strand+'\n']))
        if int(exon_nr) == 0:
            if ts_strand == '-':

                beg = str(int(beg)-40)
                end = str(int(beg)+40)

                outhandle.write('\t'.join([chrm,
                                           beg,
                                           end,
                                           ts_id,
                                           str(nov_ext),
                                           ts_strand+'\n']))
    outhandle.close()

    # 3) intersect the bed_model with the merged polyA sites
    cmd = ['intersectBed', '-wo', '-a', cuffends, '-b', merged_polyA_path]
    p = Popen(cmd, stdout=PIPE)

    novel_transcripts_with_pA = set([])
    novel_transcripts_with_PAS = set([])
    novel_transcripts_with_GoodPAS = set([])

    extended_transcripts_with_pA = set([])
    extended_transcripts_with_PAS = set([])
    extended_transcripts_with_GoodPAS = set([])

    nothing_transcripts_with_pA = set([])
    nothing_transcripts_with_PAS = set([])
    nothing_transcripts_with_GoodPAS = set([])

    genes_with_pA = set([])
    genes_with_PAS = set([])
    genes_with_GoodPAS = set([])

    goodpas = set(['AATAAA', 'ATTAAA'])
    allpas = set(['AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA',
                 'AATACA', 'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA',
                 'AATAGA'])
    # get the pa clusters represented by their coverage
    # extended could mean novel 3UTRs ! :)
    annotated = open(os.path.join(settings.here, 'junctions', 'with_pA'))
    # TODO print out the polyA reads that overlap cuffilnks exons and show the
    # buggers what it looks like in the browser.

    for line in p.stdout:
        (chrm, cbeg, cend, ts_id, nov_ext, strand,
         d, d, d, pas, covrg, d,d) = line.split()

        has_PAS = False
        has_GOOD_PAS = False
        pases = pas.split('#')
        for pa in pases:
            if pa in allpas:
                has_PAS = True
            if pa in goodpas:
                has_GOOD_PAS = True

        genes_with_pA.add(transcripts[ts_id][3])

        if has_PAS:
            genes_with_PAS.add(transcripts[ts_id][3])

            if has_GOOD_PAS:
                genes_with_GoodPAS.add(transcripts[ts_id][3])

        if int(nov_ext) == 1:
            novel_transcripts_with_pA.add(ts_id)
            if has_PAS:
                novel_transcripts_with_PAS.add(ts_id)

                if has_GOOD_PAS:
                    novel_transcripts_with_GoodPAS.add(ts_id)

        if int(nov_ext) == 2:
            extended_transcripts_with_pA.add(ts_id)
            if has_PAS:
                extended_transcripts_with_PAS.add(ts_id)

                if has_GOOD_PAS:
                    extended_transcripts_with_GoodPAS.add(ts_id)

        if int(nov_ext) == 0:
            nothing_transcripts_with_pA.add(ts_id)
            if has_PAS:
                nothing_transcripts_with_PAS.add(ts_id)

                if has_GOOD_PAS:
                    nothing_transcripts_with_GoodPAS.add(ts_id)

    #### summarize the above $
    nov_ts = len(novel_transcripts_with_pA)
    nov_PAS = len(novel_transcripts_with_PAS)
    nov_GoodPAS = len(novel_transcripts_with_GoodPAS)

    ext_ts = len(extended_transcripts_with_pA)
    ext_PAS = len(extended_transcripts_with_PAS)
    ext_GoodPAS = len(extended_transcripts_with_GoodPAS)

    not_ts = len(nothing_transcripts_with_pA)
    not_PAS = len(nothing_transcripts_with_PAS)
    not_GoodPAS = len(nothing_transcripts_with_GoodPAS)

    total_pAts = sum([nov_ts, ext_ts, not_ts])
    total_ts_PAS = sum([nov_PAS, ext_PAS, not_PAS])
    total_ts_GoodPAS = sum([nov_GoodPAS, ext_GoodPAS, not_GoodPAS])

    total_pAgenes = len(genes_with_pA)
    total_PASgenes = len(genes_with_PAS)
    total_GoodPASgenes = len(genes_with_GoodPAS)

    tot_ts = len(transcripts)
    tot_genes = len(genes)

    # TOTAL TRANSCRIPTS

    print('Total transcripts with pA at 3prime end:\t{0} ({1:.2f})'\
          .format(total_pAts, total_pAts/tot_ts))

    print('Total transcripts with PAS:\t{0} ({1:.2f})'\
          .format(total_ts_PAS, total_ts_PAS/total_pAts))

    print('Total transcripts with GOOD PAS:\t{0} ({1:.2f})'\
          .format(total_ts_GoodPAS, total_ts_GoodPAS/total_pAts))

    # NOVEL TRANSCRIPTS
    print('Novel transcripts with pA at 3prime end:\t{0} ({1:.2f})'\
          .format(nov_ts, nov_ts/total_pAts))

    print('Novel transcripts with PAS:\t{0} ({1:.2f})'\
          .format(nov_PAS, nov_PAS/nov_ts))

    print('Novel transcripts with PAS:\t{0} ({1:.2f})'\
          .format(nov_GoodPAS, nov_GoodPAS/nov_ts))

    # EXTENDED TRANSCRIPTS
    print('Extended transcripts with pA at 3prime end:\t{0} ({1:.2f})'\
          .format(ext_ts, ext_ts/total_pAts))

    print('Extended transcripts with PAS:\t{0} ({1:.2f})'\
          .format(ext_PAS, ext_PAS/ext_ts))

    print('Extended transcripts with PAS:\t{0} ({1:.2f})'\
          .format(ext_GoodPAS, ext_GoodPAS/ext_ts))

    # NOTHING TRANSCRIPTS
    print('Nothing transcripts with pA at 3prime end:\t{0} ({1:.2f})'\
          .format(not_ts, not_ts/total_pAts))

    # GENES
    print('Total genes with min 1 transcript with pA at 3prime end:\t{0} ({1})'\
          .format(total_pAgenes, total_pAgenes/tot_genes))

    print('PAS genes with min 1 transcript with pA at 3prime end:\t{0} ({1})'\
          .format(total_PASgenes, total_PASgenes/total_pAgenes))

    print('GOOD PAS genes with min 1 transcript with pA at 3prime end:\t{0} ({1})'\
          .format(total_GoodPASgenes, total_GoodPASgenes/total_pAgenes))


def merge_center(genome, settings):
    """
    Merge and center a set of bed-paths. Since these polya_counts come from the
    same sample, don't add the polyAs -- rather keep the max value
    """

    hg19 = settings.hg19_path
    genomem = {}

    for dset, dsetpath in genome.items():
        # 1) extend with 20 nt
        extended_path = dsetpath +'_extended'
        merged_path = extended_path + '_merged'

        # skip if you've done this before
        # XXX remove for debugging
        #if os.path.isfile(merged_path):
            #genomem[dset] = merged_path
            #continue

        # expand with 10
        cmd = ['slopBed', '-i', dsetpath, '-g', hg19, '-b', '10']
        p = Popen(cmd, stdout=open(extended_path, 'wb'))
        p.wait()

        # 2) merge 
        #cmd = ['mergeBed','-s', '-nms', '-scores', 'max', '-i', extended_path]
        cmd = ['mergeBed','-s', '-nms', '-scores', 'sum', '-i', extended_path]
        p = Popen(cmd, stdout=PIPE)

        out_handle = open(merged_path, 'wb')
        for line in p.stdout:
            (chrm, beg, end, pases, maxcovrg, strand) = line.split('\t')
            meanbeg = int((int(beg)+int(end))/2)
            strand = strand.rstrip()

            ps = pases.split(';')
            ps2 = [pa.split(' ') for pa in ps if pa != 'NA']
            if ps2 == []:
                pas = 'NA'
            else:
                uniq = list(set(sum(ps2, [])))
                pas = '#'.join(uniq)

            out_handle.write('\t'.join([chrm,
                                       str(meanbeg),
                                       str(meanbeg+1),
                                       pas,
                                       maxcovrg,
                                       strand+'\n']))
        out_handle.close()

        genomem[dset] = merged_path

    return genomem

def join_pAs(genome, genome_dir):

    total_pAs = {}
    for key in ['K562', 'GM12878']:
        thesefiles = [flpath for fl, flpath in genome.items() if key in fl]
        outfile = os.path.join(genome_dir, key+'_joined')

        # write the 15-nt-expanded things to file
        with open(outfile, 'wb') as outhandle:
            for td in thesefiles:
                for line in open(td, 'rb'):
                    outhandle.write(line)

        total_pAs[key] = outfile

    return total_pAs

def expand_pAs(total_pAs, genome_dir, settings):

    expanded_pAs = {}
    for key, dsetpath in total_pAs.items():
        expanded_path = os.path.join(genome_dir, dsetpath+'_expanded')
        cmd = ['slopBed', '-i', dsetpath, '-g', settings.hg19_path,
               '-b', '15']
        p = Popen(cmd, stdout=open(expanded_path, 'wb'))
        p.wait()
        expanded_pAs[key] = expanded_path

    return expanded_pAs

def merge_pAs(expanded_pAs, genome_dir):

    merged_pAs = {}
    for key, dsetpath in expanded_pAs.items():
        cmd = ['mergeBed', '-s', '-scores', 'max', '-i', dsetpath]
        p = Popen(cmd, stdout=PIPE)
        merged_path = os.path.join(genome_dir, key+'_mergedCentered')

        with open(merged_path, 'wb') as out_handle:
            for line in p.stdout:
                (chrm, beg, end, maxcovrg, strand) = line.split()
                meanbeg = int((int(beg)+int(end))/2)

                out_handle.write('\t'.join([chrm,
                                           str(meanbeg),
                                           str(meanbeg+1),
                                           '0',
                                           maxcovrg,
                                           strand+'\n']))

        merged_pAs[key] =  merged_path

    return merged_pAs

def sort_pAs(merged_pAs, genome_dir):

    sorted_pAs = {}
    for key, dsetpath in merged_pAs.items():
        sorted_path = os.path.join(genome_dir, key+'_mergedCentered_sorted')
        cmd = 'sort -k1.4,1.5 -k2,2 -n {0}'.format(dsetpath)
        p = Popen(cmd, shell=True, stdout=open(sorted_path, 'wb'))
        p.wait()

        sorted_pAs[key] = sorted_path

    return sorted_pAs

def distansiate_pAs(sorted_pAs, limit, genome_dir):

    distant_pAs = {}
    for dsetkey, dsetpath in sorted_pAs.items():

        discarders = set([]) # list of files to throw away

        sorted_handle = open(dsetpath, 'rb')

        # skip the first site: but save its key in case of throwing away
        first = sorted_handle.next()
        (chrm, beg, end) = first.split()[:3]
        key = ''.join([chrm, beg, end])

        for line in sorted_handle:

            (nchrm, nbeg, nend) = line.split()[:3]
            nkey = '_'.join([nchrm, nbeg, nend])

            if nchrm == chrm:
                dist = int(nbeg) - int(beg)
                if dist < limit:
                    discarders.add(key)
                    discarders.add(nkey)
                if dist < 0:
                    print 'error dist '

                # reset the coordinates
                (chrm, beg, end) = (nchrm, nbeg, nend)
                key = nkey

            # if new chromosome, don't compare
            else:
                pass
                # reset the coordinates
                (chrm, beg, end) = (nchrm, nbeg, nend)
                key = nkey

        # go through the file again, this time discarding the ones you marked
        # before

        distant_path = os.path.join(genome_dir, dsetkey+'_distant{0}'.format(limit))

        print dsetkey
        print len(discarders)

        disc = open(os.path.join(genome_dir, dsetkey+'discard'), 'wb')
        with open(distant_path, 'wb') as out_handle:
            for line in open(dsetpath, 'rb'):
                (chrm, beg, end) = line.split()[:3]
                key = '_'.join([chrm, beg, end])

                # if not being discarded, write
                #debug()
                if key in discarders:
                    disc.write(line)
                else:
                    out_handle.write(line)

        distant_pAs[dsetkey] = distant_path
        disc.close()

    return distant_pAs

def get_annot_polyA(genome_dir):
    class myclass(object):
        pass

    someset = myclass()
    someset.annotation_path =\
            '/users/rg/jskancke/phdproject/3UTR/Annotations/'\
            'gencode.v3c.annotation.GRCh37.gtf'
    someset.annotation_format = 'GENCODE'
    someset.min_utrlen = 200
    import annotation_parser as annpars

    # annotated polyA sites
    annot_path = os.path.join(genome_dir, 'vc3_annot_polyas')
    # get annotatd polyasites at this site, from vc3
    annpars.get_a_polyA_sites_bed(someset, annot_path)

    return annot_path


def hagen(settings, speedrun):
    """
    1) Gather up the polyA files for K562 and GM separately.
    2) Expand and merge them, keeping the PAS and coverage information.
    3.5) sort the poly(A) sites and mark all that are 200 nt within each other,
    and cut them out of the list by grep or some way
    3) Gather the regions around annotated gencode v3 sites
    4) Intersect the wide-apart poly(A) sites from K and G with the annotated
    sites. Give the anotated sites a name, an ID so you can follow them.
    5) Read the intersection of both files into a dictionary keyed by the
    annotated_site's ID. FIrst do one. Then other. If it's already found in
    other, then both are there. The dict must know about the PAS and the
    coverage for both two. Then you can run through the list and give both
    things to HAgen that he wants.
    """

    exonic = settings.only_files('3UTR-exonic')
    anti = settings.only_files('anti-3UTR-exonic')
    # skip the minus ones
    exonic = dict((k,v) for k,v in exonic.items() if ('Minus' not in k) and (
                  'K562' in k or 'GM12878' in k))

    anti = dict((k,v) for k,v in anti.items() if ('Minus' not in k) and (
                  'K562' in k or 'GM12878' in k))

    genome_dir = os.path.join(settings.here, 'hagen_stuff')
    ## join the 3utr and anti-3utr regions
    genome = join_antiexonic(exonic, anti, genome_dir)
    # TODO by not merging without sum first, you are artificially increasing the
    # count for those that are on the borders of 3utr/non-3utr.
    # it's ok, i fixed it. now re-run.

    ## Join all subset files 
    total_pAs = join_pAs(genome, genome_dir)

    ## Expand
    expanded_pAs = expand_pAs(total_pAs, genome_dir, settings)

    ## Merge
    merged_pAs = merge_pAs(expanded_pAs, genome_dir)

    ## Sort
    sorted_pAs = sort_pAs(merged_pAs, genome_dir)

    ## For both dsets, screen away those that are closer than 100/300 nt.
    limit = 300
    distant_pAs = distansiate_pAs(sorted_pAs, limit, genome_dir)

    ## 3) Retrieve the poly(A) sites from gencode version 3.
    annot_path = get_annot_polyA(genome_dir)

    distant_pAs = {'K562': os.path.join(genome_dir, 'K562_distant{0}'.format(limit)),
                   'GM12878': os.path.join(genome_dir,
                                           'GM12878_distant{0}'.format(limit))}

    annot_path = os.path.join(genome_dir, 'vc3_annot_polyas')

    hagen_path = os.path.join(genome_dir, 'hagens_ends.bed')

    # 4) Intersect the polyA sites for both cell types with the gencode

    annot_dict = {}

    # 4.1 Intersect the annotated sites with the (expanded) discoverd poly(A) sites
    for dsetkey, dsetpath in distant_pAs.items():
    #for dsetkey, dsetpath in merged_pAs.items():

        # 4.2) Expand the distant/merged sites
        expanded_dsetpath = os.path.join(genome_dir, dsetpath+'_expanded')

        cmd = ['slopBed', '-i', dsetpath, '-g', settings.hg19_path, '-b', '15']
        p = Popen(cmd, stdout=open(expanded_dsetpath, 'wb'))
        p.wait()

        cmd = ['intersectBed', '-wa', '-wb', '-a', expanded_dsetpath, '-b', hagen_path]

        #cmd = ['intersectBed','-s','-wa','-wb',
               #'-a', annot_path, '-b', expanded_dsetpath]

        found = set([])

        p = Popen(cmd, stdout=PIPE)

        for line in p.stdout:
            (achrm, abeg, aend, d, covrg, astrand, chrm, beg, end, strnd)\
                    = line.split()

            #if astrand == '+':
                #cent = abeg
            #if astrand == '-':
                #cent = aend

            #centbeg = int((int(abeg)+int(aend))/2)
            #centend = centbeg+1

            #key = '_'.join([achrm, cent, cent, astrand])
            key = '_'.join([chrm, beg, end, strnd])
            # only unique
            if key not in found:
                found.add(key)
            else:
                continue

            if key not in annot_dict:
                annot_dict[key] = [(dsetkey, int(float(covrg)))]

            else:
                annot_dict[key].append((dsetkey, int(float(covrg))))

    # 5) Parse the dict at will
    big = 1
    hagenfile = open(os.path.join(genome_dir, 'All_GM_K562'\
                                   .format(big, limit)), 'wb')
    #hagenfile1 = open(os.path.join(genome_dir, 'GMbig_Kbig{0}_dist{1}'\
                                   #.format(big, limit)), 'wb')

    hagenfile2 = open(os.path.join(genome_dir, 'GMbig_Knone{0}_dist{1}'\
                                  .format(big, limit)), 'wb')

    hagenfile3 = open(os.path.join(genome_dir, 'GMnone_Kbig{0}_dist{1}'\
                                  .format(big, limit)), 'wb')

    for annot_key, polyAsite in annot_dict.iteritems():
        if len(polyAsite) == 2:
            if polyAsite[0][0] == 'GM12878':
                gmlev = polyAsite[0][1]
                k5lev = polyAsite[1][1]
            else:
                k5lev = polyAsite[0][1]
                gmlev = polyAsite[1][1]

            if gmlev > big and k5lev > big:
                (chrm, beg, end, strand) = annot_key.split('_')
                hagfrmt1 = '\t'.join([annot_key, str(k5lev), str(gmlev)])

                #hagformt = '_'.join
                hagenfile.write(hagfrmt1 + '\n')

        if len(polyAsite) == 1:
            if polyAsite[0][0] == 'GM12878' and polyAsite[0][1] > big:
                (chrm, beg, end, strand) = annot_key.split('_')
                gmlev = polyAsite[0][1]
                hagfrmt2 = '\t'.join([annot_key, '0', str(gmlev)])
                hagenfile.write(hagfrmt2 + '\n')

            if polyAsite[0][0] == 'K562' and polyAsite[0][1] > big:
                (chrm, beg, end, strand) = annot_key.split('_')
                k5lev = polyAsite[0][1]
                hagfrmt3 = '\t'.join([annot_key, str(k5lev), '0'])
                hagenfile.write(hagfrmt3 + '\n')

    #hagenfile1.close()
    hagenfile2.close()
    hagenfile3.close()

def pas_for_hagen(settings):

    exonic = settings.only_files('3UTR-exonic')
    anti = settings.only_files('anti-3UTR-exonic')
    # skip the minus ones
    exonic = dict((k,v) for k,v in exonic.items() if ('Minus' not in k) and (
                  'K562' in k or 'GM12878' in k))

    anti = dict((k,v) for k,v in anti.items() if ('Minus' not in k) and (
                  'K562' in k or 'GM12878' in k))

    genome_dir = os.path.join(settings.here, 'hagen_stuff')
    ## join the 3utr and anti-3utr regions
    genome = join_antiexonic(exonic, anti, genome_dir)
    # TODO by not merging without sum first, you are artificially increasing the
    # count for those that are on the borders of 3utr/non-3utr.
    # it's ok, i fixed it. now re-run.

    ## Join all subset files 
    total_pAs = join_pAs(genome, genome_dir)

    ## Expand
    expanded_pAs = expand_pAs(total_pAs, genome_dir, settings)

    ## Merge
    merged_pAs = merge_pAs(expanded_pAs, genome_dir)

    ## Sort
    sorted_pAs = sort_pAs(merged_pAs, genome_dir)

    ## For both dsets, screen away those that are closer than 100/300 nt.
    limit = 300
    distant_pAs = distansiate_pAs(sorted_pAs, limit, genome_dir)
    distant_pAs = {'K562': os.path.join(genome_dir, 'K562_distant{0}'.format(limit)),
                   'GM12878': os.path.join(genome_dir,
                                           'GM12878_distant{0}'.format(limit))}

    hagen_path = os.path.join(genome_dir, 'hagens_ends.bed')

    # 4) Intersect the polyA sites for both cell types with the gencode

    annot_dict = {}
    # 4.1 Intersect the annotated sites with the (expanded) discoverd poly(A) sites
    for dsetkey, dsetpath in distant_pAs.items():
    #for dsetkey, dsetpath in merged_pAs.items():

        # 4.2) Expand the distant/merged sites
        expanded_dsetpath = os.path.join(genome_dir, dsetpath+'_expanded')

        cmd = ['slopBed', '-i', dsetpath, '-g', settings.hg19_path, '-b', '15']
        p = Popen(cmd, stdout=open(expanded_dsetpath, 'wb'))
        p.wait()

        cmd = ['intersectBed', '-wa', '-wb', '-a', expanded_dsetpath, '-b', hagen_path]

        #cmd = ['intersectBed','-s','-wa','-wb',
               #'-a', annot_path, '-b', expanded_dsetpath]

        found = set([])

        p = Popen(cmd, stdout=PIPE)

        for line in p.stdout:
            (achrm, abeg, aend, d, covrg, astrand, chrm, beg, end, strnd)\
                    = line.split()

            #if astrand == '+':
                #cent = abeg
            #if astrand == '-':
                #cent = aend

            #centbeg = int((int(abeg)+int(aend))/2)
            #centend = centbeg+1

            #key = '_'.join([achrm, cent, cent, astrand])
            key = '_'.join([chrm, beg, end, strnd])
            # only unique
            if key not in found:
                found.add(key)
            else:
                continue

            if key not in annot_dict:
                annot_dict[key] = [(dsetkey, int(float(covrg)))]

            else:
                annot_dict[key].append((dsetkey, int(float(covrg))))

    return annot_dict

def get_ts_rpkms_genc3():
    """
    """
    carrie = '/users/rg/projects/NGS/Projects/ENCODE/hg19main/'\
            'GingerasCarrie'

    thefile = '/users/rg/projects/NGS/Projects/ENCODE/hg19main/'\
            'GingerasCarrie/comparison/Trans.Expression.ENCODE.txt'

    # parse the transcript rpkm file
    thehandle = open(thefile, 'rb')
    header = thehandle.next()
    header = header.split('\t')

    gm_dirs = ['Ging003C', 'Ging004C']
    k5_dirs = ['Ging001C', 'Ging002C']

    gm_nrs = [header.index(gm_d)+1 for gm_d in gm_dirs]
    k5_nrs = [header.index(k5_d)+1 for k5_d in k5_dirs]

    gm_nr1, gm_nr2 = gm_nrs
    k5_nr1, k5_nr2 = k5_nrs

    k56 = {}
    gm1 = {}

    for line in thehandle:
        flds = line.split()
        ts_id = flds[0]

        k56[ts_id] = np.mean([float(flds[i]) for i in k5_nrs if flds[i] != '0'])
        gm1[ts_id] = np.mean([float(flds[i]) for i in gm_nrs if flds[i] != '0'])

    return gm1, k56


def hag_rpkm_refac(settings):
    """
    New version.
    """
    model = '/users/rg/jskancke/phdproject/3UTR/Annotations/'\
            'gencode.v3c.annotation.GRCh37.gtf'

    # 1) Get a gm128[ts_id] = rpkm dict (for k5 too).
    (ts_2_rpkm_GM12, ts_2_rpkm_K56) = get_ts_rpkms_genc3()

def hagens_rpkms(settings):
    """
    Hagen wants poly(A)-site specific RPKMS.
    """
    model = '/users/rg/jskancke/phdproject/3UTR/Annotations/'\
            'gencode.v3c.annotation.GRCh37.gtf'

    chr1 = True
    chr1 = False
    if chr1:
        model = '/users/rg/jskancke/phdproject/3UTR/Annotations/'\
                'gencode.v3c.annotation.GRCh37_chr1.gtf'

    an_frmt = 'GENCODE'
    import annotation_parser as annparse

    #1) Get all distant pA sites that match hagen's 
    # Pickle timesaver
    pickfile = 'PICKME'
    if not os.path.isfile(pickfile):
        pA_sites_dict = pas_for_hagen(settings)
        pickle.dump(pA_sites_dict, open(pickfile, 'wb'))
    else:
        pA_sites_dict = pickle.load(open(pickfile, 'rb'))

    #2) get dicts for all sites and all exons

    #outdir = os.path.join(settings.here, 'analysis', 'transcripts')
    #outpath = os.path.join(outdir, 'genc_transcripts.bed')
    #outhandle = open(outpath, 'rb')

    #temp_file = 'temp_exons_my_model'
    #my_handle = open(temp_file, 'wb')
    (transcripts, genes) = annparse.make_transcripts(model, an_frmt)

    end_sites = {}
    all_exons = {}
    for ts_id, ts in transcripts.iteritems():
        # all exons
        for ex in ts.exons:
            ex = list(ex)
            ex[1] += 1 # 'your' exons have 'beg'-1 compared to pipeline
            key = '_'.join([str(c) for c in ex])
            all_exons[key] = ts_id
            #my_handle.write(key+'\n')

        # only the ends
        if ts.strand == '+':
            ex = list(ts.exons[-1])
            ex[1] += 1 # 'your' exons have 'beg'-1 compared to pipeline
            end_key = '_'.join([str(c) for c in ex])
        else:
            ex = list(ts.exons[0])
            ex[1] += 1 # 'your' exons have 'beg'-1 compared to pipeline
            end_key = '_'.join([str(c) for c in ex])

        end_sites[end_key] = ts_id
    #my_handle.close()

    #3) You need a dictionary of the average rpkm for each exon in each cell
    carrie = '/users/rg/projects/NGS/Projects/ENCODE/hg19main/'\
            'GingerasCarrie'

    my_cell_lines = ('K562', 'GM12878')

    # Pickle timesaver
    pickfile = 'PICKME_2'
    if not os.path.isfile(pickfile):
        rpkm_paths = get_genc3_rpkm_paths(settings, carrie, my_cell_lines)
        rpkm_dict = get_genc3_rpkms(rpkm_paths, my_cell_lines)
        pickle.dump(rpkm_dict, open(pickfile, 'wb'))
    else:
        rpkm_dict = pickle.load(open(pickfile, 'rb'))

    #temp_file = 'temp_exons_pipeline'
    #pipeline_handle = open(temp_file, 'wb')
    #for exon in rpkm_dict.values()[0].iterkeys():
        #pipeline_handle.write(exon+'\n')

    #pipeline_handle.close()

    # 4) You need a dictionary of the RPKMS of all transcripts
    transcript_rpkms = dict((ce_l, {}) for ce_l in my_cell_lines)

    for cell_l in my_cell_lines:

        for ts_id, ts in transcripts.iteritems():
            all_present = True
            rpkmlist = []
            for ex in ts.exons:
                ex = list(ex)
                ex[1] +=1 # adjust again...
                key = '_'.join([str(c) for c in ex])

                # connect the poly(A)
                if key not in rpkm_dict[cell_l]:
                    all_present = False
                    break
                else:
                    rpkmlist.append(rpkm_dict[cell_l][key])

            # if all present, calcluate rpkms for transcript
            if all_present:
                transcript_rpkms[cell_l][ts_id] = sum(rpkmlist)

    save_dir = os.path.join(settings.here, 'hagen_stuff')
    # 5) Each of hagen-agreeing poly(A) site must know which transcript it
    # belongs to:

    # 5.1) Write out end-sites to bedfile
    endfile = 'end_sites.bed'
    endpath = os.path.join(save_dir, endfile)
    endhandle = open(endpath, 'wb')
    for end_site, ts_id in end_sites.iteritems():
        (chrm, beg, end, strand) = end_site.split('_')
        endhandle.write('\t'.join([chrm, beg, end, '0', ts_id, strand])+'\n')
    endhandle.close()

    # 5.2) Write out poly(A) sites to bedfile
    pAfile = 'polyA_sites_forexonmerge.bed'
    pApath = os.path.join(save_dir, pAfile)
    pAhandle = open(pApath, 'wb')
    for end_site, cell_line_count in pA_sites_dict.iteritems():
        (chrm, beg, end, strand) = end_site.split('_')
        pAhandle.write('\t'.join([chrm, beg, end, strand])+'\n')
    pAhandle.close()

    # 5.3) Expand the end-sites
    expandpath = add_ending(endpath, 'expanded')
    cmd1 = ['slopBed', '-b', '15', '-i', endpath, '-g', settings.hg19_path]
    p1 = Popen(cmd1, stdout = open(expandpath, 'wb'))
    p1.wait()

    # 5.4) Intersect
    # NOTE there is some kind of bug with intersectBed and the -s option.
    # Solution: run without
    #cmd2 = ['intersectBed', '-s', '-wa', '-wb', '-a', expandpath, '-b', pApath]
    cmd2 = ['intersectBed', '-wa', '-wb', '-a', expandpath, '-b', pApath]
    p2 = Popen(cmd2, stdout=PIPE)
    # 5.5) Read into a dict those poly(A) sites that intersect and the
    polyA2ts = {}
    for line in p2.stdout:
        (chra, bega, enda, d, ts_id, standa, chrb,
         begb, endb, strandb) = line.split()

        key = '_'.join([chrb, begb, endb, strandb])
        polyA2ts[key] = ts_id

    cl_rpkm_name = 'polyAsiteRPKMS_transcript'
    outpath = os.path.join(save_dir, cl_rpkm_name)
    outhandle = open(outpath, 'wb')
    keepers = 0
    missers = 0
    for polyA_site, ts_id in polyA2ts.iteritems():

        if (ts_id in transcript_rpkms['GM12878'] and
            ts_id in transcript_rpkms['K562']):
            keepers += 1
            gm_rpkm_ts = transcript_rpkms['GM12878'][ts_id]
            k5_rpkm_ts = transcript_rpkms['K562'][ts_id]

            outhandle.write('\t'.join([polyA_site, str(gm_rpkm_ts),
                                       str(k5_rpkm_ts)]) + '\n')
        else:
            missers +=1

    outhandle.close()

    print('keepers {0}'.format(keepers))
    print('missers {0}'.format(missers))

    # 5.6.1) Look up the RPKMs of the transcripts for the poly(A) sites
            #for each cell line.

    # 5.6.2) Loop up the RPKM of the terminal exon for the transcript for
            # the poly(A) sites for each cell line

def add_ending(filepath, ending):
    """
    with ending = '_temp'
    turn /var/ost.py
    into /var/ost_temp.py
    """
    mainpart, suffix = os.path.splitext(filepath)
    return mainpart+ending+suffix


def get_genc3_rpkms(rpkm_paths, my_cell_lines):
    """ Save the rpkms and average them across cell lines
    """
    rpkm_dict = dict((ce_l, {}) for ce_l in my_cell_lines)

    for cell_l, paths in rpkm_paths.items():
        for path in paths:
            cmd = ['zcat', path]
            p = Popen(cmd, stdout=PIPE)
            for line in p.stdout:
                coord, rpkm, d = line.split()
                (chrm, beg, end, strnd_info) = coord.split('_')
                if strnd_info == '-1':
                    coord = '_'.join([chrm, beg, end, '-'])
                else:
                    coord = '_'.join([chrm, beg, end, '+'])

                if coord in rpkm_dict[cell_l]:
                    rpkm_dict[cell_l][coord].append(rpkm)
                else:
                    rpkm_dict[cell_l][coord] = [rpkm]

    mean_rpkm_dict = dict((ce_l, {}) for ce_l in my_cell_lines)
    for cell_l, coord_dict in rpkm_dict.iteritems():
        for coord, rpkms in coord_dict.iteritems():
            mean_rpkm_dict[cell_l][coord] = np.mean([float(r) for r in rpkms])

    return mean_rpkm_dict

def get_genc3_rpkm_paths(settings, carrie, my_cell_lines):
    """
    Parse through gencode3 exon rpkms and put their paths in your dict
    """
    #0) make a directory -> cell line dictionary
    dir_2_cell_line = {}
    for line in open(os.path.join(carrie, 'Dataset.info.txt')):
        (dirinfo, cell_l, comp, longshort, rplica) = line.split()
        dirname = dirinfo[4:]
        dir_2_cell_line[dirname] = cell_l

    #1 collect file paths for nuclear and cytoplasmic exons
    #paths = AutoVivification()
    paths = dict((ce_l, []) for ce_l in my_cell_lines)

    for dirname, dirnames, filenames in os.walk(carrie):
        if dirname.endswith('exons'):

            main_dir = os.path.split(os.path.split(dirname)[0])[1]

            # include only those cell lines you have been using
            if main_dir not in dir_2_cell_line:
                print main_dir
                continue

            cell_line = dir_2_cell_line[main_dir]
            if cell_line not in my_cell_lines:
                continue

            compartment = main_dir[3:]

            # we're only looking at nuclear and cytoplasm
            if compartment not in ['N','C', 'WC']:
                continue

            for f in filenames:
                if f == 'All.exon.rpkm.pooled.txt.gz':

                    mypath = os.path.join(dirname, f)

                    paths[cell_line].append(mypath)

    return paths

def main():
    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    # Directory paths for figures and where the output lies
    (savedir, outputdir) = [os.path.join(here, d) for d in ('figures', 'output')]

    # Speedruns with chr1
    chr1 = False
    #chr1 = True

    # Read UTR_SETTINGS (there are two -- for two different annotations!)
    settings = Settings(os.path.join(here, 'UTR_SETTINGS'), savedir, outputdir,
                        here, chr1)

    #gencode_report(settings, speedrun=False)

    # XXX cufflinks report
    #new_cufflinks_report(settings, speedrun=False)

    # XXX same as cufflnksm but for gencode
    #cell_lines = ['All_Cell_Lines', 'GM12878', 'HEPG2', 'HUVEC', 'HeLa-S3',
                  #'K562']
    #for subset in cell_lines:
        #gencode_cufflinks_report(settings, subset)

    # Hagen's stuff
    #hagen(settings, speedrun=False)
    #hagens_rpkms(settings)
    hag_rpkm_refac(settings)

    # XXX For the length + polyA (not only) output files:
    # NOTE the difference is that with this one you have more information about
    # the UTR object (RPKM etc) and consequently it will take longer time to run
    # dsets, super_3utr = get_utrs(settings, speedrun=False)

    # for the onlypolyA files
    # you now have an additional level: the "region" (often the 3UTR)
    #pickfile = 'super_pickle'

    # TODO to make this faster and consume less memory, don't save whole
    # clusters. save only the information you need. add as necessary. when
    # you'll have more of these it will be too much.

    # Get the basic poly(A) stat file that is output with each run. It gives
    # statistics on the number of reads that fall in which strand and so forth.
    #polyAstats = get_polyA_stats(settings)

    ##Basic poly(A) read statistics for the datasetes
    #for side in ['annotated', 'opposite']:
        #polyA_summary(dsets, super_3utr, polyAstats, settings, side)

    #### RPKM distribution
    #### What is the distribution of 3UTR RPKMS? run this function and find out!
    #rpkm_dist(dsets)

    #### How does your WC = C + N model hold up? Plot relationships by RPMKM
    #wc_is_n_plus_c(dsets) # RESULT: max-min(WC,N,C)-dist: 30 +/- 66 nt
    # However, this is only valid for the rpkm-limited ones...

    #### Extend beyond distribution and number
    #### How many 3UTRs extend significantly? beyond the aTTS
    #### NOTE that, with current code, they must be checked by hand, as the
    #### extended area can overlap a genomic feature
    #beyond_aTTS(dsets)

    ##### Write a bedfile with all the polyA sites confirmed by both annotation and
    ##### 3UTR poly(A) reads (For Sarah)
    #write_verified_TTS(clusters) broken, clusters replaced with dsets

    # TODO the methods in this functions must be updated to use super3utr. It
    # should be a huge hustle.
    #polyadenylation_comparison(dsets, super_3utr, settings)

    ##### get the before/after coverage ratios for trusted epsilon ends
    #before_after_ratio(dsets)
    # RESULT these are not useful values. rna seq is too stochastic

    #utr_length_comparison(settings, dsets)

    ## The classic polyA variation distributions
    #classic_polyA_stats(settings, dsets)

    #reads(settings)

    #data_annotation_correspondence(settings)

    #UTR_processing(settings)


if __name__ == '__main__':
    main()

# DISCUSSION #
# XXX README COMING BACK FROM MALAGA

# The new poly(A) investigation: looking for small poly(A) tail remnants. It has
# been shown that human rRNA have poly(A) tails in cytoplasm and nucleus. These
# tails are likely degradation-related, as poly(A) tails assist the exosome in
# degradation. This degradation happens chiefly in the nucleus at least in
# yeast. If you find increased poladenylation in the nucleus, perhaps in the
# poly(A) minus sets (the degradation-coupled poly(A) tails are short), you
# could have found some evidence for this. The only caveat is that you have to
# watch out for splice sites. Maybe you can filter with the splice sites found
# by the gem mapper. Then show that along transcripts, poly(A) sites are
# distributed differently in the nucleus and cytocolic fractions (and chromatin
# and nucleoplasm and nuclolus). Also show that the poly(A) reads are of
# differnt type (length and composition wise).

# The general merit of the method is a WholeCell comparison to mimick Fu et al.
# 2011. Just merge your sites together and show the increase in poly(A) sites
# and also show the merger of annotated sites with annotations

# You have shown already that you have a high specificity for annotated poly(A)
# sites in annotated 3UTRs.

# Next, you must show that your gene-internal reads are not simply split-map
# reads. Simply use the gem-mapped split-map reads and do a bed-intersect.

# The reason you have more maps in the nucleus to the introns is simply because
# they don't exist in the cytoplasm! But the question is what is the nature of
# these poly(A) reads.

# evidence against split-mapped reads: In the cytoplasm you expect to have only
# split-map poly(A) errors. They will fall in in annotated exon-intron
# junctions. Simply make a region like this with annotation_parser. By comparing
# the # of poly(A) sites in the 3UTR to the number across the exon-exon
# junctions, you should get a reasonably number of the number of split-mapped
# reads.

    # THIS IS PRIMARILY FOR MY OWN ANALYSIS
    # TRYING TO SPLITMAP
    # if the directory splitmapped reads exists, and a splitmapped file
    # inside there exists as well, filter out the poly(A) reads that map to
    # splitmapped regions
    # TODO: result: you don't achieve anything.
    # For 3UTR EXONIC: before 6000/15000 were not mapping to annotated sites
    # (they were 'NA'). After, 2700/7100 were 'NA'. The ratio in both cases
    # is roughly 40%. Your screening did not improve this ratio, which is
    # odd. Maybe the splitmapping data cannot be trusted. How does this fare
    # for the intronic data?

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

# UPDATE: you damn well need to do some more research on what exactly happens to
# an mRNA as it crosses from the nucleus to the cytoplasm. You know that
# polyadentlaytion and processing happens in the nucleus. In the cytoplasm the
# mRNA is translated/degraded.

# How reproducable are the polyA reads?
# If for the same gene, for a similar before/after coverage, what is the
# variation within the number of polyA reads?

# IDEA: measure of coverage-limitation: how many polyA clusters (>5 reads) are
# found with only 1 read in other compartments/datasets?

# Q: what should we call a reliable cluster?

# Q: poly(A) cluster coverage in each 3UTR -- does it show that 1 3UTR isomer is
# preferentially used?

# Q: In the Whole Cell -- what percentage of sample is cytoplasm and what
# percentage is nucleus? I suspect that there is more cytoplasm than whole cell.
# A: Roderic confirmed. The 'Whole Cell' dataset resembles mostly the Cytoplasm.

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
# how many colors you need. Can you modify it to your needs?

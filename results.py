"""
Script for displaying and summarizing the results from utail.py.
"""

from __future__ import division
import os
import ConfigParser
import sys
from itertools import combinations as combins
from copy import deepcopy

from subprocess import Popen, PIPE

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#from matplotlib import lines

#plt.ion() # turn on the interactive mode so you can play with plots
plt.ioff() # turn off interactive mode for working undisturbed

#import csv2latex

from operator import attrgetter
from operator import itemgetter
import math

import numpy as np
#from scipy import stats

import time
import glob
import cPickle as pickle

import pybedtools

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

    def __init__(self, dsets, begend, center, coverages, covrg_sum, tail_types,
                 tail_infos, strand, nearby_PASes, PAS_distances,
                 annot_distances, PET_sup,
                 cuffL_support):

        self.beg, self.end = begend
        self.dsets = dsets
        self.coverages = coverages
        self.polyA_coordinate = center
        self.nr_support_reads = covrg_sum
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

        if any(PET_sup):
            self.PET_support = 1
        else:
            self.PET_support = 0

        if any(cuffL_support):
            self.cufflinks_support = 1
        else:
            self.cufflinks_support = 0


class Only(object):
    """ For the onlypolyA files. They don't have coverage and have fewer
    parameters than the othes.
    """

    def __init__(self, input_line):

        (chrm, beg, end, utr_ID, strand, polyA_coordinate,
         polyA_coordinate_strand, tail_info, annotated_polyA_distance, nearby_PAS,
         PAS_distance, PET_support, cuff_support, number_supporting_reads,
         nr_unique_supporting_reads, unique_reads_spread) = input_line.split('\t')

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
        self.PET_support = int(PET_support)
        self.cuff_support = cuff_support

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

    # if you had made a @staticmethod here you would not need to do
    # plot_obj=Plotter() and the plot_obj.all_lying_bar, but could instead have
    # called Plotter.all_lying_bar directly
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

            # set a grid so its easier to see
            ax.xaxis.grid(True, alpha=0.5, color='lightgrey')

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

        output_dir = os.path.join(here, 'thesis_figures')

        filename = 'All_polyA_different_compartments'
        filepath = os.path.join(output_dir, filename+'.pdf')
        fig.savefig(filepath, format='pdf')
        filepath = os.path.join(output_dir, filename+'.eps')
        fig.savefig(filepath, format='eps', papertype='A4')

    @staticmethod
    def unified_lying_bar(data_dict, regions, title, blobs, blob_order, here,
                          compartments, filename, compDsetNr):
        """
        Side-lying barplot for 2 of your bar plots. You recently joined them to
        make the code easier to maintain.
        """
        blob_order
        def blobsort(btuple):
            return blob_order.index(btuple[1])

        # The nr and names of bars in the plot
        #plot_keys = ['PAS', 'PET', 'T', 'all']
        plot_keys = ['PAS', 'PET', 'all']
        colors = {'all': 'm', 'T': 'g', 'PAS': 'b', 'PET':'g'}

        labels = {'all': 'All', 'T': 'Mapped with poly(T)',
                  'PAS': 'With downstream PAS', 'PET': 'Supported by PET'}

        sidelabels = {'Whole_Cell': 'Whole cell', 'Cytoplasm': 'Cytoplasm',
                      'Nucleus': 'Nucleus'}

        colnr = len(blobs)
        for thr in ['2+', '3+']:
            (fig, axes) = plt.subplots(3, colnr, sharex=True)

            # make enough letters to cover all the plots
            plot_names = list('IHGFEDCBA')[-colnr*3:]
            for comp_nr, comp in enumerate(compartments):
                for blob_nr, (blob, blob_name)\
                        in enumerate(sorted(blobs.items(), key = blobsort)):

                    #plotme = {'all': [], 'PAS': [], 'T': [], 'PET': []}
                    plotme = {'all': [], 'PAS': [], 'PET': []}

                    # get the height of the bars from the input
                    for region in regions:

                        for k in plot_keys:
                            plotme[k].append(data_dict[thr][comp][region][blob][k])

                    # you want to place the left foot of all at 1,2,3, etc
                    # you want to place the left foot of T at 1.25
                    # you want to place the left foot of PAS at 1.375
                    heights = {'all': 0.25, 'PET': 0.125, 'PAS': 0.125}

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

                    ax = axes[comp_nr, blob_nr]
                    rects = dict() # save the rect information

                    # put A and B letters
                    # transform=ax.transAxes means use the axes frame (0,0) to (1,1)
                    ax.text(0, 1.05, plot_names.pop(), transform=ax.transAxes,
                            fontsize=18, fontweight='bold', va='top')

                    # set a grid so it's easier to see the values
                    ax.xaxis.grid(True, ls='-', alpha=0.3, color='grey')

                    # make the actual plots
                    for pkey in plot_keys:
                        rects[pkey] = ax.barh(bottom=pos[pkey],
                                              width=plotme[pkey],
                                              height=heights[pkey],
                                              color=colors[pkey],
                                              label=labels[pkey])

                    ## print either the number or percentage
                    # XXX For now, don't do this.
                    #for pkey, rs in rects.items():
                        #for r_nr, rect in enumerate(rs):
                            #width = int(rect.get_width())
                            #xloc = width + 100
                            #yloc = rect.get_y()+rect.get_height()/2.0
                            #clr = 'black'
                            #align = 'left'

                            #if pkey == 'all':
                                #txt = width
                                #fsize=11
                            #else:
                                #divby = plotme['all'][r_nr]
                                #try:
                                    #txt = format(width/divby, '.2f')
                                #except ZeroDivisionError:
                                    #txt = '0'

                                #fsize=9
                                #yloc = yloc - 0.04

                            ## ylocation, centered at bar

                            #ax.text(xloc, yloc, txt,
                                     #horizontalalignment=align,
                                     #verticalalignment='center', color=clr,
                                     #weight='bold', fontsize=fsize)

                    # print the total number for 'all', and the percentage of
                    # 'all' for the other two
                    # specify xticks if needeed

                    # put some titles here and there
                    # get the y-ticks. they should centered
                    center = sum(heights.values())/2.0
                    yticks = np.arange(start+center, dpoints+start)

                    if blob_nr == 0:
                        slabel = sidelabels[comp]+ ' ({0} datasets)'\
                                .format(compDsetNr[blob][comp])
                        ax.set_ylabel(slabel, size=22)
                        ax.set_yticks(yticks) # set the 3utr-exonic etc
                        ax.set_yticklabels(regions, size=18) # 
                    else:
                        ax.set_yticklabels([])

                    ax.set_ylim(start-0.5, dpoints+1) # extend the view

                    if comp_nr == 0:
                        ax.set_title(blobs[blob], size=22)


            # put the legend only in the top-left corner plot
            axes[0][-1].legend(loc='upper right')

            # Set xlim (it's shared)
            xlm = ax.get_xlim()
            extra_space = 2000
            ax.set_xlim((0, xlm[1]+extra_space))

            ax.set_xticklabels([]) # remove xticks

            # set ticks for the two bottom axes
            stepsize = 5000
            xticks = np.arange(0, xlm[1]+stepsize, stepsize)
            f = lambda x: '' if x%(stepsize*2) else int(x)
            ticklabels = [f(tick) for tick in xticks]

            # set x-ticks in the two last axes
            two_last = axes[-1]
            for ax in two_last:
                ax.set_xticks(xticks)
                ax.set_xticklabels(ticklabels)

            fig.subplots_adjust(wspace=0.1)
            fig.subplots_adjust(hspace=0.2)
            #fig.suptitle(title+ 'for {0}'.format(titles[key1]), size=20)

            # different size depending on number of blobs
            if colnr == 2:
                fig.set_size_inches(16,19)
            if colnr == 3:
                fig.set_size_inches(18,19)

            #output_dir = os.path.join(here, 'Results_and_figures',
                                      #'GENCODE_report', 'Figures')
            output_dir = os.path.join(here, 'thesis_figures')

            fname = filename+'_{0}'.format(thr)
            filepath = os.path.join(output_dir, fname+'.pdf')
            fig.savefig(filepath, format='pdf')
            #filepath = os.path.join(output_dir, fname+'.eps')
            #fig.savefig(filepath, format='eps', papertype='A4')

    def AllLadderPlot(self, plotdict, settings, data_grouping, yek):

        # small color dictionary
        colors = ['m', 'r', 'b', 'g', 'k']
        cols = {}
        for indx, title in enumerate(data_grouping.keys()):
            cols[title] = colors[indx]

        (fig, axes) = plt.subplots(1,2)

        for plot_nr, (title, dictsum) in enumerate(plotdict.items()):
            ax = axes[plot_nr]

            read_counts = []
            cluster_counts = []
            PAScluster_counts = []
            PETcluster_counts = []

            # you must loop throug dependnig on lenghr
            def sorthelp(tup):
                return len(tup[0])

            for dsets, countdict in sorted(dictsum.items(), key=sorthelp):
                countdict = countdict[yek]
                # the the sum of reads from these datasets
                #x = [get_dsetreads(settings, '3UTR-exonic')[ds]
                     #for ds in dsets.split(':')]
                x = [get_dsetreads(settings, 'whole')[ds]
                     for ds in dsets.split(':')]
                read_counts.append(sum(x))

                # get all cluster counts
                cluster_counts.append(countdict['All'])
                PAScluster_counts.append(countdict['PAS'])
                PETcluster_counts.append(countdict['PET'])

            ax.plot(read_counts, cluster_counts, color=cols[title],
                          linewidth=3, label='All sites')
            ax.plot(read_counts, PAScluster_counts, ls='--', color=cols[title],
                          linewidth=3, label='Sites with PAS')
            ax.plot(read_counts, PETcluster_counts, ls=':', color=cols[title],
                          linewidth=3, label='Sites with PET')

            ax.set_xlabel('Billions of reads', size=22)
            ax.set_ylabel('Polyadenylation sites', size=22)
            #ax.set_title('Polyadenylation site discovery saturates fast', size=22)
            ax.set_title(title, size=24)

            # Sort the legends to your preference
            from matplotlib.font_manager import FontProperties
            ax.legend(loc=0, prop=FontProperties(size=15))

            # Set a grid on the y-axis
            ax.yaxis.grid(True)
            ax.xaxis.grid(True)

            xmin, xmax = ax.get_xlim()
            ymin, ymay = ax.get_ylim()

            # change the size of the ticks
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(15)

            if plot_nr == 0:
                label = 'A'
            else:
                label = 'B'

            # put A and B letters
            # transform=ax.transAxes means use the axes frame (0,0) to (1,1)
            ax.text(-0.15, 1.05, label, transform=ax.transAxes, fontsize=20,
                    fontweight='bold', va='top')

        output_dir = os.path.join(settings.here, 'thesis_figures')

        plt.tight_layout()
        fig.set_size_inches(16,8)
        filename = 'Saturation_plot_{0}'.format(yek)
        filepath = os.path.join(output_dir, filename+'.pdf')
        fig.savefig(filepath, format='pdf')
        #filepath = os.path.join(output_dir, filename+'.eps')
        #fig.savefig(filepath, format='eps', papertype='A4')

    def CyNucComp(self, data_dict, nr_comp, settings):
        """
        Make a 2X2 plot of Cyto Nucl with #clustrs and #reads falling in those
        clusters RoyalBlue is a nice color

    {'Cytoplasm': {'count': [4081, 2938, 7796],
                   'magnitude': [142051, 88724, 438443]},
     'Nucleus': {'count': [3104, 2376, 7299], 'magnitude': [75848, 45606, 282801]}}
        """

        fig, axes = plt.subplots(2,2)
        c = 'RoyalBlue'

        ylabdict = {'count': '# of clusters',
                    'magnitude': '# of poly(A) reads'}

        xpos = np.arange(1,nr_comp+1)

        for measure_nr, measure in enumerate(['count', 'magnitude']):
            for comp_nr, comp in enumerate(['Cytoplasm', 'Nucleus']):

                ax = axes[measure_nr, comp_nr]

                heights = data_dict[comp][measure]
                ax.bar(xpos, heights, color=c, width=1)
                #midpart = [' ' for x in range(len(xpos-2))]
                #ax.set_xticklabels(['beginning'] + midpart + ['end'])
                ax.set_xlim(1, nr_comp+2)

                if measure_nr == 0:
                    ax.set_title(comp)
                if comp_nr == 0:
                    ax.set_ylabel(ylabdict[measure])

        output_dir = os.path.join(settings.here, 'Results_and_figures',
                                  'GENCODE_report', 'Figures')
        filename = 'PolyA_site_usage_{0}'.format(nr_comp)
        filepath = os.path.join(output_dir, filename+'.pdf')
        fig.savefig(filepath, format='pdf')


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
            maxlines = 200

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

        super_3utr[region] = get_super_regions(dsets[region], settings,
                                               batch_key)

    return dsets, super_3utr

def get_super_regions(dsets, settings, batch_key):
    """
    Use bed-tools to cluster the regions

    Let each dna-patch know about all the cell lines and compartments where its
    poly(A) clusters are found.

    Update all the clusters in dsetswith their super-cluster key
    At the same time create a 3UTR dict. Each dna-patch has information about the
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

    cluster_file = os.path.join(settings.here, 'merge_dir', 'merged_pA_'+batch_key)
    temp_cluster_file = cluster_file + '_temp'

    out_handle = open(cluster_file, 'wb')

    line_nr = 0
    for cell_line, comp_dict in dsets.items():
        for comp, utr_dict in comp_dict.items():
            for feature_id, feature in utr_dict.iteritems():

                # store every found feature in super_features for later
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
                                     str(cls.nr_support_reads),
                                     str(cls.annotated_polyA_distance),
                                     str(cls.PET_support),
                                     str(cls.cuff_support)])
                    line_nr +=1
                    # bug at 6232 -> way too large name.
                    # cls id anti-3utr_32546
                    #if line_nr == 6232:
                        #debug()
                    # also add PAS_distance and PAS_type
                    out_handle.write('\t'.join([cls.chrm,
                                                str(cls.polyA_coordinate),
                                                str(cls.polyA_coordinate+1),
                                                name,
                                                str(0),
                                                cls.strand+'\n'
                                               ]))
    out_handle.close()

    junctions = os.path.join(settings.here, 'junctions',
                             'splice_junctions_merged.bed')

    # This is for CRG only... sigh.
    #mySubBed = '/users/rg/jskancke/programs/other/bedtools/bin/subtractBed'
    mySubBed = '/usr/local/bin/subtractBed'
    bedDir, subtrFile = os.path.split(mySubBed)
    # You need to change directory to other bedtools bin, otherwise subtract bed
    # doesnt find the bedtools main file it seems to be using
    #os.chdir(bedDir) CRG

    # XXX subtractBed fails because merged_pA_first_ladder contas a line of
    # 30mbytes at line 6231 with region=whole
    #debug()
    # i.5) cut away the exon-exon noise
    cmd = [mySubBed, '-a', cluster_file, '-b', junctions]
    p1 = Popen(cmd, stdout=open(temp_cluster_file, 'wb'), stderr=PIPE)
    p1.wait()

    # check for errors
    errors1 = [l for l in p1.stderr]
    if errors1 != []:
        debug()

    #mySlopBed = '/users/rg/jskancke/programs/other/bedtools/bin/slopBed' CRG
    mySlopBed = 'slopBed'
    # ii) expand the entries
    expanded = temp_cluster_file + '_expanded'
    cmd = [mySlopBed, '-i', temp_cluster_file, '-g', settings.hg19_path, '-b', '15']
    p2 = Popen(cmd, stdout=open(expanded, 'wb'), stderr=PIPE)
    p2.wait()

    # check for errors
    errors2 = [l for l in p2.stderr]
    if errors2 != []:
        debug()

    #myMergeBed = '/users/rg/jskancke/programs/other/bedtools/bin/mergeBed' CRG
    myMergeBed = 'mergeBed'
    # iii) merge and max scores
    cmd = [myMergeBed, '-nms', '-s', '-scores', 'sum', '-i', expanded]
    p3 = Popen(cmd, stdout=PIPE, stderr=PIPE)

    # check for errors. This doesn't work. Its it because you haven't done a
    #wait?
    #errors3 = [l for l in p3.stderr]
    #if errors3 != []:
        #debug()

    # go back to here CRG
    #os.chdir(settings.here)

    # iv) read through the centered output and update the super_cluster
    for line in p3.stdout:
        (chrm, beg, end, names, maxcovrg, strand) = line.split()

        center = int((int(beg)+int(end))/2)

        dsets = []
        f_ids = []
        tail_types = []
        tail_infos = []
        nearby_PASes = []
        PAS_distances = []
        coverages = []
        annot_distances = []
        pet_sup = []
        cuf_sup = []

        # un-split the joined names
        name_list = names.split(';')
        for names2 in name_list:
            (dset, f_id, t_type, t_info, nearby_PAS, PAS_distance, coverage,
             annot_dist, PET_support, cufflinks_support) = names2.split('%')
            dsets.append(dset)
            f_ids.append(f_id)
            tail_types.append(t_type)
            tail_infos.append(t_info)
            nearby_PASes.append(nearby_PAS)
            PAS_distances.append(PAS_distance)
            coverages.append(coverage)
            annot_distances.append(annot_dist)
            pet_sup.append(int(PET_support))
            cuf_sup.append(int(cufflinks_support))

        feature_id = set(f_ids).pop()

        polyA_coverage = sum([int(c) for c in coverages])

        # create cluster object
        cls = Super(dsets, (beg, end), center, coverages, polyA_coverage,
                    tail_types, tail_infos, strand, nearby_PASes, PAS_distances,
                    annot_distances, pet_sup, cuf_sup)

        # add the cluster to super_featuress
        super_features[feature_id].super_clusters.append(cls)


    return super_features

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

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

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

def get_clustercount(settings, region):
    """ For each dataset, get the number of total reads. The region doesn't
    matter, because the number of reads are dataset-specific.
    """

    dsetclusters = {}
    polyA_files = settings.polyAstats_files(region)
    for dset, dsetpath in polyA_files.items():

        filedict = dict((line.split('\t')[0], line.split('\t')[1])
                        for line in open(dsetpath, 'rb'))

        dsetclusters[dset] = int(filedict['Nr of clusters'].rstrip())

    return dsetclusters

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

        if cls.PET_support:
            this_dict['wPASPET']['info_dict'][keyw] += 1

            taildict = this_dict['wPASPET']['tail_lens'][keyw]
            taildict = avrg_tail(cls.tail_info, taildict)

        if cls.cufflinks_support:
            this_dict['wPASCufflinks']['info_dict'][keyw] += 1

            taildict = this_dict['wPASCufflinks']['tail_lens'][keyw]
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

    if cls.PET_support:
        this_dict['wPET']['info_dict'][keyw] += 1

        taildict = this_dict['wPET']['tail_lens'][keyw]
        taildict = avrg_tail(cls.tail_info, taildict)

    if cls.cufflinks_support:
        this_dict['Cufflinks']['info_dict'][keyw] += 1

        taildict = this_dict['Cufflinks']['tail_lens'][keyw]
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

    categories1 = ['Total clusters', 'morethan1', 'only1', 'morethan2',
                   'morethan3', 'morethan4', 'morethan10', 'morethan30']
    subcategories = ['All', 'annotated', 'wPAS', 'annotated_wPAS', 'goodPAS',
                     'bestPAS', 'wPET', 'wPASPET', 'Cufflinks', 'wPASCufflinks']

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

            bigcl['Total clusters'] = data_scooper(cls, keyw,
                                                   bigcl['Total clusters'])

            # Count clusters with 2 or more reads
            for minlim in [1,2,3,4,10,30]:
                if cls.nr_support_reads > minlim:
                    key = 'morethan{0}'.format(minlim)

                    bigcl[key] = data_scooper(cls, keyw, bigcl[key])

            # Count clusters with only 1 read
            if cls.nr_support_reads == 1:

                bigcl['only1'] = data_scooper(cls, keyw, bigcl['only1'])

    return bigcl

def super_cluster_statprinter(dsetclusters, region, thiskey, settings, filename):
    """
    Write key stats about the clusters
    """

    statdict = dsetclusters[thiskey]

    keys = ['Total clusters', 'morethan1', 'only1', 'morethan2',
            'morethan3', 'morethan4', 'morethan10', 'morethan30']

    #subkeys =  ['All', 'wPAS', 'goodPAS', 'bestPAS', 'annotated',
                #'annotated_wPAS']
    subkeys =  ['All', 'wPAS', 'goodPAS', 'annotated', 'annotated_wPAS', 'wPET',
                'wPASPET', 'Cufflinks', 'wPASCufflinks']

    datakeys = ['info_dict', 'tail_lens']

    headers = {'Total clusters': '### All clustes ###',
               'only1': '### Clusters with only 1 coverage ###',
               'morethan1': '### Clusters with 2 or more coverage ###',
               'morethan2': '### Clusters with 3 or more coverage ###',
               'morethan3': '### Clusters with 4 or more coverage ###',
               'morethan4': '### Clusters with 5 or more coverage ###',
               'morethan10': '### Clusters with 10 or more coverage ###',
               'morethan30': '### Clusters with 30 or more coverage ###'}

    subheaders = {'wPAS': 'With PAS',
                  'All': 'All',
                  'wPET': 'With PET',
                  'wPASPET': 'With PAS and PET',
                  'wPASCufflinks': 'With PAS and Cufflinks',
                  'Cufflinks': 'With Cufflinks',
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
                          'thesis_figures/thesis_stats')

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

def clusterladder(settings, speedrun, reuse_data):
    """
    The more reads, the more polyAs, up to a point.
    """

    # where you save your reused data
    plot_save_file = os.path.join(settings.here, 'thesis_figures','save_data',
                                  'clusterladder')

    region = 'whole'
    #region = '3UTR-exonic'

    #1) Make a dictionary: dataset-> nr of total reads and dataset -> nr of
    #clusters
    dsetreads = get_dsetreads(settings, region=region)
    clustercount = get_clustercount(settings, region=region)

    #2) Make super-clusters for your datasets of choice

    wc_c = [ds for ds in settings.datasets if (('Cytoplasm' in ds) or
                   ('Whole_Cell' in ds) or ('Nucleus' in ds)) and
            (not 'Minus' in ds)]

    wc_c_minus = [ds for ds in settings.datasets if (('Cytoplasm' in ds) or
                   ('Whole_Cell' in ds) or ('Nucleus' in ds)) and 'Minus' in ds]

    #data_grouping = {'Poly(A) plus': c,
                    #'Poly(A) minus': c_minus}
    data_grouping = {'Poly(A) plus': wc_c,
                     'Poly(A) minus': wc_c_minus}

    # keep a dictionary with reference to all the plots
    plotdict = {}

    if speedrun:
        max_sets = 5 # you need minimum 3. why?
        data_grouping['Poly(A) plus'] = data_grouping['Poly(A) plus'][:max_sets]
        data_grouping['Poly(A) minus'] = data_grouping['Poly(A) minus'][:max_sets]

    # dont speedrun the rest of the application!
    #speedrun = True
    speedrun = False

    # only calculate again if not reusing data
    t1 = time.time()
    if not reuse_data:
        for title, dsets in data_grouping.items():

            # sort the dsets in cell_lines by # of reads
            def mysorter(dset):
                return clustercount[dset]
            all_dsets = sorted(dsets, key=mysorter, reverse=True)

            # add more and more datasets
            subsets = [all_dsets[:end] for end in range(1, len(all_dsets)+1, 2)]
            #subsets = [all_dsets[:end] for end in range(1, len(all_dsets)+1)]

            subsetcounts = {}

            for subset in subsets:

                # Get the number of 'good' and 'all' clusters
                key = ':'.join(subset)
                batch_key = 'first_ladder'
                dsetclusters = get_dsetclusters(subset, region, settings,
                                                speedrun, batch_key)

                subsetcounts[key] = count_clusters(dsetclusters, dsetreads)

            plotdict[title] = subsetcounts

        # save plot data
        pickle_wrap(storage_file=plot_save_file, action='save', data=plotdict)

    else:
        plotdict = pickle_wrap(storage_file=plot_save_file, action='load')


    print("It took {0} seconds to merge".format(time.time()-t1))
    p = Plotter()

    for yek in ['2+', '3+']:

        p.AllLadderPlot(plotdict, settings, data_grouping, yek)


def pickle_wrap(storage_file, action, data=0):
    """
    Either save data to plot_save_file or load data from plot_save file and
    return it.
    """

    if action == 'save' and data != 0:
        pickle.dump(data, open(storage_file, 'wb'))

    if action == 'save' and data == 0:
        print('Watcha savin bro?')

    if action == 'load' and data == 0:
        output = pickle.load(open(storage_file, 'rb'))

        return output

def count_clusters(dsetclusters, dsetreads):
    """
    Parse through and count the number of clusters with +1 and number with PAS
    """

    countdict = {'2+': {
        'All': sum(dsetclusters['morethan1']['All']['info_dict'].values()),
        'PAS': sum(dsetclusters['morethan1']['wPAS']['info_dict'].values()),
        'PET': sum(dsetclusters['morethan1']['wPET']['info_dict'].values())},
                '3+': {
        'All': sum(dsetclusters['morethan2']['All']['info_dict'].values()),
        'PAS': sum(dsetclusters['morethan2']['wPAS']['info_dict'].values()),
        'PET': sum(dsetclusters['morethan2']['wPET']['info_dict'].values())}
        }

    return countdict

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
             os.path.join(settings.here,
                          'annotated_polyAsites/gencode_polyA_proper.bed'),
             'polyAdb':\
             os.path.join(settings.here, 'annotated_polyAsites/polyA_db_proper.bed')}

    region = 'whole'

    subset = all_ds
    if speedrun:
        subset = subset[:2]

    batch_key = 'venn'
    dsets, super_3utr = super_falselength(settings, region, batch_key,
                                          subset, speedrun)

    for yek in ['+2', '+3']:
        # create 'outfile' and fill up the stats dicts
        paths['whole'] = whole_tobed(super_3utr, outdir, region, yek)

        # get all paths to merged polyA files and their linenumbers
        (paths, merged_stats) = intersect_polyAs(paths, outdir, region)

        # make r-code for making the venn diagrams
        make_venn(paths, merged_stats, outdir, region, settings, yek)


def make_venn(paths, wcdict, outdir, region, settings, yek):
    """ Call upon R to summon forth the elusive Venn diagrams
    """

    outdir = os.path.join(settings.here, 'analysis', 'R_code')
    outpath = os.path.join(outdir, 'venn.r')
    outfile = open(outpath, 'wb')

    not_weighted = os.path.join(outdir, 'venn_notWeighted_{0}.pdf'.format(yek))
    weighted = os.path.join(outdir, 'venn_weighted_{0}.pdf'.format(yek))

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
        # run imagemagic to make eps?


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

    # Send in A and [B,C]
    paths = intersect_wrap(paths, outdir, extendby=10)

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

def intersect_wrap(paths, outdir, extendby=20):
    """
    1) extend each of the joiner files with some value
    2) intersect them
    3) return paths dict with link to the intersected path
    """
    #hg19 = '/home/jorgsk/work/3UTR/ext_files/hg19'
    hg19 = '/users/rg/jskancke/phdproject/3UTR/the_project/ext_files/hg19'

    # function used only here
    def extender(extendby, path, hg19, ext_path):

        out_handle = open(ext_path, 'wb')

        cmd1 = ['slopBed', '-b', str(extendby), '-i', path, '-g', hg19]
        p1 = Popen(cmd1, stdout = PIPE)

        for line in p1.stdout:
            out_handle.write(line)

        out_handle.close()

    # 1) extend all files
    for pathname, path in paths.items():

        ext_name = pathname+'_extended_{0}'.format(extendby)
        ext_path = os.path.join(outdir, ext_name)

        extender(extendby, path, hg19, ext_path)
        paths[ext_name] = ext_path

    # there are 4 intersections to be performed:
        # AIB, AIC, BIC, and AIBIC

    intersecters = ['gencode', 'polyAdb', 'whole']

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

            cmd2 = ['intersectBed', '-s', '-wa', '-wb', '-a', filea, '-b', fileb]
            p2 = Popen(cmd2, stdout=PIPE)

            out_handle = open(isect_path, 'wb')
            for line in p2.stdout:

                # re-center the intersection
                if len(line.split()) == 12:
                    (chrmA, begA, endA, nameA, valA, strndA,
                    chrmB, begB, endB, nameB, valB, strndB) = line.split()
                else:
                    debug()

                beg = min(begA, begB)
                end = max(endA, endB)

                center_dist = math.ceil((int(end)-int(beg))/2)
                center = int(int(end) - center_dist)

                # write the center of the new mergers
                out_handle.write('\t'.join([chrmA, str(center), str(center+1), '0',
                                           '0', strndA]) + '\n')
            out_handle.close()
            paths[isect_name] = isect_path

    return paths

def whole_tobed(super_3utr, outdir, region, yek):
    """ Merge the paths in 'only_these' into one file, and then run mergeBed on
    them, and finally report the centers of the mergings obtained.  Return
    'paths' with the path to the merged version.
    """

    superbed_path = os.path.join(outdir, region+'.bed')
    handle = open(superbed_path, 'wb')

    if yek == '+2':
        minLim = 1
    else:
        minLim = 2

    annotated = 0

    for utr_name, utr in super_3utr[region].iteritems():

        for cls in utr.super_clusters:

            # Write to file
            if cls.nr_support_reads>minLim or cls.annotated_polyA_distance!='NA':

                beg = cls.polyA_coordinate

                entry = '\t'.join([utr.chrm, str(beg), str(beg+1), utr.ID,
                                   str(cls.nr_support_reads), utr.strand])

                handle.write(entry + '\n')
                # write only those with PAS!
                if cls.annotated_polyA_distance!='NA':
                    annotated +=1

    handle.close()

    print 'annotated', annotated, minLim

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
    #region = 'Intergenic'
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

        batch_key = 'all_stats_cumul'
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
    #regions = ['3UTR-exonic']

    for region in regions:

        for compartments in [['Cytoplasm'], ['Nucleus'], ['Whole_Cell']]:
            for ignorers in [[], ['Minus']]:
                if ignorers == []:
                    demanders = ['Minus']
                else:
                    demanders = []

                if ignorers == []:
                    filename = '+'.join(compartments+demanders)
                else:
                    filename = '+'.join(compartments+demanders)\
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

def get_cellLineRPKM(cell_lines, compartment = 'C'):
    """
    Return a dict [cell_line][ts_id] = rpkms

    The compartment can be specified with the comp variable
    """
    carrie = '/users/rg/projects/NGS/Projects/ENCODE/hg19main/GingerasCarrie'

    #0) make a dir -> cell line dict
    dir_2_cell_line = {}
    for line in open(os.path.join(carrie, 'Dataset.info.txt')):
        (dirinfo, cell_l, comp, longshort, rplica) = line.split()
        dirname = dirinfo[4:]
        dir_2_cell_line[dirname] = cell_l

    thefile = '/users/rg/projects/NGS/Projects/ENCODE/hg19main/'\
            'GingerasCarrie/comparison/Trans.Expression.ENCODE.txt'

    #1) a cell_line 2 header
    cellLine2header = {}
    for cellLine in cell_lines:
        for dirPart, cL in dir_2_cell_line.items():
            if cL == cellLine:
                comp = dirPart[3:]
                if dirPart.startswith('0') and comp == compartment:
                    if cellLine in cellLine2header:
                        cellLine2header[cellLine].append('Ging'+dirPart)
                    else:
                        cellLine2header[cellLine] = ['Ging'+dirPart]

    # parse the transcript rpkm file
    thehandle = open(thefile, 'rb')
    header = thehandle.next()
    header = header.split('\t')

    cell2index = dict((cL, []) for cL in cell_lines)

    for cellLine, dirnames in cellLine2header.items():
        for hindex, h in enumerate(header):
            if h in dirnames:
                cell2index[cellLine].append(hindex+1)

    cL2rpkm = {}

    # for each cell line, return the RPKM of all its transcripts
    for cL in cell_lines:
        tsrpkm = {}
        for line in thehandle:
            flds = line.split()
            ts_id = flds[0]
            inds = cell2index[cL]
            meanRPKM = np.mean([float(flds[i]) for i in inds if flds[i] != '0'])
            if math.isnan(meanRPKM):
                meanRPKM = 0

            tsrpkm[ts_id] = meanRPKM

        cL2rpkm[cL] = tsrpkm

    return cL2rpkm

def get_cellLineEXONRPKM(cell_lines, tempdir, compartment='C'):
    """
    Return a dict [cell_line][ts_id] = rpkms

    The compartment can be specified with the comp variable
    """
    carrie = '/users/rg/projects/NGS/Projects/ENCODE/hg19main/GingerasCarrie'

    #0) make a dir -> cell line dict
    dir_2_cell_line = {}
    for line in open(os.path.join(carrie, 'Dataset.info.txt')):
        (dirinfo, cell_l, comp, longshort, rplica) = line.split()
        dirname = dirinfo[4:]
        dir_2_cell_line[dirname] = cell_l

    paths = dict((cl, []) for cl in cell_lines)

    # get the paths for each cell line
    t1 = time.time()
    for dirname, dirnames, filenames in os.walk(carrie):
        if dirname.endswith('exons'):

            main_dir = os.path.split(os.path.split(dirname)[0])[1]

            # include only those cell lines you have been using
            if main_dir not in dir_2_cell_line:
                print main_dir
                continue

            if dir_2_cell_line[main_dir] not in cell_lines:
                continue
            cL = dir_2_cell_line[main_dir]

            comp = main_dir[3:]

            # we're only looking at nuclear and cytoplasm
            if comp != compartment:
                continue

            for f in filenames:
                if f == 'All.exon.rpkm.pooled.txt.gz':

                    mypath = os.path.join(dirname, f)
                    paths[cL].append(mypath)

    print('time to traverse:\t{0}'.format(time.time()-t1))

    cl2rpkmPath = {}
    # get the rpkm for each cell line
    for cL, paths in paths.items():
        rpkmdict = {}

        outfile = '{0}_{1}_exon_rpkms'.format(cL, compartment)
        outpath = os.path.join(tempdir, outfile)
        outhandle = open(outpath, 'wb')

        # add the rpkms
        for path in paths:
            p = Popen(['zcat', path], stdout=PIPE)
            for line in p.stdout:
                exkey, rpkm, d = line.split()
                chrm, beg, end, rg = exkey.split('_')

                if rg == '-1':
                    strand = '-'
                else:
                    strand = '+'

                mykey = '_'.join([chrm, beg, end, strand])

                if mykey in rpkmdict:
                    rpkmdict[mykey].append(float(rpkm))
                else:
                    rpkmdict[mykey] = [float(rpkm)]

        # write the rpkms to file
        for key, rpkmlist in rpkmdict.iteritems():
            (chrm, beg, end, strand) = key.split('_')
            rpkm = np.mean([r for r in rpkmlist if r!= '0'])
            if math.isnan(rpkm):
                rpkm = 0

            outhandle.write('\t'.join([chrm, beg, end, 'e', str(rpkm),
                                       strand])+ '\n')

        outhandle.close()

        cl2rpkmPath[cL] = outpath

    return cl2rpkmPath


def rpkm_polyA_correlation(settings, speedrun):
    """
    Get the general correlation between # of poly(A) sites in a transcript and
    its RPKM. Later you can reuse this code for getting the RPKM difference in
    the nucleus and the cytoplasm.
    """

    #speedrun = True
    speedrun = False

    temp_dir = os.path.join(settings.here, 'temp_files')

    # 1) For each cell line, get a dict with the rpkm of its transcripts
    cell_lines = ('K562', 'GM12878', 'HUVEC', 'HELAS3', 'HEPG2')

    cell_lines = cell_lines[:1]

    cellLineRPKM = get_cellLineEXONRPKM(cell_lines, temp_dir, compartment='C')

    region = 'whole'

    for cL in cell_lines:
        subset = [ds for ds in settings.datasets if ('Cytoplasm' in ds) and
                       (not 'Minus' in ds) and (cL in ds)]

        #subset = subset[:1] #debugging

        batch_key = 'rpkm_link'
        dsets, super_3utr = super_falselength(settings, region,
                                              batch_key, subset,
                                              speedrun=speedrun)

        # write the poly(A) sites to dir
        polyAspath = os.path.join(temp_dir, 'polyAs_forRPKM')
        polyAshandle = open(polyAspath, 'wb')

        for piece_id, piece in super_3utr[region].iteritems():
            for cls in piece.super_clusters:

                if cls.nr_support_reads > 1: # 

                    chrm = piece.chrm
                    beg = str(cls.polyA_coordinate)
                    end = str(cls.polyA_coordinate)
                    strand = cls.strand
                    d = '0'

                polyAshandle.write('\t'.join([chrm, beg, end,d,d, strand])+'\n')
        polyAshandle.close()

        # intersect the poly(A) sites with the transcripts
        exon_path = cellLineRPKM[cL]

        cmd = ['intersectBed', '-s', '-c', '-wa', '-a', exon_path,
               '-b', polyAspath]

        # for each merged transcript, keep how many poly(A) sites and the RPKM
        # fo the merged transcripts
        # because it is merged, I expect each p(A) sit to hit once; and I expect
        rpkms = []
        cls_nr = []

        p = Popen(cmd, stdout = PIPE)
        for line in p.stdout:
            (chrm, beg, end, d, rpkm, astrnd, a_site_count) = line.split()
            if a_site_count != '0':
                rpkms.append(float(rpkm))
                cls_nr.append(int(a_site_count))

        plt.scatter(rpkms, cls_nr)
        # Now going through every single cytoplasm and whole cell alone

def barsense_counter(super_3utr, comp, frac, region, d_keys):
    """
    Must return [comp][reg][+/-][all, pas, t]
    """

    count_dict = dict((key, {'all':0, 'PAS':0, 'T':0, 'PET':0}) for key in d_keys)

    for utr_id, utr in super_3utr[region].iteritems():

        for cls in utr.super_clusters:

            for key in d_keys:
                minlim = int(key[0])

                if cls.nr_support_reads > minlim-1:

                    count_dict[key]['all'] += 1

                    # any PAS
                    if cls.nearby_PAS[0] != 'NA':
                        count_dict[key]['PAS'] += 1

                    # Get if this was an A or a T cluster
                    if cls.tail_type == 'T':
                        count_dict[key]['T'] += 1

                    if cls.PET_support:
                        count_dict[key]['PET'] += 1

    return count_dict

def side_plot(settings, speedrun, reuse_data):
    """
    Side plot of poly(A) reads in different regions. Show PAS and PET support.
    """

    # if reusing data, save to this file
    plot_save_file = os.path.join(settings.here, 'thesis_figures', 'save_data',
                                  'side_not_intersected')

    #1 Gather the data in a dictlike this [wc/n/c|/+/-]['region'] = [#cl, #cl w/pas]

    compartments = ['Whole_Cell', 'Cytoplasm', 'Nucleus']

    fractions = ['+', '-']

    #regions = ['5UTR-exonic', '5UTR-intronic', '3UTR-exonic', '3UTR-intronic',
               #'CDS-exonic', 'CDS-intronic', 'Nocoding-exonic',
               #'Noncoding-intronic', 'Intergenic']
    #regions = ['3UTR-exonic', 'CDS-exonic', 'CDS-intronic',
               #'Nocoding-exonic', 'Noncoding-intronic', 'Intergenic']
    regions = ['3UTR-exonic', 'CDS-exonic', 'CDS-intronic', 'Intergenic']

    data_dict_keys = ['2+', '3+']

    compDsetNr = {'+': {}, '-':{}}

    reuse_data = True
    #reuse_data = False

    #speedrun = True
    speedrun = False

    total = len(compartments) * len(fractions)
    current = 0

    # if not reusing data, get all the stuff anew
    data_dict = dict((key, AutoVivification()) for key in data_dict_keys)

    if reuse_data and (not os.path.isfile(plot_save_file)):
        print("Reuse data set, but no savefile found. Gathering data ...")
        reuse_data = False

    if not reuse_data:
        for comp in compartments:
            for frac in fractions:

                if frac == '+':
                    subset = [ds for ds in settings.datasets if (comp in ds) and
                              (not 'Minus' in ds)]
                if frac == '-':
                    subset = [ds for ds in settings.datasets if (comp in ds) and
                              ('Minus' in ds)]

                compDsetNr[frac][comp] = len(subset)

                if speedrun:
                    subset = subset[:1]

                speedrun = False
                for region in regions:

                    batch_key = 'side'
                    dsets, super_3utr = super_falselength(settings, region,
                                                          batch_key, subset,
                                                          speedrun=speedrun)

                    # count the number clusters with +1, of those with PAS/good_PAS
                    counted = barsense_counter(super_3utr, comp, frac, region,
                                               data_dict_keys)
                    # add for +2 and +3
                    for key in data_dict_keys:
                        data_dict[key][comp][region][frac] = counted[key]

                current += 1
                print("Finished {0} of {1}".format(current, total))

        saver = (data_dict, compDsetNr)

        pickle_wrap(storage_file=plot_save_file, action='save', data=saver)

    else:
        (data_dict, compDsetNr) = pickle_wrap(storage_file=plot_save_file, action='load')

    # Stuff for the plot
    title = 'Polyadenlyation in different regions for different'\
            ' cellular compartments'
    blobs = {'+': 'Poly(A)+', '-': 'Poly(A)-'}
    blob_order = ['Poly(A)+', 'Poly(A)-']
    filename = 'Sidebars_pA'

    Plotter.unified_lying_bar(data_dict, regions, title, blobs, blob_order,
                        settings.here, compartments, filename, compDsetNr)

def isect(settings, merged, extendby, temp_dir, comp, reg):
    """
    Expand file A and B file by 20 nt, intersect, and return
    A - intersection
    intersection
    B - intersection
    """

    # function used only here
    def extender(extendby, path, hg19, ext_path):

        out_handle = open(ext_path, 'wb')

        cmd1 = ['slopBed', '-b', str(extendby), '-i', path, '-g', hg19]
        p1 = Popen(cmd1, stdout = PIPE)

        for line in p1.stdout:
            out_handle.write(line)

        out_handle.close()

    # 1) extend all files
    expanded = {}
    for pathname, path in merged.items():

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

        aTA, aPAS, aPET = aname.split('%')
        bTA, bPAS, bPET = bname.split('%')

        pases = aPAS.split('#') + bPAS.split('#')
        pasuniq = list(set(pases))
        if len(pasuniq)>1 and 'NA' in pasuniq:
            pasuniq.remove('NA')

        PAS = '#'.join(pasuniq)
        T = aTA

        if int(aPET)+int(bPET):
            PET = 1
        else:
            PET = 0

        name = '%'.join([T, PAS, str(PET)])

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
        outpath = path + '_Sliced'
        out_handle = open(outpath, 'wb')
        cmd = ['subtractBed', '-s', '-a', path, '-b', isect_exp_path]
        p3 = Popen(cmd, stdout = PIPE)
        for line in p3.stdout:
            out_handle.write(line)

        out_handle.close()

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
                 PET_support, Cuff_support, number_supporting_reads,
                 number_unique_supporting_reads,
                 unique_reads_spread) = line.split('\t')

                # Get if it is an 'A' or a 'T' tail
                ga = [g.split('=') for g in polyA_average_composition.split(':')]
                tail_type = sorted([(float(g[1]), g[0]) for g in ga])[-1][-1]

                name = '%'.join([tail_type, nearby_PAS, PET_support])

                # write a bed format to disc
                ext_handle.write('\t'.join([chrm,
                                          polyA_coordinate,
                                          str(int(polyA_coordinate)+1),
                                          name,
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
    p = Popen(cmd, stdout=PIPE)

    # iv) center the new transcript
    finalpath = os.path.join(temp_dir,
                             '{0}_{1}_extendify_juncfree_merged{2}'.format(comp,
                                                                          reg,
                                                                          key))
    out_handle = open(finalpath, 'wb')
    for line in p.stdout:
        (chrm, beg, end, name, maxcovrg, strand) = line.split('\t')

        ta = []
        pas = []
        pet = []
        name_list = name.split(';')

        for name2 in name_list:
            (t_info, nearby_PAS, PET_support) = name2.split('%')
            ta.append(t_info)
            pet.append(int(PET_support))
            pas.append(nearby_PAS)

        # Has PET?
        if any(pet):
            pet_support = 1
        else:
            pet_support = 0

        # A or T read?
        if ta.count('T') > ta.count('A'):
            ta_type = 'T'
        else:
            ta_type = 'A'

        pure_pas = []
        for pagroup in pas:
            splitpas = pagroup.split('#')
            for pa in splitpas:
                pure_pas.append(pa)

        uniquepas = list(set(pure_pas))

        if len(uniquepas) > 1 and 'NA' in uniquepas:
            uniquepas.remove('NA')

        pases = '#'.join(uniquepas)
        name = '%'.join([ta_type, pases, str(pet_support)])

        meanbeg = int((int(beg)+int(end))/2)

        # filter on coverage
        if int(float(maxcovrg)) > min_covr:

            out_handle.write('\t'.join([chrm,
                                       str(meanbeg),
                                       str(meanbeg+1),
                                       name,
                                       maxcovrg,
                                       strand.rstrip()+'\n']))
    out_handle.close()

    return finalpath

def standard_sideplot(settings, speedrun, reuse_data):
    """
    Intersect poly(A)+ and poly(A)- for the all regions in all compartments.
    Finally merge the poly(A) sites for each region and plot them sidebarwise.
    """

    plot_save_file = os.path.join(settings.here, 'thesis_figures', 'save_data',
                                  'side_standard')

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
    min_covr = 2

    if reuse_data and (not os.path.isfile(plot_save_file)):
        print("Reuse data set, but no savefile found. Gathering data ...")
        reuse_data = False

    if not reuse_data:
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
            if speedrun:
                plus_subset = plus_subset[:2]
                minus_subset = minus_subset[:2]

            path_dict = {'plus': plus_subset, 'minus': minus_subset}
            merged = {}

            comp ='All'
            for key, paths in path_dict.items():
                merged[key] = merge_paths(settings, paths, temp_dir, key,
                                          min_covr, comp, reg)

            reg_dict[reg] = merged

        data_dict = count_these_all(reg_dict)

        pickle_wrap(storage_file=plot_save_file, action='save', data=data_dict)

    else:

        data_dict = pickle_wrap(storage_file=plot_save_file, action='load')

    p = Plotter()

    title = 'Polyadenylation in different regions'

    p.all_lying_bar(data_dict, regions, title, settings.here)

def intersection_sideplot(settings, speedrun, reuse_data):
    """
    Intersect poly(A)+ and poly(A)- for the different regions and blot
    intersection and unique ones on both sides.
    """

    # where you save your reused data
    plot_save_file = os.path.join(settings.here, 'thesis_figures','save_data',
                                  'side_intersection')

    compartments = ['Whole_Cell', 'Cytoplasm', 'Nucleus']

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
    # maybe this isn't a good idea? Maybe do it like the other function?

    temp_dir = os.path.join(settings.here, 'temp_files')
    min_covr = 2
    extendby = 20

    speedrun = True
    #speedrun = False

    reuse_data = True

    if reuse_data and (not os.path.isfile(plot_save_file)):
        print("Reuse data set, but no savefile found. Gathering data ...")
        reuse_data = False

    if not reuse_data:
        compDsetNr = {'plus_sliced': {},
                      'minus_sliced': {},
                      'intersection': {}}
        comp_dict = {}
        for comp in compartments:

            # Merge the plus and minus subsets for each region
            reg_dict = {}
            for reg in regions:

                plus_subset = [path for ds, path in dsetdict[reg].items() if 
                               (not 'Minus' in ds) and comp in ds]
                minus_subset = [path for ds, path in dsetdict[reg].items() if 
                                ('Minus' in ds) and comp in ds]
                if speedrun:
                    plus_subset = plus_subset[:2]
                    minus_subset = minus_subset[:2]

                path_dict = {'plus': plus_subset, 'minus': minus_subset}
                merged = {}

                # merge each plus and minus separately
                for key, paths in path_dict.items():
                    merged[key] = merge_paths(settings, paths, temp_dir, key,
                                              min_covr, comp, reg)

                # join plus and minus and return a dictionary
                # separated['plus'/'minus'/'intersected']
                separated = isect(settings, merged, extendby, temp_dir, comp, reg)

                reg_dict[reg] = separated

            comp_dict[comp] = reg_dict

            compDsetNr['plus_sliced'][comp] = len(plus_subset)
            compDsetNr['intersection'][comp] = len(plus_subset)
            compDsetNr['minus_sliced'][comp] = len(minus_subset)

        # save the bedfiles for the different regions. DOn't doit any more.
        #save_dir = os.path.join(settings.here, 'analysis', 'pure')
        #save_pure(comp_dict, save_dir) # note! you output the split files.

        data_dict = count_these(comp_dict)

        saver = (data_dict, compDsetNr)

        pickle_wrap(storage_file=plot_save_file, action='save', data=saver)

    else:

        # get the data from the previous run
        (data_dict, compDsetNr) = pickle_wrap(storage_file=plot_save_file,
                                              action='load')

    title = 'Polyadenlyation in different regions for different'\
            ' cellular compartments'

    blobs = {'plus_sliced': 'Poly(A)+ unique',
             'intersection': 'Common to poly(A)+/-',
             'minus_sliced': 'Poly(A)- unique'}

    filename = 'intersected_sidebars_pA'
    blob_order = ['Poly(A)+ unique', 'Common to poly(A)+/-', 'Poly(A)- unique']

    Plotter.unified_lying_bar(data_dict, regions, title, blobs, blob_order,
                              settings.here, compartments, filename, compDsetNr)


def save_pure(comp_dict, save_dir):
    """
    Save as bed-files the poly(A) sites from the pure regions
    """

    for compartment, reg_dict in comp_dict.items():
        for region, slice_dict in reg_dict.items():
            for sl, slice_path in slice_dict.items():

                outname = '_'.join([compartment, region, sl])
                outpath = os.path.join(save_dir, outname)
                outhandle = open(outpath, 'wb')

                for line in open(slice_path, 'rb'):
                    outhandle.write(line)

                outhandle.close()

def count_these(some_dict):
    """
    Count PAS and all. You know its +1 already
    """

    allpas = set(['AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA',
                 'AATACA', 'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA',
                 'AATAGA'])

    data_dict = AutoVivification()

    for comp, comp_dict in some_dict.items():
        for reg, reg_dict in comp_dict.items():
            for key, keypath in reg_dict.items():

                # get how many of the isect have PAS
                data_dict['2+'][comp][reg][key]['all'] = 0
                data_dict['2+'][comp][reg][key]['PAS'] = 0
                data_dict['2+'][comp][reg][key]['PET'] = 0
                data_dict['2+'][comp][reg][key]['T'] = 0

                data_dict['3+'][comp][reg][key]['all'] = 0
                data_dict['3+'][comp][reg][key]['PAS'] = 0
                data_dict['3+'][comp][reg][key]['PET'] = 0
                data_dict['3+'][comp][reg][key]['T'] = 0

                for line in open(keypath, 'rb'):
                    (chrm, beg, end, name, covr, strand) = line.split('\t')

                    data_dict['2+'][comp][reg][key]['all'] +=1
                    t_info, PAS, PET_support = name.split('%')

                    if int(PET_support):
                        data_dict['2+'][comp][reg][key]['PET'] +=1

                    if t_info == 'T':
                        data_dict['2+'][comp][reg][key]['T'] +=1

                    has_pas = False
                    for pa in PAS.split('#'):
                        if pa in allpas:
                            has_pas = True

                    if has_pas:
                        data_dict['2+'][comp][reg][key]['PAS'] += 1

                    if int(float(covr)) > 2:
                        data_dict['3+'][comp][reg][key]['all'] +=1

                        if t_info == 'T':
                            data_dict['3+'][comp][reg][key]['T'] +=1

                        if has_pas:
                            data_dict['3+'][comp][reg][key]['PAS'] += 1

                        if int(PET_support):
                            data_dict['3+'][comp][reg][key]['PET'] +=1

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
    # 2) For each polyA site, make a dict with how many polyA reads are gotten

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
    dsetreads = get_dsetreads(settings, region='3UTR-exonic')
    clustercount = get_clustercount(settings, region='3UTR-exonic')

    cell_lines = ['GM12878', 'HeLa-S3', 'K562']
    # cell_line color dictionary
    colors = ['m', 'r', 'b', 'g', 'k']
    cell_cols = {}
    for indx, cline in enumerate(cell_lines):
        cell_cols[cline] = colors[indx]

    #speedrun = True
    speedrun = False

    wc_c = [ds for ds in settings.datasets if ((('Cytoplasm' in ds) or
                                                ('Whole_Cell' in ds) or
                                                ('Nucleus' in ds)) and
                                               (not 'Minus' in ds))]

    wc_c_minus = [ds for ds in settings.datasets if ((('Cytoplasm' in ds) or
                                                      ('Whole_Cell' in ds)
                                                      or ('Nucleus' in ds))
                                                     and ('Minus' in ds))]

    region = 'whole'
    for cell_line in cell_lines:

        wc_cCL = [ds for ds in wc_c if cell_line in ds]
        wc_c_minusCL = [ds for ds in wc_c_minus if cell_line in ds]

        if speedrun:
            wc_cCL = wc_cCL[:2]
            wc_c_minusCL = wc_c_minusCL[:2]

        data_grouping = {'Poly(A) plus': wc_cCL,
                         'Poly(A) minus': wc_c_minusCL}

        # keep a dictionary with reference to all the plots
        plotdict = {}

        for title, dsets in data_grouping.items():

            # sort the dsets in cell_lines by # of reads
            def mysorter(dset):
                return clustercount[dset]
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

        for yek in ['2+', '3+']:
            # set figure and axes
            (fig, axes) = plt.subplots(1,2)
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
                    countdict = countdict[yek]
                    # the the sum of reads from these datasets
                    x = [get_dsetreads(settings, '3UTR-exonic')[ds] for
                         ds in dsets.split(':')]

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
            filename = 'Saturation_plot_cell_lines_{0}'.format(yek)
            filepath = os.path.join(output_dir, filename+'.pdf')
            fig.savefig(filepath, format='pdf')
            filepath = os.path.join(output_dir, filename+'.eps')
            fig.savefig(filepath, format='eps', papertype='A4')

def strand_prediction(settings):
    """
    Simple: intersect your polyA sites with the annotated ones.

    RESULTS:

        polyAdb_gencode_merged_pure.bed
        Strand specificity of annotated strands: 0.96


        gencode_polyA.bed
        Strand specificity of annotated strands: 0.97


        polyA_db_proper.bed
        Strand specificity of annotated strands: 0.95

    """
    cell_lines = ['Whole_Cell', 'Cytoplasm', 'Nucleus']
    co = cell_lines

    speedrun = True
    #speedrun = False

    # 1) get all polyA + datasets
    plus_subset = [ds for ds in settings.datasets if (not 'Minus' in ds) and
                   ((co[0] in ds) or (co[1] in ds) or (co[2] in ds))]

    if speedrun:
        plus_subset = plus_subset[:2]

    # 1.1) write each poly(A) site to file with +/- 15 and strand
    batch_key = 'strand_speff'
    region = 'whole'

    dsets, super_3utr = super_falselength(settings, region, batch_key,
                                          plus_subset, speedrun)

    outdir = os.path.join(settings.here, 'analysis', 'strand_pred')
    outfile = 'all_pAs.bed'
    outpath = os.path.join(outdir, outfile)
    outhandle = open(outpath, 'wb')

    minlim = 2

    for utr_id, utr in super_3utr[region].iteritems():
        for cls in utr.super_clusters:

            if cls.nr_support_reads > minlim: # ?

                chrm = utr.chrm
                beg = str(cls.polyA_coordinate-15)
                end = str(cls.polyA_coordinate+15)
                strand = cls.strand

                outhandle.write('\t'.join([chrm, beg, end, strand])+'\n')

    outhandle.close()

    ann_dir = '/users/rg/jskancke/phdproject/3UTR/the_project/'\
            'annotated_polyAsites'

    anns = ['polyAdb_gencode_merged_pure.bed', 'gencode_polyA.bed',
            'polyA_db_proper.bed']

    ann2title = {'polyAdb_gencode_merged_pure.bed': 'Both annotations',
                 'gencode_polyA.bed': 'GENCODE V.7',
                 'polyA_db_proper.bed': 'PolyAdb_2'}

    for_latex = {}
    title = []
    data = []
    for ann in anns:
        annotation = os.path.join(ann_dir, ann)

        cmd = ['intersectBed', '-wa', '-wb', '-a', outpath, '-b', annotation]

        total = 0
        same = 0

        p = Popen(cmd, stdout=PIPE)
        for line in p.stdout:

            total += 1

            spl = line.split()

            if len(spl) == 8:
                strand1 = line.split()[3]
                strand2 = line.split()[7]
            if len(spl) == 10:
                strand1 = line.split()[3]
                strand2 = line.split()[9]

            if strand1 == strand2:
                same += 1

        title.append(ann2title[ann])
        data.append(format(same/total, '.2f'))
        print(ann)
        print('Strand specificity of annotated strands: {0:.2f}'\
              .format(same/total))
        print('\n')

    # because there is only one line
    data = [data]

    for_latex['caption'] = 'Capturing the strand of annotated poly(A) sites'
    for_latex['label'] = 'polyA_strand_capture'
    for_latex['header'] = title
    for_latex['data'] = data
    for_latex['template_path'] = '/users/rg/jskancke/phdproject/templates/'\
                                    'simple_table.tex'

    output_dir = os.path.join(settings.here, 'Results_and_figures',
                              'GENCODE_report', 'tables')

    for_latex['savedir'] = output_dir

    # create the latex table
    csv2latex.make_latex_table(for_latex)

def intergenic_finder(settings):
    """
    TODO make a table of this for the RNA seq part in the paper.
    """
    #speedrun = True
    speedrun = False
    region = 'Intergenic'

    outhandle = open(os.path.join(settings.here, 'analysis', 'utr_extension',
                                  '3utr_extension.tsf'), 'wb')

    header = '\t'.join(['Compartment', 'Extension', 'Found', 'Same strand (\%)',
                        'PAS', 'T', 'PET', 'Cufflinks']) + '\n'
    outhandle.write(header)

    #for comp in ['Whole_Cell', 'Cytoplasm', 'Nucleus', 'All']:
    for comp in ['All']:

        subset = [ds for ds in settings.datasets if (comp in ds) and
                  ('Minus' not in ds)]

        if comp == 'All':
            subset = [ds for ds in settings.datasets if ('Minus' not in ds)]

        if speedrun:
            subset = subset[:2]

        batch_key = 'intergenicFinder'
        dsets, super_3utr = super_falselength(settings, region, batch_key, subset,
                                              speedrun)
        Fdir = '/home/jorgsk/phdproject/3UTR/extended_3UTR'
        outdir = '/home/jorgsk/phdproject/3UTR/the_project/thesis_figures/thesis_stats'

        # use all cell lines and intergenic

        outfile = 'intergenic_pAsites'
        outpath = os.path.join(outdir, outfile)
        pAouthandle = open(outpath, 'wb')

        toosmall = 1
        # Write the intergenic ones to file
        for gulp_id, gulp in super_3utr[region].iteritems():

            for cls in gulp.super_clusters:

                if cls.nr_support_reads > toosmall: # ?

                    chrm = gulp.chrm
                    beg = str(cls.polyA_coordinate)
                    end = str(cls.polyA_coordinate)
                    strand = cls.strand
                    pas = 'pasno'
                    if cls.nearby_PAS[0] != 'NA':
                        pas = 'pasyes'

                    pet = 'petno'
                    if cls.PET_support:
                        pet = 'petyes'

                    cuff = 'cuffno'
                    if cls.cufflinks_support:
                        cuff = 'cuffyes'

                    name = '%'.join([pas, pet, cuff])
                    pAouthandle.write('\t'.join([chrm, beg, end, name,
                                               cls.tail_type, strand])+'\n')

        pAouthandle.close()
        # intersect the intergenic poly(A) sites with the extended
        # transcripts that go into the intergenic region
        if comp == 'Whole_Cell':
            comp = 'Whole Cell'

        for dnr, dist in enumerate([500, 1000, 5000]):
            intergenic = os.path.join(Fdir, 'trimmed_extended_3UTRs_{0}'.\
                                     format(dist))
            if dnr != 0:
                comp = ' '
            found = {}
            cmd = ['intersectBed', '-wa', '-wb', '-a', intergenic, '-b', outpath]
            p = Popen(cmd, stdout=PIPE)

            # store the number of poly(A) sites. Select only 
            total = 0
            countdict = {'same': {'total': 0, 'PAS': 0, 'T': 0, 'PET': 0,
                                  'Cufflinks': 0},
                         'other': {'total': 0, 'PAS': 0, 'T': 0, 'PET': 0,
                                   'Cufflinks': 0}}
            for line in p.stdout:
                (chma, bega, enda, ts_id, d, stranda, chrmb, begb, endb, name,
                 t_info, strandb) = line.split()

                (pas, pet, cuff) = name.split('%')
                total +=1

                if stranda == strandb:
                    keyw = 'same'
                    identical = 1
                else:
                    keyw = 'other'
                    identical = 0

                countdict[keyw]['total'] +=1

                if pas == 'pasyes':
                    countdict[keyw]['PAS'] +=1
                if t_info == 'T':
                    countdict[keyw]['T'] +=1
                if pet == 'petyes':
                    countdict[keyw]['PET'] +=1
                if cuff == 'cuffyes':
                    countdict[keyw]['Cufflinks'] +=1

                # see if you have from both sides with identical
                checkkey = chrm+'_'+begb
                if checkkey in found:
                    found[checkkey].append(identical)
                else:
                    found[checkkey] = [identical]

            fails = 0
            for key, hits in found.iteritems():
                if len(set(hits)) != 1:
                    fails +=1

            #print('Compartment: {0} Extension: {1}'.format(comp, dist))
            print('Extension: {0}'.format(dist))
            print('Total hits: {0}'.format(total))
            #print('p(A) sites hit by extensions from both sides: {0}'\
                  #.format(fails))
            for keyw in ['same', 'other']:
                keytotal = countdict[keyw]['total']
                if total > 0:
                    ratio = format(keytotal/total, '.2f')
                elif total == 0:
                    ratio = 'NaN'
                print('\t{3} strand: {1}, ratio: {2}'\
                      .format(total, keytotal, ratio, keyw.capitalize()))

                CUF = countdict[keyw]['Cufflinks']
                PAS = countdict[keyw]['PAS']
                PET = countdict[keyw]['PET']
                T = countdict[keyw]['T']
                print("""\t{6} strand with PAS: {0} ({1:.2f}),
                     with PET: {2} ({3:.2f}),
                     with Cuffl: {4} ({5:.2f}),
                     with T: {6} ({7:.2f})"""\
                      .format(PAS, PAS/keytotal,
                              PET, PET/keytotal,
                              CUF, CUF/keytotal,
                              T, T/keytotal,
                              keyw.capitalize()))
                print('')

            ## write to a tab-delimited file for ease of re-use.
            #Extension Total #Same_strand #Same_w PAS # #Same_w T
    outhandle.close()

def write_all_pA(settings):

    co = ['Whole_Cell', 'Cytoplasm', 'Nucleus']

    subset = [ds for ds in settings.datasets if (not 'Minus' in ds) and
                       ((co[0] in ds) or (co[1] in ds) or (co[2] in ds))]
    #speedrun = True
    speedrun = False
    if speedrun:
        subset = subset[:2]

    # 1.1) write each poly(A) site to file with +/- 15 and strand
    batch_key = 'PET'
    region = 'whole'

    speedrun = False
    dsets, super_3utr = super_falselength(settings, region, batch_key,
                                          subset, speedrun)

    temp_dir = os.path.join(settings.here, 'temp_files')
    superbed_path = os.path.join(temp_dir, 'all_pA.bed')

    handle = open(superbed_path, 'wb')

    for utr_name, utr in super_3utr[region].iteritems():

        for cls in utr.super_clusters:

            beg = cls.polyA_coordinate
            if cls.PAS_distance[0] != 'NA':
                pas = 1
            else:
                pas = 0

            entry = '\t'.join([utr.chrm, str(beg-20), str(beg+20),
                               str(pas), str(cls.nr_support_reads),
                               cls.strand])
            handle.write(entry + '\n')

    handle.close()

    return superbed_path


def storePetPA(settings):
    """
    """

    petDir = '/users/rg/jskancke/phdproject/3UTR/PET'

    pet_all = os.path.join(petDir, 'all_Pet_merged_1+.bed')
    pA_all = write_all_pA(settings)

    pet_paths = {}
    pA_paths = {}

    myRange = [(1,1),(2,2),(3,5),(5,10),(10,10000000)]

    #for r in myRange:
        #pA_path = os.path.join(petDir, 'pA_range')
        #for 

    # the ranges you want to report 
    # note: put this in a separate file
    # the output without verbose is exactly like this:
    #    Z Score:  64.2057283354
    #    p-value:  0.0


    counts = np.zeros((10,10))
    PAS = np.zeros((10,10)) # keep track of the PAS % in the various regions
    paths_PET = dict((val, []) for val in myRange)
    paths_pA = dict((val, []) for val in myRange)

def pet_intersection(settings, speedrun):
    """
    Compare poly(A) sites with Tags.
    Method: for each region in each compartment, show how many poly(A) sites are
    supported by TAGs

    You might have to include the PET just like another annotation in a way in
    utail. That's a hassle. But then you could use it much more nicely; for
    example you could make the bar plot easily. I'd give that a half a day of
    coding. And a day of running.

    First I'll make a quadrant of 1,2,3-5,5-10,10+ where I separate out the
    genomic locations of these things in bed-files. Then I'll run 25 jobs in
    parallel, computing the p-value of the co-llocation of the regions. That
    will be interesting.

    Finally, do the # of poly(A) sites and the # of sequence tags correlate?
    """

    # 1) Make a 4X4 array of correlation coefficients with bc. For this you need
    # to write to bedfiles all the PAS and PET with 1,2,3-5,5-10, 10+ into files
    # and run two for-loops over them, getting the Z and P values

    (petPaths, pAPaths) = storePetPA(settings)

    debug()

    ## 1.2) Intersect the clustered poly(A) sites
    #cmd = ['intersectBed', '-s', '-wa','-wb','-a', superbed_path, '-b', pet_1plus]
    #p = Popen(cmd, stdout=PIPE)

    #pA_count = []
    #PET_count = []

    #for line in p.stdout:
        #(chrmA, begA, endA, pas, nr_reads_pA, strandA,
         #chrmB, begB, endB, P, nr_reads_PET, strandB) = line.split()
        #pA_count.append(nr_reads_pA)
        #PET_count.append(nr_reads_PET)

    # parse like this [0] != 0 and [1] >1, ==1, > 2 etc

    # How many pA? How many intersect with pet? PAS % for all PAS and for those
    # that intersect with PET. What % of 1+ p(A) have pet, 2+ p(A) have Pet? 3+
    # p(A) have pet? etc. What about the p(A) that don't intersect with PET. Why
    # not? There is a spearman correlation, but the scatter lot shows the normal
    # crazyness. Even for the 1 sin PAS and 1 PET there is an unequivocal
    # correlation. You have two totally different p/z zcores depending if you
    # use pase pair hit or region hit. Without knowing exactly what's happening
    # it's difficult to know which parameter to chose.


def pet_vs_pa(settings, speedrun):
    """
    Intersect 3+ pA sites with 3UTRs (strand specific? :) ) and compare the
    specificity and sensitifity in your own way.

    Instead of comparing against all these 3UTR exonic regions, why don't you
    compare against only those that are expressed? Would be tons better. And
    tons more work.
    """

    pA_reads2 = os.path.join(settings.here, 'temp_files', 'all_pA2+.bed')
    pA_reads3 = os.path.join(settings.here, 'temp_files', 'all_pA3+.bed')
    pA_reads4 = os.path.join(settings.here, 'temp_files', 'all_pA4+.bed')
    pA_reads5 = os.path.join(settings.here, 'temp_files', 'all_pA5+.bed')
    pet_reads10 = '/users/rg/jskancke/phdproject/3UTR/PET/all_Pet_merged10+.bed'
    pet_reads20 = '/users/rg/jskancke/phdproject/3UTR/PET/all_Pet_merged20+.bed'
    pet_reads30 = '/users/rg/jskancke/phdproject/3UTR/PET/all_Pet_merged30+.bed'
    pet_reads40 = '/users/rg/jskancke/phdproject/3UTR/PET/all_Pet_merged40+.bed'

    three_utrs = '/users/rg/jskancke/phdproject/3UTR/the_project/'\
            'source_bedfiles/3UTR_exons.bed'
    #three_utrs = '/users/rg/jskancke/phdproject/3UTR/the_project/'\
            #'source_bedfiles/3UTR-exonic_stranded.bed'

    # count 3utrs
    utr_nr = sum((1 for line in open(three_utrs, 'rb')))

    for fin in [pA_reads2, pA_reads3, pA_reads4, pA_reads5, pet_reads10,
                pet_reads20, pet_reads30, pet_reads40]:

        filename = os.path.split(fin)[1]
        nr_hit_set = set([])
        utr_hit_set = set([])

        # count pa/pet
        nr = sum((1 for line in open(fin, 'rb')))

        cmd = ['intersectBed', '-wa', '-wb', '-s', '-a', fin, '-b', three_utrs]
        #cmd = ['intersectBed', '-wa', '-wb', '-a', fin, '-b', three_utrs]
        p = Popen(cmd, stdout = PIPE)

        for line in p.stdout:
            (chrmA, begA, begB, maybePAS, counts, strandA,
            chrmB, begB, begB, ts_id, val, strandB) = line.split()

            nr_hit_set.add('_'.join([chrmA, begA, strandA]))
            utr_hit_set.add(ts_id.split('_')[0])

        # adapt ts_id to not count mulit-exons
        # count nr that intersect
        nr_hit = len(nr_hit_set)
        utr_hit = len(utr_hit_set)

        print filename
        # nr and % of pA/pet mapping
        l1 = '{0} of {1} ({2:.2f}) intersected'.format(nr_hit, nr, nr_hit/nr)

        # nr and of 3UTR covered by pA/pet
        l2 = 'Covering {0} of {1} ({2:.2f}) 3UTRs'.format(utr_hit, utr_nr,
                                                      utr_hit/utr_nr)
        print l1
        print l2
        print(' ')

def regionalizer(settings):
    """
    Show the distribution of cufflink ends, PET ends, and poly(A) ends in the
    different genomic regions. First plot shows counts of cluster, the second
    plot shows the basepairs of clusters. The second row in both figures show
    numbers normalized to the size of the genomic regions.
    """

    utrDir = '/users/rg/jskancke/phdproject/3UTR/'
    pA_ends = os.path.join(settings.here, 'temp_files', 'all_pA3+.bed')
    cuffLends = os.path.join(utrDir,'CUFF_LINKS/cufflinks_3UTR_ends_merged.bed')
    petEnds = os.path.join(utrDir, 'PET/all_Pet_merged10+.bed')

    dfiles = ['Poly(A)', 'PET', 'Cufflinks']
    dfile2path = {'Poly(A)': pA_ends, 'PET': petEnds, 'Cufflinks': cuffLends}

    #regions = ['5UTR-exonic', '5UTR-intronic', '3UTR-exonic', '3UTR-intronic',
               #'CDS-exonic', 'CDS-intronic', 'Nocoding-exonic',
               #'Noncoding-intronic', 'Intergenic']
    regions = ['3UTR-exonic', 'Intergenic', 'CDS-exonic', 'CDS-intronic',
               'Nocoding-exonic', 'Noncoding-intronic', '5UTR-exonic',
               '5UTR-intronic', '3UTR-intronic']

    # if not made, make a region-basepairs dict
    pickfile = 'region_PICKME'
    # if not pickled, calculate
    if not os.path.isfile(pickfile):
        regnorm = dict((reg, 0) for reg in regions)

    # Intersect the above with all regions
    regionDir = os.path.join(settings.here, 'genomic_regions_hg19')
    reg2path = {}
    for fpath in glob.glob(regionDir+'/*'):
        for reg in regions:
            if os.path.split(fpath)[1].startswith(reg):
                reg2path[reg] = fpath

                # if you haven't aldready made this dictionary
                if not os.path.isfile(pickfile):
                    for line in open(fpath):
                        (chrm, beg, end, d,d,d) = line.split()
                        bps = int(end)-int(beg)
                        regnorm[reg] += bps

    ## get the pcikeld file if you have aldready calcluateld
    if not os.path.isfile(pickfile):
        pickle.dump(regnorm, open(pickfile, 'wb'))
    else:
        regnorm = pickle.load(open(pickfile, 'rb'))

    data_dict = AutoVivification()

    pickfile = 'PICKME'
    # if not pickled, calculate
    if not os.path.isfile(pickfile):

        for dname, dfile in dfile2path.items():
            for reg in regions:

                cmd = ['intersectBed', '-a', dfile, '-b', reg2path[reg]]
                p = Popen(cmd, stdout=PIPE)
                data_dict['basepairsXn'][dname][reg] = 0
                data_dict['counts'][dname][reg] = 0

                for line in p.stdout:
                    (chrm, beg, end, d, n, m) = line.split()

                    data_dict['counts'][dname][reg] += 1

                    basepairs = int(end) - int(beg)
                    data_dict['basepairsXn'][dname][reg] += basepairs*int(n)

        pickle.dump(data_dict, open(pickfile, 'wb'))

    # if pickled, read from file
    else:
        data_dict = pickle.load(open(pickfile, 'rb'))

    # onw where you normalize for region size and one where you don't
    for typ_nr, typ in enumerate(['counts', 'basepairsXn']):
        dset_dict = data_dict[typ]
        fig, axes = plt.subplots(2,3)

        for regNorm_nr, regNormalized in enumerate(['', 'region_normalized']):

            for dset_nr, dset in enumerate(dfiles):
                regDict = dset_dict[dset]

                # take into account normalization for region AND
                # normalization for total basepairs in 'pA, cuff, and 
                pos = np.arange(1, len(regDict)+1)
                data = []
                names = []
                for reg_nr, reg in enumerate(regions[::-1]):
                    bpairs = regDict[reg]

                    if regNormalized:
                        data.append(bpairs/regnorm[reg])
                    else:
                        data.append(bpairs)

                    names.append(reg)

                ax = axes[regNorm_nr, dset_nr]

                if regNormalized:
                    c = 'g'
                else:
                    c = 'm'

                ax.barh(bottom=pos, width=data, height=1, color=c)

                if dset_nr == 0:
                    ax.set_yticks(pos+0.5)
                    ax.set_yticklabels(names, size=9)
                else:
                    ax.set_yticklabels([])

                ax.set_xticklabels([])

                if regNorm_nr == 0:
                    ax.set_title(dset, size=12)

                if regNorm_nr == 0 and dset_nr == 2:
                    xax2 = ax.twinx()
                    xax2.set_yticklabels([])
                    xax2.set_xticklabels([])
                    xax2.set_ylabel('Not normalized', size=8)

                if regNorm_nr == 1 and dset_nr == 2:
                    xax2 = ax.twinx()
                    xax2.set_yticklabels([])
                    xax2.set_xticklabels([])
                    xax2.set_ylabel('Normalized by region size', size=8)

        fig.subplots_adjust(left=0.23)
        fig.set_size_inches(8,4)

        if typ == 'counts':
            filename = 'Counts_of_clusters_in_regions'
        else:
            filename = 'Basepairs_of_clusters_in_regions'

        output_dir = os.path.join(settings.here, 'Results_and_figures',
                                  'GENCODE_report', 'Figures')

        filepath = os.path.join(output_dir, filename+'.pdf')
        fig.savefig(filepath, format='pdf')
        filepath = os.path.join(output_dir, filename+'.eps')
        fig.savefig(filepath, format='eps', papertype='A4')

def gencode_report(settings, speedrun):
    """ Make figures for the GENCODE report. The report is an overview of
    evidence for polyadenylation from the Gingeras RNA-seq experiments.

    Mystery: underreporting of polyA sites. I get smth like 40 percent when 75
    percent have it in the raw data.

    So far you've put numbers on several plots. Just make a table with the
    number of extended 3UTR sites and their strand, and I think that should be
    enough. Spend a few days discussing about it. In the meanwhile, you can fix
    more plots?
    """
    speedrun = True
    #speedrun = False

    # reuse plot_data from previous simulation
    # this saves precious time if you just want to change some detail in a plot

    reuse_data = True
    #reuse_data = False

    # TODO interesting for 3'utrs. Before, how long were they on average -- how
    # much were they extended by?
    # How many poly(A) sites without PAS also have PET?

    # 0.1) Core stats selected regions + genome
    #cumul_stats_printer(settings, speedrun)
    # 0.2) Core stats All regions All cell lines
    #cumul_stats_printer_all(settings, speedrun)

    # 1) nr of polyA sites obtained with increasing readnr
    # XXX seems to work. Even put letters there.
    #clusterladder(settings, speedrun, reuse_data)
    #clusterladder_cell_lines(settings)

    # 2) venn diagram of all the poly(A) sites from the whole genome for 3UTR and
    # venn_polysites(settings, speedrun)
    # Then: correct the slides with the numbers. Maybe make a table.
    # Then: cufflinks. Shit this cufflinks.
    # Maybe you have to add another 'pet' --> cufflinks data. From all the
    # cufflinks transcripts, write out the 3'end. Specifically write out the +/-
    # 50 nt around the 3'end. Put this in a BED-file and check it with your
    # poly(A) reads.
    # Then you should make two different plots --> one for 'annotated' poly(A)
    # reads and another for 'non-annotated' poly(A) reads, where you show the
    # cufflinks support.
    # ...
    # IDEA: For each compartment and for all:
        # region:
            #found in gencode (PAS %) (PET%)
            #novel found with PET (PAS% PET%)
            #novel with cufflinks (PAS% PET%)
            # the easiest thing would be if 

    # 3) plot of correlation between number of poly(A) sites expressed
    # transcript 3UTRs and the RPKM of the 3UTRs.
    #rpkm_polyA_correlation(settings, speedrun)
    # still no correlation. (and you should have kept the old code)

    # 4) Sidewise bar plots of the number of poly(A) incidents in the different
    # genomic regions, for poly(A)+ and poly(A)-
    #side_plot(settings, speedrun, reuse_data)

    # 5) Poly(A)+ pure, intersection, poly(A)-
    intersection_sideplot(settings, speedrun, reuse_data)

    # 6) All side plot!
    standard_sideplot(settings, speedrun, reuse_data)

    # 6) Cytoplasmic vs. nuclear RPKMS in 3UTR exonic exons
    #cytonuclear_rpkms_genc3(settings) # negative result
    #transcripts

    # 7) Demonstrate that you can predict the strand with 3UTR
    #strand_prediction(settings)

    # 8) Show what percentage of the intergenic ones are found within 1000 nt of
    # annotated transcript ends (and, preferably, back this up with read
    # coverage) Got results without coverage. With PET. With Cufflinks.
    #intergenic_finder(settings)

    # 9) How does the PET data fit into all this? Use both good and not-good PET
    #pet_intersection(settings, speedrun)
    # It is clear that the PET needs to be merged with more overlap. They don't
    # map to the 3' end specifically, just to the 3UTR generally.

    # 10) Compare the % of intersection of PET-3UTR and poly(A) reads
    #pet_vs_pa(settings, speedrun)

    # 11) Compare region distributio for poly(A), PET, and Cufflinks
    # actually count base-pair overlap, and report how much this is for the
    # different regions in terms of total size. You will see that PET is less
    # 3UTR-focused than p(A) data. And the same for cufflinks?
    # This could be that the annotated poly(A) sites are mostly mRNA?
    #regionalizer(settings)

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

def write_gff(settings, toosmall, minus, cell_lines, speedrun, expandby):
    """
    File must be in compliance with Sarah's file for transcription start sites.

  - field no 1: Poly(A) cluster chromosome
  - field no 2: 'encodecrg'
  - field no 3: 'TER/TTS' (transcript end region)/(transcription termination site)
  - field no 4: Poly(A) cluster start (1-based)'
  - field no 5: Poly(A) cluster end (1-based)
  - field no 6: '.'
  - field no 7: Poly(A) cluster strand
  - field no 8: '.'
  - field no 9: list of (key value1,value2,value3,...) pairs:
      * Cluster class: (Gencode/GenSeq/RnaSeqOnly)
      * Cluster center: value
      * Cluster total TPM: value
      * Expression by tpm in GM12878
      * Expression by tpm in ...
      *
      * List of Gencode TTS cluster coordinates ('.' for RnaSeqOnly)
      * list of Gencode TTS cluster genes ('.' for RnaSeqOnly)

    NOTE TO SELF:

    To be able to do the above I must modify the information I keep during
    clustering. For example, I must keep the positions of the inidividual polyA
    sites during clustering, or at least the min/max. I must also keep
    information about how many reads have been found for each cell
    line/fraction. Also obtain the number of reads for each cell line /
    fraction.

    Recalleth that they proposed another measure: poly(A) reads divided by nr of
    other reads covering this site. I should pick a spot 10 nt downstream the
    site actually (or the downstream end of the cluster). Will you have to use
    bedtools to get this number out of the bam files? it will take-ya forever.
    """

    co = cell_lines
    # 1) get all polyA- datasets
    if minus:
        subset = [ds for ds in settings.datasets if ('Minus' in ds) and
                       ((co[0] in ds) or (co[1] in ds) or (co[2] in ds))]

    # 2) or polyA+ datasets
    if not minus:
        subset = [ds for ds in settings.datasets if (not 'Minus' in ds) and
                       ((co[0] in ds) or (co[1] in ds) or (co[2] in ds))]

    # reduce the nr of datasets
    if speedrun:
        subset = subset[:2]

    # 1.1) write each poly(A) site to file +/- expandby
    batch_key = 'cuffer'
    region = 'whole'

    dsets, super_3utr = super_falselength(settings, region, batch_key,
                                          subset, speedrun)

    outdir = os.path.join(settings.here, 'genome_wide_dir')

    if speedrun:
        outfile = 'all_pAs_speedrun.gff'
    else:
        outfile = 'all_pAs.gff'

    outpath = os.path.join(outdir, outfile)
    outhandle = open(outpath, 'wb')

    for utr_id, utr in super_3utr[region].iteritems():
        for cls in utr.super_clusters:

            if cls.nr_support_reads > toosmall:

                core_fields = '\t'.join([utr.chrm,
                                           'encodeCRG',
                                           'TER/TTS',
                                           cls.beg,
                                           cls.end,
                                           '.',
                                           cls.strand,
                                           '.'])


                clusterClass = 'clusterClass: RnaSeqOnly'
                clusterCenter = 'clusterCenter: {0}'.format(cls.polyA_coordinate)
                clusterTPM = 'clusterTotalReads: {0}'.format(cls.nr_support_reads)

                dsetcount = ' '.join([' '.join(t) for t in zip(cls.dsets, cls.coverages)])

                keyvals = ' '.join([clusterClass, clusterCenter, clusterTPM,
                                   dsetcount])

                outhandle.write(core_fields + '\t' + keyvals + '\n')

    outhandle.close()

    return outpath

def merge_polyAs(settings, toosmall, minus, cell_lines, speedrun, expandby):

    co = cell_lines
    # 1) get all polyA- datasets
    if minus:
        subset = [ds for ds in settings.datasets if ('Minus' in ds) and
                       ((co[0] in ds) or (co[1] in ds) or (co[2] in ds))]

    # 2) or polyA+ datasets
    if not minus:
        subset = [ds for ds in settings.datasets if (not 'Minus' in ds) and
                       ((co[0] in ds) or (co[1] in ds) or (co[2] in ds))]

    # reduce the nr of datasets
    if speedrun:
        subset = subset[:2]

    # 1.1) write each poly(A) site to file +/- expandby
    batch_key = 'merge_all'
    region = 'whole'

    dsets, super_3utr = super_falselength(settings, region, batch_key,
                                          subset, speedrun)

    outdir = os.path.join(settings.here, 'genome_wide_dir')

    if speedrun:
        outfile = 'all_pAs_speedrun.bed'
    else:
        outfile = 'all_pAs.bed'

    if minus:
        outfile = 'minus_'+ outfile
    else:
        outfile = 'plus_'+ outfile

    outpath = os.path.join(outdir, outfile)
    outhandle = open(outpath, 'wb')

    for utr_id, utr in super_3utr[region].iteritems():
        for cls in utr.super_clusters:

            if cls.nr_support_reads > toosmall:

                chrm = utr.chrm
                beg = str(cls.polyA_coordinate-expandby)
                end = str(cls.polyA_coordinate+expandby)
                dsets = '#'.join(cls.dsets)
                pas = '#'.join(cls.nearby_PAS)
                annot_dist = str(cls.annotated_polyA_distance)
                covr = str(cls.nr_support_reads)
                strand = cls.strand

                outhandle.write('\t'.join([chrm, beg, end, annot_dist, dsets, covr, pas,
                                           strand])+'\n')

    outhandle.close()

    return outpath

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
        if float(exon_nr) == 0.5:
            debug()

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

def gencode_cufflinks_report(settings):
    """
    Same as for cufflinks, but for gencode
    """

    # 1) and 2) merge all polyA files (not poly(A) minus)
    #minus = True
    minus = False

    expandby = 0

    mincovr = 2
    compartments = ['Whole_Cell', 'Cytoplasm', 'Nucleus']

    speedrun = True
    #speedrun = False

    polyA_path = merge_polyAs(settings, mincovr, minus, compartments,
                                     speedrun, expandby)
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
    cmd = ['intersectBed', '-wa', '-wb', '-a', gencends, '-b', polyA_path]
    p = Popen(cmd, stdout=PIPE)

    goodpas = set(['AATAAA', 'ATTAAA'])

    allpas = set(['AATAAA', 'ATTAAA', 'TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA',
                 'AATACA', 'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA',
                 'AATAGA'])

    # a dictionary that counts the number of transcripts and if there is a PAS
    # or not. When a new transcript type is identified, a new dictionary is made
    # for it
    subdict = {'PAS': set([]), 'Good PAS': set([]), 'No PAS': set([]), 'Cov':[]}

    event_counter = {'all': deepcopy(subdict)}

    for line in p.stdout:
        (chrm, beg, end, ts_type, ts_id, strand, d,d,d,
                                         pas, cov, d) = line.split()

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
    output_dir = os.path.join(settings.here, 'Results_and_figures',
                              'GENCODE_report', 'csv_files')

    output_path = 'gencode_transcript_types_{0}_{1}_Minus{2}\
                          '.format(mincovr, 'all+celllines', str(minus))

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

        outkey = ' '.join(key.split('_')).capitalize()
        outcov = format(mean_cov, '.0f')

        gencsum_handle.write('{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4:.2f}\t{5}\n'\
                             .format(outkey, nr_found, pcnt_of_an, pcnt_PAS,
                                     pcnt_Good_PAS, outcov))
    gencsum_handle.close()


def get_cufflink_ends(cuff_dir, model, settings):
    """
    Return the region +/- 50 nucleotides around the 3'UTR of cufflink end points
    """
    novel_ext = os.path.join(cuff_dir, 'deNovoCuffTrs.merge64.'\
                             'filteredG7.trTypes.txt')
    # get dicts for the cufflinks model
    genes, transcripts, exons = get_gentrex(novel_ext, model)

    ending = os.path.splitext(os.path.split(model)[1])[0].split('_')[-1]
    cuffends = os.path.join(settings.here, 'source_bedfiles',
                            'cufflinks_ends_{0}.bed'.format(ending))

    outhandle = open(cuffends, 'wb')

    for ex_id, (chrm, beg, end, exon_nr, ts_id, gene_id, lst) in exons.iteritems():

        # get information from the transcript
        nov_ext, ts_exon_nr, ts_strand, gene_id = transcripts[ts_id]

        if int(exon_nr) == ts_exon_nr:
            if ts_strand == '+':

                beg = str(int(end)-50)
                end = str(int(end)+50)

                outhandle.write('\t'.join([chrm,
                                           beg,
                                           end,
                                           ts_id,
                                           str(nov_ext),
                                           ts_strand+'\n']))
        if int(exon_nr) == 0:
            if ts_strand == '-':

                beg = str(int(beg)-50)
                end = str(int(beg)+50)

                outhandle.write('\t'.join([chrm,
                                           beg,
                                           end,
                                           ts_id,
                                           str(nov_ext),
                                           ts_strand+'\n']))
    outhandle.close()

def new_cufflinks_report(settings, speedrun=False):
    """ Basic statistics on the cufflinks data.
    1) How many p(A)? How many overlap annotated?
    2) How many novel?
    3) How many fall very close to 3' ends?
    4) We verify that XXX of these genes are polyadenylated (min 1 transcript)

    1) Get poly(A)s for the whole genome for all the cell lines 
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
    expandby = 0
    toosmall = 1
    minus = False
    cell_lines = ['Whole_Cell', 'Cytoplasm', 'Nucleus']
    speedrun = True

    polyA_path = merge_polyAs(settings, toosmall, minus, cell_lines,
                                     speedrun, expandby)

    juncfree_pA_path = polyA_path + '_juncfree'

    # remove areas around exon junctions, because it causes biases
    junctions = os.path.join(settings.here, 'junctions',
                             'splice_junctions_merged.bed')

    # i.5) cut away the exon-exon noise
    cmd = ['subtractBed', '-a', polyA_path, '-b', junctions]
    p = Popen(cmd, stdout=open(juncfree_pA_path, 'wb'))
    p.wait()

    cuff_dir = '/users/rg/jskancke/phdproject/3UTR/CUFF_LINKS'
    model = cuff_dir+ '/CuffCode_CSHL_exons_sans_GENCODE7_fromAngelika.gtf'
    # NOTE doesn't work for model2. parsing error.

    cuffends = get_cufflink_ends(cuff_dir, model, settings)

    # 3) intersect the bed_model with the merged polyA sites
    cmd = ['intersectBed', '-wo', '-a', cuffends, '-b', juncfree_pA_path]
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

        cmd = ['intersectBed', '-wa', '-wb', '-a', expanded_dsetpath,
               '-b', hagen_path]

        found = set([])

        p = Popen(cmd, stdout=PIPE)

        for line in p.stdout:
            (achrm, abeg, aend, d, covrg, astrand, chrm, beg, end, strnd)\
                    = line.split()

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

def get_ts_for_hag(model, settings, save_dir):
    #1) Get all distant pA sites that match hagen's 
    # Pickle timesaver
    pickfile = 'PICKME'
    if not os.path.isfile(pickfile):
        pA_sites_dict = pas_for_hagen(settings)
        pickle.dump(pA_sites_dict, open(pickfile, 'wb'))
    else:
        pA_sites_dict = pickle.load(open(pickfile, 'rb'))

    an_frmt = 'GENCODE'
    import annotation_parser as annparse

    (transcripts, genes) = annparse.make_transcripts(model, an_frmt)

    end_sites = {}
    for ts_id, ts in transcripts.iteritems():

        #  the ends
        if ts.strand == '+':
            ex = list(ts.exons[-1])
            ex[1] += 1 # 'your' exons have 'beg'-1 compared to pipeline
            end_key = '_'.join([str(c) for c in ex])
        else:
            ex = list(ts.exons[0])
            ex[1] += 1 # 'your' exons have 'beg'-1 compared to pipeline
            end_key = '_'.join([str(c) for c in ex])

        end_sites[end_key] = ts_id

    # 2) Each of hagen-agreeing poly(A) site must know which transcript it
    # belongs to:

    # 2.1) Write out end-sites to bedfile
    endfile = 'end_sites.bed'
    endpath = os.path.join(save_dir, endfile)
    endhandle = open(endpath, 'wb')
    for end_site, ts_id in end_sites.iteritems():
        (chrm, beg, end, strand) = end_site.split('_')
        endhandle.write('\t'.join([chrm, beg, end, '0', ts_id, strand])+'\n')
    endhandle.close()

    # 2.2) Write out poly(A) sites to bedfile
    pAfile = 'polyA_sites_forexonmerge.bed'
    pApath = os.path.join(save_dir, pAfile)
    pAhandle = open(pApath, 'wb')
    for end_site, cell_line_count in pA_sites_dict.iteritems():
        (chrm, beg, end, strand) = end_site.split('_')
        pAhandle.write('\t'.join([chrm, beg, end, strand])+'\n')
    pAhandle.close()

    # 2.3) Expand the end-sites
    expandpath = add_ending(endpath, 'expanded')
    cmd1 = ['slopBed', '-b', '15', '-i', endpath, '-g', settings.hg19_path]
    p1 = Popen(cmd1, stdout = open(expandpath, 'wb'))
    p1.wait()

    # 2.4) Intersect
    # NOTE there is some kind of bug with intersectBed and the -s option.
    # Solution: run without
    cmd2 = ['intersectBed', '-wa', '-wb', '-a', expandpath, '-b', pApath]
    p2 = Popen(cmd2, stdout=PIPE)

    # 2.5) Read into a dict those poly(A) sites that intersect and the
    polyA2ts = {}

    for line in p2.stdout:
        (chra, bega, enda, d, ts_id, standa, chrb,
         begb, endb, strandb) = line.split()

        key = '_'.join([chrb, begb, endb, strandb])
        if key in polyA2ts:
            polyA2ts[key].append(ts_id)
        else:
            polyA2ts[key] = [ts_id]

    # fill up the empty polyAs -- the ones that don't have a ts
    for pA, coutns in pA_sites_dict.iteritems():
        if pA not in polyA2ts:
            polyA2ts[pA] = []

    return polyA2ts

def hag_rpkm_refac(settings):
    """
    New version.
    """
    model = '/users/rg/jskancke/phdproject/3UTR/Annotations/'\
            'gencode.v3c.annotation.GRCh37.gtf'
    save_dir = os.path.join(settings.here, 'hagen_stuff')

    #chr1 = True
    chr1 = False
    if chr1:
        model = '/users/rg/jskancke/phdproject/3UTR/Annotations/'\
                'gencode.v3c.annotation.GRCh37_chr1.gtf'

    # 1) Get a gm128[ts_id] = rpkm dict (for k5 too).
    (ts_2_rpkm_GM12, ts_2_rpkm_K56) = get_ts_rpkms_genc3()

    # 2) Write out the ends of the transcripts and intersect them with hagen's
    hag_2_ts = get_ts_for_hag(model, settings, save_dir)

    # XXX if more than 1 transcript lands there, hagen wants the sum.

    outname = 'pA_GM12_K562_transcript_RPKMS'
    outhandle = open(os.path.join(save_dir, outname), 'wb')
    for pA, ts_ids in hag_2_ts.iteritems():

        gm_rpkm = 0
        k6_rpkm = 0

        for ts_id in ts_ids:
            if ts_id in ts_2_rpkm_GM12:

                grpkm = ts_2_rpkm_GM12[ts_id]
                krpkm = ts_2_rpkm_K56[ts_id]

                if not math.isnan(grpkm):
                    gm_rpkm += grpkm

                if not math.isnan(krpkm):
                    k6_rpkm += krpkm

        outhandle.write('\t'.join([pA, str(k6_rpkm), str(gm_rpkm)]) + '\n')

    outhandle.close()

    # ends: this gives a hagen_end -> transcript link. Yes, there may be
    # multiple transcripts for each end.

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

def enders(beg, end, nr_comp):

    beg, end = int(beg), int(end)
    length = end - beg
    interval_length = math.ceil(length/nr_comp)

    checkpoints = [beg+interval_length*n for n in range(nr_comp-1)]

    return checkpoints + [end]


def EML(settings):
    """
    Simply count the number and cluster number of poly(A) site in early, medium,
    andl late 3UTR, divided by 3.

    Maybe the extension ruins this? Extend by 100!

    # Select only those that have 2/3 polyA sites and all are found in both
    # nuclear and cytoplasmic.
    """

    region = '3UTR'

    nr_comp = 3

    counters = [0 for n in range(nr_comp)]

    #speedrun = True
    speedrun = False

    data_dict = {'Cytoplasm': {'count': deepcopy(counters),
                               'magnitude': deepcopy(counters)},
                'Nucleus': {'count': deepcopy(counters),
                            'magnitude': deepcopy(counters)}}

    for comp in ['Cytoplasm', 'Nucleus']:
        subset = [ds for ds in settings.datasets if comp in ds and
                     ('Minus' not in ds)]

        #subset = subset[:2] # XXX REMOVE

        batch_key = 'EML'
        dsets, super_3utr = super_falselength(settings, region, batch_key,
                                              subset, speedrun)

        for utr_name, utr in super_3utr[region].iteritems():

            if len(utr.super_clusters) != nr_comp:
                continue

            for cls in utr.super_clusters:

                # divide the 3UTR into X pieces, and see where the poly(A) sites
                utrEnd, utrBeg = int(utr.end), int(utr.beg)

                # skip those not on the same strand
                if cls.strand != utr.strand:
                    break

                if cls.nr_support_reads < 2:
                    break

                checkpoints = np.linspace(utrBeg, utrEnd, nr_comp)[1:]

                this_pA = cls.polyA_coordinate

                # TODO something is strange, the distribution is turned on its
                # head; lowhighlow. fixaftermeeting.
                for c_nr, cpoint in enumerate(checkpoints):

                    if utr.strand == '+':
                        index = c_nr
                    if utr.strand == '-':
                        index = nr_comp - (c_nr+1)

                    if this_pA < cpoint:
                        data_dict[comp]['count'][index] += 1
                        data_dict[comp]['magnitude'][index] += cls.nr_support_reads
                        break

    p = Plotter()
    p.CyNucComp(data_dict, nr_comp, settings)

def write_all_polyAs(settings):
    """
    Simply merge all polyA sites
    """

    #minus = False
    minus = True

    expandby = 0

    covr_more_than = 1 # minimum 2 in coverage
    compartments = ['Whole_Cell', 'Cytoplasm', 'Nucleus']

    speedrun = True
    #speedrun = False
<<<<<<< HEAD

    #gff_path = write_gff(settings, covr_more_than, minus, compartments, speedrun, expandby)

    # make gff/bed files of all 3UTR ends and poly(A) ends from the Gencode
    # annotation (obtained from EST data). Merge this information with the
    # information from the poly(A) investigation.
    gencode_path = get_gencode_ends()

def get_gencode_ends():
    """
    Parse two Gencode files. One about poly(A) sites specifically, and one about
    3UTR ends in general.

    poly(A) sites:
        1) get only the polyA_site entries
        2) merge them with +/- 15 (slopbed + merge bed)

    transcripts:
        1) get only transcripts (not genes, exons, utrs)
        2) get the transcript 3' ends
        3) merge them with +/- 15
        4) write to file (keep each transcript ID and gene ID)

    Then you have to merge these two files. The output will be in terms of the
    merged clusters. It will look like the transcript 3' end file but with
    poly(A) evidence as an adjunct.

    Return the path of this final merged object


    UPDATE: this is not what Julien wants.
    """

    #### PolyA sites ####
    polyA = BedTool('/home/jorgsk/phdproject/3UTR/gencode/gencode.v13.polyAs.gtf')

    # 1) get only polyA_site entries 
    polyA_pure = BedTool(e for e in polyA if e.fields[2] == 'polyA_site')

    # 2) merge them +/- 15
    merged_polyA = polyA_pure.slop(g=chromsizes('hg19'), b=15).merge(s=True,
                                                                     n=True)

    #### Transcript 3' end-sites ####
    #ant = BedTool('/home/jorgsk/phdproject/3UTR/gencode/gencode.v13.annotation.gtf')
    ant = BedTool('/home/jorgsk/phdproject/3UTR/gencode/gencode.v13.annotation_chr1.gtf')

    # 1) get only transcripts
    ts = BedTool(e for e in ant if e.fields[2] == 'transcript')

    # 2) make a new entry from just the three prime end site
    threEnd = '/home/jorgsk/phdproject/3UTR/gencode/transcript_ends.bed'
    end_handle = open(threEnd, 'wb')

    for entry in ts:

        chrm, bs, ts, ts_beg, ts_end, dot, strand, dot, rest = entry.fields

        restSplit = rest.split()

        # extract the geneID and transcriptID
        geneID, tsID = restSplit[1], restSplit[3]
        geneID = geneID[1:-2]
        tsID = tsID[1:-2]

        name = '#'.join([geneID, tsID])

        if strand == '+':
            beg = ts_end
            end = str(int(ts_end) + 1)

        if strand == '-':
            beg = ts_beg
            end = str(int(ts_beg) + 1)

        end_handle.write('\t'.join([chrm, beg, end, name, '.', strand]) + '\n')

    end_handle.close()

    # don't merge transcripts; you need the info on each one.
    ts_ends = BedTool(threEnd)

    #5) intersect with polyA file and create a final bed-entry from them

    # will be written in genome_wide_dir as all_polyA
    merge_polyAs(settings, covr_more_than, minus, compartments, speedrun, expandby)

def main():
    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    # Directory paths for figures and where the output lies
    (savedir, outputdir) = [os.path.join(here, d) for d in ('figures', 'output')]

    # Keep houskeeping information in a settings object
    settings = Settings(os.path.join(here, 'UTR_SETTINGS'), savedir, outputdir,
                        here, chr1=False)

    #gencode_report(settings, speedrun=False)

    # early, medium, late poly(A) in cytoplasm and nucleus
    #EML(settings)

    # XXX cufflinks report
    #new_cufflinks_report(settings, speedrun=False)

    # XXX same as cufflnksm but for gencode
    #cell_lines = ['All_Cell_Lines', 'GM12878', 'HEPG2', 'HUVEC', 'HeLa-S3',
                  #'K562']
    # NOTE must fix for individual cell lines
    #gencode_cufflinks_report(settings)

    write_all_polyAs(settings)

    # Hagen's stuff
    #hagen(settings, speedrun=False) # negative results for your pA
    #hag_rpkm_refac(settings) # trying with transcript RPKMS


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


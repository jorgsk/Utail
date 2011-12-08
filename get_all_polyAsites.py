"""
This script clusters all polyA sites you find and outputs a bedfile with the
regions -60, +10 for these clustered sites.

This bedfile should later be used to obtain the uniquely mapping reads from the
rna-seq experiments that land there.

You should run gem -> bed, with the seq as name, and then do intersect-bed with
the -60, +10 polyA-sites you have found, and you must demand that the polyA-site
overlaps with say, 40%?

Now, you don't know if reads have been reverse transcribed yet. However, this
comes at a later stage.

At that stage, you should try both directions of the read to see which overlaps
with the hg19 sequence in local alignment; when you know that, you can keep it
it that strandedness. Is there a better way to do this? How about sam/bam?

Simply do gem -> bed and bed-intersect; for each dataset you end up with the
relevant reads. Later you can make cell_line specific files just for the polyA
sites of each specific cell line.

This program returns 2 things: 1) a bedfile with the regions around pA sites and
2) the hg19_sequences of those regions in 5->3
"""

import os
from results import Settings
import math

from IPython.Debugger import Tracer
debug = Tracer()

def get_pA_centers(region, settings):

    chrms = range(1,23) + ['M', 'X', 'Y']
    # for storing the pA sites
    genome = {'+': dict(('chr' + str(v), []) for v in chrms),
              '-': dict(('chr' + str(v), []) for v in chrms)}

    center_clusters = []

    # put all the polyA sites in lists unique to strand and chromosome
    for dset, pA_path in settings.only_files(region).items():

        handle = open(pA_path, 'rb')
        header = handle.next()

        for pA_entry in handle:
            (chrm, beg, end, utr_ID, strand, pA_coord, pA_coor_strand) = pA_entry.split()[:7]

            genome[pA_coor_strand][chrm].append(pA_coord)

    # sort each list and then cluster the list
    for strand, chrm_dict in genome.items():
        for chrm, pA_sites in chrm_dict.items():
            pA_sites.sort(key=int)

            # we're using chr1 in the beginning
            if chrm != 'chr1':
                continue

            # all clusters
            clusters = []

            # initial values for clustering
            mean = int(pA_sites[0])
            clustsum = mean
            clustcount = 1
            this_cluster = [mean]

            for val in pA_sites[1:]:
                ival = int(val)

                # If dist between new entry and cluster mean is < 20, keep in cluster
                if abs(ival - mean) < 24:
                    this_cluster.append(ival)

                    clustsum = clustsum + ival
                    clustcount += 1
                    mean = clustsum/clustcount

                else: # If not, start a new cluster, and save the old one
                    clusters.append(this_cluster)
                    clustsum = ival
                    clustcount = 1
                    this_cluster = [ival]
                    mean = ival

            # append the last cluster
            clusters.append(this_cluster)

            # take the mean of the subclusters
            cl_mean = [int(math.floor(sum(clus)/float(len(clus)))) for clus in clusters]

            # add each center to the pyfast dict we need for getting sequences
            for center in cl_mean:

                center_clusters.append((chrm, center, strand))

    return center_clusters


def get_pA_surroundingSeq(region):
    """
    For the given region, return the sequences surrounding (-60, +20) all
    pA_sites found for that region in the available datasets from the
    UTR_SETTINGS file
    """

    # Just for the settings 
    here = os.path.dirname(os.path.realpath(__file__))
    (savedir, outputdir) = [os.path.join(here, d) for d in ('figures', 'output')]
    # Keep houskeeping information in a settings object
    settings = Settings(os.path.join(here, 'UTR_SETTINGS'), savedir, outputdir,
                            here, chr1=False)

    center_cls = get_pA_centers(region, settings)

    pyfas_input = {}
    for chrm, center, strand in center_cls:
        key = '|'.join([chrm, str(center), strand])

        if strand == '+':
            pyfas_input[key] = (chrm, center-60, center+20, strand)

        else:
            pyfas_input[key] = (chrm, center-20, center+60, strand)

    from annotation_parser import get_seqs
    import time

    t1 = time.time()
    surroundSeqs = get_seqs(pyfas_input, settings.hg19_fasta)
    t2 = time.time() - t1

    print('Time to get surrounding seqs for {0}: {1:2f}'.format(region, t2))

    return surroundSeqs

def get_pA_surroundingBed(region):
    """
    For the given region, return a bedfile for the the regions surrounding (-60,
    +20) all pA_sites found for that region in the available datasets from the
    UTR_SETTINGS file
    """

    # Just for the settings 
    here = os.path.dirname(os.path.realpath(__file__))
    (savedir, outputdir) = [os.path.join(here, d) for d in ('figures', 'output')]
    # Keep houskeeping information in a settings object
    settings = Settings(os.path.join(here, 'UTR_SETTINGS'), savedir, outputdir,
                            here, chr1=False)

    center_cls = get_pA_centers(region, settings)

    # where you will save a bedfile with all the reads that map close to any of your polyA sites
    snp_dir = os.path.join(here, 'SNP_analysis', 'non_polyA_PAS_reads')
    if not os.path.isdir(snp_dir):
        os.mkdir(snp_dir)

    out_bed = os.path.join(snp_dir, 'surrounding_pA_{0}.bed'.format(region))
    handle = open(out_bed, 'wb')

    for chrm, center, strand in center_cls:
        if strand == '+':
            handle.write('\t'.join([chrm, str(center-60), str(center+20), strand])+'\n')

        else:
            handle.write('\t'.join([chrm, str(center-20), str(center+60), strand])+'\n')

    handle.close()

    return out_bed

def main():
    region = '3UTR-exonic'
    bedfile = get_pA_surroundingBed(region)
    seqs = get_pA_surroundingSeq(region)

    debug()

if __name__ == '__main__':
    main()


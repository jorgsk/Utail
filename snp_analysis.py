"""
A separate file for the snp analysis
"""
from __future__ import division
import os
import ConfigParser

from subprocess import Popen, PIPE

import cPickle as pickle

import annotation_parser as genome
import utail
import results
from operator import itemgetter

# For sequence alignment
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment as MultSeqAlign
from Bio import AlignIO

from scipy import stats

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

def seq_getter(cls):
    """
    Get sequences for the different cell lines for this cls.seqs object. Return
    seqs as a set to remove duplicates.
    """
    cl_seqs = {}

    # add the seqs for the different cell lines to dict
    for metaseqs in cls.seqs:
        (seqs, cellL, comp) = metaseqs.split('|')

        # error source: a-reads in a t-group?

        for seq in seqs.split('$'):

            # reverse complement if of 't' type
            if cls.tail_type == 'T':
                seq = utail.reverseComplement(seq)

            # add sequences as a set; this removes any duplicates and thus
            # minimizes errors based on pcr amplification
            if cellL in cl_seqs:
                cl_seqs[cellL].add(seq)
            else:
                cl_seqs[cellL] = set([seq])

    return cl_seqs

def get_matchmismatches(columns, count_matrix, max_len):
    """
    Helper function to reduce code in the big one ...

    Count the number of matches and mismatches in the alignment of the RNA-seq
    data. This will later be used to infer the sequencing error rate.
    """
    mismatches = 0
    matches = 0

    # idea: #of mismatches = #seqs - max_count of any nt
    for col_nr in range(columns):
        local_max = max(count_matrix[col_nr].values())

        # don't count -s in the alignment; make a new local max
        if '-' in count_matrix[col_nr]:
            misses = count_matrix[col_nr]['-']

            # if more misses than any other, see if any other has
            # more than 1 seq; otherwise skip this
            if misses == local_max:
                local_max = max([v for l, v in
                                 count_matrix[col_nr].items() if l
                                 != '-'])
            if local_max < 2:
                continue
            newmax = max_len - misses

            mismatches += newmax - local_max
            matches += newmax
        else:
            mismatches += max_len - local_max
            matches += max_len

    return matches, mismatches

def get_alignment(cl, key, seqs, hg19seq, settings):
    """
    Write seqs to file and run clustalW on them and return a biopython align
    object
    """

    infile = os.path.join(settings.here, 'SNP_analysis',
                             'temp_files', key + cl+'error_rate.fasta')
    handle = open(infile, 'wb')

    # NOTE: only accept sequences longer than 40 nt. Shorter
    # sequences will give you alignment problems.
    written = 0
    for seqnr, seq in enumerate(seqs):
        if len(seq) > 40:
            handle.write('>someseq{0}\n{1}\n'.format(seqnr, seq))
            written += 1

    handle.write('>hg19\n{0}\n'.format(hg19seq))
    handle.close()

    # skip if 0 or 1 seq written
    if written < 3:
        return False

    outfile = os.path.join(settings.here, 'SNP_analysis',
                             'temp_files', key + cl+'.aln')

    # call clustalw and get consensus with score
    cmd = 'clustalw -infile={0} -outfile={1} -quiet'.format(infile,
                                                           outfile)
    # run clustalw and send stdout (messages) to devnull
    Popen(cmd, shell=True, stdout=open(os.devnull, 'w')).wait()

    alignment = AlignIO.read(outfile, "clustal")

    return alignment

def get_sequencing_error_rate(settings, region, cell_lines, super_3utr, hg19Seqs):
    """
    For each cell line:
        for each poly(A) site with > 5 reads:
            align the reads and get the dumb consensus at 99%; for the smallest
            common sequence-stretch (CSS), count the number of X'es in the
            consensus and divide by the length of the CSS and divide by the
            number of sequences. This will be the site-specific error rate.

    Average the error rate distribution; hopefully it's gaussian; if not there
    are some outliers that should be taken care of.

    For saving the snp-info, when you find AAAAAA (6) in RNA-seq and a non-A in
    hg19, test for stats.binom_test(6,6, P(sequencing_error)). So you must save
    the 6-A and 6max info. The P(seq_err) will come later.

    What to do with AAATTT (g) ? Then both seem likely
    What about AAAATT (G) ? Check if either or just one is different from (G)
    and give p-values to both? This could be the result of sampling the same
    piece twice, which is not good behind the assumption of independent draws
    (the draws are the sequences themselves)

    Beware when the hg19 has an '-'. This means an insertion into the rnaseq
    dataset (or very unlikely a deletion missed in hg19).
    """

    # for the error rates
    matchCounter = dict((cl, {'matches_non_strict':0,
                             'mismatches': 0}) for cl in cell_lines)

    # for future p-value calculations
    snp_site_counter = dict((cl, {}) for cl in cell_lines)

    # for PAS change info
    # you must store the position information relative to hg19 too
    # TODO
    pas_changes = {}

    for utr_name, utr in super_3utr[region].iteritems():
        for cls in utr.super_clusters:

            # skip those with < 5 supported reads
            if cls.nr_support_reads < 5:
                continue

            # an ID for this polyA coordinate
            key = '_'.join([utr.chrm, str(cls.polyA_coordinate), cls.strand])

            # the sequence for hg19
            hg19seq = hg19Seqs[key]

            # the seqs for the different cell lines at this poly(A) site
            for cl, seqs in seq_getter(cls).items():

                # get a biopython alignment object for these sequences + hg19
                alignment = get_alignment(cl, key, seqs, hg19seq, settings)

                # if the alignment cound not be made (too few secs), abort
                if not alignment:
                    continue

                # Get the start and end coordinates of the alignment according
                # to the hg19 sequence
                start, stop, hg19row = get_startstop(alignment)

                ## crop the alignment according to the hg19 sequence
                crop_alignment = alignment[:, start:stop]

                # Get the aligned hg19 sequence (may contain gaps, which the
                # hg19seq doesn't)
                aln_hg19seq = crop_alignment[hg19row].seq.tostring()

                # remove the hg19 sequence from the alignment
                rnaseq_align = MultSeqAlign([a for a in crop_alignment if a.id != 'hg19'])

                # get the count_matrix for the rna-seq reads only
                count_matrix = AlignInfo.SummaryInfo(rnaseq_align).pos_specific_score_matrix()

                max_len = len(rnaseq_align[:,0])

                nr_columns = stop - start

                ##### 1) Get the # of mismatches/matches_non_strict
                non_str_matches, mismatches = get_matchmismatches(nr_columns,
                                                                  count_matrix,
                                                                 max_len)
                matchCounter[cl]['matches_non_strict'] += non_str_matches
                matchCounter[cl]['mismatches'] += mismatches

                ##### 2) Get the stats for the places with hg19 mismatch
                for col_nr in range(nr_columns):
                    col = count_matrix[col_nr]
                    # don't count -s in the alignment; make a new local max
                    max_nt, max_count = sorted(col.items(), key=itemgetter(1))[-1]

                    # if by some obscure reason - is the max, take the second
                    # highest
                    if max_nt == '-':
                        max_nt, max_count = sorted(col.items(), key=itemgetter(1))[-2]

                    if max_nt != aln_hg19seq[col_nr]:
                        col_key = key + '_col_nr'

                        if aln_hg19seq[col_nr] == '-':
                            snp_site_counter[cl][col_key] = (col[max_nt],
                                                             max_len,
                                                             'insertion')
                        else:
                            snp_site_counter[cl][col_key] = (col[max_nt],
                                                             max_len,
                                                             'substitution')
    # convert matches/mismatches into error rates
    error_rates = make_error_rates(matchCounter)

    return error_rates, snp_site_counter, pas_changes

def make_error_rates(matchCounter):
    """
    Convert matches/mismataches into an error rates
    """
    error_rates = {}
    for cl, matchdict in matchCounter.items():
        nsmatch = matchdict['matches_non_strict']
        mismatch = matchdict['mismatches']

        if nsmatch == 0:
            continue

        error_rates[cl] = {'strict': mismatch/(nsmatch - mismatch),
                           'non_strict': mismatch/(nsmatch)}

    return error_rates


def get_startstop(alignment):
    """
    Given a sequence aligment, find where the hg19 element start and stops
    """

    for rownr, seqr in enumerate(list(alignment)):
        if seqr.description == 'hg19':
            hg19gapseq = seqr.seq.tostring()

            # get the start
            start = 0
            for pos in hg19gapseq:
                if pos == '-':
                    start += 1
                else:
                    break

            # get the stop
            stop = len(hg19gapseq)
            for pos in reversed(hg19gapseq):
                if pos == '-':
                    stop -= 1
                else:
                    break

            return start, stop, rownr

def get_all_pvalues(error_rates, hg19Seqs, settings, region, cell_lines, super_3utr):
    """
    Re-align all sequences and use the sequencing error rates to give p-values
    to each nucleotide that is different in the rna-seq set compared to hg19

    When you find AAAAAA (6) in RNA-seq and a non-A in hg19, test for
    stats.binom_test(6,6, P(sequencing_error) = 0.0036)

    What to do with AAATTT (g) ? Then both seem likely
    What about AAAATT (G) ? Check if either or just one is different from (G)
    and give p-values to both? This could be the result of sampling the same
    piece twice, which is not good behind the assumption of independent draws
    (the draws are the sequences themselves)
    """

    # dict of list of tuples (pval, cell_line+polyAID)
    pvals = dict((cl, []) for cl in cell_lines)

    for utr_name, utr in super_3utr[region].iteritems():

        for cls in utr.super_clusters:

            # include only those with > 5 supported reads
            if cls.nr_support_reads < 5:
                continue

            key = '_'.join([utr.chrm, str(cls.polyA_coordinate), cls.strand])

            hg19seq = hg19Seqs[key]

            # for saving sequences for writing to file
            cl_seqs = {}

            # add the seqs for the different cell lines to dict
            for metaseqs in cls.seqs:
                (seqs, cellL, comp) = metaseqs.split('|')

                for seq in seqs.split('$'):

                    # reverse complement if of 't' type
                    if cls.tail_type == 'T':
                        seq = utail.reverseComplement(seq)

                    # add sequences as a set; this removes any duplicates and thus
                    # minimizes errors based on pcr amplification
                    if cellL in cl_seqs:
                        cl_seqs[cellL].add(seq)
                    else:
                        cl_seqs[cellL] = set([seq])

            # write the seqs to a fasta file and align
            for cl, seqs in cl_seqs.items():

                # get the probability of sequence error for this cell line

                probSeqErr = error_rates[cl]['strict']


                infile = os.path.join(settings.here, 'SNP_analysis',
                                         'temp_files', key + cl+'pvalue.fasta')
                handle = open(infile, 'wb')

                # NOTE: only accept sequences longer than 40 nt. Shorter
                # sequences will give you alignment problems.
                written = 0
                for seqnr, seq in enumerate(seqs):
                    if len(seq) > 40:
                        handle.write('>someseq{0}\n{1}\n'.format(seqnr, seq))
                        written += 1

                # finally write hg19
                handle.write('>hg19\n{0}\n'.format(hg19seq))
                handle.close()

                # skip if less than 4 written; we will have poor p-values; no
                # point in even checking
                if written < 3:
                    continue

                outfile = os.path.join(settings.here, 'SNP_analysis',
                                         'temp_files', key + cl+'.aln')

                # call clustalw and get consensus with score
                cmd = 'clustalw -infile={0} -outfile={1} -quiet'.format(infile,
                                                                       outfile)
                # run clustalw and send stdout (messages) to devnull
                Popen(cmd, shell=True, stdout=open(os.devnull, 'w')).wait()

                alignment = AlignIO.read(outfile, "clustal")

                # Get the start and end coordinates of the alignment according
                # to the hg19 sequence; restrict yourself to the alignment
                # matrix within the 40 hg19.

                start, stop, hg19row = get_startstop(alignment)

                ## crop the alignment according to the hg19 sequence
                crop_alignment = alignment[:, start:stop]

                # Get the aligned hg19 sequence (may contain gaps, which the
                # hg19seq doesn't)

                aln_hg19seq = crop_alignment[hg19row].seq.tostring()

                # remove the hg19 sequence from the alignment
                rnaseq_align = MultSeqAlign([a for a in crop_alignment if a.id != 'hg19'])

                # get the count_matrix for the rna-seq reads only
                count_matrix = AlignInfo.SummaryInfo(rnaseq_align).pos_specific_score_matrix()

                # You now have your alingment in the coordinates of the hg19
                # sequence; compare each column; if 

                max_len = len(rnaseq_align[:,0])

                nr_columns = stop - start

                # for each column, check if the max-count is different from
                # hg19; if it is, compare the p-value for its difference

                for col_nr in range(nr_columns):
                    col = count_matrix[col_nr]
                    # don't count -s in the alignment; make a new local max
                    max_nt, max_count = sorted(col.items(), key=itemgetter(1))[-1]

                    # if by some obscure reason this is the max, take the second
                    # highest
                    if max_nt == '-':
                        max_nt, max_count = sorted(col.items(), key=itemgetter(1))[-2]

                    #if '-' in col:
                        #misses = col['-']

                    # if you have a mismatch with the human genome
                    # hi-ho! This becomes a problem when there's a gap in the
                    # alignement -- you cannot anymore compare to the human
                    # genome, you must compare with the aligned version which
                    # contains the gap. But how to deal with the gap???

                    if max_nt != aln_hg19seq[col_nr]:
                        pval = stats.binom_test(col[max_nt], max_len,
                                                probSeqErr)

                        #print pval
                        #print ''
                        #print rnaseq_align[:, col_nr-3:col_nr+4]
                        #print ''
                        #print aln_hg19seq[col_nr-3:col_nr+4]

                        if aln_hg19seq[col_nr] == '-':
                            pvals[cl].append((pval, key, 'insertion'))
                        else:
                            pvals[cl].append((pval, key, 'substitution'))

    return pvals


def get_hg19_seqs(settings, super_3utr, region):

    hg19input = {}
    for utr_name, utr in super_3utr[region].iteritems():
        for cls in utr.super_clusters:

            if cls.strand == '+':
                end = cls.polyA_coordinate
                beg = end - 50
            elif cls.strand == '-':
                beg = cls.polyA_coordinate
                end = beg + 50

            key = '_'.join([utr.chrm, str(cls.polyA_coordinate), cls.strand])

            hg19input[key] = (utr.chrm, beg, end, cls.strand)

    hg19Seqs = genome.get_seqs(hg19input, settings.hg19_fasta)

    return hg19Seqs

def correct_pvalues(p_values, cell_lines, settings, super_3utr, region):
    """
    Use the Benjamini-Hochman method to correct for multiple testing

    For each cell line, sort the p-values in increasing order (low->high).
    Compare these p-values to the modified alpha: (k/n)*alpha where k is the
    rank of the pvalue and n is the total number of tests (total number of
    nucleotide positions tested. n is equal to the number of tries you did.

    Note, there is a problem here. You only do tests when the max_nt is
    different from the hg19 -> this leaves you with very few tests for the n.
    However, if you had tested all of them, n would be a lot lot bigger.
    """

    alpha = 0.05

    cor_pval = dict((cl, []) for cl in cell_lines)

    for cl, pvals in p_values.items():
        n = len(pvals)
        for rank_nr, (pval, key, mutation_type) in enumerate(sorted(pvals)):
            q = ((rank_nr+1)/n)*alpha
            if pval < q:
                cor_pval[cl].append((pval, key, mutation_type, rank_nr+1, 'significant'))

            print mutation_type
            print 'pval', pval
            print 'qval', q
            if pval < q:
                print 'significant!'
            else:
                print 'insignificant ...'

            print '-----------------------------------------------'

    return cor_pval

def get_snps_corrected(snp_stats, error_rates, cell_lines):
    """
    Calculate pvalues for all the snp sites; then order then and re-calculate
    their significance using the Benjamini-Hochman FDR
    """
    alpha = 0.05

    cor_snp = dict((cl, []) for cl in cell_lines)

    for cl, sstats in snp_stats.items():
        # skip empty ones
        if sstats == {}:
            continue
        n = len(sstats)
        probErr = error_rates[cl]['strict']

        # get the pvals
        pvals = []
        for rank_nr, (key, (max_count, N, mutation_type))\
                in enumerate(sorted(sstats.items())):

            pval = stats.binom_test(max_count, N, probErr)
            pvals.append((pval, key, mutation_type))

        # sort pvals and calculate q-s
        for rank_nr, (pval, key, mutation_type) in enumerate(sorted(pvals)):

            q = ((rank_nr+1)/n)*alpha

            if pval < q:
                cor_snp[cl].append((pval, key, mutation_type, rank_nr+1, 'significant'))

            print mutation_type
            print 'pval', pval
            print 'qval', q
            if pval < q:
                print 'significant!'
            else:
                print 'insignificant ...'

    return cor_snp

def snp_analysis(settings):
    """
    For all the poly(A) sites, align the seqs for each dataset separately and
    find either PAS or PAS_with_snips.
    """
    # work on the following cell_lines
    cell_lines = ('K562', 'GM12878', 'HUVEC', 'HeLa-S3', 'HEPG2','H1HESC',
                  'NHEK', 'NHLF', 'HSMM', 'MCF7', 'AG04450')


    # Get all PAS and all the snips for thos PAS
    #pas_snips, allpas = pas_and_pas_snips()
    #all_ds = [ds for ds in settings.datasets if ((('Cytoplasm' in ds) or
                                                #('Whole_Cell' in ds) or
                                                #('Nucleus' in ds)) and
                                               #(not 'Minus' in ds))]
    all_ds = [ds for ds in settings.datasets if (('Cytoplasm' in ds) and
                                               (not 'Minus' in ds))]
    #speedrun = True
    speedrun = False

    outdir = os.path.join(settings.here, 'SNP_analysis')

    region = '3UTR-exonic'

    subset = all_ds

    batch_key = 'SNP'

    # Get the clustered poly(A) sites
    dsets, super_3utr = results.super_falselength(settings, region, batch_key,
                                          subset, speedrun)

    # Get hg19 sequences 50bp downstream all polyA sites
    hg19Seqs = get_hg19_seqs(settings, super_3utr, region)

    # 1 round of clustalw, 3 important pieces of information
    error_rates, snp_stats, PAS_mutants = get_sequencing_error_rate(settings, region,
                                                                    cell_lines,
                                                                    super_3utr,
                                                                   hg19Seqs)

    # use the error rates to calculate P-values for the mutation sites using
    # snp_stats; simultaneously correct for multiple testing using
    # Benjamini-Hochberg's method
    corr_snps = get_snps_corrected(snp_stats, error_rates, cell_lines)

    # check the PAS-changes for snps; if they exist, way to go.
    # TODO

    # Get the pvalues for every nucleotide with poly(A) reads
    # At the same time, to save speed, get the PAS_mutant list
    p_values, pas_mutants = get_all_pvalues(error_rates, hg19Seqs,settings,
                                      region, cell_lines, super_3utr)

    # Get the corrected p-values
    corr_pval = correct_pvalues(p_values, cell_lines, settings,
                                super_3utr, region)

    # go over the PAS_mutant list with the corrected p values

    # Check for PAS in both rna-seq consensus and hg19; when found, compare with
    # the other for mismatch; and if it exists, report a pval for the
    # mismatch and what type of mismatch it is (there are 6 different)

    #pas_mutatns = get_PAS_mutatns(corr_pval)

    debug()

    for utr_name, utr in super_3utr[region].iteritems():

        for cls in utr.super_clusters:

            key = '_'.join([utr.chrm, str(cls.polyA_coordinate), cls.strand])

            # Separate the rna-seqs from the different cell_lines

            # make a temporary file with the 
            temp_file = os.path.join(settings.here, 'SNP_analysis',
                                     'temp_files', key)
            # for each dataset, align the seqs
            # TODO you don't really know how to solve this problem
            # you are thinking of aligning all the sequences, including the hg19
            # sequence when there is a PAS, so you know where the PAS fits.

            cl_seqs = {}

            # add the seqs for the different cell lines to dict
            for metaseq in cls.seqs:
                (seq, cellL, comp) = metaseq.split('|')

                # reverse complement if of 't' type
                if cls.tail_type == 'T':
                    seq = utail.reverseComplement(seq)

                if cellL in cl_seqs:
                    cl_seqs[cellL].append(seq)
                else:
                    cl_seqs[cellL] = [seq]

            # for each cell line, write the seqs to a fasta file and align; give
            # a score to each aligned base and find the PAS or snip-PAS; give a
            # score to the PAS/Snip-PAS; align with hg19 sequence to get a
            # coordinate reference to compare with the other sequences. A 'hg19'
            # coordinate.

            # IDEA: if a deletion/insertion is detected, stop for a second and
            # take a look

            for cellLn, seqs in cl_seqs.items():
                infile = os.path.join(settings.here, 'SNP_analysis',
                                         'temp_files', key + cellLn+'.fasta')
                handle = open(infile, 'wb')

                # NOTE: only accept sequences longer than 40 nt. Shorter
                # sequences will give you alignment problems.
                for seqnr, seq in enumerate(seqs):
                    if len(seq) > 40:
                        handle.write('>someseq{0}\n{1}\n'.format(seqnr, seq))

                handle.close()

                outfile = os.path.join(settings.here, 'SNP_analysis',
                                         'temp_files', key + cellLn+'.aln')

                # call clustalw and get consensus with score
                cmd = 'clustalw -infile={0} -outfile={1} -quiet'.format(infile,
                                                                       outfile)
                # run clustalw and send stdout (messages) to devnull
                Popen(cmd, shell=True, stdout=open(os.devnull, 'w')).wait()

                # We accept 1 in 7 wrong
                summary_align = AlignInfo.SummaryInfo(AlignIO.read(outfile, "clustal"))
                consensus = summary_align.dumb_consensus(threshold=0.7).tostring()

                # get PAS/PAS-snips for both rna-seq and hg19. Index them from
                # the hg19 (40 nt downstream).

                # compare the PAS/PAS-snips coordinates between rnaseq and hg19;
                # if there is a mismatch, save it to the cls object.
                # Rules
                # RNA-seq PAS -> hg19 PAS-snip == novel PAS
                # RNA-seq PAS-snip -> hg19 PAS == loss of PAS

                # call clustalw on consensus and hg19 and get hg19-coordinates
                # of the PAS/PAS-snips

                debug()

            # remember: hg19 output is always in the 5-3 direction; while if
            # type == T you must reverse transcribe the poly(A) reads.
            #               TAAATGCCAACATTTAAAATAATTTGATAAAGTCTTTGTAGCCACTCAG
            #ATACACATGAATACTTAAATGCCAACATTTAAAATAATTTGATAAAGTCTTTTTAGCCACTCAGCC

            if cls.strand == utr.strand:
                keyw = 'same'
            else:
                keyw = 'opposite'

            debug()

def main():
    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    (savedir, outputdir) = [os.path.join(here, d) for d in ('figures', 'output')]

    # Keep houskeeping information in a settings object
    settings = results.Settings(os.path.join(here, 'UTR_SETTINGS'), savedir, outputdir,
                        here, chr1=False)

    snp_analysis(settings)

if __name__ == '__main__':
    main()

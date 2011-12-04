"""
A separate file for the snp analysis
"""
from __future__ import division
import os
import re

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

class Mutant(object):
    """
    A simple class to describe a mutant PAS site for hg19 or rna-seq sequences
    """

    def __init__(global_coord, my_PAS, other_PAS, hg19_coord,
                 mutant_type, identifier):
        self.global_coord = global_coord
        self.my_PAS = my_PAS
        self.other_PAS = other_PAS
        self.hg19_coord = hg19_coord
        self.mytant_type = mutant_type
        self.identifier = identifier

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

    # skip if less than 5 written
    if written < 5:
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

def align_seqs(settings, region, cell_lines, super_3utr, hg19Seqs):
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

    # list of PAS hexamers -- needed for downstream code
    high = ['AATAAA', 'ATTAAA']
    low = ['TATAAA', 'AGTAAA', 'AAGAAA', 'AATATA', 'AATACA',
           'CATAAA', 'GATAAA', 'AATGAA', 'TTTAAA', 'ACTAAA', 'AATAGA']

    # compile regular expressions of PAS patterns
    pas_patterns = [re.compile(pas) for pas in high+low]

    # for the error rates
    matchCounter = dict((cl, {'matches_non_strict':0,
                             'mismatches': 0}) for cl in cell_lines)

    # for future p-value calculations
    snp_site_counter = dict((cl, {}) for cl in cell_lines)

    # for PAS change info
    pas_changes = dict((cl, {}) for cl in cell_lines)

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

                ##### 2) Get the stats for binomial testing for the places with hg19 mismatch
                binom_stats = get_binom_stats(key, nr_columns, count_matrix,
                                              aln_hg19seq, max_len)

                for col_key, bin_stats in binom_stats:
                    snp_site_counter[cl][col_key] = bin_stats

                ##### 3) Look for PAS-mutants to and from hg19
                pasmut = get_pasmutants(rnaseq_align, aln_hg19seq, pas_patterns,
                                        high, low)
                # only add if PAS mutants were found
                if pasmut:
                    pas_changes[cl][key] = pasmut

    # convert matches/mismatches into error rates
    error_rates = make_error_rates(matchCounter)

    return error_rates, snp_site_counter, pas_changes

def get_binom_stats(key, nr_columns, count_matrix, aln_hg19seq, max_len):
    """
    Simply collect the n and r you need for doing the binomial test:
        binom_test(r, n, P[err])

    index by the column key + cls key + column nr from hg19 seq
    """
    outp = []

    for col_nr in range(nr_columns):
        col = count_matrix[col_nr]
        # don't count -s in the alignment; make a new local max
        max_nt, max_count = sorted(col.items(), key=itemgetter(1))[-1]

        # if by some obscure reason - is the max, take the second
        # highest
        if max_nt == '-':
            max_nt, max_count = sorted(col.items(), key=itemgetter(1))[-2]

        if max_nt != aln_hg19seq[col_nr]:
            col_key = key + '_{0}'.format(col_nr)

            if aln_hg19seq[col_nr] == '-':
                outp.append((col_key, (col[max_nt], max_len, 'insertion')))
            else:
                outp.append((col_key, (col[max_nt], max_len, 'substitution')))

    return outp

def compare_PAS(seqsource, this_pas, other_hex, pos, high, low, other_seq):
    """
    Compare the PAS/hexamer. Comparing hg19 -> rna-seq or rna-seq -> hg19
    matters in only 1 case: if a '-' is found, it's always either an insertion
    or deletion in rna-seq. We assume that hg19 is the gold standard.
    """
    # check for high/low for your PAS
    if this_pas in high:
        this_pasType = 'high'
    elif this_pas in low:
        this_pasType = 'low'

    if '-' in other_hex:
        # 
        # if there is a '-' in hg19, it's actually an insertion in rna-seq
        if seqsource == 'rnaseq':
            other_mutType = 'neutral'
            this_mutType = 'insertion'

        # if there is a '-' in rna-seq, it's actually a deletion rna-seq
        if seqsource == 'hg19':
            other_mutType = 'deletion'
            this_mutType = 'neutral'

        # try to recover the PAS in hg19
        var1 = ''.join(other_seq[pos-1:pos+6].split('-'))
        var2 = ''.join(other_seq[pos:pos+7].split('-'))

        # if the recovered PAS is the same as the rna-seq pas, this
        # means that the insertion didn't actually change the PAS
        if var1 == this_pas or var2 == this_pas:
            return False, False
        else:
            # if the hg19 PAS is high or low, label it like that. it
            # means that the insertion changed the rna-seq PAS from
            # one type to another.
            if var1 in high or var2 in high:
                other_pasType = 'high'
            elif var1 in low or var2 in low:
                other_pasType = 'low'
            else:
                other_pasType = 'NaP' #Not a PAS ...

    # if there is no '-' in hg19, the comparison is easy
    else:
        other_mutType = 'neutral'
        this_mutType = 'neutral'

        if other_hex in high:
            other_pasType = 'high'
        elif other_hex in low:
            other_pasType = 'low'
        else:
            other_pasType = 'NaP'

    mut_switch = '_'.join([this_mutType, 'to', other_mutType])
    pas_switch = '_'.join([this_pasType, 'to', other_pasType])

    return mut_switch, pas_switch

def get_pasmutants(rnaseq_align, aln_hg19seq, pas_patterns, high, low):
    """
    Return a dict [hg19/rnaseq] = (this_motif, coord, other_motif, type of
    change). For example [hg19] = (AATAAA, (3,9), AAGAAA, 'high_to_low')

    The changes are of these types:

        high_to_low (canonical to noncanonical)
        low_to_high
        high_to_broken (canonical to not funcational PAS: AATGTA f.ex)
        low_to_broken
        broken_to_high
        broken_to_low

        + you must look for insertions/deletions

    if you ever want to count the # of cases you have, this can be a good place

    """
    muts = {'hg19': [], 'rnaseq': []}
    # Get the consensus sequence at 51% (you're using 'max' (>50) to determine snp)
    # NOTE changed from dumb to gap-consensus
    # NOTE you will reduce the errors if you increase treshold 
    consensus = AlignInfo.SummaryInfo(rnaseq_align).\
            gap_consensus(threshold=0.91).tostring()

    for seqsource in ['rnaseq', 'hg19']:
        if seqsource == 'rnaseq':
            this_seq = consensus
            other_seq = aln_hg19seq
        elif seqsource == 'hg19':
            this_seq = aln_hg19seq
            other_seq = consensus

        # Get all the PAS in consensus
        pastup = get_pas_and_distance(pas_patterns, this_seq)

        # skip when no PAS was found
        if pastup == ('NA', 'NA'):
            continue

        for this_pas, pos in zip(*pastup):

            # Don't process 'na' pases!
            if this_pas == 'NA':
                continue

            other_hex = other_seq[pos:pos+6]

            # don't include PAS too close to the 3' end
            if pos > 36:
                continue

            # do some checks on other_hex
            if 'X' in other_hex:
                continue

            # ignore if 3 or more gaps
            # TODO, try to reconciliate 2 gaps
            if other_hex.count('-') > 2:
                continue

            # don't do anything if the PAS is the same...
            if this_pas == other_hex:
                continue

            # ... but if it's different, you struck gold!
            # test both cases; hg19 is this_seq, or hg19 is other_seq
            # you need to test if hexamer is high, low, insertion, or deletion.
            # the last two depends on where the '-' is. If it's in hg19, then
            # there as an insertion in rnas-seq. if it's in rnaseq, there was a
            # deletion in rnaseq.
            mut_switch, pas_switch = compare_PAS(seqsource, this_pas,
                                                  other_hex, pos, high, low,
                                                  other_seq)
            print seqsource
            print seqsource +': ' + this_pas

            try:
                print other_seq[pos-1]+'|'+other_hex+'|'+other_seq[pos+7]
            except IndexError:
                print '|'+other_hex+'|'

            print mut_switch
            print pas_switch
            debug()
            #if mut_switch:
                #mutobj = Mutation()

            # oops: in 1 dataset you seem to have 2 candidates. Does this mean
            # you will have only 20 in total? Is this a lot or very little? Who
            # knows? I was hoping to see new polyA sites in known 3UTRs, maybe
            # to shorten them? solution: add more datasets?

            # XXX you don't capture when there has been 2 indels. either you
            # must ignore it or you must treat it.

            # XXX you occasionally have totally different rna-seq reads aligned
            # for a polyA site. A subsection aligns well with hg19, while the
            # others could be wrongly reverse-transcribed or result from another
            # error. If the 'wrong' sites domeniate, you'll have a PAS -> NaP
            # mutation which is wrong. These cases are rare, so maybe you want
            # to ignore them? However, it's not clear how you can detect them.
            # Maybe with some of the bio python tools.

            # occasionally the 'wrong' sequence dominates 20 times over the
            # 'right' as compared with hg19.

            # XXX not a problem, but whenever you have a high->low, high->high,
            # or low->low, you will have the same site repeated for both rna-seq
            # and hg19s point of view

            # TODO make a mutation object from the mutswitchpasswitch
            #mut_obj = 

            #muts[seqsource].append(this_pas, (pos,pos+6), other_hex, switch)

    # Only return dict if you've added something
    if muts != {'hg19': [], 'rnaseq': []}:
        return muts
    else:
        return False


def get_pas_and_distance(pas_patterns, sequence):
    """
    Go through the -40 from the polya read average. Collect PAS and distance
    as you find them. Must supply poly(A) sites relative to UTR-beg. The
    sequence must also be relative to UTR-beg (3->5 direction)
    """

    pases = []
    temp_pas = []

    for pas_exp in pas_patterns:
        temp_pas.append([(m.group(), m.start()) for m in
                        pas_exp.finditer(sequence)])

    if sum(temp_pas, []) != []:
        has_pas = sum(temp_pas, [])
        if len(has_pas) == 1:
            (has_pas, pas_dist) = ((has_pas[0][0],), (has_pas[0][1],))
        else:
            (has_pas, pas_dist) = zip(*has_pas)

        #(('GATAAA', 'AATGAA'), (22, 28))
        # This is what you return.
        pases.append((has_pas, pas_dist))

    # in case no pas are found, return 'na'
    else:
        # If both were 'NA', return it
        pases.append(('NA', 'NA'))

    # Join the pases together. if the length is 1, it means that both are
    # 'NA', 'NA'. Just return one of them.
    if len(set(pases)) == 1:
        return(pases[0])

    # If not, join the two together. A zip perhaps? If one of them is 'NA',
    # remove it, and return the other
    else:
        try:
            pases.remove(('NA', 'NA'))
            return pases[0]
        except ValueError:
            zipases = zip(*pases)
            return (sum(zipases[0], ()), sum(zipases[1], ()))


def make_error_rates(matchCounter):
    """
    Convert matches/mismataches into error rates
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
    #speedrun = True
    dsets, super_3utr = results.super_falselength(settings, region, batch_key,
                                          subset, speedrun)

    # Get hg19 sequences 50bp downstream all polyA sites
    hg19Seqs = get_hg19_seqs(settings, super_3utr, region)

    # 1 round of clustalw, 3 important pieces of information
    error_rates, snp_stats, PAS_mutants = align_seqs(settings, region,
                                                     cell_lines, super_3utr,
                                                     hg19Seqs)

    # use the error rates to calculate P-values for the mutation sites using
    # snp_stats; simultaneously correct for multiple testing using
    # Benjamini-Hochberg's method
    corr_snps = get_snps_corrected(snp_stats, error_rates, cell_lines)

    # Go though the PAS_mutants and quantify if the change you detected is
    # significant. make some basic statistics on the rate of
    # changes compared to expected and maybe compared to sequence changes in
    # non-PAS regions.

    debug()


def main():
    # The path to the directory the script is located in
    here = os.path.dirname(os.path.realpath(__file__))

    (savedir, outputdir) = [os.path.join(here, d) for d in ('figures', 'output')]

    # Keep houskeeping information in a settings object
    settings = results.Settings(os.path.join(here, 'UTR_SETTINGS'), savedir,
                                outputdir, here, chr1=False)

    snp_analysis(settings)

if __name__ == '__main__':
    main()

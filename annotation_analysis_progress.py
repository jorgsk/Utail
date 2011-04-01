"""
Create bed files from annotation file. Specifically, create bed-files for all
annotated introns and exons in 5UTR, CDS, and 3UTR, as well as all annotated TTS
sites.

From the original annotation file, I must extract the following
aTTS
3UTR_exons
3UTR_terminal_exon
3UTR_introns
CDS_exons
CDS_introns
5UTR_exons
5UTR_introns

As well, the program has a function to get only 3UTR regions.

"""
def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

# only get the debug functio prop
if run_from_ipython():
    from IPython.Debugger import Tracer
    debug = Tracer()
else:
    def debug(): pass

from operator import itemgetter
from os.path import join as path_join
import os
here = os.path.dirname(os.path.realpath(__file__))
import time
import sys
from subprocess import Popen, PIPE, call
sys.path.append(os.path.join(here,'py_modules'))
sys.path.append(os.path.join(here,'py_modules/pyfasta'))
from fasta import Fasta
import math


class Region(object):

    def __init__(self, name):
        self.name = name

        # Note exons and introns should be on the form
        # (chr1, beg, end, strand), and sorted.
        self.exons = []
        self.length = 0

class Transcript(object):

    def __init__(self, t_id, g_id, t_type, chrm, beg, end, strand):
        self.ts_id = t_id
        self.gene_id = g_id
        self.t_type = t_type
        self.chrm = chrm
        self.beg = beg
        self.end = end
        self.strand = strand

        self.three_utr = Region('3UTR')
        self.cds = Region('CDS')
        self.five_utr = Region('5UTR')
        self.exons = []

    def add_cds(self, cds):
        self.cds.exons.append(cds)

    def get_cds_introns(self):
        self.cds.exons.sort(key = itemgetter(1)) # sort by exon start coord
        chrm = (self.chrm,)*len(self.cds.exons)
        strnd = (self.strand,)*len(self.cds.exons)
        useful = sum([[str(ex[1]-1), str(ex[2]+1)] for ex in self.cds.exons],[])

        # Zip together the beginnings and the ends of the introns
        return zip(chrm, useful[1::2], useful[2::2], strnd)

    def add_exon(self, exon):
        self.exons.append(exon)

    def get_aTTS(self, with_utr=True):
        """ Return the aTTS of a transcript. If processed, also return for those
        who have not a 3UTR"""

        if with_utr:
            if self.three_utr.exons != []:
                if self.strand == '+':
                    ex = self.three_utr.exons[-1]
                    return (ex[0], ex[2], str(ex[2]+1), ex[3])
                elif self.strand == '-':
                    ex = self.three_utr.exons[0]
                    return (ex[0], ex[1], str(ex[1]+1), ex[3])

        if not with_utr:
            if self.exons != []:
                self.exons.sort(key = itemgetter(1))
                if self.strand == '+':
                    ex = self.exons[-1]
                    return (ex[0], ex[2], str(ex[2]+1), ex[3])
                elif self.strand == '-':
                    ex = self.exons[0]
                    return (ex[0], ex[1], str(ex[1]+1), ex[3])
        return 0


    def add_utr(self, utr):
        # Determine if 5' or 3' by comparing with cds_beg and strand add utr to
        # the corresponding list.
        # if no CDS was annotated, don't annotate the UTR either because you
        # don't know if it's 5' or 3'
        if self.cds.exons == []:
            return
        self.cds.exons.sort(key=itemgetter(1))
        cds_beg = self.cds.exons[0][1]
        utr_beg = utr[1]

        if self.strand == '+':
            if utr_beg <= cds_beg:
                self.five_utr.exons.append(utr)
            if utr_beg >= cds_beg:
                self.three_utr.exons.append(utr)
                self.three_utr.length = self.three_utr.length + utr[2]-utr[1]

        if self.strand == '-':
            if utr_beg <= cds_beg:
                self.three_utr.exons.append(utr)
                self.three_utr.length = self.three_utr.length + utr[2]-utr[1]
            if utr_beg >= cds_beg:
                self.five_utr.exons.append(utr)

    def get_utr_introns(self, side=3):

        # You've a problem. Once or twice you might have an intron of length 1,
        # and this would give you a bed-entry of length -1, causing a crash.

        if side == 5:
            # 5 utr
            chrm = (self.chrm,)*len(self.five_utr.exons)
            strnd = (self.strand,)*len(self.five_utr.exons)
            self.five_utr.exons.sort(key = itemgetter(1))
            useful = sum([[str(ex[1]-1), str(ex[2]+1)] for ex in
                          self.five_utr.exons], [])

            return zip(chrm, useful[1::2], useful[2::2], strnd)

        if side == 3:
            # 3 utr
            chrm = (self.chrm,)*len(self.three_utr.exons)
            strnd = (self.strand,)*len(self.three_utr.exons)
            self.three_utr.exons.sort(key = itemgetter(1))
            useful = sum([[str(ex[1]-1), str(ex[2]+1)] for ex in
                          self.three_utr.exons], [])

            return zip(chrm, useful[1::2], useful[2::2], strnd)


def make_transcripts(annotation):
    """Loop through a (GENCODE) annotation file and get the exon fields.
    Finally compose them into sets of transcripts. As well return a dictionary
    of genes with the list of transcripts belonging to that gene."""

    # Skip first 4 lines of annotation (for gencode at least.......)
    a_handle = skip_lines(annotation, 7)

    # Initialize dictionaries
    utr, cds, exons, transcripts, genes = (dict(), dict(), dict(), dict(),
                                           dict())

    # Go though annotation and classify transcripts, genes, and utr + cds exons
    for line in a_handle:
        (chrm, source, feature, beg, end, d, strand, d, d, g_id, d, t_id, d, d,
        d, d, d, d, d, t_type) = line.split()[:20]
        beg = int(beg)
        end = int(end)

        t_type = t_type.rstrip(';').strip('"') #   remove some characters
        t_id = t_id.rstrip(';').strip('"')
        g_id = g_id.rstrip(';').strip('"')

        if feature == 'transcript':
            transcripts[t_id] = Transcript(t_id, g_id, t_type, chrm, beg, end, strand)
            # The genes 
            if g_id in genes:
                genes[g_id].append(t_id)
            else:
                genes[g_id] = [t_id]

        if feature == 'UTR':
            if t_id in utr:
                utr[t_id].append((chrm, beg, end, strand))
            else:
                utr[t_id] = [(chrm, beg, end, strand)]

        if feature == 'CDS':
            if t_id in cds:
                cds[t_id].append((chrm, beg, end, strand))
            else:
                cds[t_id] = [(chrm, beg, end, strand)]

        if feature == 'exon':
            if t_id in exons:
                exons[t_id].append((chrm, beg, end, strand))
            else:
                exons[t_id] = [(chrm, beg, end, strand)]

    # Add cdses to their transcripts and order them
    # NEEDS TO BE FIRST! :S
    for (t_id, cdses) in cds.iteritems():
        for cds in cdses:
            transcripts[t_id].add_cds(cds)

    # Add utrs to their transcripts, determine 3' or 5', and order them
    for (t_id, utrs) in utr.iteritems():
        for utr in utrs:
            transcripts[t_id].add_utr(utr)

    # Add exons to their transcripts, determine 3' or 5', and order them
    for (t_id, exons) in exons.iteritems():
        for exon in exons:
            transcripts[t_id].add_exon(exon)

    return (transcripts, genes)

def skip_lines(file_path, nr):
    handle = open(file_path, 'rb')
    """ Skip 'nr' lines in the file handle and return it. If before those 'nr' a
    'chr' is encountered, call function on itself, skipping only as far as
    needed """
    for i in range(nr):
        line = handle.next()
        if line.startswith('chr'):
            handle.close()
            return skip_lines(file_path, i)
    return handle

def write_beds(transcripts, bed_dir, chr_sizes, *opts):

    (merge, extend, no_overlapping, chr1) = opts

    """" Write bed files for 3utr, cds, and 5utr introns and exons.
    Write bed files for intergenic regions and annotated transcription
    termination sites (aTTS).
    """
    output_names = ['five_utr_exons', 'five_utr_introns', 'three_utr_exons',
                   'three_utr_introns', 'cds_exons', 'cds_introns',
                   'aTTS_coding', 'aTTS_processed']

    output = dict((name, path_join(bed_dir, name)+'.bed') for name in output_names)

    # BED-files with chr, beg, end, strand.
    fiveUTR_exons = open(output['five_utr_exons'], 'wb')
    fiveUTR_introns = open(output['five_utr_introns'], 'wb')
    threeUTR_exons = open(output['three_utr_exons'], 'wb')
    threeUTR_introns = open(output['three_utr_introns'], 'wb')
    CDS_exons = open(output['cds_exons'], 'wb')
    CDS_introns = open(output['cds_introns'], 'wb')
    aTTS_coding = open(output['aTTS_coding'], 'wb')
    aTTS_processed = open(output['aTTS_processed'], 'wb')
    # intergenic will be added at the bottom inside subfunction

    for ts_id, ts in transcripts.iteritems():

        # CDS exons
        for exon in ts.cds.exons:
            beg = str(exon[1])
            end = str(exon[2])
            if exon[2]-exon[1] < 3:
                continue
            outp = '\t'.join([exon[0], beg, end, '.', '.', exon[3]])
            CDS_exons.write(outp+'\n')

        # CDS introns
        for intr in ts.get_cds_introns():
            if feature_length(intr) < 2:
                continue
            outp = '\t'.join([intr[0], intr[1], intr[2], '.', '.', intr[3]])
            CDS_introns.write(outp+'\n')

        #5 UTR introns
        for intr in ts.get_utr_introns(side=5):
            if feature_length(intr) < 2:
                continue
            outp = '\t'.join([intr[0], intr[1], intr[2], '.', '.', intr[3]])
            fiveUTR_introns.write(outp+'\n')

        #5 UTR exons
        for exon in ts.five_utr.exons:
            beg = str(exon[1])
            end = str(exon[2])
            if exon[2]-exon[1] < 3:
                continue
            outp = '\t'.join([exon[0], beg, end, '.', '.', exon[3]])
            fiveUTR_exons.write(outp+'\n')

        #3 UTR introns
        for intr in ts.get_utr_introns(side=3):
            if feature_length(intr) < 2:
                continue
            outp = '\t'.join([intr[0], intr[1], intr[2], '.', '.', intr[3]])
            threeUTR_introns.write(outp+'\n')

        #3 UTR exons
        for exon in ts.three_utr.exons:
            beg = str(exon[1])
            end = str(exon[2])
            if exon[2]-exon[1] < 3:
                continue
            outp = '\t'.join([exon[0], beg, end, '.', '.', exon[3]])
            threeUTR_exons.write(outp+'\n')

        #aTTS coding transcripts
        if ts.t_type == 'protein_coding':

            # Only get it if coding_transcript has UTR. If not, don't get it.
            aTTS = ts.get_aTTS(with_utr=True)
            if aTTS == 0:
                continue

            beg = str(aTTS[1])
            end = str(aTTS[2])
            outp = '\t'.join([aTTS[0], beg, end, '.', '.', aTTS[3]])
            aTTS_coding.write(outp+'\n')

        if ts.t_type == 'processed_transcript':
            #aTTS from last exon
            aTTS = ts.get_aTTS(with_utr=False)
            beg = str(aTTS[1])
            end = str(aTTS[2])
            outp = '\t'.join([aTTS[0], beg, end, '.', '.', aTTS[3]])
            aTTS_processed.write(outp+'\n')


    fiveUTR_exons.close()
    fiveUTR_introns.close()
    threeUTR_exons.close()
    threeUTR_introns.close()
    CDS_exons.close()
    CDS_introns.close()
    aTTS_coding.close()
    aTTS_processed.close()

    # Extend the aTTS files by 40nt
    extend_by = 40
    if extend:
        output = extend_aTTS(output, chr_sizes, extend_by)

    # Create merged versions of all output if requested
    if merge:
        output = merge_output(bed_dir, output)

    # Remove overlap between the .bed files
    if no_overlapping:
        output = remove_overlap(bed_dir, output)

    ## Create the intergenic bed as the adjoint of all the bed files in output
    output = get_intergenic(bed_dir, output, chr_sizes, chr1)

    # remove all lines in the output if the length is < 1
    # THIS IS TODO -- only if you have more overlap problems and the such.
    #output = remove_bp(output)

    return output

def get_intergenic(bed_dir, output, chr_sizes, chr1):
    """Create a .bed file for the whole genome. Then successively remove each
    region from this file. Optionally divide into 3UTR proximal/distal..."""

    genBed = get_genomeBed(chr_sizes, bed_dir, chr1)
    output['intergenic'] = genBed

    keyword= ''
    all_other = output.keys()
    all_other.remove('intergenic')
    output = subtractor(['intergenic'], all_other, keyword, output)

    return output

def get_genomeBed(chrm_sizes, bed_dir, chr1):
    genBedPath = path_join(bed_dir, 'intergenic.bed')
    genBed = open(genBedPath, 'wb')

    if chr1:
        stop = 2
    else:
        stop = 25
    for line in [l for l in open(chrm_sizes, 'rb')][1:stop]:
        (chrm, end) = line.split()
        out_plus = '\t'.join([chrm, '1', end, '.', '.', '+'])+'\n'
        out_minus = '\t'.join([chrm, '1', end, '.', '.', '-'])+'\n'
        genBed.write(out_plus)
        genBed.write(out_minus)

    genBed.close()

    return genBedPath


def feature_length(feature):
    """Return the length of a (chr, beg, end, strnd) tuple """
    return (int(feature[2])-int(feature[1]))

def remove_overlap(bed_dir, output):
    """ Make sure that 3UTR_exons has no overlap with either aTTS_coding or
    aTTS_processed. Make sure that cds_exons has no overlap with aTTS_processed.
    Make sure that aTTS_processed does not overlap with aTTS_coding.
    """
    # 1) Remove entries of aTTS_coding from aTTS_processed
    aTTS_nonc = ['aTTS_processed']
    aTTS_cod = ['aTTS_coding']

    keyword = 'wo_coding'
    output = subtractor(aTTS_nonc, aTTS_cod, keyword, output)

    # 2) For all introns and exons, remove overlap with aTTS
    introns = ['three_utr_introns', 'five_utr_introns', 'cds_introns']
    exons = ['five_utr_exons', 'cds_exons', 'three_utr_exons']
    int_ex = introns + exons

    aTTS = ['aTTS_coding', 'aTTS_processed']
    keyword = 'nov_atts'
    output = subtractor(int_ex, aTTS, keyword, output)

    # 3) For all introns, remove all exons.
    keyword = 'nov_exon'
    output = subtractor(introns, exons, keyword, output)

    # 4) For CDS exons, remove 3/5 UTR exons
    cds_exon = ['cds_exons']
    other_exons = ['five_utr_exons', 'three_utr_exons']

    keyword = 'nov_cds'
    output = subtractor(cds_exon, other_exons, keyword, output)

    return output

def subtractor(a_names, b_names, keyword, output):

    for a_name in a_names:

        count = 0
        for b_name in b_names:
            count += 1

            a_file = output[a_name]
            b_file = output[b_name]

            ending = '_'+keyword+'_'+str(count)+'.bed'
            outp =  output[a_name].rpartition('.')[0] + ending

            output[a_name] = my_subtractBed(a_file, b_file, outp)

    return output

def bed_length(bed_file):
    """Return the summed lengths of all entries in bed file"""
    size = 0
    for entry in open(bed_file, 'rb'):
        (chrm, beg, end) = entry.split()[:3]
        size += float(end)-float(beg)

    return size

def my_subtractBed(file_a, file_b, outfile):

    cmd = 'subtractBed -s -a {0} -b {1} > {2}'.\
            format(file_a, file_b, outfile)

    call(cmd, shell=True)
    return outfile

def extend_aTTS(inpot, chr_sizes, extend_by):
    """ Extend the aTTS sites and replace original file with extended """

    change_us = ['aTTS_coding', 'aTTS_processed']

    for name in change_us:
        atts_path = inpot[name]

        out_file = atts_path.rpartition('.')[0] + '_extended.bed'
        cmd = 'slopBed -i {0} -g {1} -b {2} > {3}'.\
                format(atts_path, chr_sizes, extend_by, out_file)

        call(cmd, shell=True)

        inpot[name] = out_file

    return inpot

def merge_output(bed_dir, output):
    """ Run mergeBed on all the files in output and replace original """
    # mergeBed will delete the two life-important dots!!
    # I must reintroduce them 

    for (name, bed_file) in output.items():

        out_path = bed_file.rpartition('.')[0] + '_merged.bed'
        out_file = open(out_path, 'wb')

        cmd = ['mergeBed','-s', '-i', bed_file]
        p = Popen(cmd, shell=False, stdout=PIPE, stderr=PIPE)

        for merged_entry in p.stdout:
            (chrm, beg, end, strand) = merged_entry.split()
            out_file.write('\t'.join([chrm, beg, end, '0', '0', strand]) + '\n')

        out_file.close()

        output[name] = out_path

    return output

def get_3utr_bed_all_exons(settings, outfile_path):
    """Get the regions of all 3UTRs in the annotation. Cluster them by the
    UTR-beg, and save for each cluster a 'super-UTR' that contains all the
    exons in that cluster. Save these regions to a bedfile in outfile_path.
    Make sure that the last exon in each UTR is extended by the value set in
    settings.extendby"""

    out_handle = open(outfile_path, 'wb')
    extendby = settings.extendby

    # Get transcripts and genes from annotation
    (transcripts, genes) = make_transcripts(settings.annotation_path)

    # Prune the genes transcripts for short transcripts

    transcripts = dict((ts, ts_obj) for (ts, ts_obj) in transcripts.iteritems() if
                       ts_obj.three_utr.length > settings.min_utrlen)

    # Prune the genes for short transcripts
    genes = {}
    for (ts_id, ts_obj) in transcripts.iteritems():
        gene_id = ts_obj.gene_id
        if gene_id in genes:
            genes[gene_id].append(ts_id)
        else:
            genes[gene_id] = [ts_id]

    # 1) single out the transcripts that belong to genes with only 1-exon
    # 3UTRS. Give them the old clustering treatment.
    one_exon_transcripts = {}
    # 2) For the rest, make a new clustering method
    multi_exon_transcripts = {}

    for ts_id, ts_obj in transcripts.iteritems():

        # Don't consider utrs without exons
        if ts_obj.three_utr.exons != []:

            # If has one exon and all other utrs for gene has one as well, add
            # to list of transcripts that will receive old clustering approach

            if len(ts_obj.three_utr.exons) == 1:

                keeper = True
                # see if all other transcripts of this gene also has only one exon
                for other_id in genes[ts_obj.gene_id]:
                    other_obj = transcripts[other_id]

                    if other_obj.three_utr.exons != []:
                        if len(other_obj.three_utr.exons) != 1:
                            keeper = False

                # If it has only 1 exon and so do all other utrs for this gene
                if keeper:
                    one_exon_transcripts[ts_id] = ts_obj

                # If it has only 1 exon -- but others in the family have more
                else:
                    multi_exon_transcripts[ts_id] = ts_obj

            # If it has more than 1 exon
            else:
                multi_exon_transcripts[ts_id] = ts_obj

    # NOTE you seem to end up with 4 different methods for determining the 3UTR
    # tail. I have a suggestion for determining the 3UTR tail in the
    # mega-merged ones that should leave you happy.
    # When merging, include in the 3UTR the most-downstream UTR-start and the
    # most-downstream UTR-end! This way you will never confound your results
    # with CDS exons!

    # DUDE! You have to be aware that many 3UTRs, when there are multiple of
    # them, do not coincide with CDS, but with introns. In these cases it is
    # entirely justifiable to include several UTRs from one gene. They would be
    # independend, except in the case of intron retention. You will lose a lot
    # of good 3utrs if you take such a rough approach as above.

    # TODO for the 1exon-utrs, you have a new and better way :) you will have
    # 1gene-> multiple utrs. However, it is just 1 method and it is robust.

    # Cluster and write the one-exon utrs to file. Return several utrs per
    # gene. One utr returned per cluster. UTRs don't overlap CDS.
    cluster_by_utrbeg_one_exon(one_exon_transcripts, genes, extendby, out_handle)
    debug()


    # Make mega-utrs from multi-exon buggers and write them to file
    # This option has one-gene -> multi-utrs
    #cluster_by_utrbeg_multi_exon(multi_exon_transcripts, genes, extendby, out_handle)

    # Make mega-merged-utr from genes with multiple exons
    # This option has one-gene -> one utr. Chooses for each gene the
    # downstream-most UTR-start and the downstream-most UTR-end
    merge_all_exons(multi_exon_transcripts, genes, extendby, out_handle)

    out_handle.close()

def cluster_by_utrbeg_one_exon(one_exon_transcripts, genes, extendby, out_handle):
    """Takes a set of Transcript objects whose genes should only have UTRs with
    one exon. Clusters the UTRs of the transcripts together by their
    UTR start sites. Removes or trims UTRs that overlap CDS exons for that
    gene."""

    chrms = ['chr'+str(nr) for nr in range(1,23) + ['X','Y','M']]
    tsdict = dict((chrm, {'+':[], '-':[]}) for chrm in chrms)

    for ts_id, ts_obj in one_exon_transcripts.iteritems():
        # The the chrm, beg, end, and strand of the first and last exon
        # if only one exon, they will be the same
        (chrm, beg, end, strand) = ts_obj.three_utr.exons[0]

        # We want to append the beginning of the 3UTR
        if strand == '+':
            tsdict[chrm][strand].append((beg, ts_id))
        if strand == '-':
            tsdict[chrm][strand].append((end, ts_id))

    # you will need this to count how many times the gene has been used
    gene_utrcount = dict((gene_id, 0) for gene_id in genes)

    # Cluster the reads as you go through them; when leaving a cluster write
    # the cluster to disc; but first check for rouge overlapping CDS exons
    for chrm, strand_dict in tsdict.iteritems():
        for strand, beg_and_id in strand_dict.iteritems():

            # Skip empty ones; could be that chr1 only is run
            if beg_and_id == []:
                continue

            # initial values for the loop
            clustsum = 0
            clustcount = 0
            this_cluster_begs = []
            this_cluster_ids = []

            for indx, begid in enumerate(sorted(beg_and_id)):
                val = int(begid[0])
                this_id = begid[1]

                clustsum = clustsum + val   # sum, for mean
                clustcount += 1             # n elements, for mean
                mean = clustsum/clustcount

                # If dist between new entry and cluster mean is < 20, keep in cluster
                if abs(val - mean) < 20:
                    this_cluster_begs.append(val)
                    this_cluster_ids.append(this_id)

                # If not, save the old cluster to file
                else:
                    # get all the 3UTR exons from the cluster
                    exon_list = []
                    for ts_id in this_cluster_ids:
                        for ex in one_exon_transcripts[ts_id].three_utr.exons:
                            (chrm, beg, end, strand) = ex
                            exon_list.append(((beg, end)))

                    # OK the problem seems first of all to be that the super
                    # utr is only one exon in length
                    super_utr = merge_set(exon_list)

                    # of course! These are the one-exon utrs. you cluster them,
                    # therefore you will only get the ones one 'super'-utr out.

                    # check for redundancy with the coding sequence exons of
                    gene_id = one_exon_transcripts[this_cluster_ids[0]].gene_id

                    stripped_utr = strip_set(genes[gene_id], strand,
                                             one_exon_transcripts, super_utr)

                    # getting gene_id from just one of the ts in the cluster
                    gene_utrcount[gene_id] += 1

                    exon_nr = len(super_utr)

                    for nr, exon in enumerate(super_utr):
                        beg = str(exon[0])
                        end = str(exon[1])

                        # Extend the utrs for the 3'-most exon
                        if (nr == exon_nr-1) and (strand == '+'):
                            end = int(end)
                            end += extendby
                            end = str(end)
                        if (nr == 0) and (strand == '-'):
                            beg = int(beg)
                            beg -= extendby
                            beg = str(beg)
                        e_nr = '_'+str(gene_utrcount[gene_id])+ '_'+str(nr+1)
                        out_handle.write('\t'.join([chrm, beg, end, gene_id +
                                                    e_nr, '0', strand])+'\n')

                    # Get the values for the new cluster
                    clustsum = val
                    clustcount = 1
                    this_cluster_begs = [val]
                    this_cluster_ids = [this_id]

def strip_set(ts_ids, strand, single_exon_transcripts, utr_exon_list):
    # For each cds exon, check if it falls within any of the exons in the
    # utr_exon_list. Trim or delete the utr_exon accordingly.
    # THe approach is to check the utr exons against all the cds exons. for
    # each cds exon the utr exon can be marked as deleted or trimmed. when the
    # checking round is over the new utr exon is marked as deleted (if at least
    # one deletion mark) or trimmed (the most 3' conservative one if several).
    new_exon_list = []

    # Does it make sense it works if two cds fall across a single 3UTR?

    (beg, end) = (0, 1)

    if strand == '+':
        for utr_ex in utr_exon_list:
            temp_exons = []
            for ts_id in ts_ids:
                for cds_exon in single_exon_transcripts[ts_id].cds.exons:
                    cds_ex = (cds_exon[1], cds_exon[2])

                    # If the CDS-beg lands in the UTR
                    if utr_ex[beg] < cds_ex[beg] < utr_ex[end]:

                        # If the CDS-end lands in the UTR too
                        if utr_ex[beg] < cds_ex[end] < utr_ex[end]:
                            #if 'ENST00000449969' in ts_ids:
                                #debug()
                            temp_exons.append((cds_ex[end], utr_ex[end]))
                            continue
                        # If the CDS-end extends beying the UTR, add empty
                        temp_exons.append((0,0)) # if encountered, scrap UTR

                    # If the CDS-beg is outside the UTR
                    else:
                        # If the CDS-end lands in the UTR
                        if utr_ex[beg] < cds_ex[end] < utr_ex[end]:
                            #if 'ENST00000449969' in ts_ids:
                                #debug()
                            temp_exons.append((cds_ex[end], utr_ex[end]))
                            continue
                        # If neither the CDS-beg or the CDS-end lands in the utr,
                        # add nothing. the UTR is untouched for now

            if temp_exons == []:
                new_exon_list.append(utr_ex)
            # Collect the most conservative utr in temp_exons
            else:
                # Sort and choose the last in case of + strand
                temp_exons.sort()
                new_exon_list.append(temp_exons[-1])

    if strand == '-':
        for utr_ex in utr_exon_list:
            temp_exons = []
            for ts_id in ts_ids:
                for cds_exon in single_exon_transcripts[ts_id].cds.exons:
                    cds_ex = (cds_exon[1], cds_exon[2])

                    # If the CDS-end lands in the UTR
                    if utr_ex[beg] < cds_ex[end] < utr_ex[end]:

                        # If the CDS-beg lands in the UTR too, save slice
                        if utr_ex[beg] < cds_ex[beg] < utr_ex[end]:
                            temp_exons.append((utr_ex[beg], cds_ex[beg]))
                            break
                        # If the CDS-beg started before the UTR, add empty
                        temp_exons.append((0,0)) # if encountered, scrap UTR

                    # If the CDS-end is outside the UTR
                    else:
                        # If the CDS-beg lands in the UTR, save slice
                        if utr_ex[beg] < cds_ex[beg] < utr_ex[end]:
                            temp_exons.append((utr_ex[beg], cds_ex[beg]))
                            break
                        # If neither the CDS-end or the CDS-beg lands in the utr,
                        # add nothing. the utr is untouched for now.

            if temp_exons == []:
                new_exon_list.append(utr_ex)
            # Collect the most conservative utr in temp_exons
            else:
                # Sort and choose the first in case of + strand
                # The assumption is that the 0 coordiante of the tuples
                # should be the same, so the sort will act o the 1 coordinate
                temp_exons.sort()
                new_exon_list.append(temp_exons[0])

    if new_exon_list != utr_exon_list:
        if strand == '-':
        # TODO you found a star example -- and your code seems to work. verify
        # with a few more examples.
            debug()

    return new_exon_list

def merge_all_exons(multi_exon_transcripts, genes, extendby, out_handle):
    # Go through all the multi-exon transcripts; get all the transcripts for
    # that gene; merge the 3UTR exons of all those transcripts together; write
    # to file. Add those transcripts to a set of completed transcripts. Check
    # all multi exon transcripts against this set before starting the process.

    ts_written = set()
    #genes = dict((gene, tscripts) for gene, tscripts in genes.iteritems() if )
    for ts_id in multi_exon_transcripts:

        # Skip ts if you have already processed it
        if ts_id in ts_written:
            continue

        # Get all the transcripts of this transcripts' gene
        ts_obj = multi_exon_transcripts[ts_id]
        gene_id = ts_obj.gene_id
        ts_list = genes[gene_id]

        # Get a super utr from all these transcripts
        super_utr = get_super_utr(ts_list, multi_exon_transcripts)

        # Add the transcripts to the 'dunnit' set
        for ts_id in ts_list:
            ts_written.add(ts_id)

        # Write the super-exon
        chrm = ts_obj.chrm
        strand = ts_obj.strand

        exon_nr = len(super_utr)

        # Write each exon to the bed-file
        for nr, exon in enumerate(super_utr):
            beg = exon[0]
            end = exon[1]

            # Extend the utrs for the 3'-most exon
            if (nr == exon_nr-1) and (strand == '+'):
                end += extendby
                end = str(end)
            if (nr == 0) and (strand == '-'):
                beg -= extendby
                beg = str(beg)
            e_nr = '_1_'+str(nr+1)
            out_handle.write('\t'.join([chrm, beg, end, gene_id+e_nr, '0', strand])
                             + '\n')

def cluster_by_utrbeg_multi_exon(multi_exon_transcripts, genes, extendby,
                                 out_handle):
    """Cluster multi exon transcripts by their utr-start sites. From the utrs
    in these clusters, get the union of all their exons. As a result, some of
    the resulting exons will not exist in the annotation. They are
    merger-exons. They represent the sequence potential of the 3UTR."""

    chrms = ['chr'+str(nr) for nr in range(1,23) + ['X','Y','M']]
    tsdict = dict((chrm, {'+':[],'-':[]}) for chrm in chrms)

    for ts_id, ts_obj in multi_exon_transcripts.iteritems():
        # The the chrm, beg, end, and strand of the first and last exon
        # if only one exon, they will be the same
        (chrm, first_beg, first_end, strand) = ts_obj.three_utr.exons[0]
        (chrm, last_beg, last_end, strand) = ts_obj.three_utr.exons[-1]

        # We want to append the beginning of the 3UTR
        if strand == '+':
            tsdict[chrm][strand].append((first_beg, ts_id))
        if strand == '-':
            tsdict[chrm][strand].append((last_end, ts_id))

    # go through the annotated 3UTR starts and cluster them
    gene_utrcount = dict((gene_id, 0) for gene_id in genes)

    for chrm, strand_dict in tsdict.iteritems():
        for strand, beg_and_id in strand_dict.iteritems():
            # Skip empty ones; could be that chr1 only is run
            if beg_and_id == []:
                continue
            clustsum = 0
            clustcount = 0
            this_cluster_begs = []
            this_cluster_ids = []
            for indx, begid in enumerate(sorted(beg_and_id)):
                val = int(begid[0])
                this_id = begid[1]

                clustsum = clustsum + val   # sum, for mean
                clustcount += 1             # n elements, for mean
                mean = clustsum/clustcount

                # If dist between new entry and cluster mean is < 20, keep in cluster
                if abs(val - mean) < 20:
                    this_cluster_begs.append(val)
                    this_cluster_ids.append(this_id)

                # If not, save the old cluster to file
                else:
                    super_utr = get_super_utr(this_cluster_ids, multi_exon_transcripts)
                    # getting gene_id from just one of the ts in the cluster
                    gene_id = multi_exon_transcripts[this_cluster_ids[0]].gene_id
                    gene_utrcount[gene_id] += 1

                    exon_nr = len(super_utr)

                    for nr, exon in enumerate(super_utr):
                        beg = str(exon[0])
                        end = str(exon[1])

                        # Extend the utrs for the 3'-most exon
                        if (nr == exon_nr-1) and (strand == '+'):
                            end += extendby
                        if (nr == 0) and (strand == '-'):
                            beg -= extendby
                        e_nr = '_'+str(gene_utrcount[gene_id])+ '_'+str(nr+1)
                        out_handle.write('\t'.join([chrm, beg, end, gene_id +
                                                    e_nr, '0', strand])+'\n')

                    # Get the values for the new cluster
                    clustsum = val
                    clustcount = 1
                    this_cluster_begs = [val]
                    this_cluster_ids = [this_id]

def get_unequivocal_utr(ts_ids, one_exon_transcripts):
    """Given the utrs of all the transcripts, return a merged utr that starts
    at the 3'-most start of all the utrs and ends at the 3'-most end"""

    # 1) extract a tuple-list of the utr exons
    # trim all exons at the 5' end to begin at the 3'most beginning
    # (4,8), (1,3), (1,7) - > (4,8), (4,7) (the 1,3 is removed)
    oet = one_exon_transcripts
    # get a list of exons. each row is a different ts


def merge_set(exon_list):

    exon_list.sort()

    if len(exon_list) == 1:
        return exon_list

    if len(exon_list) == 0:
        debug()

    mega_utr = []

    beg = exon_list[0][0]
    end = exon_list[0][1]

    utrlen = len(exon_list)-1

    for exnr, exon in enumerate(exon_list[1:]):

        this_beg = exon[0]
        this_end = exon[1]

        # Check if previous end falls within this exon
        if end > this_beg:
            # Check if this exon extends further than the previous one
            if this_end > end:
                end = this_end

        if this_beg > end:
            mega_utr.append((beg, end))

            beg = this_beg
            end = this_end

        ## Special case for last exon. Not part of main algorithm.
        if exnr == utrlen -1:
            mega_utr.append((beg, end))
            break

    return mega_utr

def get_a_polyA_sites_bed(settings, outfile_path):
    """Get the polyA sites (end position of last exon) of annotated 3UTRs.
    Save these positions to a bedfile in outfile_path. Cluster the polyA sites
    and return the averages of the clusters."""

    out_handle = open(outfile_path, 'wb')
    # Get transcripts from annotation
    (transcripts, genes) = make_transcripts(settings.annotation_path)

    # Only get the longest of the 3UTRs when they overlap
    chrms = ['chr' + str(nr) for nr in range(1,23) + ['X','Y','M']]
    tsdict = dict((chrm, dict((('+', []), ('-', [])))) for chrm in chrms)

    for ts_id, ts_obj in transcripts.iteritems():
        # Don't consider utrs without exons
        if ts_obj.three_utr.exons != []:

            # The the chrm, beg, end, and strand of the first and last exon
            # if only one exon, they will be the same
            (chrm, first_beg, first_end, strand) = ts_obj.three_utr.exons[0]
            (chrm, last_beg, last_end, strand) = ts_obj.three_utr.exons[-1]

            # Skip the utrs shorter than utrlen
            if last_end-first_beg < settings.min_utrlen:
                continue

            if strand == '+':
                tsdict[chrm][strand].append(last_end)
            if strand == '-':
                tsdict[chrm][strand].append(first_beg)

    # go through the annotated 3UTR ends and cluster them
    for chrm, strand_dict in tsdict.items():
        for strand, ends in strand_dict.items():
            # Skip empty ones; could be that chr1 only is run
            if ends == []:
                continue
            clustsum = 0
            clustcount = 0
            this_cluster = []
            clusters = []
            for indx, val in enumerate(sorted(ends)):
                ival = int(val)

                clustsum = clustsum + ival
                clustcount += 1
                mean = clustsum/clustcount

                # If dist between new entry and cluster mean is < 40, keep in cluster
                if abs(ival - mean) < 40:
                    this_cluster.append(ival)
                else: # If not, start a new cluster, and save the old one
                    clusters.append(this_cluster)
                    clustsum = ival
                    clustcount = 1
                    this_cluster = [ival]

            # Append the last cluster
            clusters.append(this_cluster)
            # Get the mean of the clusters with length greater than one
            cl_means = [int(math.floor(sum(clu)/len(clu))) for clu in clusters]

            # Save the means to the output file. write polyA site 'toward' ts.
            if strand == '-':
                for m in cl_means:
                    out_handle.write('\t'.join([chrm, str(m), str(m+1), '0',
                                                '0', '-']) + '\n')
            if strand == '+':
                for m in cl_means:
                    out_handle.write('\t'.join([chrm, str(m-1), str(m), '0',
                                                '0', '+']) + '\n')
    out_handle.close()

def get_seqs(utr_dict, hgfasta):
    """Use the pyfasta module to get sequences quickly from an indexed version
    of the human genome fasta file"""
    f = Fasta(hgfasta)
    seq_dict = {}
    for ts_id, ts_param in utr_dict.iteritems():
        (chrm, beg, end, strand) = ts_param
        seq_dict[ts_id] = f.sequence({'chr':chrm, 'start':beg, 'stop':end,
                                      'strand':strand}).upper()

    return seq_dict

def get_all_exons(annotation, chr1):
    if chr1: chrom1 = '_chr1'
    else: chrom1 = ''

    out_path = os.path.join(os.path.dirname(annotation),
                            'just_exons'+chrom1+'.bed')

    in_handle = open(annotation, 'rb')
    out_handle = open(out_path, 'wb')


    for line in open(annotation, 'rb'):
        # Skip first commented out lines
        if line[0]=='#':
            continue
        (chrm, source, feature, beg, end, d, strand) = line.split('\t')[:7]
        if feature == 'exon':
            out_handle.write('\t'.join([chrm, beg, end, source, d, strand]) +
                             '\n')

def main():
    t1 = time.time()

    chr1 = True

    #chr1 = False

    annotation = '/users/rg/jskancke/phdproject/3UTR/'\
            'gencode5/gencode5_annotation.gtf'
    if chr1:
        annotation = '/users/rg/jskancke/phdproject/3UTR/'\
                'gencode5/gencode5_annotation_chr1.gtf'

    ### TESTING START

    #(transcripts, genes) = make_transcripts(annotation)

    #new_transcripts = cluster_by_utrbeg(transcripts)

    # Only get the longest of the 3UTRs when they overlap
    #get_all_exons(annotation, chr1)
    #debug()


    ### TESTING OVER

    # Read annotation file to get exons and introns of 3utr, cds, and 5utr
    (transcripts, genes) = make_transcripts(annotation)

    debug()

    all_begends = [(ts.chrm, ts.beg, ts.end) for ts in transcripts.values()]
    set_begends = set(all_begends)
    print time.time() - t1
    #transcripts = []

    # Create bed-files for introns and exons in 3-5 utr and cds, as well as aTTS
    # sites and intergenic regions
    bed_dir = '/users/rg/jskancke/phdproject/3UTR/key_bedfiles'
    chr_s = '/users/rg/jskancke/phdproject/3UTR/BodyMap_3UTR/overlap/hg19genome'

    (merge, extend, no_overlapping) = (True, True, True)
    opts = (merge, extend, no_overlapping, chr1)

    write_beds(transcripts, bed_dir, chr_s, *opts)

    print time.time() - t1


if __name__ == '__main__':
    #anno= '/users/rg/jskancke/phdproject/3UTR/'\
            #'gencode5/gencode5_annotation_chr1.gtf'
    #get_3utr_bed(anno, 'hell')
    main()

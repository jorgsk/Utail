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


class Gene(object):

    def __init__(self, chrm, beg, end, strand):
        self.chrm = chrm
        self.beg = beg
        self.end = end
        self.strand = strand


class Transcript(object):

    def __init__(self, ID, t_type, chrm, beg, end, strand):
        self.ID = ID
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

        if self.strand == '-':
            if utr_beg <= cds_beg:
                self.three_utr.exons.append(utr)
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

    # Skip first 4 lines of annotation (for gencode at least.......)
    a_handle = skip_lines(annotation, 7)

    # Initialize dictionaries
    utr, cds, exons, transcripts = (dict(), dict(), dict(), dict())

    # Go though annotation and classify transcripts, genes, and utr + cds exons
    for line in a_handle:
        (chrm, source, feature, beg, end, d, strand, d, d, d, d, t_ID, d, d,
        d, d, d, d, d, t_type) = line.split()[:20]
        beg = int(beg)
        end = int(end)

        t_type = t_type.rstrip(';').strip('"') #   remove some characters
        t_ID = t_ID.rstrip(';').strip('"')

        if feature == 'transcript':
            transcripts[t_ID] = Transcript(t_ID, t_type, chrm, beg, end, strand)

        if feature == 'UTR':
            if t_ID in utr:
                utr[t_ID].append((chrm, beg, end, strand))
            else:
                utr[t_ID] = [(chrm, beg, end, strand)]

        if feature == 'CDS':
            if t_ID in cds:
                cds[t_ID].append((chrm, beg, end, strand))
            else:
                cds[t_ID] = [(chrm, beg, end, strand)]

        if feature == 'exon':
            if t_ID in exons:
                exons[t_ID].append((chrm, beg, end, strand))
            else:
                exons[t_ID] = [(chrm, beg, end, strand)]

    # Add cdses to their transcripts and order them
    for (t_ID, cdses) in cds.iteritems():
        for cds in cdses:
            transcripts[t_ID].add_cds(cds)

    # Add utrs to their transcripts, determine 3' or 5', and order them
    for (t_ID, utrs) in utr.iteritems():
        for utr in utrs:
            transcripts[t_ID].add_utr(utr)

    # Add exons to their transcripts, determine 3' or 5', and order them
    for (t_ID, exons) in exons.iteritems():
        for exon in exons:
            transcripts[t_ID].add_exon(exon)

    return transcripts

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

    for ts_ID, ts in transcripts.iteritems():

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
            out_file.write('\t'.join([chrm, beg, end, '.', '.', strand]) + '\n')

        out_file.close()

        output[name] = out_path

    return output

def cluster_by_utrbeg(transcripts):
    """Place the transcripts (that have only one utr exon) in clusters depending
    on the start site of their 3utr. Choose from each cluster only the
    transcript with the longest 3utr.
    """
    chrms = ['chr'+str(nr) for nr in range(1,23) + ['X','Y','M']]
    tsdict = dict((chrm, {'+':{(0,0):[]},'-':{(0,0):[]}}) for chrm in chrms)

    for ts_id, ts_obj in transcripts.iteritems():
        # Only deal with the ones with one utr exon
        if len(ts_obj.three_utr.exons) == 1:
            found = False # trace if it is being placed in an existing cluster

            # The the chrm, beg, end, and strand of the 3UTR
            (chrm, beg, end, strand) = ts_obj.three_utr.exons[0]

            # Clustering by beg/end depending on +/- for the strand
            if strand == '+':
                utr_start = beg
            else:
                utr_start = end

            for cluster, c_ts in tsdict[chrm][strand].items():
                # if current beg is in in the range of the cluster
                if (utr_start > cluster[0]) and (utr_start < cluster[1]):
                    # Extend cluster by +/- new utr_start if needed
                    (new_min, new_max) = (utr_start-20, utr_start+20)
                    cluster_changed = False
                    if new_min < cluster[0]:
                        n_cluster = (new_min, cluster[1])
                    if new_max > cluster[1]:
                        n_cluster = (cluster[0], new_max)
                    if cluster_changed:
                        # Get the ids, remove old, and add new
                        ids = tsdict[chrm][strand].items()
                        tsdict[chrm][strand].pop(cluster)
                        tsdict[chrm][strand][n_cluster] = ids
                        cluster = n_cluster

                    # Add ts_id to cluster (new cluster or not)
                    tsdict[chrm][strand][cluster].append(ts_id)
                    # Claim as found and break out of the for loop
                    found = True
                    break

            # If not found in a cluster, create a new cluster for the ts
            if not found:
                tsdict[chrm][strand][(utr_start-20, utr_start+20)] = [ts_id]

    # Remove the two initial, dummy clusters
    for chrm in chrms:
        for strand in ['+', '-']:
            tsdict[chrm][strand].pop((0,0))

    # Now select from the tsdict clusters only the ts with the longest utr
    new_ts = {}

    for chrm, strand_dict in tsdict.items():
        for strand, cluster_dict in strand_dict.items():
            for cluster, ids in cluster_dict.items():
                # a [(utrlength, id)] kind of list
                utr_lens =\
                [(feature_length(transcripts[ts_id].three_utr.exons[0]), ts_id)
                 for ts_id in ids]
                longest_length, longest_id = sorted(utr_lens, reverse=True)[0]
                # Choose only the ts with the longest 3UTR
                new_ts[longest_id] = transcripts[longest_id]

    return new_ts

def get_3utr_bed(annotation_path, outfile_path, settings):
    """Get the annotated regions of 3UTRs that have only one exon. Save these
    regions to a bedfile in outfile_path"""

    out_handle = open(outfile_path, 'wb')
    # Get transcripts from annotation
    transcripts = make_transcripts(annotation_path)

    # Only get the longest of the 3UTRs when they overlap
    new_transcripts = cluster_by_utrbeg(transcripts)
    transcripts = new_transcripts

    for ts_id, ts_obj in transcripts.items():
        if len(ts_obj.three_utr.exons) == 1:
            (chrm, beg, end, strand) = ts_obj.three_utr.exons[0]

            # Skip the utrs shorter than utrlen
            if end-beg < settings.min_utrlen:
                continue

            # If extendby, extend the 3UTR in the direction of extendby
            if settings.extendby:
                if strand == '+':
                    end = end + settings.extendby
                if strand == '-':
                    beg = beg - settings.extendby

            out_handle.write('\t'.join([chrm, str(beg), str(end), ts_obj.ID,
                                            '0', strand])+'\n')
    out_handle.close()

def get_polyA_sites_bed(annotation_path, outfile_path, settings):
    """Get the annotated regions of 3UTRs that have only one exon. Save these
    regions to a bedfile in outfile_path"""

    out_handle = open(outfile_path, 'wb')
    # Get transcripts from annotation
    transcripts = make_transcripts(annotation_path)

    # Only get the longest of the 3UTRs when they overlap
    chrms = ['chr' + str(nr) for nr in range(1,23) + ['X','Y','M']]
    tsdict = dict((chrm, dict((('+', []), ('-', [])))) for chrm in chrms)

    for ts_id, ts_obj in transcripts.iteritems():
        # Only deal with the ones with one utr exon
        if len(ts_obj.three_utr.exons) == 1:
            # The the chrm, beg, end, and strand of the transcripts
            (chrm, beg, end, strand) = ts_obj.three_utr.exons[0]

            # Skip the utrs shorter than utrlen
            if end-beg < settings.min_utrlen:
                continue

            if strand == '+':
                tsdict[chrm][strand].append(end)
            if strand == '-':
                tsdict[chrm][strand].append(beg)

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

    #transcripts = make_transcripts(annotation)

    #new_transcripts = cluster_by_utrbeg(transcripts)

    # Only get the longest of the 3UTRs when they overlap
    #get_all_exons(annotation, chr1)
    #debug()


    ### TESTING OVER

    # Read annotation file to get exons and introns of 3utr, cds, and 5utr
    transcripts = make_transcripts(annotation)

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

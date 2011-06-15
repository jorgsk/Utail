"""
**Module for representing GENCODE transcripts as class objects and performing
calculations on them.**
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
sys.path.append(os.path.join(here,'modules'))
sys.path.append(os.path.join(here,'modules/pyfasta'))
from pyfasta import FastaRecord
from pyfasta import Fasta
import math


class Region(object):

    def __init__(self, name):
        self.name = name

        # Note exons and introns should be on the form
        # (chr1, beg, end, strand), and sorted.
        self.exons = []
        self.length = 0

class Transcript(object):
    """
    Represents transcripts as they occur in GENCODE annotation files.

    Instances of the transcript class are initialized with the input:

    - transcript_id
    - gene_id
    - transcript type
    - chromosome
    - beg
    - end
    - strand

    These are fields from the GENCODE gft file.

    Further, it initializes by itself four attributes that will be filled as the
    annotation file is parsed

    - three_utr (holds the 3UTR exons and methods to act on them)
    - cds (as for CDS)
    - five_utr (as for 5UTR)
    - exons (to hold all exons associated with this transcript)
    """

    def __init__(self, t_id, g_id, t_type, chrm, beg, end, strand):
        """
        Take transcriptID, geneID, transcript type, chromosome, beginning and
        end coordinates, and strand. Add these to the transcript object.

        Additionally, create placeholders for 3UTR exons, CDS exons, 5UTR exons,
        and an empty list of all exons.
        """
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

        # the start_codon and stop_codon are only used for some annotation
        # formats
        self.start_codon = ''
        self.stop_codon = ''

    def add_cds(self, cds):
        """
        Add CDS exons to the CDS list of exons.
        """
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
        #Return the aTTS of a transcript. If processed, also return for those
        #who have not a 3UTR

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
        """
        Determine if UTR is 5' or 3' by comparing with cds_beg and strand and
        add to the corresponding list.

        Do nothing if no CDS was annotated
        """

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

    def __repr__(self):
        return self.ts_id

    def __str__(self):
        return "\nChrm\t{0}\nBeg\t{1}\nEnd\t{2}\nStrand\t{3}\n"\
                .format(self.chrm, self.beg, self.end, self.strand)

def gencode_reader(file_handle):

    (utr, cds, exons, transcripts, genes) = (dict(), dict(), dict(), dict(),
                                           dict())

    # Go though annotation and classify transcripts, genes, and utr + cds exons
    for line in file_handle:
        (chrm, source, feature, beg, end, d, strand, d, d, g_id, d, t_id, d, d,
        d, d, d, d, d, t_type) = line.split()[:20]

        beg = int(beg) - 1
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

    # Add cdses to their transcripts
    # NEEDS TO BE FIRST! :S
    for (t_id, cdses) in cds.iteritems():
        for cds in cdses:
            transcripts[t_id].add_cds(cds)

    # Add utrs to their transcripts
    for (t_id, utrs) in utr.iteritems():
        for utr in utrs:
            transcripts[t_id].add_utr(utr)

    # Add exons to their transcripts
    for (t_id, exons) in exons.iteritems():
        for exon in exons:
            transcripts[t_id].add_exon(exon)

    return (transcripts, genes)

def ensembl_reader(file_handle):
    """
    Different from GENCODE, ensembl doesn't have the 'transcript feature. You
    have to assemble it yourself.'
    """

    (utr, cds, exons, transcripts, genes) = (dict(), dict(), dict(), dict(),
                                           dict())

    (start_codons, stop_codons) = (dict(), dict())

    # Go though annotation and classify transcripts, genes, and utr + cds exons
    for line in file_handle:

        (chrm, t_type, feature, beg, end, d, strand,\
         d, d, g_id, d, t_id) = line.split()[:12]

        # skip the non-chromosome chrm entries (if len chrm > 3)
        if len(chrm) > 2:
            continue

        # skip non-protein coding 3UTRs
        t_type = t_type.rstrip(';').strip('"') #   remove some characters
        if t_type != 'protein_coding':
            continue

        chrm = 'chr'+chrm
        beg = int(beg) - 1
        end = int(end)

        t_id = t_id.rstrip(';').strip('"')
        g_id = g_id.rstrip(';').strip('"')

        # features: CDS, exon, stop_codon, start_codon

        if t_id not in transcripts:
            # initialize the transcript with the first featre you find
            transcripts[t_id] = Transcript(t_id, g_id, t_type, chrm, beg, end, strand)

            # The genes 
            if g_id in genes:
                genes[g_id].append(t_id)
            else:
                genes[g_id] = [t_id]

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

        if feature == 'start_codon':
            if t_id in start_codons:
                print('More than 1 start codon for 1 transxcript???!')
                #debug()
            else:
                start_codons[t_id] = (chrm, beg, end, strand)

        if feature == 'stop_codon':
            if t_id in stop_codons:
                print('More than 1 stop codon for 1 transxcript???!')
                #debug()
            else:
                stop_codons[t_id] = (chrm, beg, end, strand)

    # Add cdses to their transcripts
    # NEEDS TO BE FIRST! :S
    for (t_id, cdses) in cds.iteritems():
        for cd in cdses:
            transcripts[t_id].add_cds(cd)

    # Add stop codons to their transcripts
    for (t_id, stop_cdon) in stop_codons.iteritems():
        transcripts[t_id].stop_codon = stop_cdon

    # Add start codons to their transcripts
    for (t_id, start_cdon) in start_codons.iteritems():
        transcripts[t_id].start_codon = start_cdon

    # Add exons to their transcripts
    for (t_id, exns) in exons.iteritems():
        for exon in exns:
            transcripts[t_id].add_exon(exon)

    # Add utrs to their transcripts, determine 3' or 5'
    # For each exon, find out if it is before or after the start codon or stop
    # codon

    # Finally update the transcript beg/end from the utrs or exons
    for (t_id, ts) in transcripts.iteritems():

        # don't do anything if this transcript for some reason doesn'th ave
        # exons
        if ts.exons == []:
            continue

        # if not both start and stop codon are found, don't go further!
        if ts.start_codon == '' or ts.stop_codon == '':
            continue

        ts.exons.sort(key=itemgetter(1))

        ts.beg = ts.exons[0][1]
        ts.end = ts.exons[-1][2]

        if ts.strand == '+':
            for exon in ts.exons:

                # if exon is 5' to start codon
                if exon[2] <= ts.start_codon[2]:
                    ts.five_utr.exons.append(exon)

                # if exon surrounds start codon, slice out the 5UTR-bit
                if exon[1] < ts.start_codon[2] < exon[2]:
                    utr = (exon[0], exon[1], ts.start_codon[2], exon[3])
                    ts.five_utr.exons.append(utr)

                # if exon is 3' to stop codon
                if exon[1] >= ts.stop_codon[1]:
                    ts.three_utr.exons.append(exon)

                # if exon surrounds stop codon, slice out the 3UTR-bit
                if exon[1] < ts.stop_codon[1] < exon[2]:
                    utr = (exon[0], ts.stop_codon[1], exon[2], exon[3])
                    ts.three_utr.exons.append(utr)

        if ts.strand == '-':

            for exon in ts.exons:

                # if exon is 5' to start codon
                if exon[1] >= ts.start_codon[1]:
                    ts.five_utr.exons.append(exon)

                # if exon surrounds start codon, slice out the 5UTR-bit
                if exon[1] < ts.start_codon[1] < exon[2]:
                    utr = (exon[0], ts.start_codon[1], exon[2], exon[3])
                    ts.five_utr.exons.append(utr)

                # if exon is 3' to stop codon
                if exon[2] <= ts.stop_codon[2]:
                    ts.three_utr.exons.append(exon)

                # if exon surrounds stop codon, slice out the 3UTR-bit
                if exon[1] < ts.stop_codon[2] < exon[2]:
                    utr = (exon[0], exon[1], ts.stop_codon[2], exon[3])
                    ts.three_utr.exons.append(utr)

    return (transcripts, genes)

def make_transcripts(annotation, file_format='GENCODE'):
    """
    Loop through the (GENCODE) annotation file (the input), create instancees
    of :class:`Transcript`, and save these in a dictionary. Also create a
    dictionary, *genes*, with the transcript -> gene correspondence.

    Since we go from GFF/GTF to .BED format, it's important to subtract 1 from
    the beg-coordinate. 9 9 in GTF is 8 9 in BED.

    :returns: transcripts, genes
    """

    # Skip first 4 lines of annotation (for gencode at least.......)
    a_handle = skip_lines(annotation, 7)

    # Initialize dictionaries
    utr, cds, exons, transcripts, genes = (dict(), dict(), dict(), dict(),
                                           dict())

    if file_format == 'GENCODE':
        a_handle = skip_lines(annotation, 7)

        (transcripts, genes) = gencode_reader(a_handle)

    if file_format == 'ENSEMBL':
        a_handle = open(annotation, 'rb')

        (transcripts, genes) = ensembl_reader(a_handle)

    # Go through all the damn transcripts and their damn exons, utrs, and cds,
    # and sort them.
    for (ts_id, ts_obj) in transcripts.iteritems():
        ts_obj.five_utr.exons.sort()
        ts_obj.cds.exons.sort()
        ts_obj.three_utr.exons.sort()

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

def get_last_exons(transcripts):

    last_exons = {}

    for tsid, ts in transcripts.iteritems():

        if ts.cds.exons != []:
            gene_id = ts.gene_id

            if ts.strand == '+':
                lultim = ts.cds.exons[-1][2]

                if gene_id in last_exons:
                    if lultim > last_exons[gene_id]:
                        last_exons[gene_id] = lultim
                else:
                    last_exons[gene_id] = lultim

            if ts.strand == '-':
                lultim = ts.cds.exons[0][1]

                if gene_id in last_exons:
                    if lultim < last_exons[gene_id]:
                        last_exons[gene_id] = lultim
                else:
                    last_exons[gene_id] = lultim

    return last_exons

def get_3utr_bed_all_exons(settings, outfile_path):
    """
    Get the regions of all 3UTRs in the annotation (from *settings*). Cluster
    them by the UTR-beg, and save for each cluster a 'super-UTR' that contains
    all the exons in that cluster. Save these regions to a bedfile in
    outfile_path. Make sure that the last exon in each UTR is extended by the
    value set in *settings.extendby*

    :returns: writes 3UTR-bedfile to *outfile_path*
    """

    raw_path = os.path.splitext(outfile_path)[0] + '_raw.bed'
    raw_handle = open(raw_path, 'wb')
    extendby = settings.extendby

    # Get transcripts and genes from annotation
    if settings.chr1:
        (transcripts, genes) = make_transcripts(settings.annotation_path_chr1,
                                                settings.annotation_format)
    else:
        (transcripts, genes) = make_transcripts(settings.annotation_path,
                                               settings.annotation_format)

    all_transcripts = transcripts

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
                # Skip short utrs
                if ts_obj.three_utr.length < settings.min_utrlen:
                    continue

                keeper = True
                # see if all other transcripts of this gene also have only one
                # exon utrs (as long as the utr is long enough)
                for other_id in genes[ts_obj.gene_id]:
                    other_obj = transcripts[other_id]

                    if other_obj.three_utr.exons != []:

                        # if other object is long enough to be allowed to
                        # intervene
                        if other_obj.three_utr.length > settings.min_utrlen:
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
                # if it is long enough
                if ts_obj.three_utr.length > settings.min_utrlen:
                    multi_exon_transcripts[ts_id] = ts_obj


    # Cluster and write single-exon utrs
    one_exon_cluster_write(one_exon_transcripts, all_transcripts,  genes,
                           extendby, raw_handle)

    # Cluster and write the multi-exon utrs to file.
    cluster_by_utrbeg_multi_exon(multi_exon_transcripts, all_transcripts,
                                 genes, extendby, raw_handle)

    raw_handle.close()
    # There is a problem where the wrong utr is being extended; and we have a
    # rare 'utr within a utr' problem, but this sholdn't be a big deal.

    # TODO one way: since you're doing this utr-by-utr, maybe you can intersect
    # the extensions one by one? would take time but wtf is that now

    # TODO here you need to add a step that tries to extend all 3UTR ends by
    # 1000 (or so) nt; intersects the extension with all exons; and finally
    # trims the extensions so that no extension intersects into another genomic
    # featuer. Save the extensions in a dict. When you save to file, include
    # this extension and also save how much you extended by.
    # You need to go carefully through the code after to check all downstream
    # effects. For example, all "settings.extendby" must be modified to
    # "utr.extendby", or something like that. Each 3UTR can now have its own
    # 'extendby'. 
    # How to keep the information about extensions? You output to a .bed-file
    # and you are already using the value-field. Maybe use it for two tings?
    # value_extendby?

    # Remove utrs that have CDS exons in them and write to outfile_path
    remove_intersects_and_extend(raw_path, outfile_path, all_transcripts,
                                 settings)


def remove_intersects_and_extend(unfiltered_path, outfile_path, all_transcripts,
                                settings):
    """Remove UTR intersections with CDS exons. Remove UTR intersections with
    UTRs from other genes. If sequences should be extended, do just that.

    Check how far 3UTRs can be extended given the annotation. Extend that far.
    """

    temp_exons_path = os.path.join(os.path.dirname(outfile_path), 'temp_exons.bed')
    temp_extension_path = os.path.join(os.path.dirname(outfile_path),
                                       'temp_extension.bed')

    # Write non-3UTR exons
    temp_exons = open(temp_exons_path, 'wb')
    for ts_id, ts_obj in all_transcripts.iteritems():

        # Write all CDS exons
        if ts_obj.cds.exons != []:
            for cds in ts_obj.cds.exons:
                temp_exons.write('\t'.join([cds[0], str(cds[1]), str(cds[2]),
                                            cds[3]])+'\n')

        # And write and 5UTR exons
        if ts_obj.five_utr.exons != []:
            for f_utr in ts_obj.five_utr.exons:
                temp_exons.write('\t'.join([f_utr[0], str(f_utr[1]),
                                            str(f_utr[2]), f_utr[3]])+'\n')

    temp_exons.close()

    # Set of utrs to remove
    remove_these_utrs = set()

    # 1) Run intersect bed on bedfile with CDS exons and with itself
    # Ignore strandedness
    cmd1 = ['intersectBed', '-wb', '-a', temp_exons_path, '-b', unfiltered_path]
    cmd2 = ['intersectBed', '-wo', '-a', unfiltered_path, '-b', unfiltered_path]

    # Run the above command and loop through output
    f1 = Popen(cmd1, stdout=PIPE)
    f2 = Popen(cmd2, stdout=PIPE)

    # 2) Add utrs to-be-removed from the CDS intersection
    for line in f1.stdout:
        utr_id = line.split()[7]
        # remove id is gene_id + '_' + utr_nr
        remove_id = '_'.join(utr_id.split('_')[:2])
        remove_these_utrs.add(remove_id)

    # 3) Add utrs to-be-removed from the self intersection
    for line in f2.stdout:
        # Get the IDs of the intersecting utrs (most will be self-intersect)
        (d,d,d, utr_id_1, d,d,d,d,d, utr_id_2, d,d,d) = line.split()

        # If intersection occurred on utrs from DIFFERENT GENES, remove them
        # same-3utr intersection is accepted (because we choose the longest utr)
        # Hey, sometimes this doesn't work. There is overlap within the 3UTR.
        if utr_id_1.split('_')[0] != utr_id_2.split('_')[0]:
            remove_these_utrs.add('_'.join(utr_id_1.split('_')[:2]))
            remove_these_utrs.add('_'.join(utr_id_2.split('_')[:2]))

    # Create the output file
    outfile_handle = open(outfile_path, 'wb')

    extendby = settings.extendby

    # Get the maximum number of extensions you can get

    # NOTE the 3utrs in extend_me are ONLY those that actually intersect with
    # something. If they don't intersect, you can use the 1000 extension.
    extend_me = maximum_extension(unfiltered_path, all_transcripts, extendby,
                                  remove_these_utrs, temp_extension_path,
                                  temp_exons_path)

    # 4) Loop through the original file and write to the 'filtered' file
    with open(outfile_path, 'wb') as outfile_handle:

        if not extendby:

            for line in open(unfiltered_path, 'rb'):
                utr_id = '_'.join(line.split()[3].split('_')[:2])
                if utr_id not in remove_these_utrs:
                    outfile_handle.write(line)

        if extendby:

            for line in open(unfiltered_path, 'rb'):
                utrId = '_'.join(line.split()[3].split('_')[:2])

                # Only write the file if it's not marked for removal
                if utrId not in remove_these_utrs:

                    (chrm, beg, end, utr_ex_id, val, strand) = line.split()
                    id_split = utr_ex_id.split('_')

                    # Extend the utrs for the final exon in the utr-batch
                    if len(id_split) == 4: # it has an EXTENDBY mark

                        #only calculate when you do the extension. once is enough
                        utr_id = '_'.join(id_split[:3])

                        if  utr_id in extend_me:
                            extendby = extend_me[utr_id]

                        if strand == '+':
                            end = str(int(end) + extendby)

                        if strand == '-':
                            beg = str(int(beg) - extendby)

                        # update val to also include information about extension
                        val = val + '+{0}'.format(extendby)

                    # val = total Exons in UTR
                    # val = 3 means 3 exon in this UTR. See utr_id for which one
                    # BUT if this is an extended exon, val = val
                    # +extension_length

                    outfile_handle.write('\t'.join([chrm, beg, end, utr_id, val,
                                                   strand]) + '\n')

def maximum_extension(unfiltered_path, all_transcripts, extendby,
                      remove_these_utrs, temp_extension_path, temp_exons_path):

    """ Should you do it in a redundant manner for all UTR_IDs, and then discard
    those that don't make it? I say yes.

    Extend all 3UTRs in raw_path that have 'extendby' for them. Intersect this
    region
    """
    temp_ext = open(temp_extension_path, 'wb')

    for line in open(unfiltered_path, 'rb'):
        utrId = '_'.join(line.split()[3].split('_')[:2])

        # Skip utr_ids to be removed
        if utrId in remove_these_utrs:
            continue

        (chrm, beg, end, utr_ex_id, val, strand) = line.split()
        id_split = utr_ex_id.split('_')

        # Extend the utrs for the final exon in the utr-batch
        if len(id_split) == 4: # it has an EXTENDBY mark

            if strand == '+':
                beg = str(int(end)+1)
                end = str(int(beg) + extendby)

            if strand == '-':
                end = str(int(beg)-1)
                beg = str(int(end) - extendby)

            utr_id = '_'.join(id_split[:3])

            # Write the extension of the 3UTR
            temp_ext.write('\t'.join([chrm, beg, end, utr_id, val, strand]) +
                           '\n')

    temp_ext.close()

    # run intersectBed on temp_exons_path and temp_extension
    cmd1 = ['intersectBed', '-wb',
           '-a', temp_exons_path,
           '-b', temp_extension_path]

    # with -wb, he will write FIRST the overlapping portion of A, and THEN the
    # ORIGINAL B entry (the 1000nt extension's beg and ends)

    # Run the above command and loop through output
    f1 = Popen(cmd1, stdout=PIPE)

    extend_me = {}

    # Get the extension lenghts for CDS and 5prime exons
    for line in f1.stdout:
        (int_chrm, int_beg, int_end, int_strand) = line.split()[:4]
        (o_chrm, o_beg, o_end, o_ID, o_val, o_strand) = line.split()[4:]

        # Initialize with 1000 nt extension if o_ID is not found
        if o_ID not in extend_me:
            extend_me[o_ID] = 1000

        if o_strand == '+':
            max_extension = int(int_beg) - int(o_beg)

        if o_strand == '-':
            max_extension = int(o_end) - int(int_end)

        # get the minimum intersection
        extend_me[o_ID] = min(max_extension, extend_me[o_ID])

    return extend_me


def one_exon_cluster_write(one_exon_transcripts, all_transcripts, genes,
                           extendby, out_handle):
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

    gene_clusters = {}

    # Cluster the reads as you go through them; when leaving a cluster check
    # for rouge overlapping CDS exons; then write the cluster to disc
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

                # If not, filter UTR by checking for CDS sequences
                else:
                    # get all the 3UTR exons from the cluster
                    utr_list = []
                    for ts_id in this_cluster_ids:
                        this_utr = []
                        for ex in one_exon_transcripts[ts_id].three_utr.exons:
                            (chrm, beg, end, strand) = ex
                            this_utr.append(((beg, end)))

                        utr_list.append(this_utr)

                    gene_id = one_exon_transcripts[this_cluster_ids[0]].gene_id
                    # The the longest of the utrs in the cluster
                    print_utr = longest_utr(utr_list)
                    # Count how many times you have printed from this gene
                    gene_utrcount[gene_id] += 1

                    exon = print_utr[0]
                    beg = str(exon[0])
                    end = str(exon[1])

                    e_nr = '_'+str(gene_utrcount[gene_id])+ '_1_EXTENDME'

                    out_handle.write('\t'.join([chrm, beg, end, gene_id + e_nr,
                                                '1', strand]) + '\n')

                    # Start the variables for the new cluster
                    clustsum = val
                    clustcount = 1
                    this_cluster_begs = [val]
                    this_cluster_ids = [this_id]


def cluster_by_utrbeg_multi_exon(multi_exon_transcripts, all_transcripts,
                                 genes, extendby, out_handle):
    """Cluster multi exon transcripts by their utr-start sites."""

    chrms = ['chr'+str(nr) for nr in range(1,23) + ['X','Y','M']]
    tsdict = dict((chrm, {'+':[],'-':[]}) for chrm in chrms)

    for ts_id, ts_obj in multi_exon_transcripts.iteritems():
        # The the chrm, beg, end, and strand of the first and last exon
        # if only one exon, they will be the same
        # BUG this is not true: I got me a 3UTR whose exons aren't sorted!
        (chrm, first_beg, first_end, strand) = ts_obj.three_utr.exons[0]
        (chrm, last_beg, last_end, strand) = ts_obj.three_utr.exons[-1]

        # We want to append the beginning of the 3UTR
        if strand == '+':
            tsdict[chrm][strand].append((first_beg, ts_id))
        if strand == '-':
            tsdict[chrm][strand].append((last_end, ts_id))

    # keep count of the number of times a utr has been written to file from a
    # given gene
    gene_utrcount = dict((gene_id, 0) for gene_id in genes)

    # go through the annotated 3UTR starts and cluster them
    for chrm, strand_dict in tsdict.iteritems():
        for strand, beg_and_id in strand_dict.iteritems():
            # Skip empty ones; could be that chr1 only is run
            if beg_and_id == []:
                continue

            # Initiate cluster parameters
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
                    # Get all the exons of the cluster in a list
                    utr_list = []
                    for ts_id in this_cluster_ids:
                        this_utr = []
                        for ex in multi_exon_transcripts[ts_id].three_utr.exons:
                            (chrm, beg, end, strand) = ex
                            this_utr.append(((beg, end)))

                        utr_list.append(this_utr)

                    # get gene_id from just one of the ts in the cluster
                    gene_id = multi_exon_transcripts[this_cluster_ids[0]].gene_id
                    # Get the longest utr in the cluster
                    print_utr = longest_utr(utr_list)
                    # Count how many times you've written from this gene
                    gene_utrcount[gene_id] += 1
                    # Get the nr of exons in this utr
                    exon_nr = len(print_utr)

                    # Write each exon to file
                    for nr, exon in enumerate(print_utr):

                        beg = str(exon[0])
                        end = str(exon[1])

                        # Code for utr_exon_id gene_id + #utr from this gene +
                        # #exon in this utr _ possibly EXTENDME

                        ex_nr = str(nr+1)

                        e_nr = '_'.join([str(gene_utrcount[gene_id]), ex_nr])

                        if extendby:
                            # Mark the final utrs for extension
                            if (nr == exon_nr-1) and (strand == '+'):
                                e_nr += '_EXTENDME'

                            if (nr == 0) and (strand == '-'):
                                e_nr += '_EXTENDME'

                        out_handle.write('\t'.join([chrm, beg, end, gene_id +
                                                    '_' + e_nr, str(exon_nr) ,
                                                    strand])+'\n')

                    # Get the values for the new cluster
                    clustsum = val
                    clustcount = 1
                    this_cluster_begs = [val]
                    this_cluster_ids = [this_id]

def longest_utr(utr_list):
    new_list = []
    for utr in utr_list:
        if utr == []:
            new_list.append((0, (0,0)))
        else:
            length = sum([ex[1]-ex[0] for ex in utr])
            new_list.append((length, utr))

    new_list.sort()

    # Return the longest of the 3utrs
    return new_list[-1][-1]
    # Get length of each utr in utr_list

def get_a_polyA_sites_bed(settings, outfile_path):
    """Get the polyA sites (end position of last exon) of annotated 3UTRs.
    Save these positions to a bedfile in *outfile_path*. Cluster the polyA sites
    and return the averages of the clusters along with the number of reads part
    of the cluster.

    :returns: annotated TTS sites as bedfile in *outfile_path*
    """

    out_handle = open(outfile_path, 'wb')
    # Get transcripts from annotation
    (transcripts, genes) = make_transcripts(settings.annotation_path,
                                            settings.annotation_format)

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
    """
    Use the *pyfasta* module to get sequences quickly from an indexed version
    of the human genome fasta file.

    :returns: *utr_dict*, an exon -> sequence dictionary
    """
    f = Fasta(hgfasta, record_class=FastaRecord)
    #f = Fasta(hgfasta)
    seq_dict = {}
    for ts_id, ts_param in utr_dict.iteritems():
        (chrm, beg, end, strand) = ts_param
        if strand == '+':
            seq_dict[ts_id] = f.sequence({'chr':chrm,'start':beg+1, 'stop':end-1,
                                          'strand':strand}).upper()
        if strand == '-':
            seq_dict[ts_id] = f.sequence({'chr':chrm,'start':beg+2, 'stop':end,
                                          'strand':strand}).upper()
    return seq_dict

def main():
    t1 = time.time()

    chr1 = True

    #chr1 = False

    #annotation = '/users/rg/jskancke/phdproject/3UTR/'\
            #'gencode5/gencode5_annotation.gtf'
    annotation = '/users/rg/jskancke/phdproject/3UTR/Annotations/'\
            'Mus_musculus.NCBIM37.62.gtf'
    if chr1:
        #annotation = '/users/rg/jskancke/phdproject/3UTR/'\
                #'gencode5/gencode5_annotation_chr1.gtf'
        annotation = '/users/rg/jskancke/phdproject/3UTR/Annotations/'\
                'Mus_musculus.NCBIM37.62_chr1.gtf'

    an_frmt = 'GENCODE'
    an_frmt = 'ENSEMBL'

    #algo = get_3utr_bed_all_exons(annotation, 'somewhere')
    ### TESTING START

    (transcripts, genes) = make_transcripts(annotation, an_frmt)

    debug()

    #new_transcripts = cluster_by_utrbeg(transcripts)

    debug()


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

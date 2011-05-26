.. _introduction:

=====================
Introduction to Utail
=====================
This document gives the biological background for why this Utail was made, it
lists the dependencies of the Utail, and it shows two different ways of using
the Utail -- with or without an annotation.


On 3UTRs and RNA-seq
====================
The 3' untranslated regions (3'UTRs) of mRNAs are known to contain sequence and
structural motifs that regulate stability, localization and translation of the
mRNAs. By varying the length of 3'UTRs, the cell can thus include or exclude
regulatory regions.

High-throughput RNA-sequencing gives a high-resolution coverage of the
transcriptome. Utail uses RNA-seq reads to cover the 3UTRs of genes. It
outputs the length of the 3UTR as well as evidence of polyadenylation events if
possible.

We made Utail because we wanted a to characterize the 3UTR of transcripts from
any kind of RNA-seq experiment.

Mandatory dependencies
======================
These are required for the most basic function:

* Python 2.6.4 ++
* `bedTools <http://code.google.com/p/bedtools/>`_
* RNA-seq reads, either in bed-format or in `Gem <http://sourceforge.net/apps/mediawiki/gemlibrary/index.php?title=Gem_mapper_man_page>`_ format


Optional dependencies
=====================
For finding polyadenylation events the pipeline needs to remap reads to the
genome with the Gem mapper:

* `Gem mapper <http://sourceforge.net/apps/mediawiki/gemlibrary/index.php?title=Gem_mapper_man_page>`_

To use the mapper you have to create an index file from a fasta file of the
human genome (this takes hours) using::

    gem-do-index --reverse-complement emulate

Then you must provide a path to this index file in UTR_SETTINGS

Included dependencies
=====================
The following software packages are shipped with Utail for ease of use.
Credits go to their authors.

* `pyFasta <https://github.com/brentp/pyfasta>`_, a Python module by Brent
  Pedersen for extracting genomic sequences fast 
* `bedGraphToBigWig
  <http://130.91.8.212/GenomeBrowser/goldenPath/help/bigWig.html>`_, a program
  for converting from bedGraph to bigWig

Setting up UTR_SETTINGS
=======================
The file ``UTR_SETTINGS`` contains the necessary parameters for Utail to
run. The file is written in the Microsoft INI format, and read by the Python
`ConfigParser <http://docs.python.org/library/configparser.html>`_ module.
Briefly, the file is divided in to ``[sections]`` which have a series of ``name
= value`` entries. An example::

    [ANNOTATION]
    # The annotation should be in the gencode format
    annotation = /users/rg/jskancke/phdproject/3UTR/gencode5/gencode5_annotation.gtf

See the ``UTR_SETTINGS`` file for explanations of the settings.

Next follow some examples of ``UTR_SETTINGS`` set-ups for different scenarios.

Case1: Run with annotation, Gem-mapped reads, and poly(A)-output
======================================================================
You have mapped the RNA-seq reads using the Gem-mapper to the file
``myreads.gem.tar.gz``, and you have saved this in the directory
``/home/me/mymappedreads/``. You have downloaded the GENCODE annotation that
suits you from ftp://ftp.sanger.ac.uk/pub/gencode/. You now wish to run Utail
to get information about the 3UTRs in your experiment. You set the following
settings in ``UTR_SETTINGS``:

Supply the datasets::

    [DATASETS]
    mydir = /home/me/mymappedreads
    dset1 = %(mydir)s/myreads.gem.tar.gz

Give the path to the annotation you wish to use::

    [ANNOTATION]
    annotation = /home/me/annotations/gencode5_annotation.gtf

Specify that you want to extract poly(A) reads, and provide the path to the
gem-index files, without an extension::

    [POLYA_READS]
    polya = true 
    gem_mapper_index = /home/me/gem_index/H.sapiens.genome.hg19.main

You are using the annotation to get the 3UTR regions, so you are not supplying
a bedfile for this::

    [SUPPLIED_3UTR_BEDFILE]
    utr_bed_path = false


Case2: Run without annotation but with polyA reads
=====================================================
You have mapped the RNA-seq reads with your mapper of choice. You select the
reads that you want, and you save them in a file ``myreads.bed`` (or
``myreads.bed.gz``), in the directory ``/home/me/mymappedreads/``.

Further, you extracted the poly(A) reads from the original mapping
yourself. You have saved these as ``/home/me/mymappedreads/polyareads.bed``

Using your own scripts, you have selected some genomic regions you wish to
investigate (for example some 3UTRs of your choice). You save the regions in
.bed-format in ``/home/me/mybed/genomic_region.bed``.

You now wish to run Utail. You set the following settings in ``UTR_SETTINGS``:

Supply the datasets::

    [DATASETS]
    mydir = /home/me/mymappedreads
    dset1 = %(mydir)s/myreads.bed

Set the annotation path to 'false' since you will be annotation-independent::

    [ANNOTATION]
    annotation = false # false with small 'f'

Give the path to the poly(A) reads you have extracted yourself. Set
gem_mapper_index to false since you have already supplied the polyA reads.::

    [POLYA_READS]
    polya = /home/me/mymappedreads/polyareads.bed
    gem_mapper_index = false

Supply your own bedfile with genomic regions (f.ex some 3UTRs)::

    [SUPPLIED_3UTR_BEDFILE]
    # The path to your genomic regions of interest
    utr_bed_path = /home/me/mybed/genomic_region.bed


Case3: Run without annotation and without polyA reads
========================================================
You have mapped the RNA-seq reads with your mapper of choice. You select the
reads that you want, and you save them in a file ``myreads.bed`` (or
``myreads.bed.gz``), in the directory ``/home/me/mymappedreads/``.

Using your own scripts, you have selected some genomic regions you wish to
investigate (for example some 3UTRs of your choice). You save the regions in
.bed-format in ``/home/me/mybed/genomic_region.bed``.

You now wish to run Utail. You set the following settings in ``UTR_SETTINGS``:

Supply the datasets::

    [DATASETS]
    mydir = /home/me/mymappedreads
    dset1 = %(mydir)s/myreads.bed

Set the annotation path to 'false' since you will be annotation-independent::

    [ANNOTATION]
    annotation = false # false with small 'f'

Specify 'false' for polyA reads and gem_mapper_index::

    [POLYA_READS]
    polya = false
    gem_mapper_index = false

Supply your own bedfile with genomic regions (f.ex some 3UTRs)::

    [SUPPLIED_3UTR_BEDFILE]
    utr_bed_path = /home/me/mybed/genomic_region.bed


Remaining settings
==================

The rest of the settings generally don't change::

    [CHROMOSOME1]
    only_chr1 = false # you wish to run all chromosomes, not just chrm1

    [MIN_3UTR_LENGTH]
    min3utrlen = 200 # default setting

    [CPU_CORES]
    max_cores = default # Use default: max(cpu_count()) - 1 

    [RESTRICT_READS]
    restrict_reads = false # you want to process all reads

    [EXTEND]
    extend_by = 100 # default (recommended setting)

    [HG_FASTA]
    hg_fasta = /home/me/fastafiles/hg19.fa # 3GB file

    [BIGWIG]
    bigwig = all # You want to make bigWig files for all datasets


Running Utail
=============
Go to the directory where the Utail scripts are saved and run::

    $ python Utail.py

Depending on the size of your dataset this takes from 10 minutes to 1-2 hours.
When the script has finished running, the output is found in the folder
``output`` as two files: ``length_dset1`` and ``polya_dset``. See section XXX
for a description of the output parameters.


Length output parameters
========================

===========================  ===================================================
Parameter                     Description  
===========================  ===================================================
chrm                         The chromosome of the 3UTR 
beg                          Start-coordinate of 3UTR
end                          End-coordinate of 3UTR
3utr_extended_by             How far this 3UTR was extended
strand                       The strand (positive or negative)
utr_ID                       Format: Ensembl_ID+internal numbering of 3UTR
epsilon_coord                3UTR end-coordinate as determined by the cumulative
                             coverage
epsilon_rel_size             A value of 0.6 would mean that this 3UTR is 60% of
                             the annotated length
epsilon_downstream_covrg     Average coverage 50nt downstream the epsilon site
epsilon_upstream_covrg       Average coverage 50nt upstream the epsilon site
annotation_distance          Distance in nucleotides from annotated 3UTR end
                             (within +/- 100 nucleotides; if distance is more
                             than that, 'NA' is set)
annotation_downstream_covrg  The same as for epsilon, but for the annotated end
annotation_upstream_covrg    The same as for epsilon, but for the annotated end 
epsilon_PAS_type             Space-separated list of PASe found within 40nt
                             downstream the site. 'NA' if none
epsilon_PAS_distance         Space-separated list of distances to those PAS.
                             'NA if none
epsilon_covrg_beyond_aTTS    If the epsilon_rel_size is > epsilon itself, the
                             program attempts to get the coverage of the
                             downstream region into the extended region. This
                             value can be used to judge if the 3UTR extends
                             beyond the annotated region
3utr_RPKM                    The RPKM of the 3UTR
3utr_average_coverage        The average read coverage of the 3UTR
===========================  ===================================================

PolyA output parameters
========================

========================= =====================================================
Parameter                  Description  
========================= =====================================================
chrm                      The chromosome of the 3UTR
beg                       Start-coordinate of the 3UTR
end                       End-coordinate of the 3UTR
utr_ID                    Format: Ensembl_ID+internal numbering of 3UTR 
polyA_number              Each 3UTR can have more than one polyA cluster.
                          This is the number for this cluster for this 3UTR
strand                    Strand (positive or negative)
polyA_coordinate          The genomic coordinate of the poly(A) cluster
number_supporting_reads   The number of polyA reads supporting this cluster
coverage_50nt_downstream  Average coverage 50nt downstream the cluster
coverage_50nt_upstream    Average coverage 50nt upstream the cluster
annotated_polyA_distance  If annotated polyA end (TTS) nearby, give the
                          distance. If none, write 'NA'
nearby_PAS                Space-separated list of PAS found within 40nt
                          downstream. 'NA' if none
PAS_distance              Space-separated list of distances to those PAS. 'NA'
                          if none
3utr_RPKM                 The RPKM of the 3UTR
========================= =====================================================


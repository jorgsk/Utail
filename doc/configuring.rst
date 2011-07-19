.. _configuration:

========================
Configuring UTR_SETTINGS
========================

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

Case1: Run with annotation, with Gem-mapped reads, and with polyA reads
=========================================================================
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
==================================================
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
=====================================================
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



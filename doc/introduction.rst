.. _introduction:

============
Introduction
============
This document gives the biological background for why this program was made, it
lists the dependencies of the program, and it shows two different ways of using
the program -- with or without an annotation.


On 3UTRs and RNA-seq
====================
The 3' untranslated regions (3'UTRs) of mRNAs are known to contain sequence and
structural motifs that regulate stability, localization and translation of the
mRNAs. By varying the length of 3'UTRs, the cell can thus include or exclude
regulatory regions.

High-throughput RNA-sequencing gives a high-resolution coverage of the
transcriptome. This program uses RNA-seq reads to cover the 3UTRs of genes. It
outputs the length of the 3UTR as well as evidence of polyadenylation events if
possible.

Mandatory dependencies
======================
These are required for the most basic function:

* Python 2.6.4 ++
* `bedTools <http://code.google.com/p/bedtools/>`_
* RNA-seq reads, either in bed-format or in `Gem <http://sourceforge.net/apps/mediawiki/gemlibrary/index.php?title=Gem_mapper_man_page>`_ format


Optional dependencies
======================
For finding polyadenylation events the pipeline needs to remap reads to the
genome with the Gem mapper:

* `Gem mapper <http://sourceforge.net/apps/mediawiki/gemlibrary/index.php?title=Gem_mapper_man_page>`_

To use the mapper you have to create an index file from a fasta file of the
human genome (this takes hours) using::
    gem-do-index --reverse-complement emulate

Then you must provide a path to this index file in UTR_SETTINGS

Included dependencies
====================
The following software packages are shipped with this program for ease of use.
Credits go to their authors.

* `pyFasta <https://github.com/brentp/pyfasta>`_, a Python module by Brent
  Pedersen for extracting genomic sequences fast 
* `bedGraphToBigWig
  <http://130.91.8.212/GenomeBrowser/goldenPath/help/bigWig.html>`_, a program
  for converting from bedGraph to bigWig

Setting up UTR_SETTINGS
=======================
The file ``UTR_SETTINGS`` contains the necessary parameters for the program to
run. The file is written in the Microsoft INI format, and read by the Python
`ConfigParser <http://docs.python.org/library/configparser.html>`_ module.
Briefly, the file is divided in to ``[sections]`` which have a series of ``name
= value`` entries. An example::

    [ANNOTATION]
    # The annotation should be in the gencode format
    annotation = /users/rg/jskancke/phdproject/3UTR/gencode5/gencode5_annotation.gtf

The setting 

Example of usage - with an annotation
=====================================



Example of usage - without annotation
=====================================

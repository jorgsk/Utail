.. _installation:

============
Installation
============

* Download the latest version from github
* Install using python setup.py install
* For help, contact help@helme.con
  
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


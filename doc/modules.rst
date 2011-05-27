.. _modules:

=======
Modules
=======

The central classes and methods of Utail

=================
:mod:`utail`
=================

.. automodule:: utail

Classes
"""""""
.. autoclass:: Settings
    :members:

.. autoclass:: UTR
    :members:

.. autoclass:: FullLength
    :members:

.. autoclass:: PolyAReads
    :members:

Methods
"""""""

.. autofunction:: read_settings
.. autofunction:: get_bed_reads
.. autofunction:: zcat_wrapper
.. autofunction:: get_rpkm
.. autofunction:: output_writer
.. autofunction:: make_bigwigs

========================
:mod:`annotation_parser`
========================

.. automodule:: annotation_parser

Classes
"""""""
.. autoclass:: Transcript
    :members:

Methods
"""""""
.. autofunction:: make_transcripts
.. autofunction:: get_3utr_bed_all_exons
.. autofunction:: get_a_polyA_sites_bed
.. autofunction:: get_seqs


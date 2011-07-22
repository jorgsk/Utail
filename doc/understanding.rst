.. _understanding:

================
The output files
================

Two output files. One concerning the length of the 3UTRs, the other concerning
the polyA-read clusters of the 3UTRs. Here is an explanation of the output
variables.

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


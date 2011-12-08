"""
Read gem-mapped reads and output to bed with

chrm beg end flow_cell seq strand
"""

from sys import stdin
from re import compile as recompile

# Accept up to two mismatches. Make as set for speedup.
acceptable = set(('1:0:0', '0:1:0', '0:0:1'))
getstrand = {'R':'-', 'F':'+'}
start_re = recompile('[0-9]*')

for line in stdin:

    try:
        (ID, seq, quality, mapinfo, position) = line.split('\t')
    except:
        print 'Could not process {0}'.format(line)
        continue

    # Acceptable and poly(A) reads are mutually exclusive.
    if mapinfo in acceptable:
        # Get chromosome
        chrm, rest = position.split(':')
        # Get the strand
        strand = getstrand[rest[0]]
        # Get read beg
        beg = start_re.match(rest[1:]).group()

        flowcell = ID.split(':')[1]

        # Write to file
        print '\t'.join([chrm, beg, str(int(beg)+len(seq)), flowcell, seq, strand])

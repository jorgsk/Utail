"""Create statistics and plots from the tuning file. Will help in determining
the correct cut-off value for this dataset."""

from __future__ import division
import os
import matplotlib.pyplot as plt
from scipy import stats

# only get the debug function if run from Ipython
def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False
if run_from_ipython():
    from IPython.Debugger import Tracer
    debug = Tracer()
else:

    def debug(): pass

def main():
    here = os.path.dirname(os.path.realpath(__file__))
    filenames = ['cumul_k562_whole_cell.stat']
    for filename in filenames:
        tuning_file = os.path.join(os.path.split(here)[0], 'output', filename)
        #tuning_file = os.path.join(os.path.split(here)[0], 'output', 'test.stat')
        contents = {0:'cumul', 1:'pA_to_cumul_dist', 2:'pA_cumul', 3: 'd_stream_covr',
                    4: 'u_stream_covr', 5: 'rpkm', 6:'utr_length', 7: 'strand'}

        tuning_handle = open(tuning_file, 'rb')
        header = tuning_handle.next().split()

        # Get the stats dictionary
        data = {}

        for line in tuning_handle:
            (utr_id, default_cumul, pA_to_cumul_dist, pA_cumul, d_stream_covr,
             u_stream_covr, rpkm, utr_length, strand) = line.split()

            if float(rpkm) < 2:
                continue

            data[utr_id] = (int(default_cumul), int(pA_to_cumul_dist),
                            float(pA_cumul), float(d_stream_covr),
                            float(u_stream_covr), float(rpkm), int(utr_length),
                            strand)

        # Print the distribution of pA_cumul
        # Print the mean and the standard deviation as well
        # TODO wait for the calculation to finish. Then look at the mean and std and
        # maybe plot as well. For now, I move on!
        # AS WELL! Get a measure on how good your changes are: get the mean of the
        # distances from the cut-off to the actual polyas. This distance should
        # decrease with each iteration.

        # The relative cumulative length of the pA clusters
        pA_cumuls = [vals[2] for vals in data.itervalues()]
        (n_cumul, min_max_cumul, mean_cumul, var_cumul) = stats.describe(pA_cumuls)[:4]

        # The before/after coverage ratio of the pA clusters
        beg_aft = [vals[3]/vals[4] for vals in data.itervalues() if vals[4]!=0]
        (n_ratio, min_max_ratio, mean_ratio, var_ratio) = stats.describe(pA_cumuls)[:4]


        #box_plot(pA_cumuls)
    debug()

def box_plot(values):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.boxplot(values)
    ax.set_ylabel('Cumulative coverages of poly(A) clusters', size=20)
    ax.set_ylim(0,1.2)
    fig.show()
    pass

if __name__ == '__main__':
    main()

import os
import numpy as np
from scipy.signal import argrelextrema    #, argrelmax


def read2list(datafile):
    l = []
    with open(datafile) as f:
        while True:
            line = f.readline()
            if len(line) == 0:
                break
            l.append(line.split())
    return l


datasets = {}    # A dictionary, compound: [(x1, y1), (x2, y2), ...], x/y being exprFreq/deviationOfCalc
l = read2list('zzz-LIST')
for i in range(0, len(l)):
    expr_f = l[i][0]
    calc_f = l[i][1]
    if os.path.isfile(expr_f) and os.path.isfile(calc_f):
        compound = expr_f.split('-')[0]
        datasets[compound] = []
        expr = np.loadtxt(expr_f, usecols=range(0,2), dtype=np.float32)
        calc = np.loadtxt(calc_f, usecols=range(0,2), dtype=np.float32)
        comm_range_start = np.maximum(expr[0][0], calc[0][0])    # for common freq range
        comm_range_end = np.minimum(expr[-1][0], calc[-1][0])    # for common freq range
        # curtail to common freq range
        expr_comm = np.delete(expr, np.where(
            (expr[:,0] < comm_range_start) | (expr[:,0] > comm_range_end))[0], axis=0)
        calc_comm = np.delete(calc, np.where(
            (calc[:,0] < comm_range_start) | (calc[:,0] > comm_range_end))[0], axis=0)
        expr_max_i = argrelextrema(expr_comm[:,1], np.greater)[0]    # array of indices of local maxima
        calc_max_i = argrelextrema(calc_comm[:,1], np.greater)[0]    # array of indices of local maxima
        expr_max_freq = expr_comm[:,0][expr_max_i]    # freq of all peaks (local maxima)
        calc_max_freq = calc_comm[:,0][calc_max_i]    # freq of all peaks (local maxima)
        for freq in calc_max_freq:
            idx_match = (np.abs(expr_max_freq - freq)).argmin()
            match = expr_max_freq[idx_match]    # the expr freq that's found to match the calc freq
            datapoint = match, freq-match
            datasets[compound].append(datapoint)

print(datasets)

datasets_plain = []    # [(x1, y1), (x2, y2), ...] containing all data points of all compounds
for comp in datasets.keys():
    datasets_plain.extend(datasets[comp])

filename = "plot_all_data_points.dat"
np.savetxt(filename, np.array(datasets_plain))


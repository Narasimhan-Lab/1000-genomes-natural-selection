#!/usr/bin/env python3
""" Scan Result Compiler

Compiles the results of monogenic_scan.py into a single tab-delimited .txt file given monogenic_config.yml. The output file will saved in /config/out/path/{target}_scan_results.txt

Command-line Arguments:
    [1] num_chunks: number of segments the genome was divided into
    [2] chunk_size: number of positions in each segment
"""
import sys
import math

import yaml
import numpy as np
from scipy.stats import ncx2
import statistics

num_chunks = int(sys.argv[1])
chunk_size = int(sys.argv[2])

config = yaml.safe_load(open("monogenic_config.yml"))

epoch = config['populations']['target']
source1 = config['populations']['source1']
source2 = config['populations']['source2']

scan_prefix = config['paths']['out'] + "/"
out_path = scan_prefix + epoch + "_scan_results.txt"
scan_prefix += "result_arrs/"
pos_path = config['paths']['chrom_pos'] 

#########################################################
#########################################################

start_indices = np.arange(1, num_chunks) * chunk_size
p_vals = np.load(scan_prefix + "p_vals0.npy")
stats = np.load(scan_prefix + "stats0.npy")
t_freqs = np.load(scan_prefix + "targetFreq0.npy")
s1_freqs = np.load(scan_prefix + "source1Freq0.npy")
s2_freqs = np.load(scan_prefix + "source2Freq0.npy")
expected = np.load(scan_prefix + "expected0.npy")
expected_h0 = np.load(scan_prefix + "expectedH00.npy")

s1_n = np.load(scan_prefix + "source1SampleSize0.npy")
s2_n = np.load(scan_prefix + "source2SampleSize0.npy")
t_n = np.load(scan_prefix + "targetSampleSize0.npy")

s1_h0s = np.load(scan_prefix + "s1H00.npy")
s2_h0s = np.load(scan_prefix + "s2H00.npy")
h1_neg_logs = np.load(scan_prefix + "h1NegLog0.npy")
h0_neg_logs = np.load(scan_prefix + "h0NegLog0.npy")

pert_p = np.load(scan_prefix + "pertP0.npy")
pert_stats = np.load(scan_prefix + "pertStat0.npy")

total_sample_sizes = [len(np.load(config['paths']['source1_indices'])), len(np.load(config['paths']['source2_indices'])), len(np.load(config['paths']['target_indices']))]

file_names = ["source1Freq", "source2Freq", "targetFreq", "expected", "s1H0", "s2H0", "expectedH0", "h1NegLog", "h0NegLog", "stats", "p_vals", "pertStat", "pertP", "source1SampleSize", "source2SampleSize", "targetSampleSize"]
arr_names = [s1_freqs, s2_freqs, t_freqs, expected, s1_h0s, s2_h0s, expected_h0, h1_neg_logs, h0_neg_logs, stats, p_vals, pert_stats, pert_p, s1_n, s2_n, t_n]

# load all result arrays
i = 0
while i < num_chunks - 1:
    for j in range(len(file_names)):
        file_name = scan_prefix + file_names[j] + str(start_indices[i]) + ".npy"
        arr_names[j] = np.append(arr_names[j], np.load(file_name))
    i += 1

# calculate non-centrality parameter
lam = np.nanmedian(pert_stats)

chros_file = open(pos_path, 'r')

out_file = open(out_path, 'w')
out_file.write("CHROM\tPOSITION\tREF\tALT\tSOURCE1_FREQ\tSOURCE2_FREQ\tTARGET_FREQ\tEXPECTED\tS1_H0\tS2_H0\tEXPECTED_H0\tH1_NEG_LOG\tH0_NEG_LOG\tSTATISTIC\tP_VALUE\tPERT_STAT\tPERT_P\tNCX2_P\tS1_N\tS2_N\tT_N\n")

j = 0
while j < len(arr_names[0]):
    chrom_line = chros_file.readline()
    chrom_split = chrom_line.split()

    out_file.write(chrom_line.rstrip() + '\t')
    for k in range(13):
        out_file.write(str(arr_names[k][j]) + '\t')

    # writing the corrected perturbed statistic
    out_file.write(str(ncx2.sf(arr_names[11][j], 1, lam)))

    # write sample sizes
    for k in range(13, 16):
        out_file.write('\t' + str(arr_names[k][j]))

    out_file.write('\n')
    j += 1

chros_file.close()
out_file.close()

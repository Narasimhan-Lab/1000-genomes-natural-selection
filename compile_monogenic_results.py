#!/usr/bin/env python3
import numpy as np
import math
from scipy.stats.distributions import chi2
import statistics
import sys

#########################################################
#########################################################

num_chunks = int(sys.argv[1])
chunk_size = int(sys.argv[2])

scan_prefix = "/work2/07754/meganle/lonestar/run_v4.5/scan/EN/"
out_path = scan_prefix + "EN_scan_results.txt"
pos_path = "/work2/07754/meganle/lonestar/admixture_scans/data/alleleChrPos.txt"

#########################################################
#########################################################

start_indices = np.arange(1, num_chunks) * chunk_size
p_vals = np.load(scan_prefix + "p_vals0.npy")
stats = np.load(scan_prefix + "stats0.npy")
t_freqs = np.load(scan_prefix + "targetFreq0.npy")
s1_freqs = np.load(scan_prefix + "source1Freq0.npy")
s2_freqs = np.load(scan_prefix + "source2Freq0.npy")
expected = np.load(scan_prefix + "expected0.npy")

i = 0
while i < num_chunks- 1:
	file_name = scan_prefix + "p_vals" + str(start_indices[i]) + ".npy"
	p_vals = np.append(p_vals, np.load(file_name))

	file_name = scan_prefix + "stats" + str(start_indices[i]) + ".npy"
	stats = np.append(stats, np.load(file_name))

	file_name = scan_prefix + "targetFreq" + str(start_indices[i]) + ".npy"
	t_freqs = np.append(t_freqs, np.load(file_name))

	file_name = scan_prefix + "source1Freq" + str(start_indices[i]) + ".npy"
	s1_freqs = np.append(s1_freqs, np.load(file_name))

	file_name = scan_prefix + "source2Freq" + str(start_indices[i]) + ".npy"
	s2_freqs = np.append(s2_freqs, np.load(file_name))

	file_name = scan_prefix + "expected" + str(start_indices[i]) + ".npy"
	expected = np.append(expected, np.load(file_name))

	i += 1

chros_file = open(pos_path, 'r')
chros_file.readline()

out_file = open(out_path, 'w')
out_file.write("CHROM\tPOSITION\tREF\tALT\tSOURCE1_FREQ\tSOURCE2_FREQ\tTARGET_FREQ\tEXPECTED\tSTATISTIC\tP_VALUE\n")

j = 0
while j < len(p_vals):
	chrom_line = chros_file.readline()
	out = chrom_line.rstrip() + "\t" + str(s1_freqs[j]) + "\t" + str(s2_freqs[j]) + "\t" + str(t_freqs[j]) + "\t" + str(expected[j]) + "\t" + str(stats[j]) + "\t" + str(p_vals[j]) + "\n"
	out_file.write(out)
	
	j += 1

chros_file.close()
out_file.close()

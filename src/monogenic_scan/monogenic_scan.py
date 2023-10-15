#!/usr/bin/env python3 
""" Monogenic Scan

Performs a genome-wide scan for SNPs whose allele frequencies deviate from expectation based on genome-wide mixing proportions.

Configure paths and scan arguments in monogenic_config.yml. The output can be summarized in a .txt file using compile_monogenic_results.py

Command-line Arguments:
    [1] p: relative position of the SNP where computation should start
    [2] n: number of positions to perform calculations on
"""
import math
import sys
import os

import yaml
import tables
import numpy as np
from scipy.stats import binom
from scipy.optimize import minimize
from scipy.stats.distributions import chi2

# Load configurations
config = yaml.safe_load(open("monogenic_config.yml"))

# population info here
# e.g. target = "C", source1 = "A", source2 = "B" for the model C = A + B
target = config['populations']['target']
source1 = config['populations']['source1']
source2 = config['populations']['source2']

populations = [source1, source2, target]

# qpAdm results
alph = config['qpAdm_results']['alpha']
pop_proportions = [alph, 1 - alph]
multiplier = config['qpAdm_results']['multiplier']
prop_se = config['qpAdm_results']['standard_error']

assert alph > 0 and alph < 1

# amount to perturb mixing proportions by
eps_q = multiplier * prop_se
perturbed_props = [[pop_proportions[0]+eps_q, pop_proportions[1]-eps_q], [pop_proportions[0]-eps_q, pop_proportions[1]+eps_q]]

assert (np.sum(perturbed_props[0]) == 1)
assert (np.sum(perturbed_props[1]) == 1)

mat_prefix = "/work2/07754/meganle/lonestar/run_v19/data/"
data_prefix = "/work2/07754/meganle/lonestar/run_v21/data/"

# paths of hdf5 files with read count data
ref_path = config['paths']['ref_reads']
alt_path = config['paths']['alt_reads']

# indices of samples in hdf5 file for each population
target_npy = config['paths']['target_indices']
s1_npy = config['paths']['source1_indices']
s2_npy = config['paths']['source2_indices']
ids_list = [s1_npy, s2_npy, target_npy]

##########################################################
##########################################################

# these arguments are used to split the scan up into multiple chunks for parallel computation
# first user argument, specifies which relative position to start at (out of 1150477 positions in our case)
start_index = int(sys.argv[1])
# second user argument, number of sequential positions to compute (starting at start_index)
num_positions = int(sys.argv[2])
pos_index = 0

sequencing_error = 0.001
sample_indices = []
ref_tables = []
alt_tables = []
admixed_idx = len(populations) - 1

# calculates negative log likelihood for a single population
# theta is an array with a single entry that contains a frequency p (needs to be in this format to use scipy.minimize)
def neg_log_likelihood_single(theta, pop_index):
    p = theta[0]
    total_counts = ref_tables[pop_index][pos_index] + alt_tables[pop_index][pos_index]
    p2 = p * p * binom.pmf(ref_tables[pop_index][pos_index], total_counts, 1-sequencing_error)
    pq = 2 * p * (1-p) * binom.pmf(ref_tables[pop_index][pos_index], total_counts, 0.5)
    q2 = (1-p) * (1-p) * binom.pmf(ref_tables[pop_index][pos_index], total_counts, sequencing_error)
    res = np.add(p2, pq)
    res = np.add(q2, res)
    res = vfunc(res)

    ans = np.sum(res)
    return ans

# theta is an array with the frequencies for both source populations
def neg_log_likelihood_null(theta, alpha):
    ans = 0
    curr_pop = 0
    for p in [theta[0], theta[1], compute_expected(alpha, theta[0], theta[1])]:
        total_counts = ref_tables[curr_pop][pos_index] + alt_tables[curr_pop][pos_index]
        p2 = p * p * binom.pmf(ref_tables[curr_pop][pos_index], total_counts, 1-sequencing_error)
        pq = 2 * p * (1-p) * binom.pmf(ref_tables[curr_pop][pos_index], total_counts, 0.5)
        q2 = (1-p) * (1-p) * binom.pmf(ref_tables[curr_pop][pos_index], total_counts, sequencing_error)

        res = np.add(p2, pq)
        res = np.add(q2, res)
        res = vfunc(res)

        ans += np.sum(res)
        curr_pop += 1
    
    return ans

    
def neg_log(x):
    if (x == 0):
        return 0
    return -1 * math.log(x)

# finds MLE
def optimize_single(pop_index):
    ref_counts = np.sum(ref_tables[pop_index][pos_index])
    alt_counts = np.sum(alt_tables[pop_index][pos_index])
    if(ref_counts + alt_counts == 0):
        # reference and alternate read counts are both 0, so mark NA
        return np.nan, 0
    # initial guess is fraction of reference read counts
    theta = np.array([ref_counts/(ref_counts + alt_counts)])
    # minimize the negative log likelihood
    res = minimize(neg_log_likelihood_single, theta, args = (pop_index), method = 'SLSQP', bounds=((0.001, 0.999),), options={'disp':False})
    # number of samples with nonzero read counts
    n = 0
    for j in range(len(ref_tables[pop_index][pos_index])):
        if ref_tables[pop_index][pos_index][j] + alt_tables[pop_index][pos_index][j] > 0:
            n += 1
    return res.x[0], n

def optimize_null(alpha):
    # return nan if any population has 0 read counts
    for pop_index in range(3):
        if(np.sum(ref_tables[pop_index][pos_index] + alt_tables[pop_index][pos_index]) == 0):
            return np.nan, np.nan

    s1_ref = np.sum(ref_tables[0][pos_index])
    s1_alt = np.sum(alt_tables[0][pos_index])

    s2_ref = np.sum(ref_tables[1][pos_index])
    s2_alt = np.sum(alt_tables[1][pos_index])

    # initial guess is fraction of reference read counts
    theta = np.array([s1_ref/(s1_ref + s1_alt), s2_ref/(s2_ref + s2_alt)])

    res = minimize(neg_log_likelihood_null, theta, args = (alpha), method = 'SLSQP', bounds=((0.001, 0.999), (0.001, 0.999),), options={'disp':False})
    return res.x[0], res.x[1]

def compute_expected(alpha, source_a, source_b):
    assert alpha < 1 and alpha > 0
    res =  alpha*source_a + (1-alpha)*source_b
    assert (res < 1 and res > 0) or np.isnan(res)
    return res

def likelihood_ratio_test(mixing_freq, estimates):
    observed = estimates[admixed_idx]
    expected = compute_expected(mixing_freq[0], estimates[0], estimates[1])

    try: 
        stat_res =  2 * (neg_log_likelihood_single([expected], admixed_idx) - neg_log_likelihood_single([observed], admixed_idx))
    except ValueError:
        return np.nan, np.nan, np.nan

    return p_value, stat_res, expected

# these arrays will store the results of the scan
p_values = np.zeros(num_positions)
statistics = np.zeros(num_positions)
target_af = np.zeros(num_positions)
source1_af = np.zeros(num_positions)
source2_af = np.zeros(num_positions)
expected_af = np.zeros(num_positions)
expected_h0 = np.zeros(num_positions)
source1_sample_sizes = np.zeros(num_positions)
source2_sample_sizes = np.zeros(num_positions)
target_sample_sizes = np.zeros(num_positions)

s1_h0s = np.zeros(num_positions)
s2_h0s = np.zeros(num_positions)
h1_neg_logs = np.zeros(num_positions)
h0_neg_logs = np.zeros(num_positions)

# save the highest results after mixing proportion perturbation
pert_p = np.zeros(num_positions)
pert_stat = np.zeros(num_positions)

# create a list of arrays of read count data indices for each population
for i in range(len(ids_list)):
    indices = np.load(ids_list[i])
    sample_indices.append(indices)  

# create a list of reference read count data by extracting the relevant columns for each population
with tables.File(ref_path, 'r') as ref_file:
    r_table = ref_file.get_node('/dataset0').read(start=start_index, stop=num_positions+start_index)
    for i in range(len(sample_indices)):
        ref_tables.append(r_table[:, sample_indices[i]].astype(np.int32))

# create a list of alternate read count data by extracting the relevant columns for each population
with tables.File(alt_path, 'r') as alt_file:
    a_table = alt_file.get_node('/dataset0').read(start=start_index, stop=num_positions+start_index)
    for i in range(len(sample_indices)):
        alt_tables.append(a_table[:, sample_indices[i]].astype(np.int32))

vfunc = np.vectorize(neg_log)

# iterate over all the populations for this chunk
while(pos_index < num_positions):
    # calculate MLE of observed allele frequencies
    pos_sample_sizes = []
    frequencies = []
    for i in range(len(populations)):
        # minimize negative log-likelihood
        result, n = optimize_single(i)
        frequencies.append(result)
        pos_sample_sizes.append(n)

    source1_sample_sizes[pos_index] = pos_sample_sizes[0]
    source2_sample_sizes[pos_index] = pos_sample_sizes[1]
    target_sample_sizes[pos_index] = pos_sample_sizes[2]

    # add together neg log likelihoods for each population's MLE frequency
    h1_neg_log = 0
    for pop_index in range(len(frequencies)):
        h1_neg_log += neg_log_likelihood_single([frequencies[pop_index]], pop_index)
    
    # optimize likelihood under neutrality
    s1_h0, s2_h0 = optimize_null(pop_proportions[0])
    h0_neg_log = neg_log_likelihood_null([s1_h0, s2_h0], pop_proportions[0])

    try:
        max_stat = 2 * (h0_neg_log - h1_neg_log)
        max_p = chi2.sf(max_stat, 1)
    except:
        max_stat = np.nan
        max_p = np.nan

    s1_h0s[pos_index] = s1_h0
    s2_h0s[pos_index] = s2_h0
    h1_neg_logs[pos_index] = h1_neg_log
    h0_neg_logs[pos_index] = h0_neg_log

    source1_af[pos_index] = frequencies[0]
    source2_af[pos_index] = frequencies[1]
    target_af[pos_index] = frequencies[2]
    expected_af[pos_index] = compute_expected(pop_proportions[0], frequencies[0], frequencies[1])

    expected_h0[pos_index] = compute_expected(pop_proportions[0], s1_h0, s2_h0)

    p_values[pos_index] = max_p
    statistics[pos_index] = max_stat

    # save the results in our arrays
    if not np.isnan(max_p):
        for j in range(len(perturbed_props)):
            temp_props = perturbed_props[j]
            temp_s1_h0, temp_s2_h0 = optimize_null(temp_props[0])
            temp_h0_neg_log = neg_log_likelihood_null([temp_s1_h0, temp_s2_h0], temp_props[0])

            try:
                temp_stat = 2 * (temp_h0_neg_log - h1_neg_log)
                temp_p = chi2.sf(temp_stat, 1)
                if temp_p > max_p:
                    max_p = temp_p
                    max_stat = temp_stat
            except:
                print("Error in perturbation test")

    pert_p[pos_index] = max_p
    pert_stat[pos_index] = max_stat

    print(max_p, frequencies)

    pos_index += 1

ref_file.close()
alt_file.close()

# save all of the results arrays to output files
out_path = config['paths']['out'] + target + '/result_arrs/'
if not os.path.exists(out_path):
    os.makedirs(out_path)

np.save(out_path + "p_vals" + str(start_index), p_values)
np.save(out_path + "stats" + str(start_index), statistics)
np.save(out_path + "targetFreq" + str(start_index), target_af)
np.save(out_path + "source1Freq" + str(start_index), source1_af)
np.save(out_path + "source2Freq" + str(start_index), source2_af)
np.save(out_path + "targetSampleSize" + str(start_index), target_sample_sizes)
np.save(out_path + "source1SampleSize" + str(start_index), source1_sample_sizes)
np.save(out_path + "source2SampleSize" + str(start_index), source2_sample_sizes)
np.save(out_path + "expected" + str(start_index), expected_af)
np.save(out_path + "expectedH0" + str(start_index), expected_h0)
np.save(out_path + "s1H0" + str(start_index), s1_h0s)
np.save(out_path + "s2H0" + str(start_index), s2_h0s)
np.save(out_path + "h1NegLog" + str(start_index), h1_neg_logs)
np.save(out_path + "h0NegLog" + str(start_index), h0_neg_logs)
np.save(out_path + "pertP" + str(start_index), pert_p)
np.save(out_path + "pertStat" + str(start_index), pert_stat)

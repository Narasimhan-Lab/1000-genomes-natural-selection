#!/usr/bin/env python3 
import tables, math, time
import sys
import numpy as np
from scipy.stats import binom
from scipy.optimize import minimize
from scipy.stats.distributions import chi2

##########################################################
# population info here
# e.g. target = "C", source1 = "A", source2 = "B" for the model C = A + B
target = "EN"
source1 = "M"
source2 = "AN"

populations = [source1, source2, target]

# proportions of source 1 then source 2
# e.g. pop_proportions = [0.6, 0.4] for C = (0.6)A + (0.4)B
pop_proportions = [0.163, 0.837]

assert (np.sum(pop_proportions) == 1)

##########################################################

# paths of hdf5 files with read count data
ref_path = "/work2/07754/meganle/lonestar/run_v4.4/data/v4_v9_combined_BA_refAlleles.hdf5"
alt_path = "/work2/07754/meganle/lonestar/run_v4.4/data/v4_v9_combined_BA_altAlleles.hdf5"

# indices of samples in hdf5 file for each population
t_npy = "../../data/" + target + "_indices.npy"
s1_npy = "../../data/" + source1 + "_indices.npy"
s2_npy = "../../data/" + source2 + "_indices.npy"
ids_list = [s1_npy, s2_npy, t_npy]

# these arguments are used to split the scan up into multiple chunks for parallel computation
# first user argument, specifies which relative position to start at (out of 1150477 positions in our case)
start_index = int(sys.argv[1])
# second user argument, number of sequential positions to compute (starting at start_index)
num_positions = int(sys.argv[2])
pos_index = 0

# these arrays will store the results of the scan
p_values = np.zeros(num_positions)
statistics = np.zeros(num_positions)
target_af = np.zeros(num_positions)
source1_af = np.zeros(num_positions)
source2_af = np.zeros(num_positions)
expected_af = np.zeros(num_positions)

error = 0.001
sample_indices = []
ref_tables = []
alt_tables = []


# calculates negative log likelihood for a single population
# theta is an array with a single entry that contains a frequency p (needs to be in this format to use scipy.minimize)
def neg_log_likelihood(theta, pop_index):
	p = theta[0]
	total_counts = ref_tables[pop_index][pos_index] + alt_tables[pop_index][pos_index]
	p2 = p * p * binom.pmf(ref_tables[pop_index][pos_index], total_counts, 1-error)
	pq = 2 * p * (1-p) * binom.pmf(ref_tables[pop_index][pos_index], total_counts, 0.5)
	q2 = (1-p) * (1-p) * binom.pmf(ref_tables[pop_index][pos_index], total_counts, error)
	res = np.add(p2, pq)
	res = np.add(q2, res)
	vfunc = np.vectorize(func)
	res = vfunc(res)
	ans = np.sum(res)
	return ans
	
def func(x):
	if (x == 0):
		return 0
	return -1 * math.log(x)

# finds MLE
def optimize_likelihood(pop_index):
	if(np.sum(ref_tables[pop_index][pos_index] + alt_tables[pop_index][pos_index]) == 0):
		# reference and alternate read counts are both 0, so mark NA
		return np.nan
	theta = np.array([0.5])
	# minimize the negative log likelihood
	res = minimize(neg_log_likelihood, theta, args = (pop_index), method = 'SLSQP', bounds=((0.01, 0.99),), options={'disp':False})
	return res.x[0]

def likelihood_ratio_test(mixing_freq, estimates):
	expected = 0
	# computes expected frequency using mixture proportions
	for i in range(len(mixing_freq)):
		expected += mixing_freq[i]  * estimates[i]
	admixed_idx = len(populations) - 1
	# observed is the MLE of the target population
	observed = estimates[admixed_idx]
	statistic = 0
	try:
                # likelihood ratio test
		statistic = 2 * (neg_log_likelihood([expected], admixed_idx) - neg_log_likelihood([observed], admixed_idx))
	except ValueError:
		# experiened NAN because read counts were 0 for some population
		return np.nan, np.nan, np.nan
	p_value = chi2.sf(statistic, 1)
	return p_value, statistic, expected

print("refAlleles location:", ref_path[0: ref_path.rindex("/")])
print("altAlleles location:", alt_path[0: alt_path.rindex("/")])

# create a list of arrays of read count data indices for each population
for i in range(len(ids_list)):
	indices = np.load(ids_list[i])
	sample_indices.append(indices)	

# create a list of reference read count data by extracting the relevant columns for each population
with tables.File(ref_path, 'r') as ref_file:
	r_table = ref_file.get_node('/dataset0').read(start=start_index, stop=num_positions+start_index)
	for i in range(len(sample_indices)):
		ref_tables.append(r_table[:, sample_indices[i]])

# create a list of alternate read count data by extracting the relevant columns for each population
with tables.File(alt_path, 'r') as alt_file:
	a_table = alt_file.get_node('/dataset0').read(start=start_index, stop=num_positions+start_index)
	for i in range(len(sample_indices)):
		alt_tables.append(a_table[:, sample_indices[i]])

# iterate over all the populations for this chunk
while(pos_index < num_positions):
	# calculate frequencies
	frequencies = []
	for i in range(len(populations)):
		# minimize negative log-likelihood
		result = optimize_likelihood(i)
		frequencies.append(result)

	lrt_result, stat, exp = likelihood_ratio_test(pop_proportions, frequencies)
	# save the results in our arrays
	p_values[pos_index] = lrt_result
	statistics[pos_index] = stat
	source1_af[pos_index] = frequencies[0]
	source2_af[pos_index] = frequencies[1]
	target_af[pos_index] = frequencies[2]
	expected_af[pos_index] = exp

	pos_index += 1


ref_file.close()
alt_file.close()
# save all of the results arrays to output files
np.save("p_vals" + str(start_index), p_values)
np.save("stats" + str(start_index), statistics)
np.save("targetFreq" + str(start_index), target_af)
np.save("source1Freq" + str(start_index), source1_af)
np.save("source2Freq" + str(start_index), source2_af)
np.save("expected" + str(start_index), expected_af)


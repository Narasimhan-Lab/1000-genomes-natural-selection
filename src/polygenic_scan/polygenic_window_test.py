#!/usr/bin/env python3
""" Polygenic Window Test

Performs a scan for evidence of polygenic selection on a given GWAS trait using the monogenic selection scan statistics.

Configure paths and options polygenic_config.yml.

Command-line Arguments:
    [1] TRAIT: GWAS trait ID (e.g., A02B, A10, AA...)
    [2] P_EXP: p-value cutoff exponent for choosing significant GWAS SNPs (e.g., a `P_EXP` of 6 corresponds to using a cutoff of 1e-6)
    [3] DATASET: indicates whether the GWAS summary statistics are from Biobank Japan or UK Biobank
"""
import sys
import os
import math
import random
import copy

import yaml
import numpy as np
import glob

config = yaml.safe_load(open("polygenic_config.yml"))
mono_config = yaml.safe_load(open("../monogenic_scan/monogenic_config.yml"))

EPOCH = mono_config['populations']['target']
USE_EFFECT_SIZE = config['options']['use_effect_size']

# first argument: trait ID
TRAIT = sys.argv[1]

# second argument: p-value threshold when choosing associated SNPs
P_EXP = int(sys.argv[2])
P_CUTOFF = np.float_power(10, -1 * P_EXP)

# third argument: which biobank/dataset to use
DATASET = sys.argv[3]

if DATASET == "UKB":
    gwas_path = config['gwas_files']['path'] + DATASET + "/hum0197.v3.EUR." + TRAIT + ".v1/*.auto.txt"
    gwas_file = open(glob.glob(gwas_path)[0], 'r')
elif DATASET == "BBJ":
    gwas_path = config['paths']['gwas_path'] + DATASET + "/hum0197.v3.BBJ." + TRAIT + ".v1/*.auto.txt"
    gwas_file = open(glob.glob(gwas_path)[0], 'r')
else:
    assert 0

# quantiles computed from recombination rate map
RR_CUTOFFS = [1.42309722e-19, 1.33096160e-06, 7.99490360e-04, 1.00316488e-02, 5.29890833e-02, 1.98940513e-01, 8.59677415e-01, 9999999999]

# use global file markers/lines so we only have to pass through the files once
global gwas_line
global rr_line
global b_line

# parse GWAS file header
gwas_line = gwas_file.readline()
header = gwas_line.split()
gwas_indices = {'chr': header.index(config['gwas_files']['chr']),
                'pos': header.index(config['gwas_files']['pos']),
                'p': header.index(config['gwas_files']['p']),
                'allele1': header.index(config['gwas_files']['allele1']),
                'allele2': header.index(config['gwas_files']['allele2']),
                'beta': header.index(config['gwas_files']['beta']),
                'freq': header.index(config['gwas_files']['freq']),
                'freq_a': header.index(config['gwas_files']['freq_a'])
            }
gwas_line = gwas_file.readline()

# parse monogenic selection scan header
scan_file = open(config['selection_file']['path'], 'r')
scan_line = scan_file.readline()
header = scan_line.split()
scan_indices = {'chr': header.index(config['selection_file']['chr']),
                'pos': header.index(config['selection_file']['pos']),
                'ref': header.index(config['selection_file']['ref']),
                'alt': header.index(config['selection_file']['alt']),
                't_freq': header.index(config['selection_file']['target_freq']),
                'exp': header.index(config['selection_file']['expected_target_freq']),
                'stat': header.index(config['selection_file']['statistic'])
                }

scan_line = scan_file.readline()

# parse recombination rate map header
rr_file = open(config['rr_file']['path'], 'r')
rr_line = rr_file.readline()
header = rr_line.split()
rr_indices = {'chr': header.index(config['rr_file']['chr']),
            'start': header.index(config['rr_file']['start']),
            'end': header.index(config['rr_file']['end']),
            'r_rate': header.index(config['rr_file']['r_rate'])
            }
rr_line = rr_file.readline()

# parse b file header 
b_file = open(config['b_file']['path'], 'r')
b_line = b_file.readline()
header = b_line.split()
b_indices = {'chr': header.index(config['b_file']['chr']),
            'pos': header.index(config['b_file']['pos']),
            'bd': header.index(config['b_file']['bd']), 
            'anc': header.index(config['b_file']['anc']),
            'der': header.index(config['b_file']['der'])
            }
b_line = b_file.readline()

NUM_BINS = 8
NUM_TRIALS = config['options']['num_trials']
WINDOW_SIZE = config['options']['window_size']

print("Population:", EPOCH)
print("GWAS Path:", gwas_path)
print("Threshold:", P_CUTOFF)
print("Number of trials:", NUM_TRIALS)
print("Number of bins:", NUM_BINS)

def get_gwas_line(find_chr, find_loc):
    """ Parses GWAS file line for a given loci
    Args:
        find_chr: chromosome of location to get information for
        find_loc: position of location to get information for
    Returns:
        dictionary containing GWAS p-value, GWAS beta, ref allele, alt allele, allele frequency
    """
    global gwas_line
    while(gwas_line):
        g_split = gwas_line.split()
        g_chr = g_split[gwas_indices['chr']]
        g_loc = int(float(g_split[gwas_indices['pos']]))

        if g_chr == "X" or g_chr == "Y":
            return np.nan

        # check for match
        if g_chr == find_chr and find_loc == g_loc:
            try:
                test = float(g_split[gwas_indices['freq']])
            except:
                gwas_line = gwas_file.readline()
                return np.nan
            return {'p_val': float(g_split[gwas_indices['p']]), 
                    'ref': g_split[gwas_indices['allele1']], 
                    'alt': g_split[gwas_indices['allele2']], 
                    'beta': float(g_split[gwas_indices['beta']]), 
                    'freq': float(g_split[gwas_indices['freq']]), 
                    'freq_allele': g_split[gwas_indices['freq_a']]
                    }

        # check if we've gone too far
        if int(g_chr) > int(find_chr) or (g_chr == find_chr and g_loc > find_loc):
            return np.nan

        # keep searching 
        gwas_line = gwas_file.readline()
    
    # already reached end of file
    return np.nan

def get_recombination_rate(find_chr, find_loc):
    """ Finds recombination rate for a given location
    Args:
        find_chr: chromosome of location to get recombination rate for
        find_loc: position of location to get recombination rate for
    Returns:
        recombination rate
    """
    global rr_line
    while(rr_line):
        rr_split = rr_line.split()
        rr_chr = rr_split[rr_indices['chr']]
        rr_begin = int(rr_split[rr_indices['start']])
        rr_end = int(rr_split[rr_indices['end']])

        if rr_chr == "X":
            return np.nan

        if rr_chr == find_chr and find_loc >= rr_begin and find_loc < rr_end:
            return float(rr_split[rr_indices['r_rate']])

        if int(rr_chr) > int(find_chr) or (rr_chr == find_chr and rr_begin > find_loc):
            return np.nan

        rr_line = rr_file.readline()
    return np.nan

def get_bdecile(find_chr, find_loc):
    """ Finds B statistic decile for a given location
    Args:
        find_chr: chromosome of location to get B statistic for
        find_loc: position of location to get B statistic for
    Returns:
        B statistic decile
    """
    global b_line
    while(b_line):
        b_split = b_line.split()
        b_chr = b_split[b_indices['chr']]
        b_loc = int(b_split[b_indices['pos']])

        # reached the end of the relevant lines
        if b_chr == "23":
            return np.nan

        # check for match
        if b_chr == find_chr and find_loc == b_loc:
            return {'decile': int(b_split[b_indices['bd']]), 
                    'anc': b_split[b_indices['anc']], 
                    'der': b_split[b_indices['der']]}

        # check if we've gone too far
        if int(b_chr) > int(find_chr) or (b_chr == find_chr and b_loc > find_loc):
            return np.nan

        # keep searching
        b_line = b_file.readline()
    return np.nan

# initialize data structures and window variables
LOWEST_P = np.nan
LOWEST_STAT = np.nan
LOWEST_DAF_BIN = np.nan
LOWEST_B_VALUE = np.nan
LOWEST_RR_BIN = np.nan

bin_counts = np.zeros(shape=(NUM_BINS, 10, NUM_BINS))
other_variants = []
for i in range(NUM_BINS):
    other_variants.append([])

    for j in range(10):
        other_variants[i].append([])

        for k in range(NUM_BINS):
            other_variants[i][j].append([])

num_lowest = 0
lowest_sum = 0

start = 1
cur_chr = "1"

while(scan_line):
    split_line = scan_line.split()
    
    # parse admixture scan results file line
    chrom = split_line[scan_indices['chr']]
    loc = int(split_line[scan_indices['pos']])
    ref = split_line[scan_indices['ref']]
    alt = split_line[scan_indices['alt']]
    
    # get GWAS data for this chromosome and position, if it exists
    gwas_data = get_gwas_line(chrom, loc)

    # check that this position overlaps with GWAS
    if not isinstance(gwas_data, float) and (split_line[scan_indices['stat']] != "NA") and ((ref == gwas_data['ref'] and alt == gwas_data['alt']) or (ref == gwas_data['alt'] and alt == gwas_data['ref'])) and not math.isinf(gwas_data['beta']):

        stat = float(split_line[scan_indices['stat']])
        t_freq = float(split_line[scan_indices['t_freq']])
        t_exp = float(split_line[scan_indices['exp']])

        # get B statistic for this chromosome and position, if it exists
        bscore_data = get_bdecile(chrom, loc)

        # get recombination rate for this chromosome and position, if it exists
        r_rate = get_recombination_rate(chrom, loc)

        # verify overlap between scan file, B statistic file, and recombination rate map
        if np.isnan(r_rate) or isinstance(bscore_data, float) or isinstance(bscore_data['anc'], float) or isinstance(bscore_data['der'], float): 
            scan_line = scan_file.readline()
            continue

        # check B score file alleles against scan alleles
        if not((bscore_data['anc'] == ref and bscore_data['der'] == alt) or (bscore_data['anc'] == alt and bscore_data['der'] == ref)):
            scan_line = scan_file.readline()
            continue

        # assign the derived allele frequency
        if bscore_data['der'] == gwas_data['freq_allele']:
            daf = gwas_data['freq']
        elif bscore_data['anc'] == gwas_data['freq_allele']: 
            daf = 1 - gwas_data['freq']
        else:
            assert 0
        
        # calculate the derived allele frequency bin
        daf_bin = (NUM_BINS - 1) if daf == 1 else int(daf // (1 / NUM_BINS))

        # assign the direction of selection from scan results
        if t_freq > t_exp:
            direction = 1
        elif t_freq < t_exp:
            direction = -1
        else:
            direction = 0

        # get trait-decreasing allele from effect size
        inc_allele = gwas_data['ref'] if gwas_data['beta'] < 0 else gwas_data['alt'] 

        # if the alternate allele is the trait-increasing allele, reverse the sign of polarized statistic
        if alt == inc_allele:
            direction *= -1
        else:
            assert ref == inc_allele

        if USE_EFFECT_SIZE:
            # compute the polarized statistic using magnitude of beta
            polarized = abs(stat) * direction * abs(gwas_data['beta']) 
        else:
            # compute the polarized statistic only using hte direction of beta
            if gwas_data['beta'] == 0:
                polarized = 0
            else:
                polarized = abs(stat) * direction
            
        # iterate to find correct recombination rate bin
        rr_bin = 0
        while r_rate > RR_CUTOFFS[rr_bin]:
            rr_bin += 1

        # check if we're in a new window
        if (cur_chr != chrom or (cur_chr == chrom and loc >= start + WINDOW_SIZE)):
            # save current lowest variant information and reset 
            if not np.isnan(LOWEST_P):
                bin_counts[LOWEST_DAF_BIN][LOWEST_B_VALUE][LOWEST_RR_BIN] += 1
                lowest_sum += LOWEST_STAT
                num_lowest += 1
                
                # reset the lowest variant information
                LOWEST_P = np.nan
                LOWEST_STAT = np.nan
                LOWEST_DAF_BIN = np.nan
                LOWEST_B_VALUE = np.nan
                LOWEST_RR_BIN = np.nan
                LOWEST_CHR = np.nan
                LOWEST_POS = np.nan

            # jump to the next chromosome if necessary
            if cur_chr != chrom:
                cur_chr = chrom 
                start = 1
            # move the window forward until the current variant's location is inside the window
            while(loc >= start + WINDOW_SIZE):
                start += WINDOW_SIZE

        assert (cur_chr == chrom and loc < start + WINDOW_SIZE)

        # check if the current variant is lower than current lowest and the p value is lower than the threshold in this window
        if (np.isnan(LOWEST_P) or gwas_data['p_val'] < LOWEST_P) and gwas_data['p_val'] <= P_CUTOFF:
            if not np.isnan(LOWEST_P):
                # bin the current lowest, since it's no longer the lowest in this window
                other_variants[LOWEST_DAF_BIN][LOWEST_B_VALUE][LOWEST_RR_BIN].append(LOWEST_STAT) 
            # save the current variant as the current lowest in this window
            LOWEST_P = gwas_data['p_val']
            LOWEST_STAT = polarized
            LOWEST_DAF_BIN = daf_bin
            LOWEST_B_VALUE = bscore_data['decile']
            LOWEST_RR_BIN = rr_bin
            LOWEST_CHR = chrom
            LOWEST_POS = loc
        else:
            # this variant is not lower than the current lowest, so bin it
            other_variants[daf_bin][bscore_data['decile']][rr_bin].append(polarized)

    scan_line = scan_file.readline()

# check the last current lowest statistic
if not np.isnan(LOWEST_P):
    lowest_sum += LOWEST_STAT
    bin_counts[LOWEST_DAF_BIN][LOWEST_B_VALUE][LOWEST_RR_BIN] += 1
    num_lowest += 1

if num_lowest != 0:
    lowest_avg = lowest_sum / num_lowest 
    print("Number of lowest variants:", num_lowest)
    print("Lowest average:", lowest_avg)

    lower_count = 0
    higher_count = 0
    equal_count = 0

    for n in range(NUM_TRIALS):
        trial_sum = 0
        trial_valid = 0
        for i in range(NUM_BINS):
            for j in range(10):
                for k in range(NUM_BINS):
                    # get the number of lowest variants in this bin
                    num_sample = int(bin_counts[i][j][k])
                    if num_sample != 0:
                        # sample the number of lowest variants in this bin from the list of non-lowest variants in this bin
                        sampled = random.sample(other_variants[i][j][k], num_sample)
                        trial_sum += np.sum(sampled)
                        trial_valid += num_sample
        assert trial_valid == num_lowest

        # calculate average polarized statistic for this trial
        trial_avg = trial_sum / trial_valid
        if trial_avg < lowest_avg:
            lower_count += 1
        elif trial_avg > lowest_avg:
            higher_count += 1
        else:
            equal_count += 1

    print("Number of trials lower:", lower_count)       
    print("Number of trials higher:", higher_count)     
    print("Number of trials equal:", equal_count)       

    out_path = config['out_path'] + "trial_results/" 
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    
    # save the results to a .npy file to be auto compiled into a spreadsheet
    np.save(out_path + TRAIT + "_" + str(P_EXP) + "_" + EPOCH + ".npy", [num_lowest, lower_count, higher_count])

rr_file.close()
scan_file.close()
b_file.close()
gwas_file.close()

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
    gwas_path = config['paths']['gwas_path'] + DATASET + "/hum0197.v3.EUR." + TRAIT + ".v1/*.auto.txt"
    print(gwas_path)
    gwas_file = open(glob.glob(gwas_path)[0], 'r')
elif DATASET == "BBJ":
    gwas_path = config['paths']['gwas_path'] + DATASET + "/hum0197.v3.BBJ." + TRAIT + ".v1/*.auto.txt"
    gwas_file = open(glob.glob(gwas_path)[0], 'r')
else:
    assert 0

# quantiles computed from LD score file
LD_CUTOFFS = [1.42309722e-19, 1.33096160e-06, 7.99490360e-04, 1.00316488e-02, 5.29890833e-02, 1.98940513e-01, 8.59677415e-01, 9999999999]

ld_file = open(config['paths']['ld_file'], 'r')
scan_file = open(config['paths']['selection_scan_path'], 'r')
b_file = open(config['paths']['b_file'], 'r')

NUM_BINS = config['options']['num_bins']
NUM_TRIALS = config['options']['num_trials']
WINDOW_SIZE = config['options']['window_size']

# BBJ headers
COL_HEADERS = [["CHR", "POS", "p.value", "Allele1", "Allele2", "BETA", "AF_Allele2_UKB", "Allele2"], ["CHR", "BP", "P_BOLT_LMM_INF", "ALLELE0", "ALLELE1", "BETA", "A1FREQ", "ALLELE1"], ["CHR", "POS", "p.value", "Allele1", "Allele2", "BETA", "AF_Allele2", "Allele2"]]

COL_CHR = 0
COL_POS = 1
COL_P = 2
COL_REF = 3
COL_ALT = 4
COL_BETA = 5
COL_FREQ = 6
COL_FREQA = 7

print("Population:", EPOCH)
print("GWAS Path:", gwas_path)
print("Threshold:", P_CUTOFF)
print("Number of trials:", NUM_TRIALS)
print("Number of bins:", NUM_BINS)

invalid_count = 0

def get_gwas_line(find_chr, find_loc):
    """ Parses GWAS file line for a given loci
    Args:
        find_chr: chromosome of location to get information for
        find_loc: position of location to get information for
    Returns:
        dictionary containing GWAS p-value, GWAS beta, ref allele, alt allele, allele frequency
    """
    global invalid_count
    global gwas_line
    while(gwas_line):
        g_split = gwas_line.split()
        g_chr = g_split[COL_INDICES[COL_CHR]]
        g_loc = int(float(g_split[COL_INDICES[COL_POS]]))

        if g_chr == "X" or g_chr == "Y":
            return np.nan

        # check for match
        if g_chr == find_chr and find_loc == g_loc:
            try:
                test = float(g_split[COL_INDICES[COL_FREQ]])
            except:
                gwas_line = gwas_file.readline()
                return np.nan
            return {'p_val': float(g_split[COL_INDICES[COL_P]]), 'ref': g_split[COL_INDICES[COL_REF]], 'alt':g_split[COL_INDICES[COL_ALT]], 'beta':float(g_split[COL_INDICES[COL_BETA]]), 'freq':float(g_split[COL_INDICES[COL_FREQ]]), 'freq_allele':g_split[COL_INDICES[COL_FREQA]]}

        # check if we've gone too far
        if int(g_chr) > int(find_chr) or (g_chr == find_chr and g_loc > find_loc):
            return np.nan

        # keep searching 
        gwas_line = gwas_file.readline()
    
    # already reached end of file
    return np.nan

def get_ld(find_chr, find_loc):
    """ Finds LD score for a given location
    Args:
        find_chr: chromosome of location to get LD score for
        find_loc: position of location to get LD score for
    Returns:
        LD score
    """
    global ld_line
    while(ld_line):
        ld_split = ld_line.split()
        ld_chr = ld_split[0]
        ld_begin = int(ld_split[1])
        ld_end = int(ld_split[2])

        if ld_chr == "X":
            return np.nan

        if ld_chr == find_chr and find_loc >= ld_begin and find_loc < ld_end:
            return float(ld_split[3])

        if int(ld_chr) > int(find_chr) or (ld_chr == find_chr and ld_begin > find_loc):
            return np.nan

        ld_line = ld_file.readline()
    return np.nan

def get_bscore(find_chr, find_loc):
    """ Finds B statistic for a given location
    Args:
        find_chr: chromosome of location to get B statistic for
        find_loc: position of location to get B statistic for
    Returns:
        B statistic
    """
    global b_line   
    while(b_line):
        b_split = b_line.split()
        b_chr = b_split[1]
        b_loc = int(b_split[3])

        b_ref = b_split[4]
        b_alt = b_split[5]

        anc_der = b_split[12]

        b_anc = np.nan
        b_der = np.nan

        if anc_der != "." and anc_der != ".,.":
            ref_code = int(anc_der.split(",")[0])
            alt_code = int(anc_der.split(",")[1])
            if ref_code == 1 and alt_code == 0:
                b_anc = b_alt
                b_der = b_ref
            elif ref_code == 0 and alt_code == 1:
                b_anc = b_ref
                b_der = b_alt
            else:
                assert 0

        # reached the end of the relevant lines
        if b_chr == "23":
            return np.nan

        # check for match
        if b_chr == find_chr and find_loc == b_loc:
            if b_split[6] == "." or b_split[6][-1:] == "-":
                return np.nan
            return {'score':int(b_split[6][-1:]), 'anc':b_anc, 'der':b_der}

        # check if we've gone too far
        if int(b_chr) > int(find_chr) or (b_chr == find_chr and b_loc > find_loc):
            return np.nan

        # keep searching
        b_line = b_file.readline()
    return np.nan

def get_column_indices(header):
    """ Returns column indices for values in GWAS files (so that different GWAS file formats can be submitted in the same batch of jobs) 
    Args:
        header line
    Returns:
        list of indices corresponding to relevant columns
    """
    h_split = header.split()
    indices = np.zeros(len(COL_HEADERS[0]), dtype=int)

    found = False
    for head_type in COL_HEADERS:
        try:
            j = 0
            for col in head_type:
                indices[j] = h_split.index(col)
                j += 1
            found = True
        except:
            continue
    assert found
    return indices

# use global files/lines so we only have to pass through the files once
global gwas_line
global ld_line
global b_line
global scan_line

scan_line = scan_file.readline()
header = scan_line.split()

# column indices in the admixture scan result file
COL_SCAN_CHR = header.index("CHROM")
COL_SCAN_POS = header.index("POSITION")
COL_SCAN_REF = header.index("REF")
COL_SCAN_ALT = header.index("ALT")
COL_SCAN_T_FREQ = header.index("TARGET_FREQ")
COL_SCAN_T_EXP = header.index("EXPECTED")
COL_SCAN_STAT = header.index("correctedStat")

scan_line = scan_file.readline()

gwas_line = gwas_file.readline()
COL_INDICES = get_column_indices(gwas_line)
gwas_line = gwas_file.readline()


ld_line = ld_file.readline()
ld_line = ld_file.readline()

b_line = b_file.readline()
# skip over header lines
while (b_line[0] == '#'):
    b_line = b_file.readline()

# initialize data structures and window variables
LOWEST_P = np.nan
LOWEST_STAT = np.nan
LOWEST_DAF_BIN = np.nan
LOWEST_B_VALUE = np.nan
LOWEST_LD_BIN = np.nan

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
    chrom = split_line[COL_SCAN_CHR]
    loc = int(split_line[COL_SCAN_POS])
    ref = split_line[COL_SCAN_REF]
    alt = split_line[COL_SCAN_ALT]
    
    # get GWAS data for this chromosome and position, if it exists
    gwas_data = get_gwas_line(chrom, loc)

    # check that this position overlaps with GWAS
    if not isinstance(gwas_data, float) and (split_line[COL_SCAN_STAT] != "NA") and ((ref == gwas_data['ref'] and alt == gwas_data['alt']) or (ref == gwas_data['alt'] and alt == gwas_data['ref'])) and not math.isinf(gwas_data['beta']):

        stat = float(split_line[COL_SCAN_STAT])
        t_freq = float(split_line[COL_SCAN_T_FREQ])
        t_exp = float(split_line[COL_SCAN_T_EXP])

        # get B statistic for this chromosome and position, if it exists
        bscore_data = get_bscore(chrom, loc)

        # get LD score for this chromosome and position, if it exists
        ld_score = get_ld(chrom, loc)

        # verify overlap between scan file, B statistic file, and LD score file
        if np.isnan(ld_score) or isinstance(bscore_data, float) or isinstance(bscore_data['anc'], float) or isinstance(bscore_data['der'], float): 
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
            
                # iterate to find correct LD bin
        ld_bin = 0
        while ld_score > LD_CUTOFFS[ld_bin]:
            ld_bin += 1

        # check if we're in a new window
        if (cur_chr != chrom or (cur_chr == chrom and loc >= start + WINDOW_SIZE)):
            # save current lowest variant information and reset 
            if not np.isnan(LOWEST_P):
                bin_counts[LOWEST_DAF_BIN][LOWEST_B_VALUE][LOWEST_LD_BIN] += 1
                lowest_sum += LOWEST_STAT
                num_lowest += 1
                
                # reset the lowest variant information
                LOWEST_P = np.nan
                LOWEST_STAT = np.nan
                LOWEST_DAF_BIN = np.nan
                LOWEST_B_VALUE = np.nan
                LOWEST_LD_BIN = np.nan
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
                other_variants[LOWEST_DAF_BIN][LOWEST_B_VALUE][LOWEST_LD_BIN].append(LOWEST_STAT) 

            # save the current variant as the current lowest in this window
            LOWEST_P = gwas_data['p_val']
            LOWEST_STAT = polarized
            LOWEST_DAF_BIN = daf_bin
            LOWEST_B_VALUE = bscore_data['score']
            LOWEST_LD_BIN = ld_bin
            LOWEST_CHR = chrom
            LOWEST_POS = loc
        else:
            # this variant is not lower than the current lowest, so bin it
            other_variants[daf_bin][bscore_data['score']][ld_bin].append(polarized)

    scan_line = scan_file.readline()

# check the last current lowest statistic
if not np.isnan(LOWEST_P):
    lowest_sum += LOWEST_STAT
    bin_counts[LOWEST_DAF_BIN][LOWEST_B_VALUE][LOWEST_LD_BIN] += 1
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

    out_path = config['paths']['out'] + "trial_results/" 
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    
    # save the results to a .npy file to be auto compiled into a spreadsheet
    np.save(out_path + TRAIT + "_" + str(P_EXP) + "_" + EPOCH + ".npy", [num_lowest, lower_count, higher_count])

ld_file.close()
scan_file.close()
b_file.close()
gwas_file.close()

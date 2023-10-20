#!/usr/bin/env python3
""" Combine Polygenic Results

Compiles results from multiple polygenic_window_test.py runs. Uses the settings in polygenic_config.yml and ../monogenic_scan/monogenic_config.yml. 

Command-line Arguments:
    [1] DATASET: indicates whether the GWAS summary statistics are from Biobank Japan, UK Biobanki, or within-sib consortium
    [2] THRESHOLDS: list of p-value threshold exponents, cannot have spaces (e.g., [2,4,6,8])
    [3] REPORT_THRESH: one p-value threshold exponent to use to report final significant traits

For example, the command 'combine_polygenic_results.py UKB [2,4,6,8] 6' will search for runs that used the UK Biobank and p-value thresholds of 1e-2, 1e-4, 1e-6, and 1e-8. The trait list output file will report significant results that used the 1e-6 threshold.

Outputs:
    /prefix/in/polygenic_config/{DATASET}_full_results.txt: all trials results
    /prefix/in/polygenic_config/{DATASET}_sig_trials.txt: all significant trial results
    /prefix/in/polygenic_config/{DATASET}_traits_list.txt: all significant trial results using a threshold of REPORT_THRESH
"""
import yaml
import numpy as np
import sys

config = yaml.safe_load(open("polygenic_config.yml"))
mono_config = yaml.safe_load(open("../monogenic_scan/monogenic_config.yml"))
prefix = config['out_path']

# defaults to one population (specified in monogenic_config.yml), but can add more manually here
POPULATIONS = [mono_config['populations']['target']]
DATASET = sys.argv[1]
THRESHOLDS = sys.argv[2][1:-1].split(',')
REPORT_THRESH = sys.argv[3]


# list of trait IDs for each biobank
if DATASET == "UKB":
    traits = ["AA", "Abo", "AC", "AD", "AG", "AIH", "Ang", "ApA", "AP", "ARF", "AR", "As", "ATR", "BC", "BD", "BP", "Bro", "Br", "BtC", "BT", "Car", "Cat", "CA_Un", "CB", "CC", "CD", "CeAn", "CeC", "CF", "CG", "CHB", "CHC", "CHF", "ChG", "Cho", "ChP", "ChS", "Cir", "CL", "COPD", "Co", "CRF", "CSOM", "CS", "CTS", "CVD", "Cys", "DC", "Dep", "Dip", "DN", "Dys", "EcP", "EC", "EM", "EnC", "Ep", "EV", "FA", "FNP", "FS", "GC", "GD", "GERD", "Gla", "Goi", "GP", "GU", "HAV", "HC", "HD", "Her", "HI", "HL", "Hype", "Hypo", "IDA", "IH", "ILD", "Ile", "Ins", "IN", "Iri", "IS", "ITP", "JRA", "Kel", "LiC", "LuC", "Mas", "MA", "Men", "MG", "MI", "ML", "Myo", "NB", "NMI", "NP", "NS", "OC", "OvC", "PaC", "PAD", "PA", "PD", "PeD", "Per", "PE", "PF", "PKD", "PLC", "Pl", "PM", "Pneu", "Pne", "Pn", "Pol", "PrC", "PsV", "PT", "Pye", "RA", "RD", "RF", "RW", "SAP", "Sar", "SAS", "Sch", "SCS", "SD", "SH", "SkC", "SLE", "SS", "T1D", "T2D", "Tat", "ThC", "Ton", "TS", "Typ", "UAP", "UC", "UF", "UP", "Uro", "Urt", "Uve", "Var", "VA", "Zos"]
elif DATASET == "BBJ":
    traits = ["A02B", "A10", "AA", "Abo", "AC", "AD", "AG", "AIH", "Alb", "ALP", "ALT", "Ang", "ApA", "AP", "ARF", "AR", "AST", "As", "ATR", "B01A", "BAS", "BC", "BD", "BMI", "BP", "Bro", "Br", "BtC", "BT", "BUN", "BW", "C01D", "C02", "C03", "C07", "C08", "C09", "C10AA", "Car", "Cat", "CA_Un", "Ca", "CB", "CC", "CD", "CeAn", "CeC", "CF", "CG", "CHB", "CHC", "CHF", "ChG", "Cho", "ChP", "ChS", "Cir", "CL", "COPD", "Co", "CRF", "CRP", "CSOM", "CS", "CTS", "CVD", "Cys", "DBP", "DC", "Dep", "Dip", "DN", "Dys", "EcP", "EC", "EM", "EnC", "EOS", "Ep", "EV", "FA", "FNP", "FS", "GBP", "GC", "GD", "GERD", "GGT", "Gla", "Glucose", "Goi", "GU", "H03A", "HAV", "HbA1c", "Hb", "HC", "HDLC", "HD", "Hei", "Her", "HI", "HL", "Ht", "Hype", "Hypo", "IDA", "IH", "ILD", "Ile", "Ins", "IN", "Iri", "IS", "ITP", "JRA", "Kel", "L04", "LDLC", "LiC", "LuC", "LYM", "M01A", "M05B", "MAP", "Mas", "MA", "MCHC", "MCH", "MCV", "Men", "MG", "MI", "ML", "MON", "Myo", "N02A", "N02BA", "N02BE", "N02C", "N06A", "NB", "NEU", "NMI", "NP", "NS", "OC", "OvC", "PaC", "PAD", "PA", "PD", "PeD", "Per", "PE", "PF", "PKD", "PLC", "PLT", "Pl", "PM", "Pneu", "Pne", "Pn", "Pol", "PP", "PrC", "PsV", "PT", "Pye", "R03A", "R03BA", "R06A", "RA", "RBC", "RD", "RF", "RP", "RW", "S01E", "SAP", "Sar", "SAS", "SBP", "Sch", "sCr", "SCS", "SD", "SH", "SkC", "SLE", "SS", "T1D", "T2D", "Tat", "TBil", "TC", "TG", "ThC", "Ton", "TP", "TS", "Typ", "UAP", "UA", "UC", "UF", "UP", "Uro", "Urt", "Uve", "Var", "VA", "WBC", "Zos"]
elif DATASET == "sib":
    traits = np.arange(4813, 4861)
else:
    assert 0

full_file = open(prefix + DATASET + "_full_results.txt",'w')
sig_file = open(prefix + DATASET + "_sig_trials.txt", 'w')
trait_file = open(prefix + DATASET + "_traits_list.txt", 'w')

# write headers
full_file.write("Population\tTrait ID\tThreshold\tNumber SNPs\tNumber Lowest\tNumber Higher\tPercent Lower\tPercent Higher\n")
sig_file.write("Population\tTrait ID\tThreshold\tNumber SNPs\tNumber Lowest\tNumber Higher\tPercent Lower\tPercent Higher\tDirection\n")
trait_file.write("Population\tTrait ID\tNumber SNPs\tDirection\tPercent Lower\tPercent Higher\n")

i = 0
while(i < len(traits)):
    trait = str(traits[i])
    for pop in POPULATIONS:
        for thresh in THRESHOLDS:
            try:
                results = np.load(prefix + "/trial_results/" + DATASET + "/" + trait + "_" + thresh + "_" + pop + ".npy") 
            except:
                continue

            n_snps = results[0]
            n_low = results[1]
            n_high = results[2]
            p_low = round(n_low / config['options']['num_trials'] * 100, 2)
            p_high = round(n_high / config['options']['num_trials'] * 100, 2)

            full_file.write(pop + "\t" + trait + "\t1e-" + thresh + "\t" + str(n_snps) + "\t" + str(n_low) + "\t" + str(n_high) + "\t" + str(p_low) + "\t" + str(p_high) + "\n")    
            if (p_low <= 2.5):
                sig_file.write(pop + "\t" + trait + "\t1e-" + thresh + "\t" + str(n_snps) + "\t" + str(n_low) + "\t" + str(n_high) + "\t" + str(p_low) + "\t" + str(p_high) + "\tDecrease\n") 
                if thresh == REPORT_THRESH:
                    trait_file.write(pop + "\t" + trait + "\t" + str(n_snps) + "\tDecrease\t" + str(p_low) + "\t" + str(p_high) + "\n")
            if (p_high <= 2.5):
                sig_file.write(pop + "\t" + trait + "\t1e-" + thresh + "\t" + str(n_snps) + "\t" + str(n_low) + "\t" + str(n_high) + "\t" + str(p_low) + "\t" + str(p_high) + "\tIncrease\n") 
                if thresh == REPORT_THRESH:
                    trait_file.write(pop + "\t" + trait + "\t" + str(n_snps) + "\tIncrease\t" + str(p_low) + "\t" + str(p_high) + "\n")

    i += 1  

full_file.close()
sig_file.close()
trait_file.close()

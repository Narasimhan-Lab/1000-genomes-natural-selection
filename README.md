# 1,000 ancient genomes uncover 10,000 years of natural selection in Europe
The bioRxiv preprint can be found here: https://biorxiv.org/content/10.1101/2022.08.24.505188

The corresponding authors can be contacted at vagheesh@utexas.edu, reich@genetics.med.harvard.edu, and arbelharpak@utexas.edu.

## Monogenic Scan (src/monogenic_scan)

### Running the selection scan

The `monogenic_scan.py` script performs a genome-wide scan for SNPs whose allele frequencies deviate from expectation based on genome-wide mixing proportions. `monogenic_config.yml` is set up for our scan in the Neolithic as an example, and there are a couple of variables that need to be modified to run a different model:

* `populations`: population labels
* `qpAdm_results`: qpAdm output
    * `alpha`: mixing proportion estimate
    * `standard_error`: standard error (jackknife)
    * `multiplier`: amount to multiply the standard error by during perturbation tests
* `paths`: data and output paths
    * `ref_reads`, `alt_reads`: .hdf5 files with read count matrices
    * `target_indices`, `s1_indices`, `s2_indices`: NumPy arrays that contain the matrix indices of the samples for a given population (e.g., for a population P composed of three samples that were in columns 0, 14, and 23 of the read count matrices, the .npy file for that population would be the array [0, 14, 23])
    * `chrom_pos`: tab-delimited file (with no header) that specifies chromosome, position, ref, and alt for each position being scanned (used for result summary)

The script is set up for a two-way admixture model but can be extended by adding source populations to the `populations` and `ids_list` lists and modifying `monogenic_config.yml`.

After the header has been modified, the script can be run with
```
./monogenic_scan.py p n
```
where `p` is the relative position of the SNP where computation should start (with zero-indexing),  and `n` is the number positions to perform calculations on (so that the genome can be split up into segments for parallel computation). For example, `./monogenic_scan.py 5000 1000` would scan 1000 sequential positions starting from the 5001st row in the read count matrices.

The ```divide_jobs.sh``` script is an example of how to automatically divide the genome to be scanned as parallel segments.

The script will save the following information for each position into NumPy files: estimated allele frequencies for each population, expected allele frequency of the target population, number of samples with non-zero read count data for each population, uncorrected chi-squared statistic for the likelihood ratio test, uncorrected p-value of the likelihood ratio test, and intermediate calculations from the likelihood ratio test. In the previous example of running `./monogenic_scan.py 5000 1000`, the p-values from the likelihood ratio test would be saved as a 1000-length array in the file `p_vals5000.npy`, the chi-squared statistics as `stats5000.npy`, and so on.

### Compiling scan results

The `compile_monogenic_results.py` file will compile the output NumPy files from the selection scan into a .txt file. The script can be run with
```
./compile_monogenic_results.py n k
```
where `n` is the number of segments the genome was divided into (e.g. using ```divide_jobs.sh```, and `k` is the number of positions in each segment.

## Polygenic Scan

### Running the scan

The GWAS summary statistics used with this script can be found at https://humandbs.biosciencedbc.jp/files/hum0197.org/. The script is run as follows:
```
./polygenic_window_test.py TRAIT P_EXP DATASET
```
with the following parameter descriptions:
* `TRAIT`: GWAS trait ID (e.g., A02B, A10, AA...)
* `P_EXP`: p-value cutoff exponent for choosing significant GWAS SNPs (e.g., a `P_EXP` of 6 corresponds to using a cutoff of 1e-6)
* `DATASET`: UKB or BBJ, indicates whether the GWAS summary statistics are from Biobank Japan or UK Biobank (if you use a different dataset, you will have to edit the paths and headers in `polygenic_window_test.py`)

The following variables will also need to be adjusted in polygenic_config.yml:
* `out_path`: directory to save the polygenic scan results
* `gwas_files`: 
    * `path`: path to directory containing GWAS summary statistics, where `gwas_file` contains a `UKB` directory and a `BBJ` directory that contain the unzipped summary statistic files. You will have to edit the paths in `polygenic_window_test.py` to use different datasets
    * column names for chromosome position, GWAS p-value, effect size, alleles, and allele frequency
* `selection_file`: file with monogenic selection scan statistics (i.e., compiled monogenic scan results after correction for genomic inflation)
    * `path`: path to file with statistics (must be sorted by location)
    * column names for chromosome, position, ref/alt alleles, target population frequency, expected target population allele frequency, and selection statistic (should be NA for invalid results)
* `rr_file`: 
    * `path`: path to file with recombination rate map (must be sorted by location)
    *  column names for chromosome, start and end of region, and recombination rate
* `b_file`: 
    * `path`: path to B statistic deciles file (must be sorted by location)
    * column names for chromosome, position, B statistic decile (i.e., 0, 1, ..., 9), and ancestral/derived alleles
* options: hyperparameters
    * `use_effect_size`: whether to use the magnitude of the GWAS beta when calculating polygenic statistic
    * `num_trials`: number of trials used to generate null distribution
    * `window_size`: size of window (in bp) when picking most significant SNP

### Compiling scan results

Each run of `polygenic_window_test.py` will print the results of the scan. You can also use `combine_polygenic_results.py` to generate a summary of multiple runs into tab-delimited files. The script is run with
```
./combine_polygenic_results.py DATASET THRESHOLDS REPORT_THRESH
```
with the following parameter descriptions:
* `DATASET`: UKB or BBJ, indicating whether the GWAS summary statistics are from UK Biobank or Biobank Japan
* `THRESHOLDS`: list of p-value threshold exponents to include, cannot have spaces (e.g., [2,4,6,8])
* `REPORT_THRESH`: one p-value threshold exponent to use to report final significant traits

For example, the command `combine_polygenic_results.py UKB [2,4,6,8] 6` will search for runs that used the UK Biobank and p-value thresholds of 1e-2, 1e-4, 1e-6, and 1e-8. The trait list output file will report significant results that used the 1e-6 threshold.

The combined results will be reported in the following files:
* `/prefix/in_polygenic_config/{DATASET}_full_results.txt`: all trials results
* `/prefix/in_polygenic_config/{DATASET}_sig_trials.txt`: all significant trial results
* `/prefix/in_polygenic_config/DATASET}_traits_list.txt`: all significant trial results using a threshold of REPORT_THRESH

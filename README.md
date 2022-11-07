# 1,000 ancient genomes uncover 10,000 years of natural selection in Europe
The bioRxiv preprint can be found here: https://biorxiv.org/content/10.1101/2022.08.24.505188
The corresponding authors can be contacted at vagheesh@utexas.edu, reich@genetics.med.harvard.edu, and arbelharpak@utexas.edu.

## Monogenic Scan

### Running the selection scan

The `monogenic_scan.py` script performs a genome-wide scan for SNPs whose allele frequencies deviate from expectation based on genome-wide mixing proportions. The script is set up for our scan in the Neolithic as an example, and there are a couple of variables in the header that need to be modified to run a different model:

* `target`, `source1`, `source2`: population labels
* `pop_proportions`: list of mixing proportions, where the first entry is the mixing proportion for source1, and the second entry is the mixing proportion for source2
* `ref_path`, `alt_path`: paths to .hdf5 files containing the matrices of reference and alternate allele read counts, where rows correspond to positions and columns correspond to samples
* `target_npy, s1_npy, s2_npy`: paths to numpy arrays/files that contain the matrix indices of the samples for a given population (e.g., for a population P composed of three samples that were in columns 0, 14, and 23 of the read count matrices, the .npy file for that population would be the array [0, 14, 23])

The script is set up for a two-way admixture model but can be extended by adding source populations to the `populations` and `ids_list` lists. 

After the header has been modified, the script can be run with
```
./monogenic_scan.py p n
```
where `p` is the relative position of the SNP to start computation at (with zero-indexing),  and `n` is the number positions to perform calculations on (so that the genome can be split up into segments for parallel computation). For example, `./monogenic_scan.py 5000 1000` would scan 1000 sequential positions starting from the 5001st row in the read count matrices.

The script will save the following information for each position into numpy files: estimated allele frequencies for each population, expected allele frequency of the target population, uncorrected chi-squared statistic for the likelihood ratio test, and uncorrected p-value of the likelihood ratio test. In the previous example of running `./monogenic_scan.py 5000 1000`, the p-values from the likelihood ratio test would be saved as a 1000-length array in the file `p_vals5000.npy`, the chi-squared statistics as `stats5000.npy`, and so on.

### Compiling scan results

The `compile_monogenic_results.py` file will compile the output numpy files from the selection scan into a .txt file. The following paths should be edited in the header:

* `scan_prefix`: path to the directory that contains the .npy output files
* `out_path`: desired path for the output file
* `pos_path`: path to tab-delimited file containing chromosome, position, reference allele, and alternate allele information for each position, where the kth row of the file should correspond to the kth rows of the read count matrices. If this file does not have a header, either add a filler line to the beginning of the file or delete line 52 of the script. 

The script can then be run with
```
./compile_monogenic_results.py n k
```
where `n` is the number of segments the genome was divided into, and `k` is the number of positions in each segment.

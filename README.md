# minSNPs_nanopore
This documents the analysis performed in using minSNPs to analyse Nanopore sequence data.
The analysis were performed on a HPC using Slurm scheduler, and the scripts may need to be adjusted accordingly when ran elsewhere. 

# Testing with minSNPs
## Identifying the SNPs for the new PubMLST samples
**Steps:**
1. Use SNIPPY to call SNPs for the new sample (see `Scripts/SNIPPY_RUN.sh`).
2. Extract all the SNPs defined in megaalignment from SNIPPY result (see `Scripts/extract_megaalignment_snps.R`), output of final alignment `mega_with_newdata.fasta` is in figshare [see here](https://figshare.com/s/464f38a92cde2fb067cc).
3. Assign CC metadata to new samples based on most similar sample in megaalignment (see `Scripts/assign_cc_meta.R`), see `Results/new_sample_most_similar_mega.csv` for the result containing the most similar sample in megaalignment and the CC to be assigned, `Results/disagreement.csv` for the samples with CC metadata from PubMLST that is different from the CC to be assigned.
4. Neighbour Joining tree is created with MEGA with `mega_with_newdata.fasta` as the input and all default parameters, output: `Result/mega_with_newdata_NJT.nwk`. see [here](https://microreact.org/project/nPEW3sbjQL3tD6EhMwNynj-minsnps-and-nanopore-analysis) for the interactive tree.

## Sampling 400 random SNPs
All the SNPs are scrambled and the first 400 is taken, see `Scripts/random_snps_selection.R`.

## Search string generation
For
- SNPs, see `Scripts/generate_snp_search_sequence.R`.
- gene sequences, see `Scripts/generate_gene_search_sequence.R`.

## Major lineage assignment

### Testing with lab generated Nanopore data
**Steps:**
1. Scan lab generated Nanopore data `Scripts/scan_lab_nanopore.R`.
2. Transform most similar isolate to most likely CC and aggregated all results for different number of SNPs or reads used : `Scripts/aggregate_lab_result.R`.

### Testing with simulated long-read data
**Steps:**
1. Simulate long read with pbsim2 (see `Scripts/simulated_long_read_generation.R`), see `Data/BIGSdb_3343897_1178826571_39602.csv` for the list of data download from pubMLST and `Results/sim_n_reads.csv` output for number of reads generated.
2. Assign most likely CC for tested samples based on SNPs distance: `Scripts/assign_cc_meta.R`
3. Scan simulated Nanopore data `Scripts/scan_simulated_read.R`.
4. Transform most similar isolate to most likely CC and aggregated all results for different number of SNPs or reads used : `Scripts/aggregate_simulation_result.R`.

## Gene detection
### Testing with lab generated Nanopore data
1. Scan lab generated Nanopore data `Scripts/scan_lab_gene.R`.
2. Transform most similar isolate to most likely CC and aggregated all results for different number of SNPs or reads used : `Scripts/aggregate_lab_gene.R`.

### Testing with simulated data
1. Scan simulated Nanopore data `Scripts/scan_simulated_gene.R`.
2. Transform most similar isolate to most likely CC and aggregated all results for different number of SNPs or reads used : `Scripts/aggregate_simulation_gene.R`.


# Comparison
## Comparison with Krocus
- Simulated data:
    - Script: `Scripts/krocus_comparison.R`
    - Result: `Results/krocus_full_result.csv` and  `Results/krocus_summary.csv`
- Lab generated data:
    - Script: `Scripts/lab_krocus_comparison.R`
    - Result: `Results/lab_krocus_full_result.csv` and `Results/lab_krocus_summary.csv`

## Comparison with Sketchy
- Simulated data:
    - Script: `Scripts/sketchy_comparison.R`
    - Result: `Results/sketchy_gene_summary.csv` and `Results/sketchy_summary.csv`
- Lab generated data:
    - Script: `Scripts/lab_sketchy_comparison.R`
    - Result: `Results/lab_sketchy_gene_summary.csv` and `Results/lab_sketchy_summary.csv`
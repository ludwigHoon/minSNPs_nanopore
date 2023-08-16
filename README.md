# minSNPs_nanopore
This documents the analysis performed in using minSNPs to analyse Nanopore sequence data.
The analysis were performed on a HPC using Slurm scheduler, and the scripts may need to be adjusted accordingly when ran elsewhere. 

# Steps
1. Simulate long read with pbsim2 (see `Scripts/simulated_long_read_generation.R`), see `Data/BIGSdb_3343897_1178826571_39602.csv` for the list of data download from pubMLST and `Results/sim_n_reads.csv` output for number of reads generated.
2. Use SNIPPY to call SNPs for the new sample (see `Scripts/SNIPPY_RUN.sh`).
3. Extract all the SNPs defined in megaalignment from SNIPPY result (see `Scripts/extract_megaalignment_snps.R`), output of final alignment is in figshare [see here](https://figshare.com/s/464f38a92cde2fb067cc).
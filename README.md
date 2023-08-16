# minSNPs_nanopore
This documents the analysis performed in using minSNPs to analyse Nanopore sequence data.
The analysis were performed on a HPC using Slurm scheduler, and the scripts may need to be adjusted accordingly when ran elsewhere. 

# Steps
1. Simulated long read with pbsim2 (see `Scripts/simulated_long_read_generation.R`), see `Data/BIGSdb_3343897_1178826571_39602.csv` for the list of data download from pubMLST and `Results/sim_n_reads.csv` output for number of reads generated.

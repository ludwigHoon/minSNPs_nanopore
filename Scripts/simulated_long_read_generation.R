## This scripts uses pbsim2 to simulate long-read data for testing
## It expects pbsim2 to be available in system path. 
## Pbsim2 can be built and installed from: https://github.com/yukiteruono/pbsim2 
## Parameter used:
## - seed: 42
## - hmm_model: R95 (FIC-HMM of R95 chemistry)
## - difference-ratio: 23:31:46 (substitution:insertion:deletion), as per recommendation for Nanopore simulation
## - depth: 100
## - accuracy mean: 0.9
PubMLST_export_list <- "BIGSdb_3343897_1178826571_39602.csv"
library(BiocParallel)
library(data.table)
all_data <- read.csv(PubMLST_export_list)

sim_creation_result <- bplapply(seq_len(nrow(all_data)), function(dat, all_data){
    id <- all_data[dat,"id"]
    fas_file <- list.files(pattern = paste0("^", id, ".*fas$"))
    if (length(fas_file) == 0){
        warning(paste0("No fasta file for ", id))
    }
    command <- sprintf("pbsim --seed 42 --hmm_model ~/pbsim2/data/R95.model --difference-ratio 23:31:46 --depth 100 --accuracy-mean 0.9 --prefix %s %s", paste0("S", id), fas_file)
    system(command, wait = TRUE, intern=FALSE)
    unlink(paste0("S", id, "_0001.", c("ref", "maf")))
    output_file <- paste0("S", id, "_0001.fastq")
    n_reads <- system(paste0("wc -l ", output_file), wait = TRUE, intern=TRUE)
    return(data.frame(id = id, fas_file = fas_file, n_reads = n_reads))
}, all_data = all_data, BPPARAM = BatchtoolsParam(workers = 391, cluster="slurm", template="~/slurm_bc_template.tmpl",
    resources=list(walltime=60*60*24*5, ncpus=1)))

n_reads <- rbindlist(sim_creation_result)
n_reads$n_reads <- sapply(strsplit(n_reads$n_reads, split = " "), `[`, 1)
n_reads$n_reads <- as.numeric(n_reads$n_reads)/4
write.csv(n_reads, "sim_n_reads.csv", row.names = FALSE)
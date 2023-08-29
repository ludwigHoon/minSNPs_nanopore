# PREV_RESULT_FILE can be downloaded from:
# SSTRING_MINSNPS_200 
# SSTRING_RANDOM_400
# SLURM_TEMPLATE is a slurm bash template, for details, look at: https://rdrr.io/bioc/BiocParallel/f/inst/doc/Introduction_To_BiocParallel.pdf

###### These files need to be updated to reflect the actual directory ####
PREV_RESULT_FILE <- "sup_5_6_mega_hd_2_steps.txt"
SSTRING_MINSNPS_200 <- "search_sequence_200_non_unique_SNPs.csv"
SSTRING_RANDOM_400 <- "search_sequence_400_random_SNPs.csv"
SLURM_TEMPLATE <- "slurm_bc_template.tmpl"
######

library(BiocParallel)
bp <- BatchtoolsParam(workers = 391, cluster="slurm", template=SLURM_TEMPLATE,
    resources=list(walltime=60*60*24*5, ncpus=4))

sim_files <- list.files(pattern = ".*\\.fastq$")
res <- bplapply(sim_files, function(file, SSTRING_RANDOM_400){
    search_table2 <- read.csv(SSTRING_RANDOM_400)
    search_table2 <- search_table2[search_table2$strand == "+", ]
    search_table2$SNP <- sapply(strsplit(search_table2$id, split = "_"), `[`, 1)
    library(minSNPs) 
    library(BiocParallel)

    # Look for 100, 200, 300, 400 SNPs
    SNPs <- unique(search_table2$SNP)
    for (max_snps in c(7,13,18,23,27, 100, 200, 300, 400)){
        search_table <- search_table2[search_table2$SNP %in% SNPs[1:max_snps], ]
        previous_result <- NULL
        for (read_i in c(0, 1, 2, 3, 4, 5, 6, 7)){

            temp_result <- search_from_fastq_reads(fastq_file = file, search_tables = search_table, 
                skip_n_reads = (1000 * read_i), progress = TRUE, max_n_reads = 1000,
                quality_offset = 33, output_temp_result = FALSE, temp_result_folder = "./temp_results",
                simplify_id = TRUE, output_read_length = TRUE, bp = BiocParallel::MulticoreParam(workers = 8)
                )
            combined_result <- combine_fastq_search_result(temp_result, search_table, previous_result, bp = MulticoreParam(workers = 2))
            previous_result <- combined_result
            inferred_result <- infer_from_combined(combined_result, search_table, 2800000)
            tt <- inferred_result$result
            tt$estimated_coverage <- rep(inferred_result$estimated_coverage, nrow(tt))
            tt$reads_used <- rep( (1000*(read_i+1)), nrow(tt))
            write.csv(tt, paste0("random_400/sim_", file, "_", max_snps, "_", read_i, ".csv"), row.names = FALSE)
            }
    }
}, SSTRING_RANDOM_400 =SSTRING_RANDOM_400, BPPARAM=bp)

res <- bplapply(sim_files, function(file, SSTRING_MINSNPS_200){
    search_table2 <- read.csv(SSTRING_MINSNPS_200)
    search_table2 <- search_table2[search_table2$strand == "+", ]
    search_table2$SNP <- sapply(strsplit(search_table2$id, split = "_"), `[`, 1)
    library(minSNPs)
    library(BiocParallel)

    snps_set <- process_result_file(PREV_RESULT_FILE)
    # n_snps =  7  13  18  23  27  99 199 300 385
    for (n_sets in c(1, 2, 3, 4, 5, 30, 79, 147, 200)){
        search_table <- search_table2[search_table2$SNP %in% unique(unlist(snps_set[c(1:n_sets)])), ]
        previous_result <- NULL
        for (read_i in c(0, 1, 2, 3, 4, 5, 6, 7)){
            temp_result <- search_from_fastq_reads(fastq_file = file, search_tables = search_table, 
                skip_n_reads = (1000 * read_i), progress = TRUE, max_n_reads = 1000,
                quality_offset = 33, output_temp_result = FALSE, temp_result_folder = "./temp_results",
                simplify_id = TRUE, output_read_length = TRUE, bp = BiocParallel::MulticoreParam(workers = 8)
                )
            combined_result <- combine_fastq_search_result(temp_result, search_table, previous_result, bp = MulticoreParam(workers = 2))
            previous_result <- combined_result
            inferred_result <- infer_from_combined(combined_result, search_table, 2800000)
            tt <- inferred_result$result
            tt$estimated_coverage <- rep(inferred_result$estimated_coverage, nrow(tt))
            tt$reads_used <- rep( (1000*(read_i+1)), nrow(tt))

            write.csv(tt, paste0("non_random_200/sim_", file, "_", n_sets, "_", read_i, ".csv"), row.names = FALSE)
            }
    }
}, SSTRING_MINSNPS_200= SSTRING_MINSNPS_200, BPPARAM=bp)

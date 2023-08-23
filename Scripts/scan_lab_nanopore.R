##
##
##
##
PREV_RESULT_FILE <- "sup_5_6_mega_hd_2_steps.txt" ##
SSTRING_MINSNPS_200 <- "search_sequence_200_non_unique_SNPs.csv"
SSTRING_RANDOM_400 <- "search_sequence_400_random_SNPs.csv"
SLURM_TEMPLATE <- "slurm_bc_template.tmpl" 
##


library(BiocParallel)
bp <- BatchtoolsParam(workers = 24, cluster="slurm", template=SLURM_TEMPLATE,
    resources=list(walltime=60*60*24*5, ncpus=4))

for (call_type in c("fast", "hac", "sup")){
    sim_files <- list.files(pattern = ".*.fastq.gz$", path = call_type, full.names = TRUE)

    res <- bplapply(sim_files, function(file, SSTRING_RANDOM_400, call_type){
        search_table2 <- read.csv(SSTRING_RANDOM_400)
        search_table2 <- search_table2[search_table2$strand == "+", ]
        search_table2$SNP <- sapply(strsplit(search_table2$id, split = "_"), `[`, 1)
        library(minSNPs) 
        library(BiocParallel)
        library(tictoc)
        id <- gsub(".fastq.gz", "", strsplit(file, split = "_")[[1]][2])
        
        # Look for 100, 200, 300, 400 SNPs
        SNPs <- unique(search_table2$SNP)
        for (max_snps in c(100, 200, 300, 400)){
            search_table <- search_table2[search_table2$SNP %in% SNPs[1:max_snps], ]
            previous_result <- NULL
            for (read_i in c(0, 1, 2, 3, 4, 5, 6, 7)){
                tic(paste0("USING ", max_snps, " SNPs set, ", (1000*(read_i+1)), " reads"))
                temp_result <- search_from_fastq_reads(fastq_file = file, search_tables = search_table, 
                    skip_n_reads = (1000 * read_i), progress = TRUE, max_n_reads = 1000,
                    quality_offset = 33, output_temp_result = FALSE, temp_result_folder = "./temp_results",
                    simplify_id = TRUE, output_read_length = TRUE, bp = BiocParallel::MulticoreParam(workers = 8)
                    )
                if (is.null(temp_result)){
                    break
                }
                combined_result <- combine_fastq_search_result(temp_result, search_table, previous_result, bp = MulticoreParam(workers = 2))
                previous_result <- combined_result
                inferred_result <- infer_from_combined(combined_result, search_table, 2800000)
                tt <- inferred_result$result
                tt$estimated_coverage <- rep(inferred_result$estimated_coverage, nrow(tt))
                tt$max_reads_used <- rep( (1000*(read_i+1)), nrow(tt))
                toc()
                write.csv(tt, paste0("results/random_400/", call_type, "_", id, "_", max_snps, "_", read_i, ".csv"), row.names = FALSE)
                }
        }
    }, call_type = call_type, SSTRING_RANDOM_400 = SSTRING_RANDOM_400, BPPARAM=bp)

    res <- bplapply(sim_files, function(file, SSTRING_MINSNPS_200, PREV_RESULT_FILE, call_type){
        search_table2 <- read.csv(SSTRING_MINSNPS_200)
        search_table2 <- search_table2[search_table2$strand == "+", ]
        search_table2$SNP <- sapply(strsplit(search_table2$id, split = "_"), `[`, 1)
        library(minSNPs)
        library(BiocParallel)
        library(tictoc)
        id <- gsub(".fastq.gz", "", strsplit(file, split = "_")[[1]][2])
        snps_set <- process_result_file(PREV_RESULT_FILE)
        # n_snps =  7  13  18  23  27  99 199 300 385
        for (n_sets in c(1, 2, 3, 4, 5, 30, 79, 147, 200)){
            search_table <- search_table2[search_table2$SNP %in% unique(unlist(snps_set[c(1:n_sets)])), ]
            previous_result <- NULL
            for (read_i in c(0, 1, 2, 3, 4, 5, 6, 7)){
                tic(paste0("USING ", n_sets, " SNPs set, ", (1000*(read_i+1)), " reads"))
                temp_result <- search_from_fastq_reads(fastq_file = file, search_tables = search_table, 
                    skip_n_reads = (1000 * read_i), progress = TRUE, max_n_reads = 1000,
                    quality_offset = 33, output_temp_result = FALSE, temp_result_folder = "./temp_results",
                    simplify_id = TRUE, output_read_length = TRUE, bp = BiocParallel::MulticoreParam(workers = 8)
                    )
                if (is.null(temp_result)){
                    break
                }
                combined_result <- combine_fastq_search_result(temp_result, search_table, previous_result, bp = MulticoreParam(workers = 2))
                previous_result <- combined_result
                inferred_result <- infer_from_combined(combined_result, search_table, 2800000)
                tt <- inferred_result$result
                tt$estimated_coverage <- rep(inferred_result$estimated_coverage, nrow(tt))
                tt$max_reads_used <- rep( (1000*(read_i+1)), nrow(tt))
                toc()
                write.csv(tt, paste0("results/non_random_200/", call_type, "_", id, "_", n_sets, "_", read_i, ".csv"), row.names = FALSE)
                }
        }
    }, call_type = call_type, SSTRING_MINSNPS_200 = SSTRING_MINSNPS_200, PREV_RESULT_FILE = PREV_RESULT_FILE, BPPARAM=bp)
}
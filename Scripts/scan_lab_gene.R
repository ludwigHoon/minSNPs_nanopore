## Scan simulated reads for gene
GENE_SEARCH_SEQUENCE <- "gene_search_sequence.csv"
SLURM_TEMPLATE <- "slurm_bc_template.tmpl"

library(BiocParallel)
bp <- BatchtoolsParam(workers = 24, cluster="slurm", template=SLURM_TEMPLATE,
    resources=list(walltime=60*60*24*5, ncpus=2))

for (call_type in c("fast", "hac", "sup")){
    sim_files <- list.files(pattern = ".*.fastq.gz$", path = call_type, full.names = TRUE)
    res <- bplapply(sim_files, function(file, GENE_SEARCH_SEQUENCE, call_type){
        library(minSNPs)
        library(BiocParallel)
        set.seed(42)
        id <- gsub(".fastq.gz", "", strsplit(file, split = "_")[[1]][2])
        search_table2 <- read.csv(GENE_SEARCH_SEQUENCE)
        search_table2 <- search_table2[search_table2$strand == "+", ]
        search_table2$gene <- sapply(strsplit(search_table2$id, split = "_"), `[`, 1)
        
        for (max_sequence in c(100, 200, 300, 400, 500, 600, 700, 800, 900, "all")){
            s_table <- list()
            for (gene in unique(search_table2$gene)){
                partial_table <- search_table2[search_table2$gene == gene, ]
                if (max_sequence == "all"){
                    n_sample <- nrow(partial_table)
                } else{
                    n_sample <- max_sequence
                }
                sampled <- sample(seq_len(nrow(partial_table)), n_sample)
                s_table[[gene]] <- partial_table[sampled, ]
            }
            search_table <- do.call(rbind, s_table)
            rownames(search_table) <- NULL
            previous_result <- NULL
            for (read_i in c(0, 1, 2, 3, 4, 5, 6, 7)){
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
                write.csv(tt, paste0("gene/sim_", call_type, "_", id, "_", max_sequence, "_", read_i, ".csv"), row.names = FALSE)
            }
        }
    }, call_type = call_type, GENE_SEARCH_SEQUENCE = GENE_SEARCH_SEQUENCE, BPPARAM=bp)
}
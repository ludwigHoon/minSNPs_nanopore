###################################################### Parameters requiring adjustment ######################################################
# @title *S. aureus* CC inference & gene detection from FASTQ
###################################################### Parameters requiring adjustment ######################################################
# @title Parameters
# @markdown Number of SNP sets to use
n_sets <- "200" # @param [5, 10, 30, 50, 100, 150, 200]
n_sets <- as.numeric(n_sets)
# @markdown Number of sequence to use per gene
n_seq_gene <- "all" # @param [100, 200, "all"]
if (n_seq_gene != "all"){
  n_seq_gene <- as.numeric(n_seq_gene)
}
# @markdown Proportion of gene search sequences that must be found to a gene positive
threshold <- 0.6 #@param [0.1, 0.15, 0.2, 0.6, 0.7, 0.8, 0.9 ]
# @markdown Number of reads in fastq to use
# @markdown -1 means use all the reads in fastq file
n_reads <- "4000" # @param [1000, 2000, 3000, 4000, -1]
n_reads <- as.numeric(n_reads)
# @markdown Fastq file name
fastq_file <- "nanopore_seq.fastq.gz" # @param {type:"string"}
# @markdown Where to save the CC inference result
output_CC_RES <- "sa_CC_result.csv" # @param {type:"string"}
# @markdown Where to save the gene detection result
output_gene_RES <- "sa_GENE_result.csv" # @param {type:"string"}
# @markdown Where to save the estimated coverage result
estimated_coverage <- "sa_estimated_coverage.txt" # @param {type:"string"}




############################################## Only touch below if you know what you are doing ##############################################
# @title *S. aureus* CC inference & gene detection from FASTQ
############################################## Only touch below if you know what you are doing ##############################################
if (!file.exists("gene_seq.csv")){
  print("Downloading gene reference library")
  download.file("https://raw.githubusercontent.com/ludwigHoon/minSNPs_nanopore/refs/heads/main/Results/gene_search_sequence.csv", "gene_seq.csv")
}
if (!file.exists("snp_seq.csv")){
  print("Downloading major lineage reference library")
  download.file("https://raw.githubusercontent.com/ludwigHoon/minSNPs_nanopore/refs/heads/main/Results/CC_search_sequence_200_non_unique_SNPs.csv", "snp_seq.csv")
}
if (!file.exists("snp_sets.csv")){
  download.file("https://raw.githubusercontent.com/ludwigHoon/minSNPs_nanopore/refs/heads/main/Data/sup_5_7_mega_hd_2_steps.txt", "snp_sets.csv")
}

if (!require("minSNPs", quietly = TRUE)){
  if (!require("BiocParallel", quietly = TRUE)){
    if (!require("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    BiocManager::install("BiocParallel")
  }
  install.packages("minSNPs")
}
# Libraries
library(minSNPs)
library(BiocParallel)

SNP_search_table <- read.csv("snp_seq.csv")
SNP_search_table <- SNP_search_table[SNP_search_table$strand == "+", ]
SNP_search_table$SNP <- sapply(strsplit(SNP_search_table$id, split = "_"), `[`, 1)
snps_set <- process_result_file("snp_sets.csv")
search_table <- SNP_search_table[SNP_search_table$SNP %in% unique(unlist(snps_set[c(1:n_sets)])), ]


GENE_search_table <- read.csv("gene_seq.csv")
GENE_search_table <- GENE_search_table[GENE_search_table$strand == "+", ]
GENE_search_table$gene <- sapply(strsplit(GENE_search_table$id, split = "_"), `[`, 1)
s_table <- list()
for (gene in unique(GENE_search_table$gene)){
    partial_table <- GENE_search_table[GENE_search_table$gene == gene, ]
    if (n_seq_gene == "all"){
        n_sample <- nrow(partial_table)
    } else{
        n_sample <- n_seq_gene
    }
    sampled <- sample(seq_len(nrow(partial_table)), n_sample)
    s_table[[gene]] <- partial_table[sampled, ]
}
search_table2 <- do.call(rbind, s_table)
rownames(search_table2) <- NULL

final_search_table <- rbind(search_table[,1:6], search_table2[,1:6])
previous_result <- NULL

# gene_result <- search_from_fastq_reads(fastq_file = fastq_file, search_tables = search_table2,
#       skip_n_reads = 0, progress = TRUE, max_n_reads = n_reads,
#       quality_offset = 33, output_temp_result = FALSE, temp_result_folder = "./temp_results",
#       simplify_id = TRUE, output_read_length = TRUE, bp = BiocParallel::MulticoreParam(workers = 2)
#       )

# snp_result <- search_from_fastq_reads(fastq_file = fastq_file, search_tables = search_table,
#       skip_n_reads = 0, progress = TRUE, max_n_reads = n_reads,
#       quality_offset = 33, output_temp_result = FALSE, temp_result_folder = "./temp_results",
#       simplify_id = TRUE, output_read_length = TRUE, bp = BiocParallel::MulticoreParam(workers = 2)
#       )

full_result <- search_from_fastq_reads(fastq_file = fastq_file, search_tables = final_search_table,
      skip_n_reads = 0, progress = TRUE, max_n_reads = n_reads,
      quality_offset = 33, output_temp_result = FALSE, temp_result_folder = "./temp_results",
      simplify_id = TRUE, output_read_length = TRUE, bp = BiocParallel::MulticoreParam(workers = 2)
      )

body(infer_from_combined)[[4]] <- substitute(if (!inherits(combined_result, "combined_fastq_search_result")) {
    if (!all(c("result", "read_length") %in% names(combined_result))) {
        stop("combined_result must be a combined_fastq_search_result object")
    }
    if (!inherits(combined_result$result, "data.frame")) {
        stop("combined_result must be a combined_fastq_search_result object")
    }
    else {
        class(combined_result) <- "combined_fastq_search_result"
    }
})

process_snp_result.old <- process_snp_result
process_snp_result.new <- process_snp_result
body(process_snp_result.new)[[21]] <- substitute(
    if (length(read_count) == 0) {
        result_df <- data.frame(type = rep("SNP", 1), rank = c(1), result = NA, reads_count = NA,
        proportion_matched = NA,
        pass_filter = NA,
        proportion_scheme_found = NA,
        details = paste0("No SNPs found, removed due to conflict"))
    } else {
        result_df <- data.frame(type = rep("SNP", length(read_count)), rank = seq_len(max(1, length(read_count))), result = names(cc_result), reads_count = unname(r_count),
        proportion_matched = c(prop),
        pass_filter = (
            prop >= 0.8 &
            read_count >= 10 * prop * length(unique(snp_id)) # average read-depth of the reads with respective SNPs = 10
        ),
        proportion_scheme_found = length(unique(snp_id)) / length(unique(searched_snps_id)),
        details = paste0("SNPs found: ", paste0(snp_id, collapse = ","), "\n")
    )
    })
process_kmer_result.new <- function(partial_result, search_table, min_match_per_read = 1, ...) {
    result <- list()
    kmer_only <- partial_result
    data.table::setDF(kmer_only)
    for (gene in unique(search_table[search_table$type == "KMER","result"])) {
        total_searched <- nrow(search_table[search_table$type == "KMER" & search_table$result == gene,])

        # Discard reads with less than min_match_per_read
        all_reads <- unlist(strsplit(unlist(kmer_only[kmer_only$result == gene, "reads"]), split = ";|,"))
        if (is.null(all_reads)){
          r_count <- 0
          prop <- 0.0
        } else{
          n_match_in_read <- table(all_reads)
          filtered_reads <- names(n_match_in_read[n_match_in_read >= min_match_per_read])
          filtered_index <- apply(as.matrix(sapply(filtered_reads, grepl, x= kmer_only$reads)), 1, sum) > 0
          kmer_only2 <- kmer_only[kmer_only$result == gene & filtered_index, ]

          prop <- length(unique(kmer_only2[kmer_only2$result == gene, "sequence"])) / total_searched
          r_count <- length(unique(unlist(strsplit(unlist(kmer_only2[kmer_only2$result == gene, "reads"]), split = ";|,"))))
        }

        result[[gene]] <- data.frame(type = "KMER", rank = 1, result = gene, reads_count = r_count,
            proportion_matched = prop, pass_filter = prop >= threshold, proportion_scheme_found = prop, details = NA)
    }
    result_df <- data.table::rbindlist(result)
    return(
        list(result = result_df))
}

MinSNPs_process_methods <- list( #nolint
          "SNP" = process_snp_result.new,
          "KMER" = process_kmer_result.new
      )

combined_full_result <- combine_fastq_search_result(full_result, final_search_table, NULL, bp = MulticoreParam(workers = 2))
inferred_full_result <- infer_from_combined(combined_full_result, final_search_table, 2800000)

# combined_gene_result <- combine_fastq_search_result(gene_result, search_table2, NULL, bp = MulticoreParam(workers = 2))
# inferred_gene_result <- infer_from_combined(combined_gene_result, search_table2, 2800000)

# combined_snp_result <- combine_fastq_search_result(snp_result, search_table, NULL, bp = MulticoreParam(workers = 2))
# inferred_snp_result <- infer_from_combined(combined_snp_result, search_table, 2800000)

write.csv(inferred_full_result$result[type == "SNP"], output_CC_RES, row.names = FALSE)
write.csv(inferred_full_result$result[type == "KMER"], output_gene_RES, row.names = FALSE)
cat(inferred_full_result$estimated_coverage, file = estimated_coverage)

print("ALL DONE, result saved to file")

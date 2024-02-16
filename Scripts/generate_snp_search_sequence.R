library(minSNPs)
library(data.table)

PREV_RESULT_FILE <- "sup_5_6_mega_hd_2_steps.txt"

ref_meta <- read.csv("ref_OUTF1.csv")
snp_matrix <- read_fasta("OUTF1.fasta")
ref_seq <- read_fasta("Mu50.fasta")

random_snps <- as.numeric(readLines("random_400_SNPs.txt"))
minsnps_snps <- unique(unlist(process_result_file(PREV_RESULT_FILE)))


random_snp_search_string <- generate_snp_search_string(
    random_snps, ref_meta, ref_seq,
            snp_matrix, 7, 7, position_type = "fasta",
            extend_length = TRUE, fasta_name_as_result = TRUE)
minsnps_snp_search_string <- generate_snp_search_string(
    minsnps_snps, ref_meta, ref_seq,
            snp_matrix, 7, 7, position_type = "fasta",
            extend_length = TRUE, fasta_name_as_result = TRUE)

fwrite(random_snp_search_string, "search_sequence_400_random_SNPs.csv", row.names = FALSE)
fwrite(minsnps_snp_search_string, "search_sequence_200_non_unique_SNPs.csv", row.names = FALSE)
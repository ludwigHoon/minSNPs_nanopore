# Random SNPs generation

# Megaalignment contains 164335 SNPs
all_SNPs <- c(1:164335)
set.seed(42) # for reproducibliity
reordered_SNPs <- sample(all_SNPs, length(all_SNPs))
writeLines(paste0(reordered_SNPs[1:400], collapse = "\n"), "random_400_SNPs.txt")
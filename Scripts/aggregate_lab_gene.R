# Aggregate gene result for lab generated data
library(data.table)
library(BiocParallel)

summarised_to_meta <- fread("Nanopore_run_meta.csv")
summarised_to_meta$Barcode.ID <- sprintf("%02d", summarised_to_meta$Barcode.ID)
summarised_to_meta$MRSA <- ifelse(summarised_to_meta$MRSA == "+", TRUE, FALSE)
summarised_to_meta$PVL <- ifelse(summarised_to_meta$PVL == "+", TRUE, FALSE)

all_gene_result <- list.files(path = "gene", pattern = "^sim.*.csv", full.names = TRUE)
result <- bplapply(all_gene_result, function(file){
    details <- strsplit(basename(file), split = "_")
    basecall_method <- sapply(details, `[`, 2)
    id <- sapply(details, `[`, 3)
    n_search_sequence <- sapply(details, `[`, 4)

    sim_res_ind <- fread(file)
    sim_res_ind$id <- id
    rel_row <- summarised_to_meta[eval(summarised_to_meta$Barcode.ID == id)]
    rel_row$estimated_coverage <- unique(sim_res_ind$estimated_coverage)
    rel_row$max_reads_used <- unique(sim_res_ind$max_reads_used)
    mecA_result <- sim_res_ind[result == "mecA"]
    lukS_result <- sim_res_ind[result == "lukS-PV"]
    lukF_result <- sim_res_ind[result == "lukF-PV"]
    rel_row$mecA_rcount <- mecA_result$reads_count
    rel_row$mecA_prop <- mecA_result$proportion_matched
    rel_row$lukS_rcount <- lukS_result$reads_count
    rel_row$lukS_prop <- lukS_result$proportion_matched
    rel_row$lukF_rcount <- lukF_result$reads_count
    rel_row$lukF_prop <- lukF_result$proportion_matched
    rel_row$basecall_method <- basecall_method
    rel_row$n_sequence_used <- n_search_sequence
    return(rel_row)
}, BPPARAM = MulticoreParam(workers = 4, progress = TRUE))

all_result <- rbindlist(result)

fwrite(all_result, "lab_gene_result.csv", row.names = FALSE)

summary <- list()
for (threshold in c(0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)){
    summary[[as.character(threshold)]] <- all_result[,
    list(mean_estimated_coverage = mean(estimated_coverage),
        mecA_TP = length(which(mecA_prop >= threshold & MRSA == TRUE)),
        mecA_FP = length(which(mecA_prop >= threshold & MRSA == FALSE)),
        mecA_TN = length(which(mecA_prop < threshold & MRSA == FALSE)),
        mecA_FN = length(which(mecA_prop < threshold & MRSA == TRUE)),
        lukF_TP = length(which(lukF_prop >= threshold & PVL == TRUE)),
        lukF_FP = length(which(lukF_prop >= threshold & PVL == FALSE)),
        lukF_TN = length(which(lukF_prop < threshold & PVL == FALSE)),
        lukF_FN = length(which(lukF_prop < threshold & PVL == TRUE)),
        lukS_TP = length(which(lukS_prop >= threshold & PVL == TRUE)),
        lukS_FP = length(which(lukS_prop >= threshold & PVL == FALSE)),
        lukS_TN = length(which(lukS_prop < threshold & PVL == FALSE)),
        lukS_FN = length(which(lukS_prop < threshold & PVL == TRUE))
        ),
    by= list(basecall_method, n_sequence_used, max_reads_used)][,
    list(
        threshold = threshold,
        mecA_TP, mecA_FP, mecA_TN, mecA_FN,
        lukF_TP, lukF_FP, lukF_TN, lukF_FN,
        lukS_TP, lukS_FP, lukS_TN, lukS_FN,
        mecA_sensitivity = mecA_TP/(mecA_TP + mecA_FN),
        mecA_specificity = mecA_TN/(mecA_TN + mecA_FP),
        lukF_sensitivity = lukF_TP/(lukF_TP + lukF_FN),
        lukF_specificity = lukF_TN/(lukF_TN + lukF_FP),
        lukS_sensitivity = lukS_TP/(lukS_TP + lukS_FN),
        lukS_specificity = lukS_TN/(lukS_TN + lukS_FP),
        mean_estimated_coverage,
        basecall_method,
        n_sequence_used, max_reads_used
        )]
}

all_summary <- rbindlist(summary)
fwrite(all_summary, "lab_gene_summary.csv", row.names = FALSE)
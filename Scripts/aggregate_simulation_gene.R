# aggregate simulated result for gene
#
# The presence/absence of the interested genes are first tested with AMRFINDER
# for a in `ls *.fas`
# do
# amrfinder --plus -n $a -O Staphylococcus_aureus > $a.got
# done
#
# The result is then aggregated with this script
library(data.table)
library(BiocParallel)

all_amrfinder_results <- list.files(pattern = "*.got")
seq_name <- gsub(".fas.got", "", all_amrfinder_results)

result <- bplapply(seq_len(length(all_amrfinder_results)), function(i, all_amrfinder_results, seq_name){
    interested_columns <- c("Gene symbol", "% Coverage of reference sequence", "% Identity to reference sequence")
    temp_result <- fread(all_amrfinder_results[i], select = interested_columns)
    interested_subset <- temp_result[`Gene symbol` %in% c("mecA", "lukF-PV", "lukS-PV")]
    interested_subset$seq_name <- seq_name[i]
    return(interested_subset)
}, all_amrfinder_results = all_amrfinder_results, seq_name = seq_name, BPPARAM = MulticoreParam(workers = 4))

all_detected_from_amrfinder <- rbindlist(result)
fwrite(all_detected_from_amrfinder, "all_detected_from_amrfinder.csv", row.names = FALSE)

summary <- bplapply(seq_name, function(name, result){
    partial_result <- result[seq_name == name]
    return(data.frame(seq_name = name, mecA = any(partial_result$`Gene symbol` == "mecA"),
            lukF_PV = any(partial_result$`Gene symbol` == "lukF-PV"),
            lukS_PV = any(partial_result$`Gene symbol` == "lukS-PV"), stringsAsFactors = FALSE)
    )
}, result = all_detected_from_amrfinder)

summarised_to_meta <- rbindlist(summary)
fwrite(summarised_to_meta, "gene_summarised_to_meta.csv", row.names = FALSE)


#######################################################################
# The result of the scanned reads are then aggregated with this script
#######################################################################

summarised_to_meta <- fread("gene_summarised_to_meta.csv")
summarised_to_meta$id <- sapply(strsplit(summarised_to_meta$seq_name, split = "_"), `[`, 1)
result <- list()

all_simulated_result <- list.files(pattern = "^sim.*.csv")
result <- bplapply(all_simulated_result, function(file){
    details <- strsplit(basename(file), split = "_")
    id <- gsub("S", "", sapply(details, `[`, 2))
    n_search_sequence <- sapply(details, `[`, 4)

    sim_res_ind <- fread(file)
    sim_res_ind$id <- id
    rel_row <- summarised_to_meta[eval(summarised_to_meta$id == id)]
    rel_row$estimated_coverage <- unique(sim_res_ind$estimated_coverage)
    rel_row$reads_used <- unique(sim_res_ind$reads_used)
    mecA_result <- sim_res_ind[result == "mecA"]
    lukS_result <- sim_res_ind[result == "lukS-PV"]
    lukF_result <- sim_res_ind[result == "lukF-PV"]
    rel_row$mecA_rcount <- mecA_result$reads_count
    rel_row$mecA_prop <- mecA_result$proportion_matched
    rel_row$lukS_rcount <- lukS_result$reads_count
    rel_row$lukS_prop <- lukS_result$proportion_matched
    rel_row$lukF_rcount <- lukF_result$reads_count
    rel_row$lukF_prop <- lukF_result$proportion_matched

    rel_row$n_sequence_used <- n_search_sequence
    return(rel_row)
}, BPPARAM = MulticoreParam(workers = 32))

all_result <- rbindlist(result)

all_result[n_sequence_used == "all" & reads_used == 8000,
    list(min_mecA = min(mecA_prop), max_mecA = max(mecA_prop), mean_mecA=mean(mecA_prop)),
    by= mecA]
all_result[n_sequence_used == "all" & reads_used == 8000,
    list(min_lukS = min(lukS_prop), max_lukS = max(lukS_prop), mean_lukS=mean(lukS_prop)),
    by= lukS_PV]
all_result[n_sequence_used == "all" & reads_used == 8000,
    list(min_lukF = min(lukF_prop), max_lukF = max(lukF_prop), mean_lukF=mean(lukF_prop)),
    by= lukF_PV]

all_result[n_sequence_used == "all" & reads_used == 8000,
    list(min_mecA = min(mecA_rcount), max_mecA = max(mecA_rcount), mean_mecA=mean(mecA_rcount)),
    by= mecA]
all_result[n_sequence_used == "all" & reads_used == 8000,
    list(min_lukS = min(lukS_rcount), max_lukS = max(lukS_rcount), mean_lukS=mean(lukS_rcount)),
    by= lukS_PV]
all_result[n_sequence_used == "all" & reads_used == 8000,
    list(min_lukF = min(lukF_rcount), max_lukF = max(lukF_rcount), mean_lukF=mean(lukF_rcount)),
    by= lukF_PV]

fwrite(all_result, "simulated_gene_result.csv", row.names = FALSE)

summary <- all_result[,
    list(mean_estimated_coverage = mean(estimated_coverage),
        mecA_TP = length(which(mecA_prop >= 0.8 & mecA == TRUE)),
        mecA_FP = length(which(mecA_prop >= 0.8 & mecA == FALSE)),
        mecA_TN = length(which(mecA_prop < 0.8 & mecA == FALSE)),
        mecA_FN = length(which(mecA_prop < 0.8 & mecA == TRUE)),
        lukF_TP = length(which(lukF_prop >= 0.8 & lukF_PV == TRUE)),
        lukF_FP = length(which(lukF_prop >= 0.8 & lukF_PV == FALSE)),
        lukF_TN = length(which(lukF_prop < 0.8 & lukF_PV == FALSE)),
        lukF_FN = length(which(lukF_prop < 0.8 & lukF_PV == TRUE)),
        lukS_TP = length(which(lukS_prop >= 0.8 & lukS_PV == TRUE)),
        lukS_FP = length(which(lukS_prop >= 0.8 & lukS_PV == FALSE)),
        lukS_TN = length(which(lukS_prop < 0.8 & lukS_PV == FALSE)),
        lukS_FN = length(which(lukS_prop < 0.8 & lukS_PV == TRUE))
        ),
    by= list(n_sequence_used, reads_used)][,
    list(mecA_sensitivity = mecA_TP/(mecA_TP + mecA_FN),
        mecA_specificity = mecA_TN/(mecA_TN + mecA_FP),
        lukF_sensitivity = lukF_TP/(lukF_TP + lukF_FN),
        lukF_specificity = lukF_TN/(lukF_TN + lukF_FP),
        lukS_sensitivity = lukS_TP/(lukS_TP + lukS_FN),
        lukS_specificity = lukS_TN/(lukS_TN + lukS_FP),
        mean_estimated_coverage, 
        n_sequence_used, reads_used
        )]
fwrite(summary, "simulated_gene_summary.csv", row.names = FALSE)
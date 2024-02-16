### From most similar isolate in alignment to most likely CC based on SNP's differences
###
### FILLED_META_FILE: the metadata file in Results
###
###

FILLED_META_FILE <- "all_meta_filled.csv"

library(data.table)
library(BiocParallel)
meta <- fread(FILLED_META_FILE)[is.na(pubmlst_id)][, c("Isolate.name", "MLST.CC")]

fin_res <- bplapply(c("random_400", "non_random_200"), function(dir, meta){
    library(BiocParallel)
    library(data.table)
    
    result <- list()

    all_simulated_result <- list.files(path = dir, pattern = "^sim_S.*.csv", full.names = TRUE)
    result <- bplapply(all_simulated_result, function(file){
        details <- strsplit(basename(file), split = "_")
        id <- gsub("S", "", sapply(details, `[`, 2))
        if (dir == "random_400"){
            snps_used <- sapply(details, `[`, 4)
        } else{
            n_sets <- c(1, 2, 3, 4, 5, 30, 79, 147, 200)
            n_snps_set <- c(7, 13, 18, 23, 27, 99, 199, 300, 385)
            snps_used <-  sapply(details, function(x) {n_snps_set[which(n_sets == as.numeric(x[4]))]})
        }
        sim_res_ind <- fread(file)
        temp_result <- merge(sim_res_ind, meta, by.x = "result", by.y = "Isolate.name")
        temp_result <- temp_result[temp_result[, .I[which.max(proportion_matched)], by=list(MLST.CC)]$V1]
        temp_result$fin_score <- temp_result$proportion_matched * temp_result$proportion_scheme_found
        temp_result <- temp_result[order(fin_score, decreasing = TRUE), list(MLST.CC, fin_score, proportion_matched, proportion_scheme_found, estimated_coverage, reads_used, reads_count)]
        temp_result$snps_used <- snps_used
        temp_result$snps_type <- dir
        temp_result$pubmlst_id <- id
        return(temp_result)
    }, BPPARAM = MulticoreParam(workers = 32))
    return(rbindlist(result))
}, meta = meta, BPPARAM = BatchtoolsParam(workers = 2, cluster="slurm", template="~/slurm_bc_template.tmpl",
        resources=list(walltime=60*60*24*5, ncpus=8)))
all_result <- rbindlist(fin_res)

colnames(all_result)[which(colnames(all_result) == "MLST.CC")] <- "predicted_MLST.CC"
fwrite(all_result, "simulation_aggregated_result.csv", row.names = FALSE)

to_collapse <- all_result[all_result[, .I[fin_score == max(fin_score)], by = list(pubmlst_id, snps_type, snps_used, reads_used)]$V1]
collapsed_top_result <- to_collapse[,.(predicted_MLST.CC = paste0(predicted_MLST.CC, collapse = ","), reads_w_result = paste0(unique(reads_count), collapse = ",")), by = list(pubmlst_id, snps_type, snps_used, reads_used, fin_score, proportion_matched, proportion_scheme_found, estimated_coverage)]

newdata_meta <- fread(FILLED_META_FILE)[!is.na(pubmlst_id)][, c("pubmlst_id", "MLST.CC", "CC_filled", "Isolate.name")]
newdata_meta$pubmlst_id <- as.character(newdata_meta$pubmlst_id)
collapsed_top_result$pubmlst_id <- as.character(collapsed_top_result$pubmlst_id)
RES <- merge(collapsed_top_result, newdata_meta, by = "pubmlst_id")
RES$predicted_mlst_is_single <- sapply(RES$predicted_MLST.CC, function(x) {!grepl(",", x)})
RES$predicted_mlst_is_correct <- sapply(seq_len(nrow(RES)), function(i){
    if (RES$predicted_mlst_is_single[i] == TRUE){
        return(RES$predicted_MLST.CC[i] == RES$CC_filled[i])
    } else{
        return(any(unlist(strsplit(RES$predicted_MLST.CC[i], ",")) == RES$CC_filled[i]))
    }
})
fwrite(RES, "simulation_individual_results.csv", row.names = FALSE)

summary <- RES[,.(n_correct = length(which(predicted_mlst_is_correct)),
                n_single = length(which(predicted_mlst_is_single)),
                n_single_correct = length(which(predicted_mlst_is_single & predicted_mlst_is_correct))),
            by = list(reads_used, snps_used, snps_type)]
fwrite(summary, "simulation_result_summary.csv", row.names = FALSE)

# Estimated coverage stat
RES[snps_type == "random_400" & snps_used == 400 ][,.(reads_used, estimated_coverage)][,.(average_estimated_coverage = mean(estimated_coverage)), by = list(reads_used)]

## Also extract the predicted MLST using all 385 SNPs with 8000 reads
predicted_385_8000 <- RES[snps_used == 385 & reads_used == 8000, c("Isolate.name", "predicted_MLST.CC", "fin_score", "predicted_mlst_is_correct", "predicted_mlst_is_single")]
fwrite(predicted_385_8000, "simulated_385_8000.csv", row.names = FALSE)

## Differences in confidence score between correctly & incorrectly predicted samples
predicted_385_8000[,list(min_fscore = min(fin_score), max_fscore = max(fin_score), mean_fscore = mean(fin_score)), list(predicted_mlst_is_correct)]
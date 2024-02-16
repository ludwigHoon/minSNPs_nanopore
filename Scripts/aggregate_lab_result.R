# 
FILLED_META_FILE <- "all_meta_filled.csv"
library(data.table)
library(BiocParallel)
meta <- fread(FILLED_META_FILE)[is.na(pubmlst_id)][, c("Isolate.name", "MLST.CC")]

all_result_all_ctype <- list()
for (call_type in c("fast", "hac", "sup")){
    fin_res <- bplapply(c("random_400", "non_random_200"), function(dir, meta, call_type){
        library(BiocParallel)
        library(data.table)
        
        result <- list()

        all_lab_result <- list.files(path = paste0("results/", dir), pattern = paste0("^", call_type, ".*.csv"), full.names = TRUE)
        result <- bplapply(all_lab_result, function(file){
            details <- strsplit(basename(file), split = "_")
            barcode_id <- as.numeric(sapply(details, `[`, 2))
            if (dir == "random_400"){
                snps_used <- sapply(details, `[`, 3)
            } else{
                n_sets <- c(1, 2, 3, 4, 5, 30, 79, 147, 200)
                n_snps_set <- c(7, 13, 18, 23, 27, 99, 199, 300, 385)
                snps_used <-  sapply(details, function(x) {n_snps_set[which(n_sets == as.numeric(x[3]))]})
            }
            sim_res_ind <- fread(file)
            temp_result <- merge(sim_res_ind, meta, by.x = "result", by.y = "Isolate.name")
            temp_result <- temp_result[temp_result[, .I[which.max(proportion_matched)], by=list(MLST.CC)]$V1]
            temp_result$fin_score <- temp_result$proportion_matched * temp_result$proportion_scheme_found
            temp_result <- temp_result[order(fin_score, decreasing = TRUE), list(MLST.CC, fin_score, proportion_matched, proportion_scheme_found, estimated_coverage, max_reads_used, reads_count)]
            temp_result$snps_used <- snps_used
            temp_result$snps_type <- dir
            temp_result$barcode_id <- barcode_id
            temp_result$call_type <- call_type
            return(temp_result)
        }, BPPARAM = MulticoreParam(workers = 8))
        return(rbindlist(result))
    }, call_type = call_type, meta = meta, BPPARAM = BatchtoolsParam(workers = 2, cluster="slurm", template="~/slurm_bc_template.tmpl",
            resources=list(walltime=60*60*24*5, ncpus=4)))
    all_result_all_ctype[[call_type]] <- rbindlist(fin_res)
}
all_result <- rbindlist(all_result_all_ctype)
colnames(all_result)[which(colnames(all_result) == "MLST.CC")] <- "predicted_MLST.CC"
fwrite(all_result, "lab_aggregated_result.csv", row.names = FALSE)

to_collapse <- all_result[all_result[, .I[fin_score == max(fin_score)], by = list(barcode_id, snps_type, snps_used, max_reads_used, call_type)]$V1]
collapsed_top_result <- to_collapse[,.(predicted_MLST.CC = paste0(predicted_MLST.CC, collapse = ","), reads_w_result = paste0(reads_count, collapse = ",")), by = list(barcode_id, snps_type, snps_used, max_reads_used, call_type, fin_score, proportion_matched, proportion_scheme_found, estimated_coverage)]
collapsed_top_result$barcode_id <- as.numeric(collapsed_top_result$barcode_id)

nanopore_run_meta <- fread("Nanopore_run_meta.csv")
nanopore_run_meta$CC <- paste0("CC", nanopore_run_meta$CC)
nanopore_run_meta$Barcode.ID <- as.numeric(nanopore_run_meta$Barcode.ID)

RES <- merge(collapsed_top_result, nanopore_run_meta, by.x = "barcode_id", by.y = "Barcode.ID")
RES$predicted_mlst_is_single <- sapply(RES$predicted_MLST.CC, function(x) {!grepl(",", x)})
RES$predicted_mlst_is_correct <- sapply(seq_len(nrow(RES)), function(i){
    if (RES$predicted_mlst_is_single[i] == TRUE){
        return(RES$predicted_MLST.CC[i] == RES$CC[i])
    } else{
        return(any(unlist(strsplit(RES$predicted_MLST.CC[i], ",")) == RES$CC[i]))
    }
})

# for samples with less reads, use result from last reads
barcode_with_less_reads <- names(which(table(RES$barcode_id) < 432))
for (barc in barcode_with_less_reads){
    for (ctype in c("fast", "hac", "sup")){
        for (snps_type in c("random_400", "non_random_200")){
            for (snps_used in RES[barcode_id == barc & call_type == ctype & snps_type == snps_type]$snps_used){
                last_reads <- max(RES[barcode_id == barc & call_type == ctype & snps_type == snps_type & snps_used == snps_used]$max_reads_used)
                last_result <- RES[
                    barcode_id == barc &
                    call_type == ctype &
                    snps_type == snps_type &
                    snps_used == snps_used &
                    max_reads_used == last_reads]
                while(last_reads != 8000){
                    dup_last_result <- last_result
                    dup_last_result$max_reads_used <- last_reads + 1000
                    RES <- rbind(RES, dup_last_result)
                    last_reads <- last_reads + 1000
                }
            }
        }
    }
}
RES[order(snps_type, snps_used, call_type, barcode_id, max_reads_used)]
fwrite(RES, "lab_individual_results.csv", row.names = FALSE)

# Estimated coverage stat
RES[snps_type == "random_400" & snps_used == 400 & (! barcode_id %in% c(33, 73))][,.(max_reads_used, call_type, estimated_coverage)][,.(average_estimated_coverage = mean(estimated_coverage)), by = list(call_type,
 max_reads_used)]

summary <- RES[,.(n_correct = length(which(predicted_mlst_is_correct)),
                n_single = length(which(predicted_mlst_is_single)),
                n_single_correct = length(which(predicted_mlst_is_single & predicted_mlst_is_correct))),
            by = list(max_reads_used, snps_used, snps_type, call_type)]
fwrite(summary, "lab_result_summary.csv", row.names = FALSE)
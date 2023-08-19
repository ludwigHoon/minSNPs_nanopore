library(data.table)
fin_res <- bplapply(c("random_400", "non_random_200"), function(dir){
    library(BiocParallel)
    library(data.table)
    meta <- fread("../step_1_search_string/meta-phylogeny-profiles.csv")[,c(1,2)]
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
}, BPPARAM = BatchtoolsParam(workers = 391, cluster="slurm", template="~/slurm_bc_template.tmpl",
        resources=list(walltime=60*60*24*5, ncpus=8)))
all <- rbindlist(fin_res)
colnames(all)[which(colnames(all) == "MLST.CC")] <- "predicted_MLST.CC"
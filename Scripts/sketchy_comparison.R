## compare with sketchy

FILLED_META_FILE <- "all_meta_filled.csv"
SKETCHY_PATH <- "~/.conda/envs/medaka/envs/sketchy/bin/sketchy"
SLURM_TEMPLATE <- "~/slurm_bc_template.tmpl"
SKETCHY_DEFAULT_REFERENCE_FILE <- "~/sketchy/saureus_default/saureus_s1000_k16_full.msh"
SKETCHY_DEFAULT_GENOME_FILE <- "~/sketchy/saureus_default/genotypes_full.tsv"
SKETCHY_LARGE_REFERENCE_FILE <- "~/sketchy/saureus_large/saureus_s10000_k16_full.msh"
SKETCHY_LARGE_GENOME_FILE <- "~/sketchy/saureus_large/genotypes_full.tsv"

library(data.table)
library(BiocParallel)

all_simulated <- list.files(pattern = ".*.fastq")

### Default reference file (1000 hashes)
default_resolution_result <- bplapply(all_simulated, function(file, SKETCHY_PATH, SKETCHY_REFERENCE_FILE, SKETCHY_GENOTYPE_FILE){
    library(data.table)
    result_for_file <- list()
    for (read_n in c(1,2,3,4,5,6,7,8)){
        n_line <- read_n*4*1000
        command <- paste0("sed -n 1,", n_line, "p ", file, " |")
        command <- paste0(command, SKETCHY_PATH, " predict -r ", SKETCHY_REFERENCE_FILE, " -g ", SKETCHY_GENOTYPE_FILE, " -t 5 -H")
        result_for_file[[read_n]] <- fread(text = system(command, intern = TRUE))
    }
    result_for_file <- rbindlist(result_for_file)
    result_for_file$pubmlst_id <- gsub("S", "", gsub("_0001.fastq", "", file))
    result_for_file$resolution <- "1k"
    return(result_for_file)
}, SKETCHY_PATH= SKETCHY_PATH, SKETCHY_REFERENCE_FILE = SKETCHY_DEFAULT_REFERENCE_FILE,
    SKETCHY_GENOTYPE_FILE = SKETCHY_DEFAULT_GENOME_FILE,
    BPPARAM = BatchtoolsParam(workers = 100, progress = TRUE, cluster="slurm", template=SLURM_TEMPLATE,
        resources=list(walltime=60*60*24*5, ncpus=1))
)
agg_default_result <- rbindlist(default_resolution_result)
fwrite(agg_default_result, "sketchy_agg_result_1k.csv", row.names = FALSE)


### High resolution reference file (10,000 hashes)
high_resolution_result <- bplapply(all_simulated, function(file, SKETCHY_PATH, SKETCHY_REFERENCE_FILE, SKETCHY_GENOTYPE_FILE){
    library(data.table)
    result_for_file <- list()
    for (read_n in c(1,2,3,4,5,6,7,8)){
        n_line <- read_n*4*1000
        command <- paste0("sed -n 1,", n_line, "p ", file, " |")
        command <- paste0(command, SKETCHY_PATH, " predict -r ", SKETCHY_REFERENCE_FILE, " -g ", SKETCHY_GENOTYPE_FILE, " -t 5 -H")
        result_for_file[[read_n]] <- fread(text = system(command, intern = TRUE))
    }
    result_for_file <- rbindlist(result_for_file)
    result_for_file$pubmlst_id <- gsub("S", "", gsub("_0001.fastq", "", file))
    result_for_file$resolution <- "10k"
    return(result_for_file)
}, SKETCHY_PATH= SKETCHY_PATH, SKETCHY_REFERENCE_FILE = SKETCHY_LARGE_REFERENCE_FILE,
    SKETCHY_GENOTYPE_FILE = SKETCHY_LARGE_GENOME_FILE,
    BPPARAM = BatchtoolsParam(workers = 10, progress = TRUE, cluster="slurm", template=SLURM_TEMPLATE,
        resources=list(walltime=60*60*24*5, ncpus=1))
)
agg_high_result <- rbindlist(high_resolution_result)
fwrite(agg_high_result, "sketchy_agg_result_10k.csv", row.names = FALSE)

### COMBINED ANALYSIS
combined_result <- rbindlist(list(agg_default_result, agg_high_result))
combined_result$mlst <- as.numeric(gsub("ST", "", combined_result$mlst))

mlst_profiles <- fread("mlst_profiles.txt")
mlst_profiles <- mlst_profiles[, c("ST", "clonal_complex")]

fin_result <- merge(combined_result, mlst_profiles, by.x = "mlst", by.y = "ST")
fin_result[clonal_complex == ""]$clonal_complex <- paste0("CC", fin_result[clonal_complex == ""]$mlst)


to_collapse <- fin_result[fin_result[, .I[shared_hashes == max(shared_hashes)], by = list(pubmlst_id, resolution, reads)]$V1]
collapsed_top_result <- to_collapse[,.(predicted_MLST.CC = paste0(unique(clonal_complex), collapse = ","), predicted_MLST.ST = paste0(unique(mlst), collapse = ","), predicted_mecA = paste0(unique(meca), collapse = ","), predicted_pvl = paste0(unique(pvl), collapse = ","), shared_hashes = paste0(unique(shared_hashes), collapse = ",")), by = list(pubmlst_id, resolution, reads)]

all_filled_meta <- fread(FILLED_META_FILE)[!is.na(pubmlst_id)][, c("pubmlst_id", "MLST.CC", "CC_filled")]

all_filled_meta$pubmlst_id <- as.character(all_filled_meta$pubmlst_id)
collapsed_top_result$pubmlst_id <- as.character(collapsed_top_result$pubmlst_id)

merged_result <- merge(collapsed_top_result, all_filled_meta, by = "pubmlst_id")
merged_result$predicted_mlst_is_single <- sapply(merged_result$predicted_MLST.CC, function(x) {!grepl(",", x)})
merged_result$predicted_cc_is_correct <- sapply(seq_len(nrow(merged_result)), function(i){
    if (merged_result$predicted_mlst_is_single[i] == TRUE){
        return(merged_result$predicted_MLST.CC[i] %in% strsplit(merged_result$CC_filled[i], split = "/")[[1]])
    } else{
        return(any(unlist(strsplit(merged_result$predicted_MLST.CC[i], ",")) %in% strsplit(merged_result$CC_filled[i], split = "/")[[1]] ))
    }
})
merged_result$predicted_st_is_correct <- sapply(seq_len(nrow(merged_result)), function(i){
    return(any(paste0("CC", unlist(strsplit(merged_result$predicted_MLST.ST[i], ","))) %in% strsplit(merged_result$CC_filled[i], split = "/")[[1]] ))
})
merged_result$predicted_mlst_is_correct <- merged_result$predicted_st_is_correct | merged_result$predicted_cc_is_correct

### MAJOR LINEAGE ASSIGNMENT COMPARISON
fwrite(merged_result, "sketchy_merged_result.csv", row.names = FALSE)
summary <- merged_result[,.(n_correct = length(which(predicted_mlst_is_correct)),
                n_single = length(which(predicted_mlst_is_single)),
                n_single_correct = length(which(predicted_mlst_is_single & predicted_mlst_is_correct))),
            by = list(resolution, reads)]
fwrite(summary, "sketchy_summary.csv", row.names = FALSE)

### GENE DETECTION COMPARISON
summarised_to_meta <- fread("gene_summarised_to_meta.csv")
summarised_to_meta$pubmlst_id <- sapply(strsplit(summarised_to_meta$seq_name, split = "_"), `[`, 1)

combined_result$pubmlst_id <- as.character(combined_result$pubmlst_id)
temp <- merge(combined_result, summarised_to_meta[,-c("seq_name")], by = "pubmlst_id")

temp <- temp[,
    list(mecA, lukF_PV, lukS_PV,
        mecA_pred = paste0(meca, collapse = ","),
        pvl_pred = paste0(pvl, collapse = ","),
        shared_hashes = paste0(shared_hashes, collapse = ",")),
    by = list(resolution, reads, pubmlst_id)
]

temp$multiple_meca <- sapply(temp$mecA_pred, function(x){
    return(length(unique(strsplit(x, split = ",")[[1]])) != 1)}
)
temp$multiple_pvl <- sapply(temp$pvl_pred, function(x){
    pvl_stats <- unique(strsplit(x, split = ",")[[1]])
    pvl_stats <- gsub("\\*", "-", pvl_stats)
    return(length(unique(pvl_stats)) != 1)}
)

temp$mecA_correct <- sapply(seq_len(nrow(temp)), function(i){
    MRSA_in_closest <- grepl("MRSA", temp$mecA_pred[i])
    return( temp$mecA[i] == MRSA_in_closest)
})

temp$pvl_correct <- sapply(seq_len(nrow(temp)), function(i){
    PVL_in_closest <- grepl("PVL+", temp$pvl_pred[i])
    return( temp$lukF_PV[i] == PVL_in_closest)
})

temp$mecA_sensitive <- sapply(seq_len(nrow(temp)), function(i){
    MRSA_in_closest <- grepl("MRSA", temp$mecA_pred[i])
    return( MRSA_in_closest )
})

temp$pvl_sensitive <- sapply(seq_len(nrow(temp)), function(i){
    PVL_in_closest <- grepl("PVL+", temp$pvl_pred[i])
    return( PVL_in_closest )
})

# SKIP
temp$mecA_correct_major <- unlist(bplapply(seq_len(nrow(temp)), function(i){
    mecA_count <- table(strsplit(temp$mecA_pred[i], split = ",")[[1]])
    predicted <- names(which(mecA_count == max(mecA_count)))
    if (predicted == "MRSA"){
        MRSA_in_closest <- TRUE
    } else{
        MRSA_in_closest <- FALSE
    }
    return( temp$mecA[i] == MRSA_in_closest)
}, BPPARAM=MulticoreParam(workers = 4, progress = TRUE)))

temp$pvl_correct_major <- unlist(bplapply(seq_len(nrow(temp)), function(i){
    pvl_stats <- strsplit(temp$pvl_pred, split = ",")[[1]]
    pvl_stats <- gsub("\\*", "-", pvl_stats)
    pvl_count <- table(pvl_stats)
    predicted <- names(which(pvl_count == max(pvl_count)))
    
    if(predicted == "PVL+"){
        PVL_in_closest <- TRUE
    } else{
        PVL_in_closest <- FALSE
    }
    return( temp$lukF_PV[i] == PVL_in_closest)
}, BPPARAM=MulticoreParam(workers = 16, progress = TRUE)))

# END SKIP

temp$mecA_pred_major <- unlist(bplapply(seq_len(nrow(temp)), function(i){
    mecA_count <- table(strsplit(temp$mecA_pred[i], split = ",")[[1]])
    predicted <- names(which(mecA_count == max(mecA_count)))
    if (predicted == "MRSA"){
        MRSA_in_closest <- TRUE
    } else{
        MRSA_in_closest <- FALSE
    }
    return( MRSA_in_closest )
}, BPPARAM=MulticoreParam(workers = 4, progress = TRUE)))

temp$pvl_pred_major <- unlist(bplapply(seq_len(nrow(temp)), function(i){
    pvl_stats <- strsplit(temp$pvl_pred, split = ",")[[1]]
    pvl_stats <- gsub("\\*", "-", pvl_stats)
    pvl_count <- table(pvl_stats)
    predicted <- names(which(pvl_count == max(pvl_count)))
    
    if(predicted == "PVL+"){
        PVL_in_closest <- TRUE
    } else{
        PVL_in_closest <- FALSE
    }
    return( PVL_in_closest )
}, BPPARAM=MulticoreParam(workers = 16, progress = TRUE)))

temp$mecA_pred_first <- sapply(seq_len(nrow(temp)), function(i){
    res <- strsplit(temp$mecA_pred[i], split = ",")[[1]][1]
    if (res == "MRSA"){
        res <- TRUE
    } else{
        res <- FALSE
    }
    return( res )
})

temp$pvl_pred_first <- sapply(seq_len(nrow(temp)), function(i){
    res <- strsplit(temp$pvl_pred[i], split = ",")[[1]][1]
    res <- gsub("\\*", "-", res)
    if (res == "PVL+"){
        res <- TRUE
    } else{
        res <- FALSE
    }
    return( res )
})

fwrite(temp, "sketchy_gene_comparison.csv", row.names = FALSE)

summary <- temp[,
    list(
        mecA_major_TP = length(which(mecA_pred_major & mecA)),
        mecA_major_TN = length(which(!mecA_pred_major & !mecA)),
        mecA_major_FP = length(which(mecA_pred_major & !mecA)),
        mecA_major_FN = length(which(!mecA_pred_major & mecA)),
        pvl_major_TP = length(which(pvl_pred_major & lukF_PV)),
        pvl_major_TN = length(which(!pvl_pred_major & !lukF_PV)),
        pvl_major_FP = length(which(pvl_pred_major & !lukF_PV)),
        pvl_major_FN = length(which(!pvl_pred_major & lukF_PV)),
        mecA_sensitive_TP = length(which(mecA_sensitive & mecA)),
        mecA_sensitive_TN = length(which(!mecA_sensitive & !mecA)),
        mecA_sensitive_FP = length(which(mecA_sensitive & !mecA)),
        mecA_sensitive_FN = length(which(!mecA_sensitive & mecA)),
        pvl_sensitive_TP = length(which(pvl_sensitive & lukF_PV)),
        pvl_sensitive_TN = length(which(!pvl_sensitive & !lukF_PV)),
        pvl_sensitive_FP = length(which(pvl_sensitive & !lukF_PV)),
        pvl_sensitive_FN = length(which(!pvl_sensitive & lukF_PV)),
        mecA_first_TP = length(which(mecA_pred_first & mecA)),
        mecA_first_TN = length(which(!mecA_pred_first & !mecA)),
        mecA_first_FP = length(which(mecA_pred_first & !mecA)),
        mecA_first_FN = length(which(!mecA_pred_first & mecA)),
        pvl_first_TP = length(which(pvl_pred_first & lukF_PV)),
        pvl_first_TN = length(which(!pvl_pred_first & !lukF_PV)),
        pvl_first_FP = length(which(pvl_pred_first & !lukF_PV)),
        pvl_first_FN = length(which(!pvl_pred_first & lukF_PV))
    ),
    by = list(resolution, reads)
][,
    list(
        resolution, reads,
        mecA_major_TP, mecA_major_TN, mecA_major_FP, mecA_major_FN,
        pvl_major_TP, pvl_major_TN, pvl_major_FP, pvl_major_FN,
        mecA_sensitive_TP, mecA_sensitive_TN, mecA_sensitive_FP, mecA_sensitive_FN,
        pvl_sensitive_TP, pvl_sensitive_TN, pvl_sensitive_FP, pvl_sensitive_FN,
        mecA_first_TP, mecA_first_TN, mecA_first_FP, mecA_first_FN,
        pvl_first_TP, pvl_first_TN, pvl_first_FP, pvl_first_FN,
        mecA_major_sensitivity = mecA_major_TP / (mecA_major_TP + mecA_major_FN),
        mecA_major_specificity = mecA_major_TN / (mecA_major_TN + mecA_major_FP),
        pvl_major_sensitivity = pvl_major_TP / (pvl_major_TP + pvl_major_FN),
        pvl_major_specificity = pvl_major_TN / (pvl_major_TN + pvl_major_FP),
        mecA_sensitive_sensitivity = mecA_sensitive_TP / (mecA_sensitive_TP + mecA_sensitive_FN),
        mecA_sensitive_specificity = mecA_sensitive_TN / (mecA_sensitive_TN + mecA_sensitive_FP),
        pvl_sensitive_sensitivity = pvl_sensitive_TP / (pvl_sensitive_TP + pvl_sensitive_FN),
        pvl_sensitive_specificity = pvl_sensitive_TN / (pvl_sensitive_TN + pvl_sensitive_FP),
        mecA_first_sensitivity = mecA_first_TP / (mecA_first_TP + mecA_first_FN),
        mecA_first_specificity = mecA_first_TN / (mecA_first_TN + mecA_first_FP),
        pvl_first_sensitivity = pvl_first_TP / (pvl_first_TP + pvl_first_FN),
        pvl_first_specificity = pvl_first_TN / (pvl_first_TN + pvl_first_FP)
    )
]

fwrite(summary, "sketchy_gene_summary.csv", row.names = FALSE)
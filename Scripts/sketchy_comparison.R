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

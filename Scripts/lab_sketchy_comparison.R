## compare with sketchy

FILLED_META_FILE <- "all_meta_filled.csv"
SKETCHY_PATH <- "~/.conda/envs/medaka/envs/sketchy/bin/sketchy"
SLURM_TEMPLATE <- "~/slurm_bc_template.tmpl"
SKETCHY_DEFAULT_REFERENCE_FILE <- "~/sketchy/saureus_default/saureus_s1000_k16_full.msh"
SKETCHY_DEFAULT_GENOME_FILE <- "~/sketchy/saureus_default/genotypes_full.tsv"
SKETCHY_LARGE_REFERENCE_FILE <- "~/sketchy/saureus_large/saureus_s10000_k16_full.msh"
SKETCHY_LARGE_GENOME_FILE <- "~/sketchy/saureus_large/genotypes_full.tsv"
NANOPORE_RUN_META <- "Nanopore_run_meta.csv"

library(data.table)
library(BiocParallel)

default_resolution_result <- list()
high_resolution_result <- list()
for (basecall_method in c("fast", "hac", "sup")){
    all_reads <- list.files(path = basecall_method, pattern = ".*.fastq.gz", full.names = TRUE)

    ### Default reference file (1000 hashes)
    temp <- bplapply(all_reads, function(file, SKETCHY_PATH, SKETCHY_REFERENCE_FILE, SKETCHY_GENOTYPE_FILE, basecall_method){
        library(data.table)
        result_for_file <- list()
        for (read_n in c(1,2,3,4,5,6,7,8)){
            n_line <- read_n*4*1000
            command <- paste0("zcat ", file, " | sed -n 1,", n_line, "p |")
            command <- paste0(command, SKETCHY_PATH, " predict -r ", SKETCHY_REFERENCE_FILE, " -g ", SKETCHY_GENOTYPE_FILE, " -t 5 --consensus -H")
            result_for_file[[read_n]] <- fread(text = system(command, intern = TRUE))
        }
        result_for_file <- rbindlist(result_for_file)
        result_for_file$barcode <- strsplit(gsub(".fastq.gz", "", file), split = "_")[[1]][2]
        result_for_file$resolution <- "1k"
        return(result_for_file)
    }, SKETCHY_PATH= SKETCHY_PATH, SKETCHY_REFERENCE_FILE = SKETCHY_DEFAULT_REFERENCE_FILE,
        SKETCHY_GENOTYPE_FILE = SKETCHY_DEFAULT_GENOME_FILE, basecall_method = basecall_method,
        BPPARAM = BatchtoolsParam(workers = 100, progress = TRUE, cluster="slurm", template=SLURM_TEMPLATE,
            resources=list(walltime=60*60*24*5, ncpus=1))
    )
    default_resolution_result[[basecall_method]] <- rbindlist(temp)
    default_resolution_result[[basecall_method]]$basecall_method <- basecall_method

    ### High resolution reference file (10,000 hashes)
    temp <- bplapply(all_reads, function(file, SKETCHY_PATH, SKETCHY_REFERENCE_FILE, SKETCHY_GENOTYPE_FILE, basecall_method){
        library(data.table)
        result_for_file <- list()
        for (read_n in c(1,2,3,4,5,6,7,8)){
            n_line <- read_n*4*1000
            command <- paste0("zcat ", file, " | sed -n 1,", n_line, "p |")
            command <- paste0(command, SKETCHY_PATH, " predict -r ", SKETCHY_REFERENCE_FILE, " -g ", SKETCHY_GENOTYPE_FILE, " -t 5 --consensus -H")
            result_for_file[[read_n]] <- fread(text = system(command, intern = TRUE))
        }
        result_for_file <- rbindlist(result_for_file)
        result_for_file$barcode <- strsplit(gsub(".fastq.gz", "", file), split = "_")[[1]][2]
        result_for_file$resolution <- "10k"
        return(result_for_file)
    }, SKETCHY_PATH= SKETCHY_PATH, SKETCHY_REFERENCE_FILE = SKETCHY_LARGE_REFERENCE_FILE,
        SKETCHY_GENOTYPE_FILE = SKETCHY_LARGE_GENOME_FILE, basecall_method = basecall_method,
        BPPARAM = BatchtoolsParam(workers = 10, progress = TRUE, cluster="slurm", template=SLURM_TEMPLATE,
            resources=list(walltime=60*60*24*5, ncpus=1))
    )
    high_resolution_result[[basecall_method]] <- rbindlist(temp)
    high_resolution_result[[basecall_method]]$basecall_method <- basecall_method 
}

agg_default_result <- rbindlist(default_resolution_result)
fwrite(agg_default_result, "lab_sketchy_agg_result_1k.csv", row.names = FALSE)

agg_high_result <- rbindlist(high_resolution_result)
fwrite(agg_high_result, "lab_sketchy_agg_result_10k.csv", row.names = FALSE)

### COMBINED ANALYSIS
combined_result <- rbindlist(list(agg_default_result, agg_high_result))
combined_result$mlst <- as.numeric(gsub("ST", "", combined_result$mlst))

mlst_profiles <- fread("mlst_profiles.txt")
mlst_profiles <- mlst_profiles[, c("ST", "clonal_complex")]

fin_result <- merge(combined_result, mlst_profiles, by.x = "mlst", by.y = "ST")
fin_result[clonal_complex == ""]$clonal_complex <- paste0("CC", fin_result[clonal_complex == ""]$mlst)

fin_result <- fin_result[, .(predicted_MLST.CC = paste0(unique(clonal_complex), collapse = ","), predicted_MLST.ST = paste0(unique(mlst), collapse = ","), predicted_mecA = paste0(unique(meca), collapse = ","), predicted_pvl = paste0(unique(pvl), collapse = ","), shared_hashes = paste0(unique(shared_hashes), collapse = ",")), by = list(barcode, resolution, basecall_method, reads)]

#to_collapse <- fin_result[fin_result[, .I[shared_hashes == max(shared_hashes)], by = list(barcode, resolution, basecall_method, reads)]$V1]
#collapsed_top_result <- to_collapse[,.(predicted_MLST.CC = paste0(unique(clonal_complex), collapse = ","), predicted_MLST.ST = paste0(unique(mlst), collapse = ","), predicted_mecA = paste0(unique(meca), collapse = ","), predicted_pvl = paste0(unique(pvl), collapse = ","), shared_hashes = paste0(unique(shared_hashes), collapse = ",")), by = list(barcode, resolution, basecall_method, reads)]

nanopore_run_meta <- fread(NANOPORE_RUN_META)
all_filled_meta <- fread(FILLED_META_FILE)
merged_info <- merge(nanopore_run_meta, all_filled_meta)
merged_info$barcode <- merged_info$Barcode.ID
all_filled_meta <- merged_info[, c("barcode", "MLST.CC", "CC_filled")]

all_filled_meta$barcode <- as.character(all_filled_meta$barcode)
fin_result$barcode <- as.character(fin_result$barcode)

merged_result <- merge(fin_result, all_filled_meta, by = "barcode")
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

for (bc in c(33, 73)){
    for (res in c("1k", "10k")){
        for (ctype in c("fast", "hac", "sup")){
            max_read <- max(merged_result[
                barcode == bc &
                resolution == res &
                basecall_method == ctype]$reads)
            last_result <- merged_result[
                barcode == bc &
                resolution == res &
                basecall_method == ctype &
                reads == max_read]
            rounded <- floor(max_read / 1000)*1000
            while(rounded != 8000){
                dup_last_result <- last_result
                dup_last_result$reads <- rounded + 1000
                merged_result <- rbind(merged_result, dup_last_result)
                rounded <- rounded + 1000
            }
        }
    }
}

merged_result <- merged_result[
                reads %in% seq(1000,8000, 1000)]
merged_result <- fread("lab_sketchy_merged_result.csv")
for (bc in unique(merged_result$barcode)){
    print(nrow(merged_result[barcode == bc]))
}
### MAJOR LINEAGE ASSIGNMENT COMPARISON
fwrite(merged_result, "lab_sketchy_merged_result.csv", row.names = FALSE)
summary <- merged_result[,.(n_correct = length(which(predicted_mlst_is_correct)),
                n_single = length(which(predicted_mlst_is_single)),
                n_single_correct = length(which(predicted_mlst_is_single & predicted_mlst_is_correct))),
            by = list(resolution, basecall_method, reads)][order(resolution, basecall_method, reads)]
fwrite(summary, "lab_sketchy_summary.csv", row.names = FALSE)

### GENE DETECTION COMPARISON
summarised_to_meta <- merged_info[,list(barcode, MRSA, PVL)]
summarised_to_meta$mecA <- ifelse(summarised_to_meta$MRSA == "+", TRUE, FALSE)
summarised_to_meta$lukF_PV <- ifelse(summarised_to_meta$PVL == "+", TRUE, FALSE)
summarised_to_meta$lukS_PV <- ifelse(summarised_to_meta$PVL == "+", TRUE, FALSE)
summarised_to_meta$barcode <- sprintf("%02d", as.numeric(summarised_to_meta$barcode))
merged_result$barcode <- sprintf("%02d", as.numeric(merged_result$barcode))
temp <- merge(merged_result, summarised_to_meta, by = "barcode")

temp$mecA_pred <- FALSE
temp[predicted_mecA == "MRSA"]$mecA_pred <- TRUE
temp$pvl_pred <- TRUE
temp[predicted_pvl == "PVL-"]$pvl_pred <- FALSE
fwrite(temp, "lab_sketchy_gene_comparison.csv", row.names = FALSE)

summary <- temp[,
    .(
        mecA_TP = length(which(mecA_pred & mecA)),
        mecA_TN = length(which(!mecA_pred & !mecA)),
        mecA_FP = length(which(mecA_pred & !mecA)),
        mecA_FN = length(which(!mecA_pred & mecA)),
        pvl_TP = length(which(pvl_pred & lukF_PV)),
        pvl_TN = length(which(!pvl_pred & !lukF_PV)),
        pvl_FP = length(which(pvl_pred & !lukF_PV)),
        pvl_FN = length(which(!pvl_pred & lukF_PV))
    ),
    by = list(resolution, basecall_method, reads)
][, .(resolution, basecall_method, reads, mecA_TP, mecA_TN, mecA_FP, mecA_FN, pvl_TP, pvl_TN, pvl_FP, pvl_FN,
    mecA_sensitivity = mecA_TP / (mecA_TP + mecA_FN),
    mecA_specificity = mecA_TN / (mecA_TN + mecA_FP),
    pvl_sensitivity = pvl_TP / (pvl_TP + pvl_FN),
    pvl_specificity = pvl_TN / (pvl_TN + pvl_FP)
)][order(resolution, basecall_method, reads)]

fwrite(summary, "lab_sketchy_gene_summary.csv", row.names = FALSE)
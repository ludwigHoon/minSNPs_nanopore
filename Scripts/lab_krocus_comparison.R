## Krocus comparison for lab generated data
#
# MLST file required are downloaded with Krocus before the start of the following script with:
# krocus_database_downloader  --species "Staphylococcus aureus"
#

FILLED_META_FILE <- "all_meta_filled.csv"
SLURM_TEMPLATE <- "slurm_bc_template.tmpl"
MLST_DATABASE <- "~/mlst_files"
KROCUS_PATH <- "/home/lhoon/.conda/envs/medaka/envs/sketchy/bin/krocus"
ST_PROFILE <- "/home/lhoon/mlst_files/profile.txt"
NANOPORE_RUN_META <- "Nanopore_run_meta.csv"

library(BiocParallel)
library(data.table)

bp <- BatchtoolsParam(workers = 15, cluster="slurm", template=SLURM_TEMPLATE,
    resources=list(walltime=60*60*24*5, ncpus=1), progress = TRUE)
krocus_result <- list()
for(basecall_method in c("fast", "hac", "sup")){
    all_reads <- list.files(path = basecall_method, pattern = ".*.fastq.gz", full.names = TRUE)
    temp <- bplapply(all_reads, function(file, KROCUS_PATH, MLST_DATABASE, basecall_method){
        library(data.table)
        result_for_file <- list()
        for (n_reads in c(1,2,3,4,5,6,7,8,9,10)){
            n_line <- n_reads*4*1000
            temp_file_name <- file.path(basecall_method, paste0("temp_", basename(file)))
            command <- paste0("sed -n 1,", n_line, "p ", file, " > ", temp_file_name)
            system(command)
            command <- paste0(KROCUS_PATH, " --print_interval 80000 ", MLST_DATABASE, " ", temp_file_name)
            result <- system(command, intern=TRUE)
            unlink(temp_file_name)
            processed_result <- data.frame(barcode =
                strsplit(gsub(".fastq.gz", "", file), split = "_")[[1]][2],
                n_reads = (n_reads * 1000))
            t_result <- unlist(strsplit(result, split = "\t"))
            processed_result$st <- t_result[1]
            processed_result$coverage <- t_result[2]
            for (allele in c("arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL")){
                processed_result[[allele]] <- paste0("",gsub(paste0(allele, "\\("), "", gsub(")", "", t_result[grep(allele, t_result)])))
            }
            result_for_file[[n_reads]] <- processed_result
        }
        result <- rbindlist(result_for_file)
        return(result)
    }, KROCUS_PATH = KROCUS_PATH, MLST_DATABASE = MLST_DATABASE, basecall_method = basecall_method, BPPARAM = bp)
    krocus_result[[basecall_method]] <- rbindlist(temp)
    krocus_result[[basecall_method]]$basecall_method <- basecall_method
}

krocus_result <- rbindlist(krocus_result)
fwrite(krocus_result, "lab_krocus_full_result.csv", row.names = FALSE)

### Checking if Krocus assign sample to same major lineage
mlst_profiles <- fread("mlst_profiles.txt")
mlst_profiles <- mlst_profiles[, c("ST", "clonal_complex")]

krocus_result$st <- as.character(krocus_result$st)
mlst_profiles$ST <- as.character(mlst_profiles$ST)
fin_result <- merge(krocus_result, mlst_profiles, by.x = "st", by.y = "ST", all.x = TRUE)

nanopore_run_meta <- fread(NANOPORE_RUN_META)
all_filled_meta <- fread(FILLED_META_FILE)
merged_info <- merge(nanopore_run_meta, all_filled_meta)
merged_info$barcode <- merged_info$Barcode.ID
all_filled_meta <- merged_info[, c("barcode", "MLST.CC", "CC_filled")]

all_filled_meta$barcode <- as.character(all_filled_meta$barcode)
fin_result$barcode <- as.character(fin_result$barcode)

merged_result <- merge(fin_result, all_filled_meta, by = "barcode")
merged_result$note <- ""
### Handling special case where ST is not defined by identifying the closest ST/CC
special_handling <- merged_result[st=="ND"]
arcC <- gsub("\\*", "", special_handling$arcC)
aroE <- gsub("\\*", "", special_handling$aroE)
glpF <- gsub("\\*", "", special_handling$glpF)
gmk <- gsub("\\*", "", special_handling$gmk)
pta <- gsub("\\*", "", special_handling$pta)
tpi <- gsub("\\*", "", special_handling$tpi)
yqiL <- gsub("\\*", "", special_handling$yqiL)

ST_profile <- fread(ST_PROFILE)
ST_profile$q_string <- paste(ST_profile$arcC, ST_profile$aroE, ST_profile$glpF, ST_profile$gmk, ST_profile$pta, ST_profile$tpi, ST_profile$yqiL, sep = "_")

pbar <- txtProgressBar(max = nrow(special_handling), style = 3)
for (i in seq_len(nrow(special_handling))){
    setTxtProgressBar(pbar, i)

    q_string <- paste(arcC[i], aroE[i], glpF[i], gmk[i], pta[i], tpi[i], yqiL[i], sep = "_")
    
    slv_q_strings <- sapply(c(1:7), function(x){
        alleles <- unlist(strsplit(q_string, split = "_"))
        alleles[x] <- ".*"
        return(paste0(alleles, collapse = "_"))
    })
    SLV_result <- ST_profile[grepl(paste(slv_q_strings, collapse="|"), ST_profile$q_string)][clonal_complex != ""]
    T_CC <- table(SLV_result$clonal_complex)
    potential_SLV_nSTs <- paste0(T_CC, collapse = ",")
    potential_SLV_CCs <- paste0(names(T_CC), collapse = ",")

    if (potential_SLV_CCs != "") {
        merged_result[st=="ND"][i]$clonal_complex <- potential_SLV_CCs
        merged_result[st=="ND"][i]$note <- paste0("taken from SLV with ", potential_SLV_nSTs, " STs")
    }
    
    dlv_q_strings <- unlist(sapply(c(1:6), function(x){
        results <- rep("", 7-x)
        for (y in seq_len(7-x)){
            alleles <- unlist(strsplit(q_string, split = "_"))
            alleles[x] <- ".*"
            alleles[x+y] <- ".*"
            results[y] <- paste0(alleles, collapse = "_")
        }
        return(results)
    }))
    DLV_result <- ST_profile[grepl(paste(dlv_q_strings, collapse="|"), ST_profile$q_string)][clonal_complex != ""]
    T_CC <- table(DLV_result$clonal_complex)
    potential_DLV_nSTs <- paste0(T_CC, collapse = ",")
    potential_DLV_CCs <- paste0(names(T_CC), collapse = ",")

    if (potential_DLV_CCs != "") {
        merged_result[st=="ND"][i]$clonal_complex <- potential_DLV_CCs
        merged_result[st=="ND"][i]$note <- paste0("taken from DLV with ", potential_DLV_nSTs, " STs")
    } else {
        merged_result[st=="ND"][i]$clonal_complex <- NA
        merged_result[st=="ND"][i]$note <- "Tested up to DLV, none found"
    }
}
close(pbar)

merged_result$predicted_cc_is_correct <- sapply(seq_len(nrow(merged_result)), function(i){
    return(
        any(
            strsplit(merged_result$clonal_complex[i], split = ",")[[1]] %in%
            strsplit(merged_result$CC_filled[i], split = "/")[[1]]
        ) || 
        (
            paste0("CC", merged_result$st[i]) %in%
            strsplit(merged_result$CC_filled[i], split = "/")[[1]]
        ) ||
        merged_result$MLST.CC[i] == merged_result$clonal_complex[i]
    )
})
merged_result$predicted_mlst_is_single <- sapply(merged_result$clonal_complex, function(x) {!grepl(",", x)})
merged_result[is.na(clonal_complex)]$predicted_mlst_is_single <- NA

fwrite(merged_result, "lab_krocus_merged_result.csv", row.names = FALSE)

summary <- merged_result[,.(n_correct = length(which(predicted_cc_is_correct)),
              n_undefined = length(which(is.na(clonal_complex))),
              n_single = length(which(predicted_mlst_is_single)),
              n_single_correct = length(which(predicted_mlst_is_single & predicted_cc_is_correct))
              ),
            by = list(basecall_method, n_reads)]
fwrite(summary, "lab_krocus_summary.csv", row.names = FALSE)
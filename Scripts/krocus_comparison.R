## Krocus comparison
#
# MLST file required are downloaded with Krocus before the start of the following script with:
# krocus_database_downloader  --species "Staphylococcus aureus"
#

FILLED_META_FILE <- "all_meta_filled.csv"
SLURM_TEMPLATE <- "slurm_bc_template.tmpl"
MLST_DATABASE <- "~/mlst_files"
KROCUS_PATH <- "/home/lhoon/.conda/envs/medaka/envs/sketchy/bin/krocus"
ST_PROFILE <- "/home/lhoon/mlst_files/profile.txt"

library(BiocParallel)
library(data.table)


bp <- BatchtoolsParam(workers = 15, cluster="slurm", template=SLURM_TEMPLATE,
    resources=list(walltime=60*60*24*5, ncpus=1), progress = TRUE)


all_simulated <- list.files(pattern = ".*.fastq")
krocus_result <- bplapply(all_simulated, function(file, KROCUS_PATH, MLST_DATABASE){
    library(data.table)
    result_for_file <- list()
    for (n_reads in c(1,2,3,4,5,6,7,8)){
        n_line <- n_reads*4*1000
        temp_file_name <- paste0("temp_", file)
        command <- paste0("sed -n 1,", n_line, "p ", file, " > ", temp_file_name)
        system(command)
        command <- paste0(KROCUS_PATH, " --print_interval 80000 ", MLST_DATABASE, " ", temp_file_name)
        result <- system(command, intern=TRUE)
        unlink(temp_file_name)
        processed_a.frame(pubmlst_id = gsub("S", "", gsub("_0001.fastq", "", file)), n_reads = (n_reads * 1000))
        t_result <- unlist(strresult <- datsplit(result, split = "\t"))
        processed_result$st <- t_result[1]
        processed_result$coverage <- t_result[2]
        for (allele in c("arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL")){
            processed_result[[allele]] <- paste0("",gsub(paste0(allele, "\\("), "", gsub(")", "", t_result[grep(allele, t_result)])))
        }
        result_for_file[[n_reads]] <- processed_result
    }
    result <- rbindlist(result_for_file)
    return(result)
}, KROCUS_PATH = KROCUS_PATH, MLST_DATABASE = MLST_DATABASE, BPPARAM = bp)

krocus_result <- rbindlist(krocus_result)
fwrite(krocus_result, "krocus_full_result.csv", row.names = FALSE)

### Checking if Krocus assign sample to same major lineage
mlst_profiles <- fread("mlst_profiles.txt")
mlst_profiles <- mlst_profiles[, c("ST", "clonal_complex")]

krocus_result$st <- as.character(krocus_result$st)
mlst_profiles$ST <- as.character(mlst_profiles$ST)
fin_result <- merge(krocus_result, mlst_profiles, by.x = "st", by.y = "ST", all.x = TRUE)

all_filled_meta <- fread(FILLED_META_FILE)[!is.na(pubmlst_id)][, c("pubmlst_id", "MLST.CC", "CC_filled")]

all_filled_meta$pubmlst_id <- as.character(all_filled_meta$pubmlst_id)
fin_result$pubmlst_id <- as.character(fin_result$pubmlst_id)

merged_result <- merge(fin_result, all_filled_meta, by = "pubmlst_id")
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
    
    dlv_q_strings <- sapply(c(1:6), function(x){
        alleles <- unlist(strsplit(q_string, split = "_"))
        alleles[x] <- ".*"
        alleles[x+1] <- ".*"
        return(paste0(alleles, collapse = "_"))
    })
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

summary <- merged_result[,.(n_correct = length(which(predicted_cc_is_correct)),
              n_undefined = length(which(is.na(clonal_complex))),
              n_single = length(which(predicted_mlst_is_single)),
              n_single_correct = length(which(predicted_mlst_is_single & predicted_cc_is_correct))
              ),
            by = list(n_reads)]
fwrite(summary, "krocus_summary.csv", row.names = FALSE)
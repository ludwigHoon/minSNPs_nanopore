## Assign CC metadata to new samples based on the most similar (judged by SNP differences) sample in megaalignment 
## Metadata.csv is based on: https://microreact.org/project/uzDospiZeBSLfy4A8d1Uts-mega-alignment-phylogeny-and-minsnps-profiles
## but has CC71 changed to CC97 (following MLST convention)
library(BiocParallel)
library(minSNPs)

all_seq_matrix <- read_fasta("mega_with_newdata.fasta")
megaalignment_meta <- read.csv("Metadata.csv")
megaalignment_meta <- megaalignment_meta[megaalignment_meta$Isolate.name %in% names(all_seq_matrix),]
megaalignment_meta$MLST.CC <- gsub("(Unknown \\(|\\))", "", megaalignment_meta$MLST.CC)
newdata_meta <- read.csv("BIGSdb_3343897_1178826571_39602.csv")
newdata_seq_name <- names(all_seq_matrix)[which(!names(all_seq_matrix) %in% megaalignment_meta$Isolate.name)]
newdata_id_name <- data.frame(pubmlst_id = as.numeric(sapply(strsplit(newdata_seq_name, split = "_"), `[`, 1)),
    Isolate.name = newdata_seq_name)
newdata_meta <- merge(newdata_meta, newdata_id_name, by.x = "id", by.y = "pubmlst_id")

res <- bplapply(newdata_seq_name, function(x){
    ori_seq <- megaalignment_meta$Isolate.name

    genetic_distance <- sapply(all_seq_matrix[ori_seq], function(y){
        length(which(y != all_seq_matrix[[ x ]]))
    })
    most_similar <- names(genetic_distance)[which(genetic_distance == min(genetic_distance))]
    distance <- min(genetic_distance)
    return(data.frame(new_sample = x,
            most_similar = paste(most_similar, collapse = ","),
            distance = distance,
            cc = paste(unique(megaalignment_meta[megaalignment_meta$Isolate.name %in% most_similar,]$MLST.CC), collapse = ",")
        )
    )
}, BPPARAM = MulticoreParam(workers = 16, progress = TRUE))
library(data.table)
predicted_meta <- rbindlist(res)
write.csv(predicted_meta, "new_sample_most_similar_mega.csv", row.names = FALSE)
meta_new <- merge(predicted_meta, newdata_meta, by.x = "new_sample", by.y = "Isolate.name")
meta_new <- meta_new[,list(Isolate.name = new_sample, pubmlst_id = id, MLST.CC= clonal_complex..MLST., CC_filled = cc)]
setDT(megaalignment_meta)
megaalignment_meta$pubmlst_id <- NA
megaalignment_meta$CC_filled <- megaalignment_meta$MLST.CC
DISAGREEMENT <- merge(meta_new[meta_new$CC_filled != meta_new$MLST.CC & meta_new$MLST.CC != ""], predicted_meta, by.x = "Isolate.name", by.y = "new_sample")
write.csv(DISAGREEMENT, "disagreement.csv", row.names = FALSE)

res2 <- bplapply(newdata_seq_name, function(x){

    invalid_missing_SNPs <- length(
        which(
            !toupper(all_seq_matrix[[ x ]]) %in% c("A", "T", "C", "G")
            )
        )
    total <- length(all_seq_matrix[[ x ]]) 
    return(data.frame(Isolate.name = x,
            n_unknown = invalid_missing_SNPs,
            prop_unknown = invalid_missing_SNPs/total
        )
    )
}, BPPARAM = MulticoreParam(workers = 16, progress = TRUE))
missing_unknown_report <- rbindlist(res2)
fwrite(missing_unknown_report, "missing_unknown_report.csv", row.names = FALSE)

megaalignment_meta$n_unknown <- 0
megaalignment_meta$prop_unknown <- 0
meta_new <- merge(meta_new, missing_unknown_report, by = "Isolate.name")

FINAL_META <- rbind(megaalignment_meta, meta_new)
FINAL_META$has_disagreement <- FALSE
FINAL_META[Isolate.name %in% DISAGREEMENT$Isolate.name,]$has_disagreement <- TRUE
write.csv(FINAL_META, "all_meta_filled.csv", row.names = FALSE)
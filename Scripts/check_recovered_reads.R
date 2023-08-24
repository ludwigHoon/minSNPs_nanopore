library(BiocParallel)
library(data.table)
call_type <- c("fast", "hac", "sup")

results <- list()
for (type in call_type){
    files <- list.files(path = type, pattern = ".*.fastq.gz", full.names = TRUE)
    print(files)
    result <- bplapply(files, function(file){
        cmd <- paste("zcat", file, "| wc -l")
        output <- as.numeric(system(cmd, intern = TRUE))
        return(data.frame(call_type = type, barcode_id = gsub("combined_", "", gsub(".fastq.gz", "", basename(file))), output = output / 4))
    }, BPPARAM = MulticoreParam(workers = 8, progress = TRUE))
    results[[type]] <- rbindlist(result)
}
comb_res <- rbindlist(results)
fwrite(comb_res, "recovered_reads.csv", row.names = FALSE)
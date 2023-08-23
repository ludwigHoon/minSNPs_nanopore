## This script extracts the 164,335 SNPs defined in megaalignment from the new samples
## OUTF1.fasta can be obtained from https://figshare.com/articles/online_resource/Megaalignment_Supp_5_/19582885/1?file=36721473
## ref_OUTF1.csv can be obtained from https://figshare.com/articles/online_resource/Megaalignment_Supp_5_/19582885/1?file=36721476
## This makes use of SNIPPY's .aligned.fa output, which is a version of the reference but with - at position with depth=0 and N for 0 < depth < --mincov (does not have variants)
##
library(BiocParallel)
library(minSNPs)
interested_dirs <- gsub(".fas", "", list.files(pattern = "^.*.fas$"))
interested_files <- paste0(interested_dirs, "/snps.aligned.fa")
ref_meta <- read.csv("ref_OUTF1.csv")
all_seqs <- bplapply(interested_files, function(file){
    ind_seq <- read_fasta(file)[[1]]
    subset <- ind_seq[ref_meta$genome_position]
    return(subset)
}, BPPARAM = MulticoreParam(workers = 8, progress = TRUE))
names(all_seqs) <- interested_dirs
write_fasta(all_seqs, "mega_with_newdata.fasta")
system("cat OUTF1.fasta >> mega_with_newdata.fasta")
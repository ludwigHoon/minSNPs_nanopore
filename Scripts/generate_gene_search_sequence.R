## Generating search sequences for gene
# The sequence and the gene annotations can be found in: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000144955.2/
# Gene information are extracted from refseq prior and put in a csv table

library(data.table)
library(minSNPs)
library(BiocParallel)

geneInfo <- fread("gene_info.csv")
sample <- read_fasta("JKD6159.fasta")[[1]]

extract_relevant_sequences_from_genome <- function(start, end, genome, rev, name){
    sequence <- paste(genome[start:end], collapse = "")
    if (rev){
        sequence <- reverse_complement(sequence)
    }

    return(data.frame(name=name, sequence=sequence))
}

search_table <- list()
for (ri in seq_len(nrow(geneInfo))){
    start <- geneInfo[ri]$start
    end <- geneInfo[ri]$end
    rev <- geneInfo[ri]$rev
    gene <- geneInfo[ri]$gene
    seq <- extract_relevant_sequences_from_genome(start, end, sample, rev, gene)$sequence
    search_table[[ri]] <- generate_kmer_search_string(gene_seq = seq, k = 15, id_prefix = gene, bp = MulticoreParam(workers = 2))
}

final_result <- rbindlist(search_table)
repeated_sequence <- final_result[which(duplicated(final_result$sequence))]$sequence

final_result <- final_result[!sequence %in% repeated_sequence]
fwrite(final_result, "gene_search_sequence.csv", row.names = FALSE)
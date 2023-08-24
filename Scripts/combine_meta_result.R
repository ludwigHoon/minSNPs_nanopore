## Combine meta-data for microreact
library(data.table)
meta <- fread("all_meta_filled.csv")
ST_info <- fread("BIGSdb_3343897_1178826571_39602.csv")[,c(1:9)]
sim_result <- fread("simulated_385_8000.csv")

meta <- merge(meta, ST_info, by.x = "pubmlst_id", by.y = "id", all.x = TRUE)
meta <- merge(meta, sim_result, by = "Isolate.name", all.x = TRUE)
fwrite(meta, "all_meta_with_sim_result.csv", row.names = FALSE)
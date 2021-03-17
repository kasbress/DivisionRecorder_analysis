library(metacell)
library(dplyr)
library(tidyverse)

## initiate the metacell database
scdb_init("./metacell_db/", force_reinit=T) 


## Import metacell object and mat object
mc <- scdb_mc("MSC_mc_big2_col")
clean_mat <- scdb_mat("MSC_clean")

## Get Hashtag information and MC identity per cellcode
meta_anot <- as.data.frame(mc@mc)
hash_ids <- clean_mat@cell_metadata["hash.ID"]
hash_ids$hash.ID <- as.character(hash_ids$hash.ID)
test <- merge(meta_anot, hash_ids, all = T, by = "row.names")
test <- test[!test$Row.names %in% as.character(clean_mat@ignore_cells),]

## count cells, normalize and get ratios
cell_counts <- as.data.frame.matrix(table(test$`mc@mc`, test$hash.ID))
cell_counts_norm <- as.data.frame(apply(X = cell_counts, MARGIN = 2, function(x) (x/sum(x))*1000))
cell_counts_norm$lab_ratio <- rowSums(cell_counts_norm[,c(2,4,6)])/rowSums(cell_counts_norm[,c(1,3,5)])
lab_ratio <- cell_counts_norm$lab_ratio
names(lab_ratio) <- c(1:22)


## Save cellcounts and labeling ratio's for downstream analysis
saveRDS(lab_ratio, "./labeling_ratios2.rds")
saveRDS(cell_counts, "./cellcounts2.rds")


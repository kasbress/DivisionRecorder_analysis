library(metacell)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)

setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_memory/")

## import data
scdb_init("./metacell_db/", force_reinit=T) 
mc <- scdb_mc("DivRecMEM_MC")
lfp <- log2(mc@mc_fp)


### Effector and multipotency genes
eff <-c("Tbx21", "Prdm1", "Id2", "Gzmb", "Klrg1", "Cx3cr1", "Prf1", "Lgals1","Ifng")
mult <- c("Eomes", "Il7r","Sell", "Tcf7", "Bcl2", "Ccr7","Myb", "Bach2")
genes <- c(eff, mult)

## Focus on genes and remove terminal MCs
to.plot <- as.matrix(lfp[genes,-c(17,19,20,21)])

## make the heatmap
pdf("./metacell_figs/phenotype_clustering/heatmap_eff_vs_multWOterm.pdf", height = 5, width = 10)
pheatmap(mat = to.plot,  scale = "row", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "correlation", breaks = seq(-1.5,1.5,by=0.03),
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), cutree_rows = 2, cutree_cols = 2)
dev.off()


## Focus on genes and remove terminal MCs
to.plot <- as.matrix(lfp[genes,])

## make the heatmap
pdf("./metacell_figs/phenotype_clustering/heatmap_eff_vs_multWterm.pdf", height = 5, width = 10)
pheatmap(mat = to.plot,  scale = "row", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", breaks = seq(-1.5,1.5,by=0.03),
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), cutree_rows = 2, cutree_cols = 2)
dev.off()


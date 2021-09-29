library(tidyverse)
library(metacell)
library(pheatmap)

# Set working directory
setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_memory/")

# Initiate metacell database
scdb_init("./metacell_db/", force_reinit=T) 

# import data
mc <- scdb_mc("DivRecMEM_MC_fin")
lfp <- log2(mc@mc_fp)

# For multipotency related genes
select.genes <- c("Myb", "Tcf7","Bcl11b", "Lef1","Ets1", "Runx1", "Klf3", "Id3", "Zeb1", "Eomes",'Bcl6', "Foxo1")
pdf(file = "./metacell_figs/Heatmap_multipo.pdf", width = 3, height = length(select.genes)/3, title = "multipo" )
as.data.frame(lfp) %>% 
  dplyr::select(one_of(c('2','14','11','18','8','6'))) %>% 
  rownames_to_column("genes") %>%
  dplyr::filter(genes %in% select.genes) %>% 
  column_to_rownames("genes") %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",breaks = seq(-1.5,1.5,by=0.03), clustering_distance_rows = "euclidean", 
           border_color = NA, cutree_cols = 2, show_rownames = T, cluster_cols = F)
dev.off()


# For survival related genes
select.genes <- c('Gimap4', 'Gimap5', 'Gimap5', 'Gimap7', 'Gimap6', 
                  'Bcl2', 'Bcl11a', 'Traf5',  "Mcl1", "Il7r", "Dedd", 
                  'Birc7','Xiap','Birc6','Birc2','Birc3',"Espl1", "Traf2")
pdf(file = "./metacell_figs/Heatmap_multipo.pdf", width = 3, height = length(select.genes)/3, title = "multipo" )
as.data.frame(lfp) %>% 
  dplyr::select(one_of(c('2','14','11','18','8','6'))) %>% 
  rownames_to_column("genes") %>%
  dplyr::filter(genes %in% select.genes) %>% 
  column_to_rownames("genes") %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",breaks = seq(-1.5,1.5,by=0.03), clustering_distance_rows = "euclidean", 
           border_color = NA, cutree_cols = 2, show_rownames = T, cluster_cols = F)
dev.off()


# For cytotoxicity related genes
select.genes <- c('Gzma', 'Gzmb','Gzmm', 'Lamp1', 'Lamp2', 
                  'Raet1e', 'Prf1', 'Ctsc', 'Ctsh', 'Ifng',
                  'Ctsd',"Nkg7", "Ctsw", 'Ctsc',"Ctsl","Slamf7", "Cst7")
pdf(file = "./metacell_figs/Heatmap_multipo.pdf", width = 3, height = length(select.genes)/3, title = "multipo" )
as.data.frame(lfp) %>% 
  dplyr::select(one_of(c('2','14','11','18','8','6'))) %>% 
  rownames_to_column("genes") %>%
  dplyr::filter(genes %in% select.genes) %>% 
  column_to_rownames("genes") %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",breaks = seq(-1.5,1.5,by=0.03), clustering_distance_rows = "euclidean", 
           border_color = NA, cutree_cols = 2, show_rownames = T, cluster_cols = F)
dev.off()


# For inhibitory function related genes
select.genes <- c("Lag3", "Pdcd1", 'Havcr2', 'Cd244', 'Ctla4', 'Btla', 'Cd160', 'Tox')
pdf(file = "./metacell_figs/Heatmap_multipo.pdf", width = 3, height = length(select.genes)/3, title = "multipo" )
as.data.frame(lfp) %>% 
  dplyr::select(one_of(c('2','14','11','18','8','6'))) %>% 
  rownames_to_column("genes") %>%
  dplyr::filter(genes %in% select.genes) %>% 
  column_to_rownames("genes") %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",breaks = seq(-1.5,1.5,by=0.03), clustering_distance_rows = "euclidean", 
           border_color = NA, cutree_cols = 2, show_rownames = T, cluster_cols = F)
dev.off()

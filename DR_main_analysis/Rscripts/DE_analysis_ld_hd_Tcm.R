library(metacell)
library(Seurat)
library(limma)
library(tidyverse)
library(ggrepel)
library(gghighlight)
library(pheatmap)
library(fgsea)
library(msigdbr)
library(ggpubr)

# Set working directory
setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_memory/")

# Initiate metacell database
scdb_init("./metacell_db/", force_reinit=T) 

# Import metacell object
mc <- scdb_mc("DivRecMEM_MC")

# Import seurat object
seurat.memory <- readRDS("./misc_data/seurat_object_all_exps.rds")

# create factor containing identity info
mc@mc %>% 
  as.factor %>% 
  fct_collapse(mult.high = c('18','8','6'),
               mult.low = c('2','11','14')) -> mult.highVlow

# Add information as metadata
seurat.memory <- AddMetaData(object = seurat.memory, metadata = mult.highVlow, col.name = "mult.highVlow")

# Set Identity for DE-analysis
Idents(seurat.memory) <- "mult.highVlow"

# Perform DE-analysis
marks <- FindMarkers(object = seurat.memory , ident.1 = 'mult.low', ident.2 = 'mult.high', 
                     logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1, slot = "counts")

# Save file
write_rds(x = marks, file = "./misc_data/markers_top3VSbot3.rds")


################### Code below to generate the volcano plot

marks <- read_rds(file = "./misc_data/markers_top3VSbot3.rds")

# Prepare data for plotting
marks %>% 
  rownames_to_column("genes") %>% 
  # add significance thresholds
  mutate(sig = case_when((avg_log2FC < -0.1 | avg_log2FC > 0.1) & p_val_adj < 0.05 ~ "sig", 
                         TRUE ~ "ns") )%>% 
  mutate(sig = factor(sig, levels = c("sig", "ns"))) %>% 
  # select a couple of informative genes for highlighting and labeling
  mutate(to.label = genes %in% c("Il7r", "Tcf7", "Ltb", "Klf3", "Vim", "Lgals1", 
                                 "S100a4", "Nkg7","S100a6","Crip1", "Gzmb", "Ccr2")) -> to_plot

# Make the plot
ggplot(to_plot, aes( y= log10(p_val) , x=avg_log2FC, color = sig, label = genes) ) + 
  geom_point()+
  scale_y_reverse()+
  scale_color_manual(values=c("red", "darkgrey"))+
  geom_text_repel(data = subset(to_plot, to.label == TRUE ),box.padding = 1.4, max.overlaps = 25 )+
  labs(title = "")+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = log10(1.578786e-06), linetype = "dotted")+
  geom_vline(xintercept = c(0.1, -0.1), linetype = "dotted")+
  xlim(-.75,.4)
ggsave(filename = "./metacell_figs/DE_TOP3vsBOT3_grey.pdf", device = "pdf", width = 4.2, height = 6, useDingbats = F)



################### Code below to generate immune gene heatmaps

# Import immune related genes
immune.genes <- read_rds("./misc_data/immune_genes.rds")

# Select significantly up and downregulated immune genes
marks %>% 
  rownames_to_column("genes") %>% 
  filter(p_val < 0.05 & (avg_log2FC) > 0.05) %>% 
  filter(genes %in% immune.genes) %>% 
  slice_max(avg_log2FC, n = 100) %>% 
  pull(genes) -> gs.up
marks %>% 
  rownames_to_column("genes") %>% 
  filter(p_val < 0.05 & (avg_log2FC) < -0.05) %>% 
  filter(genes %in% immune.genes) %>%
  slice_min(avg_log2FC, n = 100) %>% 
  pull(genes) -> gs.down

# Make log2 enrichment table
lfp <- as.data.frame(log2(mc@mc_fp))

# Plot the heatmap containing genes UP in low-division Tcm
pdf("./metacell_figs/DE_analysis/heatmap_Up_genes.pdf", height = 5, width = 4)
as.data.frame(lfp) %>% 
  dplyr::select(one_of(c('2','14','11','18','8','6'))) %>% 
  rownames_to_column("genes") %>%
  dplyr::filter(genes %in% gs.up) %>% 
  column_to_rownames("genes") %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",breaks = seq(-1.5,1.5,by=0.03), clustering_distance_rows = "euclidean", 
           border_color = NA, cutree_cols = 2, show_rownames = T, cluster_cols = F)
dev.off()

# Plot the heatmap containing genes DOWN in low-division Tcm
pdf("./metacell_figs/DE_analysis/heatmap_Down_genesNew.pdf", height = 5, width = 4)
as.data.frame(lfp) %>% 
  dplyr::select(one_of(c('2','14','11','18','8','6'))) %>% 
  rownames_to_column("genes") %>%
  dplyr::filter(genes %in% gs.down) %>% 
  column_to_rownames("genes") %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",breaks = seq(-1.5,1.5,by=0.03), clustering_distance_rows = "euclidean", 
           border_color = NA, cutree_cols = 2, show_rownames = T, cluster_cols = F)
dev.off()


################### Code below to perform and plot the GSEA

# Get C7 pathways
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "C7"))

# convert pathways names + genes to named list
pathways <- split(pathways.hallmark[, 5], pathways.hallmark[, 3])

# use significant genes for gene-ranking
marks %>% 
  rownames_to_column("genes") %>% 
  filter(p_val < 0.05) %>%  
  dplyr::select(genes, avg_log2FC) %>% deframe %>% 
  sort(decreasing = T) -> stats

# Perform GSEA
fgseaRes <- fgsea(pathways=pathways, stats=stats, minSize = 10)

# filter on Kaech and goldrath pathways and edit names for plotting
fgseaRes %>% 
  filter(padj < 0.05)  %>% 
  filter(grepl("KAECH|GOLDRATH", pathway)) %>% 
  mutate(pathway = gsub("GOLDRATH", "GR", pathway))%>% 
  mutate(pathway = gsub("KAECH", "KA", pathway))%>% 
  mutate(pathway = gsub("_CD8_TCELL", "", pathway))%>% 
  mutate(pathway = gsub("_", " ", pathway)) %>% 
  mutate(pathway = fct_reorder(as.factor(pathway), NES)) %>% 
# Make the plot
ggplot( aes(x = NES, y = pathway))+
  geom_point(aes(size = size, color  = padj))+
  grids(axis = "y", linetype = "dashed")+
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted")
ggsave(filename = "./metacell_figs/DE_analysis/Pathways_C7.pdf",device =  "pdf", height = 4, width = 5)


# Make example plot Goldrath pathway
pw <-  "GOLDRATH_NAIVE_VS_MEMORY_CD8_TCELL_DN"
nes <- fgseaRes %>%  filter(pathway == pw) %>%  pull(NES) %>% round(digits = 2)
plotEnrichment(pathway = pathways[[pw]], stats)+
  ggtitle(pw)+
  geom_text(aes(label = paste0(c("NES = ", as.character(nes)), collapse = ""), x = Inf, y = Inf  ), vjust = "inward", hjust = "inward")
ggsave(filename = "./metacell_figs/GSEA_GOLDRATH_EFF.pdf", device = "pdf", width = 5, height = 3 ,useDingbats = F)

# Make example plot Kaech pathway
pw <-  "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_DN"
nes <- fgseaRes %>%  filter(pathway == pw) %>%  pull(NES) %>% round(digits = 2)
plotEnrichment(pathway = pathways[[pw]], stats)+
  ggtitle(pw)+
  geom_text(aes(label = paste0(c("NES = ", as.character(nes)), collapse = ""), x = Inf, y = Inf  ), vjust = "inward", hjust = "inward")
ggsave(filename = "./metacell_figs/GSEA_KAECH_EFF.pdf", device = "pdf", width = 5, height = 3 ,useDingbats = F)





# Get Hallmark pathways
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "H"))

# convert pathways names + genes to named list
pathways <- split(pathways.hallmark[, 5], pathways.hallmark[, 3])

# use significant genes for gene-ranking
marks %>% 
  rownames_to_column("genes") %>% 
  filter(p_val < 0.05) %>%  
  dplyr::select(genes, avg_log2FC) %>% deframe %>% 
  sort(decreasing = T) -> stats

# Perform GSEA
fgseaRes <- fgsea(pathways=pathways, stats=stats, minSize = 14 )

# plot significant pathways and edit names
fgseaRes %>% 
  filter(padj < 0.05)  %>% 
  mutate(pathway = gsub("HALLMARK", "", pathway))%>% 
  mutate(pathway = gsub("_", " ", pathway)) %>% 
  mutate(pathway = fct_reorder(as.factor(pathway), NES)) %>% 
# Make the plot
ggplot( aes(x = NES, y = pathway))+
  geom_point(aes(size = size, color = padj))+
  grids(axis = "y", linetype = "dashed")+
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted")
ggsave(filename = "./metacell_figs/DE_analysis/Pathways_H.pdf",device =  "pdf", height = 4, width = 5)


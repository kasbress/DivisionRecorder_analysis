library(tidyverse)
library(metacell)
library(ggpubr)
library(rstatix)

# Set working directory
setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_memory/")

# Import quiescence gene sets
gs.pos <- read.table(file = "./misc_data/Quiescence_pos_sig.txt")
gs.pos <- as.vector(gs.pos$V1)
gs.neg <- read.table(file = "./misc_data/Quiescence_neg_sig.txt")
gs.neg <- as.vector(gs.neg$V1)

# Initiate metacell database
scdb_init("./metacell_db/", force_reinit=T) 

# import data
mc <- scdb_mc("DivRecMEM_MC")
lfp <- log2(mc@mc_fp)


# Calculate score positive gene-set
lfp %>% 
  as.data.frame %>% 
  rownames_to_column("genes") %>% 
  filter(genes %in% gs.pos) %>% 
  dplyr::select(-genes) %>% 
  colSums() -> pos.score

# Calculate score negative gene-set
lfp %>% 
  as.data.frame %>% 
  rownames_to_column("genes") %>% 
  filter(genes %in% gs.neg) %>% 
  dplyr::select(-genes) %>% 
  colSums() -> neg.score

# Combine into dataframe
left_join(
  enframe(pos.score, name = "MCs", value = "Pos.Score"),
  enframe(neg.score, name = "MCs", value = "Neg.Score")) %>% 
  # Add phenotype cluster info as factor, and calculate combined score
  mutate(cluster = factor(fct_collapse(.$MCs,
                                       eff = c('23','15','22','8','6','10','4','5','7','18'),
                                       mult = c('9','1','12','13','14','16','11','2','3'),
                                       term = c('21','19','17','20')), 
                          levels = c("term", "eff", "mult")  ),
         Qscore = pos.score - neg.score) %>% 
  # Remove Tem 
  filter(!(cluster == "term")) %>% 
  mutate(cluster = fct_drop(.$cluster)) -> to.plot
  
# Make the waterfall plot
p1 <- ggplot(to.plot, aes(x = fct_reorder(MCs, Qscore), y = Qscore, fill = cluster ))+
  geom_bar(stat = "identity", color = "black", size = .4)+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("#4575B4" ,"#8AB8D7"))+
  scale_y_continuous(limits = c(-2, 1.7))

# Calculate statistics for boxplot
stat.test <- to.plot %>% 
  t_test( Qscore ~ cluster) %>% 
  add_xy_position(x = "cluster") 

# Make the boxplot
p2 <- ggplot(to.plot, aes(x = cluster, y = Qscore, color = cluster ))+
  geom_boxplot(outlier.shape = NA, size = 8)+
  scale_color_manual(values = c("#4575B4" ,"#8AB8D7"))+
  geom_jitter(color = "black", width = .2)+
  stat_pvalue_manual(data = stat.test, label = "p", label.size = 4, hide.ns = T, tip.length = 0.02)+
  scale_y_continuous(limits = c(-2, 1.7))

# Plot them together 
ggarrange(plotlist = list(p1,p2), ncol =  2, widths = c(2.2,1))
ggsave("./metacell_figs/Quiescence/Barchart_Qscore.pdf", device = "pdf", width = 8, height = 4)

library(tidyverse)
library(metacell)
library(rstatix)
library(ggpubr)

# Set working directory
setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_memory/")

# Initiate metacell database
scdb_init("./metacell_db/", force_reinit=T)

# Import data
mc <- scdb_mc("DivRecMEM_MC")
lfp <- log2(mc@mc_fp)

# Set genes to plot
select.genes <- c("Klrg1", "Cx3cr1", "Gzma", "Zeb2" , "Prdm1", "Tbx21", "Bcl2", "Cd27", "Xcl1", "Ccr7", "Il7r", "Sell")

# Prepare data for plotting
 as.data.frame(lfp) %>% 
  rownames_to_column("genes") %>% 
  dplyr::filter(genes %in% select.genes) %>% 
  mutate(genes = factor(genes, levels = select.genes)) %>% 
  pivot_longer(cols = -genes, names_to = "metacell", values_to =  "lfp") %>% 
  mutate(cluster = factor(fct_collapse(.$metacell,
                                       eff = c('23','15','22','8','6','10','4','5','7','18'),
                                       mult = c('9','1','12','13','14','16','11','2','3'),
                                       term = c('21','19','17','20')), levels = c("term", "eff", "mult")  ) ) -> to.plot
  
# Perform statistical test and save as object
to.plot %>%
  group_by(genes) %>%
  tukey_hsd(lfp ~ cluster)%>%
  add_y_position(scales = "free_y", fun = "max", step.increase = 0.1) -> stat.test
stat.test
  

# Make the plots
ggplot(to.plot, aes(x = cluster, y = lfp, color = cluster ))+
  geom_boxplot(outlier.shape = NA)+
  scale_color_manual(values = c("#C1272D","#4575B4" ,"#8AB8D7"))+
  geom_jitter(color = "black", width = .2, size = .8)+
  facet_wrap(~genes, scales = "free", nrow = 2)+
  stat_pvalue_manual(data = stat.test, label = "p.adj", label.size = 3, hide.ns = T, tip.length = 0.02) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  theme(strip.background = element_blank())
  
ggsave("./metacell_figs/phenotype_clustering/boxplots_lfp.pdf", device = "pdf", width = 8, height = 4)

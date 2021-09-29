library(tidyverse)
library(metacell)

# Set working directory
setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_memory/")

## import data
scdb_init("./metacell_db/", force_reinit=T) 
mc <- scdb_mc("DivRecMEM_MC")
clean_mat <- scdb_mat("DivRecMEM_clean")

# Get metacell annotations and sample hashtags and merge to dataframe
clean_mat@cell_metadata %>% 
  rownames_to_column("cellcode") %>% 
  dplyr::select(hash.ID, cellcode) %>% 
  mutate(hash.ID = as.character(hash.ID)) %>% 
  right_join(enframe(mc@mc, "cellcode", "MC")) %>% 
  mutate(MC = factor(MC))-> df


# select hash.tags from relevant experiments
hash.tags <- grep(pattern = "DR[4|5|6|7]\\.GFP|Naive", x = unique(df$hash.ID), value = T)

## Get labeling ratios
counts.table <- df %>%
  # first filter unused cells and for HTO1-4, remove NAs
  dplyr::filter(!(cellcode %in% as.character(clean_mat@ignore_cells))) %>%
  dplyr::filter(hash.ID %in% hash.tags & !is.na(MC)) %>%
  # count MCs by HTO
  dplyr::count(MC, hash.ID) %>%
  # perform the normalization within hash.IDs
  group_by(hash.ID)%>%
  mutate(normalized.count = (n/sum(n))*1000 ) %>%
  ungroup() %>% 
  # Add Tcm cluster information
  mutate(cluster = fct_collapse(.$MC,
                                eff = c('23','15','22','8','6','10','4','5','7','18'),
                                mult = c('9','1','12','13','14','16','11','2','3'),
                                term = c('21','19','17','20'))  ) %>% 
  # Add general phenotype info
  mutate(phenotype = fct_collapse(.$MC,
                                TCM = c('23','15','22','8','6','10','4','5','7','18','9','1','12','13','14','16','11','2','3'),
                                TEM = c('21','19','17','20'))  ) %>% 
  # Add sample origin info
  mutate(origin = fct_collapse(.$hash.ID,
                              Naive = c('Naive1', 'Naive2', 'Naive3', 'Naive4'),
                              Pre = c('DR4.GFP','DR5.GFP','DR6.GFP','DR7.GFP'))  ) 


# import some colors
cols <- c(brewer.pal(8, "Blues")[5:8], brewer.pal(8, "Reds")[5:8])

# Get p values
stat.test <-counts.table %>% 
  group_by(hash.ID, phenotype, origin) %>%  
  summarise(summed.count = sum(normalized.count)) %>% 
  group_by(phenotype) %>% 
  t_test(summed.count ~ origin ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  add_xy_position(fun = "mean_sd", x = "phenotype", dodge = 0.9) %>% 
  mutate(p.adj = round(p.adj, 4))

# Make barplot
counts.table %>% 
  group_by(hash.ID, phenotype, origin) %>% 
  summarise(summed.count = sum(normalized.count)) %>% 
ggplot( aes(x = phenotype, y = summed.count) )+
  geom_bar( stat = "summary", fun = "mean", position = "dodge", color = "black", aes(fill = origin))+
  scale_fill_manual(values= cols[c(3,7)] )+
  stat_summary(geom = "errorbar", fun = mean,
               fun.min = function(x) mean(x), 
               fun.max = function(x) mean(x) + sd(x), 
               width = 0.4, position = position_dodge(width = 0.9), aes(group = origin) )+
  labs(x = "Phenotype", y = "Normalized cell count")+
  geom_jitter( position = position_dodge(width = 0.9) ,aes(group = origin))+
  stat_pvalue_manual(data = stat.test,  label = "p.adj", tip.length = 0.01,hide.ns = TRUE, label.size = 4 )

ggsave("./metacell_figs/comparison_naiveVSpre/GroupedBarchart_TCM_TEM.pdf", width = 4, height = 4)

# Get p values - Tcm only
stat.test <- counts.table %>% 
  dplyr::filter(phenotype == "TCM") %>% 
  group_by(hash.ID) %>% 
  mutate(norm.norm.counts = (normalized.count/sum(normalized.count))*1000 ) %>% 
  mutate(MC = factor(MC, levels = unique(MC))) %>% 
  group_by(MC) %>% 
  t_test(norm.norm.counts ~ origin ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(fun = "mean_sd", x = "MC", dodge = 0.9) %>% 
  mutate(p.adj = round(p.adj, 2))

# Make barplot - Tcm only
counts.table %>% 
  dplyr::filter(phenotype == "TCM") %>% 
  group_by(hash.ID) %>% 
  mutate(norm.norm.counts = (normalized.count/sum(normalized.count))*1000 ) %>% 
  ggplot( aes(x = factor(MC, levels = unique(MC)), y = norm.norm.counts) )+
  geom_bar( stat = "summary", fun = "mean", position = "dodge", color = "black", aes(fill = origin))+
  scale_fill_manual(values= cols[c(3,7)] )+
  stat_summary(geom = "errorbar", fun = mean,
               fun.min = function(x) mean(x), 
               fun.max = function(x) mean(x) + sd(x), 
               width = 0.4, position = position_dodge(width = 0.9), aes(group = origin) )+
  labs(x = "Metacell", y = "Normalized cell count")+
  geom_jitter( position = position_dodge(width = 0.9), aes(group = origin))+
  stat_pvalue_manual(data = stat.test,  label = "p.adj", tip.length = 0.01,hide.ns = F, label.size = 4 )

ggsave("./metacell_figs/comparison_naiveVSpre/GroupedBarchart_TCM_MCs.pdf", width = 10, height = 4)

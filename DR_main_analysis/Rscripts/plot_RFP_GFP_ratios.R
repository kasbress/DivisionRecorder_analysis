library(tidyverse)
library(metacell)
library(ggpubr)
library(rstatix)

# Set working directory
setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_memory/")

# Initiate metacell database
scdb_init("./metacell_db/", force_reinit=T) 

# Import data
mc <- scdb_mc("DivRecMEM_MC")
clean_mat <- scdb_mat("DivRecMEM_clean")


# Get metacell annotations and sample hashtags and merge to dataframe
clean_mat@cell_metadata %>% 
  rownames_to_column("cellcode") %>% 
  dplyr::select(hash.ID, cellcode) %>% 
  mutate(hash.ID = as.character(hash.ID)) %>% 
  right_join(enframe(mc@mc, "cellcode", "MC")) %>% 
  mutate(MC = factor(MC))-> df

## select hash.tags from relevant experiments
hash.tags <- grep(pattern = "GFP|RFP", x = unique(df$hash.ID), value = T)



################################################################
##### USED THIS FOR PLOTTING

## Get labeling ratios
lab.ratio <- df %>%
  # first filter unused cells and for HTO1-4, remove NAs
  filter(!(cellcode %in% as.character(clean_mat@ignore_cells))) %>%
  filter(hash.ID %in% hash.tags & !is.na(MC)) %>%
  # count MCs by HTO
  dplyr::count(MC, hash.ID)%>%
  # Add in ratio's of DR-RFP samples as correction factor | compensate for cell count inequalities
  mutate(div.by = case_when( hash.ID %in% grep(pattern = "DR[1|2|3]|DR45|DR67", x = hash.ID, value = T) ~ 1, 
                             hash.ID == "DR5.GFP" ~ 3.4,
                             hash.ID == "DR6.GFP" ~ 3.77, 
                             TRUE ~ 1  ) ) %>%
  # correct and re-name GFP samples of the same day
  mutate(cor.count = n * div.by)%>%
  mutate(hash.ID = fct_collapse(.$hash.ID,DR45.GFP = c("DR4.GFP","DR5.GFP"),DR67.GFP = c("DR6.GFP","DR7.GFP")  ))%>%
  # pool the GFP samples of the same day
  group_by(MC, hash.ID)%>%
  summarise(cor.count = sum(cor.count))%>%
  # perform the normalization within hash.IDs
  group_by(hash.ID)%>%
  mutate(normalized.count = (cor.count/sum(cor.count))*1000 )%>%
  mutate(DivRec = ifelse(grepl( "RFP",hash.ID), "RFP", "GFP" )) %>% 
  # add together per MC and GFP/RFP
  group_by(MC, DivRec)%>%
  summarise(normalized.count = sum(normalized.count)) %>% 
  # Divide RFP by GFP per metacell
  group_by(MC)%>%
  summarise(lab.ratio = sum(normalized.count[DivRec == "RFP"])/normalized.count[DivRec == "GFP"])  %>% 
  ungroup %>% 
  mutate(cluster = factor(fct_collapse(.$MC,
                                eff = c('23','15','22','8','6','10','4','5','7','18'),
                                mult = c('9','1','12','13','14','16','11','2','3'),
                                term = c('21','19','17','20')), levels = c("term", "eff", "mult")  )  )




# first make the waterfall plot
p1 <- ggplot(lab.ratio, aes(x = fct_reorder(MC, lab.ratio), y = log2(lab.ratio), fill = cluster ))+
  geom_bar(stat = "identity", color = "black", size = .4)+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("#C1272D","#4575B4" ,"#8AB8D7"))+
  scale_y_continuous(limits = c(-1.5, 0.6))


# calculate statistics for the boxplot
stat.test <- lab.ratio %>% 
  tukey_hsd( lab.ratio ~ cluster) %>% 
  add_xy_position(x = "cluster") %>% 
  mutate(y.position = log2(y.position))

# Make the boxplot
p2 <- ggplot(lab.ratio, aes(x = cluster, y = log2(lab.ratio), color = cluster ))+
  geom_boxplot(outlier.shape = NA)+
  scale_color_manual(values = c("#C1272D","#4575B4" ,"#8AB8D7"))+
  geom_jitter(color = "black", width = .2)+
  stat_pvalue_manual(data = stat.test, label = "p.adj", label.size = 4, hide.ns = T, tip.length = 0.02)+
  scale_y_continuous(limits = c(-1.5, 0.6))

# Plot them together
ggarrange(plotlist = list(p1,p2), ncol =  2, widths = c(1.8,1))

ggsave("./metacell_figs/labeling_ratios_AllWterm.pdf", device = "pdf", width = 8, height = 4)

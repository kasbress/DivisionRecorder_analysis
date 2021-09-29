library(tidyverse)
library(metacell)
library(tidytext)

# Set working directory
setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_memory/")

# Initiate metacell database
scdb_init("./metacell_db/", force_reinit=T)

# Import data
mc <- scdb_mc("DivRecMEM_MC_fin")
lfp <- log2(mc@mc_fp)

# Import immune-related genes
immune.genes <- read_rds("./misc_data/immune_genes.rds")

## Prepare data for plotting
as.data.frame(lfp) %>% 
  select(one_of(c('2','11','14', '6','8','18'))) %>% #select MCs for plotting
  rownames_to_column("genes") %>%
  filter(genes %in% immune.genes) %>%  #either all genes, or immune genes
  pivot_longer(cols = -genes, names_to = "MC", values_to =  "lfp") %>%   #switch to long
  mutate(MC = factor(MC, levels =c('2','11','14','13', '6','8','23','18'))) %>% #set ordering
  group_by(MC) %>% 
  filter(dense_rank(lfp) <= 5 | dense_rank(desc(lfp)) <= 5) %>% #get top and bottom 5
  ungroup() %>% 
  mutate(genes = reorder_within(genes, lfp, MC) ) %>% #reorder genes by lfp, for each MC
### pipe into plot
ggplot(aes(x = genes, y = lfp, fill=lfp > 0))+
  geom_bar(stat="identity")+ 
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("blue", "red"))+
  scale_x_reordered()+
  coord_flip()+
  facet_wrap(~MC, scales = "free", nrow = 2)

# Save the plot
ggsave(filename = "./metacell_figs/Waterfall_Marks.pdf", device = "pdf", width = 5, height = 4.5 ,useDingbats = F)


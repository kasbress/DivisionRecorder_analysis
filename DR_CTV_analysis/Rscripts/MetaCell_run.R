library(metacell)
library(Seurat)
library(tidyverse)
library(viridis)
library(SingleCellExperiment)


setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_CTV_experiment/")


### START METACELL
if(!dir.exists("./metacell_db/")) dir.create("./metacell_db/") # going to make new folder, in case it doe not exists - DATA
scdb_init("./metacell_db/", force_reinit=T) # tell R that this is a "workdirectory"  - RUN THIS EVERYTIME YOU START

if(!dir.exists("./metacell_figs")) dir.create("./metacell_figs") # new folder for FIGURES
scfigs_init("./metacell_figs/") 

seurat.memory <- readRDS("./misc_data/seurat_object_CTVexp.rds")
## Generate metacell mat object
sce = as.SingleCellExperiment(seurat.memory)
mat = scm_import_sce_to_mat(sce)
scdb_add_mat(id = "DivRecCTV",  mat = mat)

## filter small cells, Adjust where nessecary
mcell_plot_umis_per_cell("DivRecCTV",min_umis_cutoff = 2000)  
mcell_mat_ignore_small_cells("DivRecCTV", "DivRecCTV", 2000)


## Clean-up
remove(sce)
remove(seurat.memory)


## Get Mito genes and lateral genes
mito_genes <- grep(pattern = "^mt-", x = mat@genes, value = TRUE)

test <- rep(1,length(grep(pattern = "^Hsp|^Trac|Trdc|Trav|^Rpl|^Rps|Mki67|Top2a|Cks1b|E2f1|Prc1", x = mat@genes, value = TRUE)))
names(test) <- grep(pattern = "^Hsp|^Trac|Trdc|Trav|^Rpl|^Rps|Mki67|Top2a|Cks1b|E2f1|Prc1", x = mat@genes, value = TRUE)
scdb_add_gset("lateral",gset_new_gset(test, "lateral"))

## Get Mitofractions
uc = Matrix::colSums(mat@mat)
mito_f = Matrix::colSums(mat@mat[mito_genes, ]) / uc

## Plot mito_fractions, adjust treshold where nessecary
plot(mito_f)
abline(h = 0.12, col="red", lwd=3, lty=2)
text(1000, 0.3, "0.12", col = "red", cex = 2)



## CLEAN UP MAT OBJECT
mcell_mat_ignore_cells(new_mat_id = "DivRecCTV_clean", mat_id = "DivRecCTV", 
                       union(mat@ignore_cells, names(uc)[mito_f >= 0.12 ]))
clean_mat = scdb_mat("DivRecCTV_clean")
mcell_mat_ignore_cells("DivRecCTV_clean", "DivRecCTV_clean", 
                       union(clean_mat@ignore_cells, names(uc)[clean_mat@cell_metadata$HTO_classification.global != "Singlet" ]))
clean_mat = scdb_mat("DivRecCTV_clean")

clean_mat@cell_metadata %>% 
  rownames_to_column %>% 
  filter(!(rowname %in% clean_mat@ignore_cells)) %>% 
  pull(nFeature_RNA) %>% 
  plot
abline(h = 2800, col="red", lwd=3, lty=2)

clean_mat@cell_metadata %>% 
  rownames_to_column  %>% 
  filter(!(rowname %in% clean_mat@ignore_cells)) %>% 
  filter(nFeature_RNA >2800) %>% 
  pull(rowname) -> to.ignore

mcell_mat_ignore_cells("DivRecCTV_clean", "DivRecCTV_clean", 
                       union(clean_mat@ignore_cells,  to.ignore))
clean_mat = scdb_mat("DivRecCTV_clean")


clean_mat@cell_metadata %>% rownames_to_column %>% filter(!(rowname %in% clean_mat@ignore_cells)) %>% pull(nFeature_RNA) %>% plot 
clean_mat@cell_metadata %>% rownames_to_column %>% filter(!(rowname %in% clean_mat@ignore_cells)) %>% pull(hash.ID) %>% table 


## make gstat
mcell_add_gene_stat(mat_id = "DivRecCTV_clean", gstat_id = "DivRecCTV_gs", force = T)
gstat <- scdb_gstat("DivRecCTV_gs")
dim(gstat)


# generate feats_gset and plot stats
mcell_gset_filter_varmean(gstat_id = "DivRecCTV_gs", gset_id = "DivRecCTV_feats", T_vm=0.06, force_new=T)
mcell_gset_filter_cov(gstat_id = "DivRecCTV_gs", gset_id = "DivRecCTV_feats", T_tot=80, T_top3=2)


feats_gset <- scdb_gset("DivRecCTV_feats")

length(names(feats_gset@gene_set))
mcell_plot_gstats(gstat_id = "DivRecCTV_gs", "DivRecCTV_feats")



# filter feats gset for lateral genes 
lateral_gset = scdb_gset("lateral")
feats_gset = gset_new_restrict_gset(feats_gset, lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
length(names(feats_gset@gene_set))
scdb_add_gset(id = "DivRecCTV_feats_filt", gset = feats_gset)

# Run MetaCell
mcell_add_cgraph_from_mat_bknn(mat_id = "DivRecCTV_clean", gset_id = "DivRecCTV_feats_filt", graph_id = "DivRecCTV_graph", K=150, dsamp=T)
mcell_coclust_from_graph_resamp(coc_id = "DivRecCTV_coc400", graph_id = "DivRecCTV_graph", min_mc_size=100, p_resamp=0.75, n_resamp=400)
mcell_mc_from_coclust_balanced(mc_id = "DivRecCTV_MC", coc_id =  "DivRecCTV_coc400", "DivRecCTV_clean", K=80, min_mc_size=80, alpha=1.8)

# Run MetaCell
mcell_add_cgraph_from_mat_bknn(mat_id = "DivRecCTV_clean", gset_id = "DivRecMEM_feats", graph_id = "DivRecCTV_graph2", K=150, dsamp=T)
mcell_coclust_from_graph_resamp(coc_id = "DivRecCTV_coc400_2", graph_id = "DivRecCTV_graph2", min_mc_size=100, p_resamp=0.75, n_resamp=400)
mcell_mc_from_coclust_balanced(mc_id = "DivRecCTV_MC2", coc_id =  "DivRecCTV_coc400_2", "DivRecCTV_clean", K=70, min_mc_size=80, alpha=2)


## Add colors
mc = scdb_mc("DivRecCTV_MC")
length(names(mc@mc))
length((mc@annots))
table(mc@mc)
mc@colors <- viridis(length(mc@annots))
scdb_add_mc("DivRecCTV_MC",mc)

## 2d graph 
mc <- scdb_mc("DivRecCTV_MC")
mcell_mc2d_force_knn(mc2d_id = "DivRecCTV_MC" ,mc_id =  "DivRecCTV_MC", "DivRecCTV_graph", ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",800, "metacell")
tgconfig::set_param("mcell_mc2d_width",800, "metacell")
mcell_mc2d_plot("DivRecCTV_MC")


lfp <- log2(mc@mc_fp)

imm.genes <- readRDS("./misc_data/immune_genes.rds")
lfp.p <- na.omit(lfp[imm.genes[imm.genes %in% row.names(lfp)],])

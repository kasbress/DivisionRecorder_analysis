library(metacell)
library(Seurat)
library(tidyverse)
library(viridis)
library(SingleCellExperiment)

# Set working directory
setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_memory/")

# Initiate metacell database and figure repository
if(!dir.exists("./metacell_db/")) dir.create("./metacell_db/") 
scdb_init("./metacell_db/", force_reinit=T) 

if(!dir.exists("./metacell_figs")) dir.create("./metacell_figs") 
scfigs_init("./metacell_figs") 

# Import Seurat object containing all memory data
seurat.memory <- readRDS("./misc_data/seurat_object_all_exps.rds")

# Generate metacell mat object and add to database
sce = as.SingleCellExperiment(seurat.memory)
mat = scm_import_sce_to_mat(sce)
scdb_add_mat(id = "DivRecMEM",  mat = mat)

# filter small cells, Adjust where nessecary
mcell_plot_umis_per_cell("DivRecMEM",min_umis_cutoff = 1500)  
mcell_mat_ignore_small_cells("DivRecMEM", "DivRecMEM", 1500)

# Clean-up
remove(sce)
remove(seurat.memory)


# Get Mito genes
mito_genes <- grep(pattern = "^mt-", x = mat@genes, value = TRUE)

# Get Mitofractions
uc = Matrix::colSums(mat@mat)
mito_f = Matrix::colSums(mat@mat[mito_genes, ]) / uc

# Plot mito_fractions, adjust threshold where necessary
plot(mito_f)
abline(h = 0.12, col="red", lwd=3, lty=2)
text(1000, 0.2, "0.12", col = "red", cex = 2)

## Clean up mat object
# First remove cells with a high mito fraction
mcell_mat_ignore_cells(new_mat_id = "DivRecMEM_clean", mat_id = "DivRecMEM", 
                       union(mat@ignore_cells, names(uc)[mito_f >= 0.12 ]))
clean_mat = scdb_mat("DivRecMEM_clean")
# Then trash all non-singlets
mcell_mat_ignore_cells("DivRecMEM_clean", "DivRecMEM_clean", 
                       union(clean_mat@ignore_cells, names(uc)[clean_mat@cell_metadata$HTO_classification.global != "Singlet" ]))
clean_mat = scdb_mat("DivRecMEM_clean")


# Check data for doublet cells ect. 
clean_mat@cell_metadata %>% 
  rownames_to_column %>% 
  filter(!(rowname %in% clean_mat@ignore_cells)) %>% 
  pivot_longer(cols = c(nCount_RNA,nFeature_RNA), names_to = "metric", values_to = "value") %>%
  ggplot(aes((value), fill = experiment, color = experiment))+
  geom_density(alpha = 0.1)+
  facet_wrap(~metric, scales = "free")


### Remove large cells
clean_mat@cell_metadata %>% rownames_to_column %>% filter(experiment == "Mem.exp.1") %>% filter(nFeature_RNA >3000) %>% pull(rowname) -> to.ignore
clean_mat@cell_metadata %>% rownames_to_column %>% filter(experiment == "Mem.exp.2.1") %>% filter(nFeature_RNA >2500) %>% pull(rowname) %>% union(to.ignore) -> to.ignore
clean_mat@cell_metadata %>% rownames_to_column %>% filter(experiment == "Mem.exp.2.2") %>% filter(nFeature_RNA >3000) %>% pull(rowname) %>% union(to.ignore) -> to.ignore
mcell_mat_ignore_cells("DivRecMEM_clean", "DivRecMEM_clean", 
                       union(clean_mat@ignore_cells,  to.ignore))
clean_mat = scdb_mat("DivRecMEM_clean")

### Remove small cells
clean_mat@cell_metadata %>% rownames_to_column %>% filter(experiment == "Mem.exp.1") %>% filter(nFeature_RNA < 1200) %>% pull(rowname) -> to.ignore
clean_mat@cell_metadata %>% rownames_to_column %>% filter(experiment == "Mem.exp.2.1") %>% filter(nFeature_RNA < 800) %>% pull(rowname) %>% union(to.ignore) -> to.ignore
clean_mat@cell_metadata %>% rownames_to_column %>% filter(experiment == "Mem.exp.2.2") %>% filter(nFeature_RNA < 1000) %>% pull(rowname) %>% union(to.ignore) -> to.ignore
mcell_mat_ignore_cells("DivRecMEM_clean", "DivRecMEM_clean", 
                       union(clean_mat@ignore_cells,  to.ignore))
clean_mat = scdb_mat("DivRecMEM_clean")


# Check again 
clean_mat@cell_metadata %>% 
  rownames_to_column %>% 
  filter(!(rowname %in% clean_mat@ignore_cells)) %>% 
  pivot_longer(cols = c(nCount_RNA,nFeature_RNA), names_to = "metric", values_to = "value") %>%
  ggplot(aes((value), fill = experiment, color = experiment))+
  geom_density(alpha = 0.1)+
  facet_wrap(~metric, scales = "free")


### Make gstat object for first gene-gene correlations
mcell_add_gene_stat(mat_id = "DivRecMEM_clean", gstat_id = "DivRecMEM_gs", force = T)
gstat <- scdb_gstat("DivRecMEM_gs")
dim(gstat)

# generate feats_gset 
mcell_gset_filter_varmean(gstat_id = "DivRecMEM_gs", gset_id = "DivRecMEM_feats_all", T_vm=0.12, force_new=T)
mcell_gset_filter_cov(gstat_id = "DivRecMEM_gs", gset_id = "DivRecMEM_feats_all", T_tot=80, T_top3=2)

# Import feature gene-set
feats_gset <- scdb_gset("DivRecMEM_feats_all")


# Create correlation matrix
gene.anchors <- names(feats_gset@gene_set)
mcell_mat_rpt_cor_anchors(mat_id = "DivRecMEM_clean", gene_anchors = gene.anchors, cor_thresh = 0.1, gene_anti = c(),
                          tab_fn = "metacell_figs/modules_gmods.txt", sz_cor_thresh = 0.2)


# Read correlation matrix and extract focus genes
gcor.mat <- read.table("metacell_figs/modules_gmods.txt", header = T)
foc.genes <- apply(gcor.mat[,-1], 1, which.max)

# Add correlated genes as geneset to database
gset <- gset_new_gset(sets = foc.genes, desc = "Diff Expressed Genes")
scdb_add_gset("All_corr_diff_genes", gset)

# cluster correlated genes and identify modules, split into 20 clusters
mcell_mat_ignore_genes(new_mat_id = "DivRecMEM_clean_diff", mat_id = "DivRecMEM_clean",
                       ig_genes = names(foc.genes), reverse = T)
mcell_gset_split_by_dsmat(gset_id = "All_corr_diff_genes" , "DivRecMEM_clean_diff", K = 20)
mcell_plot_gset_cor_mats(gset_id = "All_corr_diff_genes", scmat_id = "DivRecMEM_clean_diff")

## 4= Ribosomal
## 7 = immune genes (Macrophage)
## 8 = IFN genes
## 9 = Replication/Stress?
## 10 = Cell Cycle
## 11 = PD-1, CTLA4
## 12 = Ribosomal
## 13 = Macrophage? Monocyte?
## 14 = Ribosomal
## 16 = Large cluster, class II, B cells?
## 17 = Hemoglobin
## 18 = Myeloid?
## 20 = mitochondrial



# Import gene-set
feats_gset <- scdb_gset("All_corr_diff_genes")

## Select gene-clusters to mask
to.lateral <- feats_gset@gene_set[feats_gset@gene_set %in% c(4, 8,10,12,14,20)]
gset <- gset_new_gset(sets = to.lateral, desc = "cell cycle, stress, ribosomal, mito")
scdb_add_gset("lateral", gset)


# generate feats_gset and plot stats
mcell_gset_filter_varmean(gstat_id = "DivRecMEM_gs", gset_id = "DivRecMEM_feats_preCluster", T_vm=0.12, force_new=T)
mcell_gset_filter_cov(gstat_id = "DivRecMEM_gs", gset_id = "DivRecMEM_feats_preCluster", T_tot=80, T_top3=2)

# Load gene-set
feats_gset <- scdb_gset("DivRecMEM_feats_preCluster")


# filter feats gset for lateral genes (mask genes)
lateral_gset = scdb_gset("lateral")
feats_gset = gset_new_restrict_gset(feats_gset, lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
length(names(feats_gset@gene_set))
scdb_add_gset(id = "DivRecMEM_feats_preCluster", gset = feats_gset)

## Perform MetaCell clustering
mcell_add_cgraph_from_mat_bknn(mat_id = "DivRecMEM_clean", gset_id = "DivRecMEM_feats_preCluster", graph_id = "DivRecMEM_graph_preClust", K=100, dsamp=T)
mcell_coclust_from_graph_resamp(coc_id = "DivRecMEM_coc_preClust", graph_id = "DivRecMEM_graph_preClust", min_mc_size=80, p_resamp=0.75, n_resamp=250)
mcell_mc_from_coclust_balanced(mc_id = "DivRecMEM_MC_preClust", coc_id =  "DivRecMEM_coc_preClust", "DivRecMEM_clean", K=100, min_mc_size=30, alpha=2)

mcell_mc_split_filt(new_mc_id = "DivRecMEM_MCfilt", mc_id = "DivRecMEM_MC", mat_id = "DivRecMEM_clean_masked", T_lfc = 3)


# load MC object and add color
mc = scdb_mc("DivRecMEM_MC_preClust")
length(names(mc@mc))
length((mc@annots))
table(mc@mc)
mc@colors <- viridis(length(mc@annots))
scdb_add_mc("DivRecMEM_MC_preClust",mc)

# Make 2D plot
mcell_mc2d_force_knn(mc2d_id = "DivRecMEM_MC_preClust" ,mc_id =  "DivRecMEM_MC_preClust", "DivRecMEM_graph_preClust", ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",800, "metacell")
tgconfig::set_param("mcell_mc2d_width",800, "metacell")
mcell_mc2d_plot("DivRecMEM_MC_preClust")

# Inspect genes for potential strong batchy-ness
lfp <- log2(mc@mc_fp)



#### NEW CLUSTER


# Initiate metacell database and figure repository
scdb_init("./metacell_db_fin/", force_reinit=T) 
scfigs_init("./metacell_figs_fin") 

# Import gene-gene correlation clusters
feats_gset <- scdb_gset("All_corr_diff_genes")


# Create lateral geneset
# Here I included gene-gene clusters such as cell-division related, IFN signature, and Ribosomal
to.lateral <- feats_gset@gene_set[feats_gset@gene_set %in% c(4,8,9,10,12,14,20)]

# Genes below were extremely batchy between the two experiments
# I mask them from the analysis to reduce to some degree batch effects during the clustering
to.lateral[["Itga4"]] <- 3
to.lateral[["Xist"]] <- 3
to.lateral[["Gm42418"]] <- 3
to.lateral[["Ptprc"]] <- 3
to.lateral[["Malat1"]] <- 3

# Import the feature gene-list obtained within the first experiment
feats_gset <- scdb_gset("MSC_featsBIG")

# Sample some genes I don;t want to include in the clustering (TCR chains, Ribosomal, Heat-shock)
test <- rep(1,length(grep(pattern = "^Hsp|^Eif|Trdc|Trav|^Rpl|^Rps", x = names(feats_gset@gene_set), value = TRUE)))
names(test) <- grep(pattern = "^Hsp|^Eif|Trdc|Trav|^Rpl|^Rps", x = names(feats_gset@gene_set), value = TRUE)

# Add these to the lateral set
to.lateral <- c(test, to.lateral)
gset <- gset_new_gset(sets = to.lateral, desc = "cell cycle, stress, ribosomal, mito")

# Write out the lateral set
scdb_add_gset("lateral", gset)


# filter feats gset for lateral genes 
lateral_gset <- scdb_gset("lateral")
feats_gset = gset_new_restrict_gset(feats_gset, lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
length(names(feats_gset@gene_set))

# Add the final feature gene-set to the metacell database
scdb_add_gset(id = "DivRecMEM_feats", gset = feats_gset)


# MetaCell generation 
mcell_add_cgraph_from_mat_bknn(mat_id = "DivRecMEM_clean", gset_id = "DivRecMEM_feats", graph_id = "DivRecMEM_graph", K=150, dsamp=T)
mcell_coclust_from_graph_resamp(coc_id = "DivRecMEM_coc", graph_id = "DivRecMEM_graph", min_mc_size=100, p_resamp=0.75, n_resamp=350)
mcell_mc_from_coclust_balanced(mc_id = "DivRecMEM_MC", coc_id =  "DivRecMEM_coc", "DivRecMEM_clean", K=70, min_mc_size=600, alpha=2)



# Load in MetaCell object to add colors
mc = scdb_mc("DivRecMEM_MC")
length(names(mc@mc))
length((mc@annots))
table(mc@mc)
mc@colors <- viridis(length(mc@annots))
scdb_add_mc("DivRecMEM_MC",mc)


# Make 2D projection
mcell_mc2d_force_knn(mc2d_id = "DivRecMEM_MC" ,mc_id =  "DivRecMEM_MC", "DivRecMEM_graph", ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",800, "metacell")
tgconfig::set_param("mcell_mc2d_width",800, "metacell")
mcell_mc2d_plot("DivRecMEM_MC")

mcell_mc2d_force_knn(mc2d_id = "DivRecMEM_MC2" ,mc_id =  "DivRecMEM_MC2", "DivRecMEM_graph2", ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",800, "metacell")
tgconfig::set_param("mcell_mc2d_width",800, "metacell")
mcell_mc2d_plot("DivRecMEM_MC2")

mc1 = scdb_mc("DivRecMEM_MC")
mc2 = scdb_mc("DivRecMEM_MC2")
lfp1 <- log2(mc1@mc_fp)
lfp2 <- log2(mc2@mc_fp)

lfp <- log2(mc@mc_fp)

sort(lfp["Lyz2",])


mcell_mc2d_plot("DivRecMEM_MCfilt_col")

mcell_mc2d_plot_by_factor(mc2d_id = "DivRecMEM_MC", mat_id = "DivRecMEM_clean", meta_field = "experiment",single_plot = T,ncols = 3)
mcell_mc2d_plot_by_factor(mc2d_id = "DivRecMEM_MC", mat_id = "DivRecMEM_clean", meta_field = "hash.ID",single_plot = T,ncols = 3)




plt = function(gene1, gene2, lfp, colors) 
{
  plot(lfp[gene1, ], lfp[gene2, ], pch=21, cex=3, bg=colors, xlab=gene1, ylab=gene2)
  text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))
  
}

plt("lab_ratio1", "lab_ratio2", t(as.matrix(cell_counts_norm)), colors = mc@colors)


plt("Cd3e", "Cd8a", as.matrix(lfp), colors = mc@colors)
plt("Lgals1", "Ccr2", as.matrix(lfp), colors = mc@colors)
plt("Serpina3g", "S100a4", as.matrix(lfp), colors = mc@colors)
plt("Tim3", "Cd160", as.matrix(lfp), colors = mc@colors)
plt("S100a4", "Lgals1", as.matrix(lfp), colors = mc@colors)


plt("Rps29", "Lgals1", as.matrix(lfp), colors = mc@colors)

plt("Rpl38", "Lgals1", as.matrix(lfp), colors = mc@colors)

















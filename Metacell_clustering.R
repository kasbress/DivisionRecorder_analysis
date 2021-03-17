library(metacell)
library(Seurat)
library(viridis)

## load seurat object - imported variable name is countsfil
load("./Seuratdata.Rda")

## make output dirs for MetaCell
if(!dir.exists("./metacell_db/")) dir.create("./metacell_db/")
scdb_init("./metacell_db/", force_reinit=T)

if(!dir.exists("./metacell_figs")) dir.create("./metacell_figs")
scfigs_init("./metacell_figs") 


## Generate metacell mat object and add ADT data to metadata
sce = as.SingleCellExperiment(countsfil)
mat = scm_import_sce_to_mat(sce)

## Get ADT counts and add to mat object metadata
ADT_counts <-  data.frame(row.names = colnames(countsfil@assays$ADT@counts), 
                          cd27 = countsfil@assays$ADT@counts[1,], 
                          klrg1 = countsfil@assays$ADT@counts[2,])
mat@cell_metadata <- merge(mat@cell_metadata, ADT_counts, by = 0, all = T)
row.names(mat@cell_metadata) <- mat@cell_metadata[,1]
mat@cell_metadata <- mat@cell_metadata[,-1]

## add mat object to database
scdb_add_mat(id = "MSC",  mat = mat)

## Clean-up
remove(sce)
remove(countsfil)
remove(ADT_counts)

## filter small cells (less than 2000 UMIs)
mcell_plot_umis_per_cell("MSC")  
mcell_mat_ignore_small_cells("MSC", "MSC", 2000)


## Get Mito genes and lateral genes
mito_genes <- grep(pattern = "^mt-", x = mat@genes, value = TRUE)
bad_genes <- grep(pattern = "^Ig|^mt|Malat1", x = mat@genes, value = TRUE)
test <- rep(1,length(grep(pattern = "^Hsp|^Trac|Trdc|Trav|^Rpl|^Rps|Mki67|Top2a|Cks1b|E2f1|Prc1", x = mat@genes, value = TRUE)))
names(test) <- grep(pattern = "^Hsp|^Trac|Trdc|Trav|^Rpl|^Rps", x = mat@genes, value = TRUE)
scdb_add_gset("lateral",gset_new_gset(test, "lateral"))

## Get Mitofractions
uc = Matrix::colSums(mat@mat)
mito_f = Matrix::colSums(mat@mat[mito_genes, ]) / uc

## Clean up mat object
## Remove Ig, mt genes and Malat1
mcell_mat_ignore_genes("MSC_clean", "MSC", bad_genes, reverse=F)
clean_mat = scdb_mat("MSC_clean")
## Remove genes with mito fraction > 0.1
mcell_mat_ignore_cells("MSC_clean", "MSC_clean", 
                       union(clean_mat@ignore_cells, names(uc)[mito_f >= 0.1 ]))
clean_mat = scdb_mat("MSC_clean")
## Keep only Singlets
mcell_mat_ignore_cells("MSC_clean", "MSC_clean", 
                       union(clean_mat@ignore_cells, names(uc)[clean_mat@cell_metadata$HTO_classification.global != "Singlet" ]))
clean_mat = scdb_mat("MSC_clean")
## Remove cells with more than 3000 genes (probable doublets)
mcell_mat_ignore_cells("MSC_clean", "MSC_clean", 
                       union(clean_mat@ignore_cells, names(uc)[clean_mat@cell_metadata$nFeature_RNA >= 3000 ]))
clean_mat = scdb_mat("MSC_clean")


## Cleanup
remove(mat)
remove(mito_f)
remove(uc)

### Make gstat object
mcell_add_gene_stat(mat_id = "MSC_clean", gstat_id = "MSC_gs", force = T)
gstat = scdb_gstat("MSC_gs")
dim(gstat)

# calculate feature genes
mcell_gset_filter_varmean(gstat_id = "MSC_gs", gset_id = "MSC_featsBIG", T_vm=0.025, force_new=T)
mcell_gset_filter_cov(gstat_id = "MSC_gs", gset_id = "MSC_featsBIG", T_tot=80, T_top3=2)

feats_gset <- scdb_gset("MSC_featsBIG")
mcell_plot_gstats(gstat_id = "MSC_gs", "MSC_featsBIG")

# filter feats gset for lateral genes
lateral_gset_id = "lateral"
lateral_gset = scdb_gset("lateral")
feats_gset = gset_new_restrict_gset(feats_gset, lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
scdb_add_gset(id = "MSC_feats_filt", gset = feats_gset)

## MetaCell generation
mcell_add_cgraph_from_mat_bknn(mat_id = "MSC_clean", gset_id = "MSC_feats_filt", graph_id = "MSC_graph_big2", K=150, dsamp=T)
mcell_coclust_from_graph_resamp(coc_id = "MSC_coc400big2", graph_id = "MSC_graph_big2", min_mc_size=100, p_resamp=0.75, n_resamp=400)
mcell_mc_from_coclust_balanced(mc_id = "MC_big2", coc_id =  "MSC_coc400big2", "MSC_clean", K=60, min_mc_size=120, alpha=2)

## Reload MC object, add colors and plot 2d projection
mc = scdb_mc("MC_big2")
mc@colors <- viridis(length(mc@annots))
scdb_add_mc("MSC_mc_big2_col",mc)
mc = scdb_mc("MSC_mc_big2_col")
mcell_mc2d_force_knn(mc2d_id = "MSC_mc_big2_col" ,mc_id =  "MSC_mc_big2_col", "MSC_graph", ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",800, "metacell")
tgconfig::set_param("mcell_mc2d_width",800, "metacell")
mcell_mc2d_plot("MSC_mc_big2_col")

## get log2enrichment values for all genes for each of the MetaCells
lfp <- log2(mc@mc_fp)

## Save lfp table
saveRDS(as.data.frame(lfp), "./misc_data/lfp_table2.rds")
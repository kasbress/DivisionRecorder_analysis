library(metacell)
library(tidyverse)
library(viridis)

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

library(metacell)

# Used this script to get a broad selection of genes that have been reported to have important functions in lymphocyte biology
# As I just grep everything based on gene symbol, we'll also get family members of genes for which function is less well known


# Set working directory
setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/Analysis_memory/")

## import data
scdb_init("./metacell_db/", force_reinit=T) 
mc <- scdb_mc("DivRecMEM_MC")
lfp <- log2(mc@mc_fp)



T_cell_marks <- grep(pattern = "Sell$|Ltb|Tox|Ccr\\d+$|Cxcr\\d+$|Cd\\d+|Xcl1|Ifng|Prf1|Gzm|^Il\\d+|^Il\\d+r|Tnf|^Jak\\d+$|Itg|Ly\\d+|Tlr\\d+|Klr", 
                     x = row.names(lfp), value = T)
more <- grep(pattern = "Havcr2|Lag3|Cx3cr1|Pdcd1|^Nr4a|Ctla|Havcr2|Serpin|P2rx|Traf\\d+$|Gimap|Bcl\\d+|Lgals\\d+$|Bcl2|Cd\\d+$|Fbxo\\d+$|Lat$|Kdelr|Prr7|Tgf|Fas$|Fasl$|S100a|Clec", 
             x = row.names(lfp), value = T)
TF <- grep(pattern = "Id\\d+|Tcf\\d+|Foxo\\d+|Dnmt1|Ezh2|Bach2|Myb|Zeb\\d+|Prdm\\d+$|Klf\\d+|Sox\\d+$|Cd\\d+$|Foxp\\d+|Gata\\d+|Socs\\d+|^Stat\\d+$|Eomes|Tbx21|Icos|Jun|Runx", 
           x = row.names(lfp), value = T)
add <- grep(pattern = "S1pr|^Jak|^Tlr", 
            x = row.names(lfp), value = T)
imm_genes <- unique(c(T_cell_marks, more, TF, add))

write_rds(imm_genes, "./misc_data/immune_genes.rds")

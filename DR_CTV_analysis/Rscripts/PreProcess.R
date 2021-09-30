library(metacell)
library(Seurat)
library(viridis)
library(tidyverse)

setwd("/DATA/users/k.bresser/DivisionRecorder_scRNAseq_full/")

# For output from CellRanger >= 3.0 with multiple data types
data_dir <- "./CTV_experiment/"
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
rawdata <- Read10X(data.dir = data_dir)

seurat_object <- CreateSeuratObject(counts = rawdata$`Gene Expression` )
seurat_object[['HTO']] = CreateAssayObject(counts = rawdata$`Custom`)
seurat_object[['ADT']] = CreateAssayObject(counts = rawdata$`Antibody Capture`)

# Assign Hashtags
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
seurat_object <- NormalizeData(seurat_object, assay = "HTO", normalization.method = "CLR")
seurat_object <- HTODemux(seurat_object, assay = "HTO", positive.quantile = .99)



# Global classification results
table(seurat_object$HTO_classification.global)
table(seurat_object$hash.ID)

# Group cells based on the max HTO signal
Idents(seurat_object) <- "HTO_maxID"
RidgePlot(seurat_object, assay = "HTO", features = rownames(seurat_object[["HTO"]]), ncol = 2)
ggsave(filename = "./Analysis_CTV_experiment/seurat_plots/Ridgeplots.pdf", device = "pdf", width = 6, height = 5)

Idents(seurat_object) <- "HTO_classification.global"
FeatureScatter(seurat_object, feature1 = "MSC15-1", feature2 = "MSC15-2")
# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(seurat_object, assay = "HTO", ncells = 5000)
ggsave(filename = "./Analysis_CTV_experiment/seurat_plots/Heatmap_HTO.png", device = "png", width = 4, height = 3)

Idents(seurat_object) <- "HTO_classification.global"
VlnPlot(seurat_object, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
ggsave(filename = "./Analysis_CTV_experiment/seurat_plots/RNAcountClassification.png", device = "png", width = 5, height = 3)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = seurat_object, assay = "HTO"))))
# Calculate tSNE embeddings with a distance matrix ||| Will take long
seurat_object <- RunTSNE(seurat_object, distance.matrix = hto.dist.mtx, perplexity = 100)

DimPlot(seurat_object)
ggsave(filename = "./Analysis_CTV_experiment/seurat_plots/HTO_tSNE.png", device = "png", width = 6, height = 6)
DimPlot(seurat_object,group.by = "hash.ID")
ggsave(filename = "./Analysis_CTV_experiment/seurat_plots/HTO_tSNE_hashID.png", device = "png", width = 6, height = 6)
DimPlot(seurat_object,group.by = "HTO_maxID")
ggsave(filename = "./Analysis_CTV_experiment/seurat_plots/HTO_tSNE_HTO_maxID.png", device = "png", width = 6, height = 6)

rm(hto.dist.mtx)

# Save Seurat object
saveRDS(object = seurat_object, file = "./Analysis_CTV_experiment/misc_data/seurat_object_CTVexp.rds")


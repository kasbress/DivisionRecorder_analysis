library(Seurat)
library(tidyverse)

# Set working directory
setwd("/PATH/TO/DATA/FOLDER/")

# First we import the data from the 2nd experiment (harvest day 1)
data_dir <- "./memory_day1/"
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
rawdata <- Read10X(data.dir = data_dir)

# Create the Seurat object and add HTO and ADT data
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

# Save as separate datafile
write_rds(seurat_object, file = "./Analysis_memory/misc_data/seurat_object_day1.rds")


##### Now we do the same for the data from the 2nd experiment (harvest day 2)
data_dir <- "./memory_day2/"
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
rawdata <- Read10X(data.dir = data_dir)

# Create the Seurat object and add HTO and ADT data
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

# Save as separate datafile
write_rds(seurat_object, file = "./Analysis_memory/misc_data/seurat_object_day2.rds")


################ Fuse the above Seurat objects with the Experiment 1 Seurat object
setwd("./Analysis_memory/")

# Import experiment 2
seurat_day1 <- read_rds("./misc_data/seurat_object_day1.rds")
seurat_day2 <-read_rds("./misc_data/seurat_object_day2.rds")

# Merge the objects
seurat.memory <- merge(seurat_day1, y = seurat_day2, project = "seurat.memory")

# Recode the sample hashtags for easier reference
seurat.memory@meta.data$hash.ID %>% 
  as.factor() %>%
  fct_recode(DR4.GFP = "MSC15-1", 
             DR5.GFP = "MSC15-2",
             DR45.RFP = "MSC15-3",
             DR6.GFP = "MTO7", 
             DR7.GFP = "MTO8",
             DR67.RFP = "MTO3",
             Naive1 = "MSC15-4", 
             Naive2 = "MSC15-5",
             Naive3 = "MTO4", 
             Naive4 = "MTO5") -> seurat.memory@meta.data$hash.ID 

# Check
table(seurat.memory@meta.data$hash.ID)

# Add experiment IDs
seurat.memory@meta.data$experiment <- ifelse(seurat.memory@meta.data$hash.ID %in% c("DR4.GFP", "DR5.GFP", "DR45.GFP", "Naive1", "Naive2"),
                                             "Mem.exp.2.1", "Mem.exp.2.2")

# Load Seurat object containing data from experiment 1
seurat_old <- read_rds(file = "./misc_data/seurat_object_all_memory.rds")

# Make sure to get all the right cells
seurat_old@meta.data %>% 
  rownames_to_column("cellcode") %>% 
  filter(experiment == "Mem.exp.1") %>% 
  pull(cellcode) -> exp.1.cells
seurat_old <- subset(seurat_old, cells = exp.1.cells)

write_rds(seurat_old, "./misc_data/seurat_object_Exp1.rds")

##Combine seurat objects and save
seurat.memory <- merge(seurat.memory, y = seurat_old, project = "seurat.memory")


# Recode the sample hashtags for easier reference
seurat.memory@meta.data$hash.ID %>% 
  as.factor() %>%
  fct_recode(M1.GFP= "DR1.GFP", 
             M1.RFP = "DR1.RFP",
             M2.GFP = "DR2.GFP",
             M2.RFP = "DR2.RFP", 
             M3.GFP = "DR3.GFP",
             M3.RFP = "DR3.RFP",
             M4.GFP = "DR4.GFP", 
             `M4+5.RFP` = "DR45.RFP",
             M5.GFP = "DR5.GFP", 
             M6.GFP = "DR6.GFP",
             M7.GFP = "DR7.GFP", 
             `M6+7.RFP` = "DR67.RFP",
             M8 = "Naive1",
             M9 = "Naive2",
             M10 = "Naive3", 
             M11 = "Naive4") -> seurat.memory@meta.data$mouseID 

seurat.memory@meta.data %>% 
  rownames_to_column("cellcode") %>% 
  as_tibble() %>% 
  select(cellcode, hash.ID, mouseID, experiment, contains("HTO"), contains("_RNA")) -> metadata

write_tsv(metadata, "./misc_data/Metadata_table.tsv")

write_rds(seurat.memory, file = "./misc_data/seurat_object_all_exps.rds")


############################### To generate metadata table
library(metacell)

scdb_init("PATH/TO/MCdb")

mc <- scdb_mc("DivRecMEM_MC")

mc@mc %>% 
  enframe("cellcode", "MetaCell") %>% 
  mutate(MetaCell = factor(MetaCell)) %>% 
  full_join(metadata) %>% 
  mutate(`Tm cluster` = fct_collapse(MetaCell,
                                `Tcm(eff)` = c('23','15','22','8','6','10','4','5','7','18'),
                                `Tcm(mult)` = c('9','1','12','13','14','16','11','2','3'),
                                Tem = c('21','19','17','20')),
         phenotype = fct_collapse(MetaCell,
                                  Tcm = c('23','15','22','8','6','10','4','5','7','18','9','1','12','13','14','16','11','2','3'),
                                  Tem = c('21','19','17','20')), 
         origin = fct_collapse(hash.ID,
                               unmodified = c('Naive1', 'Naive2', 'Naive3', 'Naive4'),
                               DR.modified = c('DR4.GFP','DR5.GFP','DR6.GFP','DR7.GFP'))  )  -> metadata.full


write_tsv(metadata.full, "./misc_data/Metadata_table_DivRec.tsv")

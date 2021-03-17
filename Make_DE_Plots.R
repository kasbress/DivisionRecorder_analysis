## Perform DE analysis between MetaCells


library(metacell)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(cowplot)



## SET WD AND INITIATE MC_DB
scdb_init("./metacell_db/", force_reinit=T) 

## LOAD MC  OBJECT
mc = scdb_mc("MSC_mc_big2_col")

## GENERATE DATAFRAME WITH GROUP ANNOTATIONS (ROWNAMES ARE CELLCODES)
MCs <- mc@mc
MCs <- data.frame( ident = as.character(MCs), 
                     row.names = names(MCs) )

## CLEANUP AND LOAD SEURAT OBJECT
remove(mc)
load("./Seuratdata.Rda")


## ADD METADATA
countsfil <- AddMetaData(object = countsfil, metadata = MCs, col.name = "MC")

## NORMALIZE UMI-COUNTS, AND RE-SCALE
countsfil <- NormalizeData(countsfil)
all.genes <- rownames(countsfil)
countsfil <- ScaleData(countsfil, features = all.genes)

## USE THE IDENTS FUNCTION TO SET THE METADATA COLUMN TO USE IN COMPARISONS
Idents(countsfil) <- "MC"


## Find markers with possion test. Ran the piece below for 13vs20, 13vs19 and 13vs17.
markers <- FindMarkers(object = countsfil, ident.1 = "13", ident.2 = "20",
                       logfc.threshold = 0, test.use = "poisson")


### GET AVERAGE UMI COUNTS
umi_data <- as.data.frame(countsfil@assays$RNA@data)
cells <- rownames(subset(MCs, ident == "13"))
LO_UMI <- umi_data[,cells]
LO_UMI$Average <- rowMeans(LO_UMI)
LO_UMI <- data.frame(Average = LO_UMI$Average, row.names = row.names(LO_UMI))

### PREPARE FOR GGPLOT
to_plot <- merge(x = LO_UMI, y = markers, by = 0, all.x = F, all.y = T)
to_plot$sig <- ifelse( (to_plot$p_val < 0.05 & to_plot$avg_logFC > 0.2) | (to_plot$p_val < 0.05 & to_plot$avg_logFC < -0.2), TRUE, FALSE)

### PLOT GRAPH
ggplot(to_plot, aes( y= avg_logFC , x=Average, color = sig, label = Row.names) ) + 
  geom_point()+
  scale_color_manual(values=c("black", "red"))+
  geom_text_repel(data = subset(to_plot, sig == TRUE ),box.padding = 1 )+
  labs(title = "MC13 v MC20")+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  ylim(-1,1)

ggsave(filename = "./MC13vMC20.pdf",device = "pdf", width = 6, height = 4,useDingbats=FALSE )

write.table(x = markers, file = "./misc_data/MC13vMC20.txt", sep = "\t")



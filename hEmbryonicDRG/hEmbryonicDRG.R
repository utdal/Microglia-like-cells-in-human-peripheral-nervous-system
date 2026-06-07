
library(Seurat)
library(ggplot2)
library(data.table)
library(tidyverse)

##### Load data ####

#Download counts
counts_fpath <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE245310&format=file&file=GSE245310%5Fhuman%5FDRG%5Fgene%5Fexpr%2Etxt%2Egz'
counts_tmp <- tempfile(fileext = ".txt.gz")
options(timeout = 600)
download.file(counts_fpath, counts_tmp, mode = "wb", method = "curl",
              extra = "-L --retry 3")
counts <- fread(counts_tmp, header = TRUE, sep = "\t") |>
  as.data.frame() |>
  (\(df) { rownames(df) <- df[[1]]; df[, -1] })()

#Download and metadata
meta_fpath <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE245310&format=file&file=GSE245310%5Fhuman%5FDRG%5Fmeta%2Etxt%2Egz'
meta_tmp <- tempfile(fileext = ".txt.gz")
download.file(meta_fpath, meta_tmp, mode = "wb")
metadata <- read.table(gzfile(meta_tmp), header = TRUE, sep = "\t", row.names = 1)
metadata <- column_to_rownames( remove_rownames(metadata), var = 'cell')
colnames(counts) <- rownames(metadata)

#Make srat obj
hEmbDRG <- CreateSeuratObject( counts, meta.data = metadata)
hEmbDRG <- NormalizeData(hEmbDRG)
table(hEmbDRG$celltype)
rm(counts, metadata)
gc()

#take a look at Mac genes
DotPlot(hEmbDRG, group.by = 'celltype', features = c('MRC1', 'CSF1R', 'CD163', 'ITGAX', 'CD68', 'AIF1', 'P2RY12', 'CX3CR1' ))



##### Get macrophages/MLC and cluster ####

#get immune cells
Idents(hEmbDRG) <- 'celltype'
macs <- subset( hEmbDRG, idents = 'macrophage' )
macs[["RNA"]] <- split(macs[["RNA"]], f = macs$sampleID)
macs <- NormalizeData(macs)
macs <- FindVariableFeatures( macs )
macs <- ScaleData( macs )
macs <- RunPCA( macs )
macs <- IntegrateLayers(macs, method = HarmonyIntegration,  
                        orig.reduction = 'pca', new.reduction = 'harmony')
macs <- RunUMAP(macs, reduction = 'harmony', 
                    dims = 1:15, reduction.name = 'umap.macs', return.model = TRUE)
macs <- FindNeighbors(macs, reduction = 'harmony', dims = 1:15)
macs <- FindClusters(macs, cluster.name = 'subclusters', resolution = 1)

#Take a look
DimPlot( macs, reduction = 'umap.macs', group.by = 'subclusters', size = 3, label = T, repel = T)
DimPlot( macs, reduction = 'umap.macs', group.by = 'sampleID', size = 3, label = T, repel = T)
FeaturePlot( macs, reduction = 'umap.macs', label = T, repel = T, 
         features = 'MRC1' )
FeaturePlot( macs, reduction = 'umap.macs', label = T, repel = T, 
             features = 'P2RY12' )

#Save macs object
saveRDS(macs, 'macs.rds')



#### DOT PLOT — MAC SUBTYPES ####

genes <- c( 'CSF1R', 'AIF1', 'CD68', 'ITGAM', "CLEC10A", "CD1D", "FCER1A", "CLEC9A", "XCR1",
            "CD300E", "LYZ", "S100A9", "C1QA", "SERPINA1",
            "MARCO", "GPNMB", "CD209", "LILRB5", "MRC1", "SIGLEC1",
            "SPP1", 'ADGRG1', 'PADI2', 'ABCC4', 'WNT5A', 'P2RY12', 
            'PLAC8','TMEM119', 'CD83', 'CCL3', 'NR4A1', 
            'IL1B', 'NFKB2', 'ICAM1', 'TNFAIP3',
            'CENPM', 'MKI67', 'FANCA', 'DNAJB1', 'HSPA1A', 'HSPA1B'
)
DotPlot(macs, scale = F,
                    features = genes[genes %in% rownames(macs)],
                    group.by = "subclusters") +
  scale_color_gradientn(
    colours = c("#FFF9C4", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C", "#72190E")
  ) +
  #scale_y_discrete(limits = label_order) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_text(hjust = 0))
ggsave('myeloid_markers.svg', plot = plot_mac, width = 13, height = 4)



#### LABEL CLUSTERS AND REMOVE CONTAMINATION ####

cluster_labels <- c(
  '0' = 'MLC(Mature)',
  '1' = 'MNP5',
  '2' = 'MNP2',
  '3' = 'MNP3',
  '4' = 'MNP8',
  '5' = 'MNP4',
  '6' = 'MLC(Pre)',
  '7' = 'MNP7',
  '8' = 'MNP6',
  '9' = 'cDC',
  '10' = 'MNP1'
)

macs$Name <- unname(cluster_labels[as.character(macs$subclusters)])

#UMAP
label_order <- c(
  'MLC(Pre)',
  'MLC(Mature)',
  'MNP1',
  'MNP2',
  'MNP3',
  'MNP4',
  'MNP5',
  'MNP6',
  'cDC',
  'MNP7',
  'MNP8'
)
cbfColors <- c(
  'MLC(Pre)'    = "#882E72",
  'MLC(Mature)' = "#1965B0",
  'MNP1' = "#BBBBBB",
  'MNP2' = "#72190E",
  'MNP3' = "#DC050C",
  'MNP4' = "#E8601C",
  'MNP5' = "#F1932D",
  'MNP6' = "#F6C141",
  'cDC'  = "#F7EE55",
  'MNP7' = "#CAE0AB",
  'MNP8' = "#4EB265"
) 
macs$Name <- factor(macs$Name, levels = label_order)
DimPlot(macs, reduction = 'umap.macs', group.by = 'Name',
        cols = cbfColors, size = 3, alpha = 0.5) +
  labs(title = 'Mononuclear Phagocytes', x = 'UMAP1', y = 'UMAP2')
ggsave( 'mnp_umap.svg', width = 6, height = 4)

#save
saveRDS( macs, 'macs.rds')


#### FINAL DOT PLOT ####

genes <- c( 'CSF1R', 'ITGAM', "CLEC10A", "CLEC9A",
             "LILRB5", "MRC1", "SIGLEC1",
            "SPP1", 'ADGRG1', 'PADI2', 'ABCC4', 'WNT5A', 'P2RY12', 
            'CX3CR1','TMEM119', 'CD83', 'CCL3', 
            'MKI67', 'FANCA'
)

plot_mac <- DotPlot(macs, scale = F,
                    features = genes[genes %in% rownames(macs)],
                    group.by = "Name") +
  scale_color_gradientn(
    colours = c("#FFF9C4", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C", "#72190E")
  ) +
  #scale_y_discrete(limits = label_order) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL) +
  coord_flip() +
  theme(axis.text.y = element_text(hjust = 0))

print(plot_mac)
ggsave('myeloid_markers_final.svg', plot = plot_mac, width = 6, height = 5)



#### CELL COMPOSITION PER CLUSTER ####

cbfColors2 <- c("#7BAFDE", "#1965B0", "#4EB265", "#CAE0AB", 
                "#F7EE55", "#F6C141", "#E8601C", "#DC050C", "#72190E")

time_order <- c('GW7', 'GW8', 'GW9', 'GW10', 'GW12', 'GW14', 'GW15', 'GW17', 'GW21')
comp_df <- as.data.frame(table(Name = macs$Name, Time = macs$time))
comp_df$Time <- factor(comp_df$Time, levels = time_order)

ggplot(comp_df, aes(x = Time, y = Freq, fill = Name)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = cbfColors) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = 'Proportion of cells', fill = 'Time') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('composition_by_time.svg', width = 4, height = 3)



#### MACS COUNTS TABLE ####

library(openxlsx)

t <- table(macs$time, macs$Name)
t_df <- as.data.frame.matrix(t)
t_df$Total <- rowSums(t_df)
t_df <- rbind(t_df, Total = colSums(t_df))

write.xlsx(t_df, "macs_table.xlsx", rowNames = TRUE)

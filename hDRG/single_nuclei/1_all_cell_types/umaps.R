#Libraries
library(SingleCellExperiment)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(openxlsx)

#Color palette
cbfColors <- c(
  "#882E72",
  "#1965B0",
  "#5289C7",
  "#7BAFDE",
  "#4EB265",
  "#90C987",
  "#CAE0AB",
  "#F7EE55",
  "#F6C141",
  "#F1932D",
  "#E8601C",
  "#DC050C",
  "#72190E",
  "#BBBBBB"
)

#Check
table(hDRG$HighName)

#Neurons
Idents(hDRG) <- 'HighName'
neurons <- subset( hDRG, idents = 'Neurons' )
NeuronsFinal.umap <- readRDS("C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/hDRG_snFlex/pipeline/v2.6/output/6_Subcluster/Neurons/round7/NeuronsFinal.umap.rds")
neuron_cells <- Cells(neurons)
umap_coords <- Embeddings(NeuronsFinal.umap)[rownames(Embeddings(NeuronsFinal.umap)) %in% neuron_cells, ]
neurons[["umap"]] <- CreateDimReducObject(
  embeddings = umap_coords,
  key = "UMAP_",
  assay = DefaultAssay(neurons)
)
DimPlot( neurons, reduction = 'umap', group.by = 'Name_label', 
         alpha = 0.1, raster = F,
         cols = cbfColors[c(1,3,4,5,7,8,9,11,13)] ) +
  labs(title = 'hDRG Neurons', x = 'UMAP1', y = 'UMAP2')
ggsave( 'neu_umap.svg', width = 6, height = 5)
ggsave( 'neu_umap.png', width = 6, height = 5, dpi = 600)
neurons <- NormalizeData(neurons)
FeaturePlot( neurons, reduction = 'umap', features = c('MRGPRE', "GRN", "PDIA2", "AGRN"), 
         alpha = 0.1, raster = F) +
  labs(x = 'UMAP1', y = 'UMAP2')
ggsave( 'neu_genes.svg', width = 12, height = 10)
ggsave( 'neu_genes.png', width = 12, height = 10, dpi = 600)
neurons[['RNA']]$data <- NULL
saveRDS(neurons, 'neurons.rds')
saveRDS(neurons[['umap']], 'neurons.umap.rds')
rm(neurons, hDRG, NeuronsFinal.umap, umap_coords)
gc()

#SGC
sgc_cells <- rownames( hDRG@meta.data )[ hDRG@meta.data$HighName == 'SatGCs']
satgc <- subset( hDRG, cells = sgc_cells)
satgc[['RNA']] <- split( satgc[['RNA']], f = satgc$sample)
satgc <- NormalizeData(satgc)
sgc_hvgs <- read.table("C:/Users/NoorT/Box/gataca/personal_analysis/Khadijah/my_projects/hDRG_snFlex/pipeline/v2.6/output/6_Subcluster/SatGC/process_files/hvgs.txt")
sgc_hvgs <- sgc_hvgs$x
satgc <- ScaleData(satgc, features = sgc_hvgs)
satgc <- RunPCA(satgc, features = sgc_hvgs, npcs = 100)
satgc <- IntegrateLayers(satgc, method = HarmonyIntegration, features = sgc_hvgs,
                         orig.reduction = 'pca', new.reduction = 'harmony')
satgc <- RunUMAP(satgc, reduction = 'harmony', dims = 1:5  , return.model = TRUE)
satgc$Name <- ifelse(satgc$Name=='SGC.GABRB1', 'SGC.PILRB', satgc$Name)
DimPlot( satgc, reduction = 'umap', group.by = 'Name', 
         alpha = 0.1, raster = F,
         cols = cbfColors[c(5,13,2,12)] ) +
  labs(title = 'hDRG SGCs', x = 'UMAP1', y = 'UMAP2')
ggsave( 'sgc_umap.svg', width = 6, height = 5)
ggsave( 'sgc_umap.png', width = 6, height = 5, dpi = 600)
saveRDS(satgc[['umap']], 'sgc.umap.rds')
rm(satgc, satgc.meta)
gc()



#NMSC
nmsc_cells <- rownames( hDRG_meta )[ hDRG_meta$HighName == 'NMSchwann']
nmschwann <- readRDS("C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/hDRG_snFlex/pipeline/v2.6/output/6_Subcluster/NMSchwann/nmschwann.rds")
nmschwann <- subset( nmschwann, cells = nmsc_cells)
Idents(nmschwann) <- 'Name'
nmschwann <- subset( nmschwann, idents = c("NMSC.NRXN1.COL9A3"), invert = T )
nmschwann$Name <- dplyr::recode(nmschwann$Name,
                                 'NMSC.NRXN1.VWA1'           = 'NMSC.NRXN1',
                                 'NMSC.NRXN1.VWA1.SERPINE2'  = 'NMSC.NRXN1',
                                 'NMSC.NRXN1.MDGA1'          = 'NMSC.NRXN1',
                                 'NMSC.NRXN1.SPOCD1'         = 'NMSC.NRXN1',
                                 'NMSC.PXDN.ADGRV1'          = 'NMSC.PXDN' ,
                                 'NMSC.PXDN.NTRK3'           = 'NMSC.PXDN' )
nmschwann <- RunUMAP( nmschwann, reduction = 'harmony', reduction.name = 'umap2', dims = 1:5 )
DimPlot(nmschwann, reduction = 'umap2', group.by = 'Name', alpha = 0.3)
FeaturePlot( nmschwann, reduction = 'umap2', features = c('nCount_RNA'), max.cutoff = 2000, alpha = 0.3, order = T)
DimPlot( nmschwann, reduction = 'umap2', group.by = 'Name', 
         alpha = 0.1, raster = F,
         cols = cbfColors[c(2,12,5,13)] ) +
  labs(title = 'hDRG NMSchwann', x = 'UMAP1', y = 'UMAP2')
ggsave( 'nmsc_umap.svg', width = 6, height = 5)
ggsave( 'nmsc_umap.png', width = 6, height = 5, dpi = 600)
saveRDS(nmschwann[['umap2']], 'nmschwann.umap.rds')
rm(nmschwann)
gc()



#MSC
Idents(hDRG) <- 'HighName'
msc <- subset( hDRG, idents = 'MSchwann')
msc_cells <- colnames(msc)
mschwann <- readRDS("C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/hDRG_snFlex/pipeline/v2.6/output/6_Subcluster/MSchwann/mschwann.rds")
harmony.reduc <- mschwann[['harmony']]
harmony.reduc <- Embeddings(harmony.reduc)[msc_cells, ]
msc[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony.reduc,
  key = "harmony_",
  assay = DefaultAssay(msc)
)
msc <- RunUMAP( msc, reduction = 'harmony', reduction.name = 'umap2', dims = 1:30 )
DimPlot(msc, reduction = 'umap2', group.by = 'Name', alpha = 0.3)
msc <- NormalizeData(msc)
FeaturePlot( msc, reduction = 'umap2', features = c('BZW2'), order = T)
msc$Name <- dplyr::recode(msc$Name,
                                'MSC.BZW2.EMID1'  = 'MSC.BZW2',
                                'MSC.BZW2.LRP1B'  = 'MSC.BZW2',
                                'MSC.MAL'         = 'MSC.IL17B' )
DimPlot( msc, reduction = 'umap2', group.by = 'Name', 
         alpha = 0.3, raster = F,
         cols = cbfColors[c(1,3,5,9,12)] ) +
  labs(title = 'hDRG MSchwann', x = 'UMAP1', y = 'UMAP2')
ggsave( 'msc_umap.svg', width = 6, height = 5)
ggsave( 'msc_umap.png', width = 6, height = 5, dpi = 600)
saveRDS(msc[['umap2']], 'mschwann.umap.rds')
rm(mschwann, msc, harmony.reduc)
gc()



#Granulocytes
Idents(hDRG) <- 'HighName'
granulo <- subset( hDRG, idents = 'Granulocytes' )
granulo_cells <- colnames(granulo)
harmony.full <- readRDS("C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/hDRG_snFlex/pipeline/v2.6/output/5_SketchIntegration/harmony.full.rds")
harmony.sub <- Embeddings(harmony.full)[granulo_cells, ]
granulo[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony.sub,
  key = "harmony_",
  assay = DefaultAssay(granulo)
)
granulo <- RunUMAP( granulo, reduction = 'harmony', reduction.name = 'umap2', dims = 1:50 )
DimPlot(granulo, reduction = 'umap2', group.by = 'Name', alpha = 0.3)
granulo <- NormalizeData(granulo)
FeaturePlot( granulo, reduction = 'umap2', features = c('HDC'), order = T)
DimPlot( granulo, reduction = 'umap2', group.by = 'Name', 
         alpha = 0.3, raster = F,
         cols = cbfColors[c(2,12)] ) +
  labs(title = 'hDRG Granulocytes', x = 'UMAP1', y = 'UMAP2')
ggsave( 'granulo_umap.svg', width = 6, height = 5)
ggsave( 'granulo_umap.png', width = 6, height = 5, dpi = 600)
saveRDS(granulo[['umap2']], 'granulo.umap.rds')
rm(granulo, harmony.sub)
gc()



#Lymphocytes
Idents(hDRG) <- 'HighName'
lympho <- subset( hDRG, idents = 'Lymphocytes' )
lympho_cells <- colnames(lympho)
immune <- readRDS("C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/hDRG_snFlex/pipeline/v2.6/output/6_Subcluster/Immune/immune.rds")
harmony <- immune[['harmony']]
harmony <- Embeddings(harmony)[lympho_cells, ]
lympho[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony,
  key = "harmony_",
  assay = DefaultAssay(lympho)
)
lympho <- RunUMAP( lympho, reduction = 'harmony', reduction.name = 'umap2', dims = 1:100 )
DimPlot(lympho, reduction = 'umap2', group.by = 'Name', alpha = 0.3)
lympho <- NormalizeData(lympho)
FeaturePlot( lympho, reduction = 'umap2', features = c('KLIRG'), order = T)
lympho$Name <- dplyr::recode(lympho$Name,
                                 'Bcell.IGHG'    = 'Bcell',
                                 'Bcell.IGHM'    = 'Bcell',
                                 'Bcell.NR4A2'   = 'Bcell',
                                 'Plasma.IGHG.CADPS2' = 'PlasmaCells',
                                 'Plasma.IGHG.CIITA'  = 'PlasmaCells',
                                 'Plasma.IGHG.ITM2C'  = 'PlasmaCells',
                                 'Plasma.IGHG.SEC24A' = 'PlasmaCells',
                                 'Plasma.IGHM.IGHG1'  = 'PlasmaCells',
                                 'Plasma.IGHM.PLEKHA2'= 'PlasmaCells',
                                 'Tcell.LTK.1'  = 'Tcell.LTK',
                                 'Tcell.LTK.2'  = 'Tcell.LTK',
                                 'Lymph.Prolif'  = 'Lymph.Mitotic')
DimPlot( lympho, reduction = 'umap2', group.by = 'Name', 
         alpha = 0.3, raster = F,
         cols = cbfColors[c(1,3,5,8,10,12,13)] ) +
  labs(title = 'hDRG Lymphocytes', x = 'UMAP1', y = 'UMAP2')
ggsave( 'lympho_umap.svg', width = 6, height = 5)
ggsave( 'lympho_umap.png', width = 6, height = 5, dpi = 600)
saveRDS(lympho[['umap2']], 'lympho.umap.rds')
rm(lympho, harmony, immune, immune.cleaned.umap)
gc()



#Vascular
vasc_cells <- colnames( hDRG )[ hDRG$HighName %in% c('Endothelial', 'Mural')]
vascular <- readRDS("C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/hDRG_snFlex/pipeline/v2.6/output/6_Subcluster/Vascular/vascular.rds")
vascular <- subset( vascular, cells = vasc_cells)
vascular <- RunUMAP( vascular, reduction = 'harmony', reduction.name = 'umap2', dims = 1:50 )
DimPlot( vascular, reduction = 'umap2', group.by = 'Name', 
         alpha = 0.1, raster = F,
         cols = cbfColors[c(1,2,4,5,7:13)] ) +
  labs(title = 'hDRG Vascular Cells', x = 'UMAP1', y = 'UMAP2')
ggsave( 'vasc_umap.svg', width = 6, height = 5)
ggsave( 'vasc_umap.png', width = 6, height = 5, dpi = 600)
saveRDS(vascular[['umap2']], 'vascular.umap.rds')
rm(vascular)
gc()



#Fibroblasts
fibro_cells <- colnames( hDRG )[ hDRG$HighName == 'Fibroblasts']
fibroblasts <- readRDS("C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/hDRG_snFlex/pipeline/v2.6/output/6_Subcluster/Fibroblasts/fibroblasts.rds")
fibroblasts <- subset( fibroblasts, cells = fibro_cells)
fibroblasts <- RunUMAP( fibroblasts, reduction = 'harmony', reduction.name = 'umap2', dims = 1:50 )
fibroblasts$Name <- dplyr::recode(fibroblasts$Name,
                               'Fibro.Endo.PI16.H'    = 'Fibro.Endo.PI16',
                               'Fibro.Endo.PI16.L'    = 'Fibro.Endo.PI16')
DimPlot( fibroblasts, reduction = 'umap2', group.by = 'Name', 
         alpha = 0.1, raster = F,
         cols = cbfColors[c(1,3,5,8,10,12,13)] ) +
  labs(title = 'hDRG Fibroblasts', x = 'UMAP1', y = 'UMAP2')
ggsave( 'fibro_umap.svg', width = 6, height = 5)
ggsave( 'fibro_umap.png', width = 6, height = 5, dpi = 600)
saveRDS(fibroblasts[['umap2']], 'fibroblasts.umap.rds')
rm(fibroblasts)
gc()









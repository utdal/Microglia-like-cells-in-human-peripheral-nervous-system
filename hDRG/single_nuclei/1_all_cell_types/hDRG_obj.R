
#### Setup ####

#Libraries
library(SingleCellExperiment)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(tidyverse)

#Working dir
setwd('C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/microglia_like_cells/analysis/single_nuc_genes')

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



#### Prep all cell types ####

#Load full data
hDRG <- readRDS("C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/hDRG_snFlex/pipeline/v2.6/output/6_Subcluster/Whole/post_filt/merged.object.rds")
finalumap <- readRDS("C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/hDRG_snFlex/pipeline/v2.6/output/6_Subcluster/Whole/post_filt/finalumap.rds")

#Published samples only = 92 samples
samples <- c(
  "AJHA487", "AJHF199", "AJFE007", "AJL5162A", "AJL5162B",
  "AKIE115", "AJIK138", "AJHX078A", "AJGJ324", "AJCR033",
  "AJHC333", "AJCE093", "AJE4179", "AJIG066", "AJKN140",
  "AJLL013", "146T9R", "155T8L", "153T6L", "156T8L",
  "9AL", "16AL", "19", "21", "22",
  "24", "25", "28", "27AR", "32AR",
  "43AR", "45AL", "153T6R", "150T6L", "159T2R",
  "161T11L", "AIAG131", "AJFA324", "AJH5198", "AJHX078B",
  "AJHX209", "AKAA358A", "AKHO306", "AKK3142", "ALEA056",
  "AJGS325", "AJHR406", "AJI1393", "AJJA309", "AKCG403",
  "AKHZ124", "AKJC406", "AKK3408", "AIIK226", "AJK1119",
  "AKHD241", "AKIQ317", "AKIY142", "ALAQ111", "ALAW101",
  "ALDR082", "AIEB188", "AIEJ414", "AIEQ242", "AJKT138",
  "AKA1235", "AKAA358B", "AKGU264", "ALAD176", "ALBT459",
  "36AL", "38AR", "40AR", "50AL", "168T5R",
  "169T3R", "171T7L", "171T7R", "AKBF175", "AKBE204",
  "AJJB371", "AJEX220", "AJLA352", "AJJJ231", "AJD2236",
  "AJGB347", "AJIW220", "ALE3489", "AMDD111", "AKK1131",
  "ALDR398", "ALLI056"
)

#Subset to published
samples <- samples[ samples %in% hDRG$sample ]
Idents(hDRG) <- 'sample'
hDRG <- subset( hDRG, idents = samples)
unique(hDRG$sample) #92 samples

#Add umap
umap_coords <- finalumap@cell.embeddings
common_cells <- intersect(rownames(umap_coords), colnames(hDRG))
umap_sub <- umap_coords[common_cells, ]
hDRG[["umap"]] <- CreateDimReducObject(
  embeddings = umap_sub,
  key        = "umap_",           # matches the original key
  assay      = "RNA",
  global     = TRUE
)

#Collapse cell types for clean umap
hDRG$HighName <- sub("\\..*", "", hDRG$Name)
hDRG$HighName <- dplyr::case_when(
  
  # Neurons
  hDRG$HighName %in% c(
    "Neu","CNP","APEP","ALTMR","CLTMR","CThermo",
    "CPEP","PDIA2","AProp"
  ) ~ "Neurons",
  
  # Satellite glial cells
  hDRG$HighName %in% c("SGC","Glia") ~ "SatGCs",
  
  # Schwann cells
  hDRG$HighName == "NMSC" ~ "NMSchwann",
  hDRG$HighName == "MSC"  ~ "MSchwann",
  
  # Mononuclear phagocytes
  hDRG$HighName %in% c(
    "Mac","Mono","Monocytes","DendCells","MLC"
  ) ~ "MNPs",
  
  # Lymphocytes
  hDRG$HighName %in% c(
    "Tcell","Bcell","Plasma","NKcells","Lymph"
  ) ~ "Lymphocytes",
  
  # Fibroblasts
  hDRG$HighName == "Fibro" ~ "Fibroblasts",
  
  # Granulocytes
  hDRG$HighName %in% c("Granulocytes","MastCells") ~ "Granulocytes",
  
  # Adipocytes
  hDRG$HighName == "Adipo" ~ "Adipocytes",
  
  # Endothelial
  hDRG$HighName == "EC" ~ "Endothelial",
  
  # Mural cells
  hDRG$HighName == "MC" ~ "Mural",
  
  # SSC
  hDRG$HighName == "SSC" ~ "SSC",
  
  TRUE ~ "Other"
)

#Clear low qc neurons and fix name_label
Idents(hDRG) <- 'Name'
hDRG <- subset( hDRG, idents = c("Neu.LowCounts", "Neu.FADS2"), invert = T )
Idents(hDRG) <- 'Name_label'
map_to_category <- function(name) {
  case_when(
    grepl("ALTMR", name) ~ "ALTMR",
    grepl("AProp", name) ~ "AProp",
    grepl("APEP", name)  ~ "APEP",
    grepl("CLTMR", name) ~ "CLTMR",
    grepl("CPEP", name)  ~ "CPEP",
    grepl("CNP.MRGPRX1", name)   ~ "CNP.MRGPRX1",
    grepl("CNP.SST", name)   ~ "CNP.SST",
    grepl("CThermo", name) ~ "CThermo",
    grepl("ATF3", name)  ~ "ATF3",
    TRUE ~ name
  )
}
hDRG$Name_label <- map_to_category(as.character(Idents(hDRG)))
table(hDRG$Name_label)

#Clear low qc Schwann
Idents(hDRG) <- 'Name'
hDRG <- subset( hDRG, idents = c("NMSC.NRXN1.COL9A3"), invert = T )
hDRG <- subset( hDRG, idents = c("MSC.COL18A1"), invert = T )

#Adjust labels for MNPs
hDRG$Name_label <- dplyr::recode(hDRG$Name_label,
                                 'MLC.TMEM119'     = 'MLC(H)',
                                 'MLC.CD83'        = 'MLC(DA)',
                                 'Monocytes'       = 'Mono(Early)',
                                 'Mono.Mac.C1Q'    = 'Mono(Diff.)',
                                 'Mono.Mac.GPNMB'  = 'MDM(GPNMB)',
                                 'Mac.CD209.H'     = 'MDM(Early)',
                                 'Mac.CD209.L'     = 'MDM(Mature)',
                                 'DendCells'       = 'DendriticCells',
                                 'Mac.Stress'      = 'HeatShock',
                                 'MLC.Stress'      = 'MLC(NFKB2)',
                                 'Mac.MLC.Prolif'  = 'MNP.Mitotic'
)

#Fix other cell type names
Idents(hDRG) <- 'Name_label'
map_to_category <- function(name) {
  case_when(
    grepl("Adipo", name) ~ "Adipocytes",
    grepl("Plasma", name) ~ "PlasmaCells",
    grepl("Bcell", name)  ~ "Bcells",
    grepl("SGC.GABRB1", name) ~ "SGC.PILRB",
    grepl("NMSC.NRXN1", name) ~ "NMSC.NRXN1",
    grepl("NMSC.PXDN", name)  ~ "NMSC.PXDN",
    grepl("MSC.DOCK3", name)   ~ "MSC.DOCK3",
    grepl("MSC.COL18A1", name) ~ "MSC.DOCK3",
    grepl("MSC.ME1", name) ~ "MSC.ME1",
    grepl("MSC.IL17B", name) ~ "MSC.IL17B",
    grepl("MSC.COL18A1", name) ~ "MSC.IL17B",
    grepl("MSC.MAL", name)  ~ "MSC.IL17B",
    grepl("MSC.BZW2", name)  ~ "MSC.BZW2",
    grepl("Fibro.Endo", name)  ~ "Fibro.Endo",
    grepl("Tcell.LTK", name)  ~ "Tcell.LTK",
    grepl("Lymph.Prolif", name)  ~ "Lymph.Mitotic",
    TRUE ~ name
  )
}
hDRG$Name_label <- map_to_category(as.character(Idents(hDRG)))
hDRG$Name_label <- ifelse( hDRG$Name == 'Fibro.Epi.NFKB2', 'Fibro.Epi.NFKB2', hDRG$Name_label )
table(hDRG$Name_label)

#Rm sample specific clusters
hDRG <- subset(hDRG, subset = HighName != "SSC")

#Rm extra meta cols
colnames( hDRG@meta.data )
hDRG@meta.data <- hDRG@meta.data[ , c(1:4,6,12,16,17)]

#Save obj
saveRDS( hDRG, 'hDRG.rds' )
saveRDS( hDRG@meta.data, 'hDRG_meta.rds' )
gc()


#### Plots ####

#UMAP plot
cbfColors <- c(
  'Adipocytes'   = "#72190E",
  'Endothelial'  = "#DC050C",
  'Fibroblasts'  = "#E8601C",
  'Granulocytes'   = "#F6C141",
  'Lymphocytes'   = "#F7EE55",
  'MNPs'  = "#4EB265",
  'MSchwann'     = "#CAE0AB",
  'Mural'    = "#7BAFDE",
  'Neurons'      = "#1965B0",
  'NMSchwann'       = "#882E72",
  'SatGCs'   = "#BBBBBB"
)
DimPlot( hDRG, reduction = 'umap', group.by = 'HighName', 
         alpha = 0.1, cols = cbfColors, raster = F) +
  labs(title = 'Human DRG', x = 'UMAP1', y = 'UMAP2')
ggsave( 'all_umap.svg', width = 10, height = 8)



#### Ref for deconvolution ####

#Down-sample
set.seed(123)
keep_cells <- hDRG@meta.data %>% rownames_to_column("barcodes") %>%
  select(barcodes, sample, Name_label) %>% group_by(sample) %>%
  mutate(sample_total = n()) %>% ungroup() %>% filter(sample_total > 1) #doesnt take cells that are the only one in a sample (avoid errors downstream)
keep_cells <- keep_cells %>% group_by(sample, Name_label) %>%
  mutate(.rand = runif(n())) %>% arrange(.rand, .by_group = TRUE) %>%
  slice_head(n = 30) %>% ungroup() %>% pull(barcodes) #selects 20 cells per cluster per sample
ref.76K <- subset( hDRG, cells = keep_cells )

#Save obj
#ref.76K <- JoinLayers( ref.76K ) #not needed, already merged
saveRDS(ref.76K, 'ref.76K.rds')

#Save inputs for cell2loc
counts <- GetAssayData(ref.76K, slot = "counts")
meta.60k <- ref.76K@meta.data
write.csv(as.matrix(counts), "ref76k_counts.csv")
write.csv(meta.60k, "ref76k_meta.csv")

#Save counts as h5
library(Matrix)
library(hdf5r)
h5file <- H5File$new("ref76k_counts.h5", mode = "w")
h5file[["x"]] <- counts@x        # values
h5file[["i"]] <- counts@i        # row indices
h5file[["p"]] <- counts@p        # column pointers
h5file[["rownames"]] <- rownames(counts)
h5file[["colnames"]] <- colnames(counts)
h5file[["dims"]] <- dim(counts)
h5file$close_all()

# #Load h5 in python w:
# import h5py
# import scipy.sparse as sp
# import pandas as pd
# import numpy as np
# 
# with h5py.File("input_sn_hDRG/ref76k_counts.h5", "r") as f:
#   x    = f["x"][:]
# i    = f["i"][:]
# p    = f["p"][:]
# rows = f["rownames"][:].astype(str)
# cols = f["colnames"][:].astype(str)
# dims = f["dims"][:]
# 
# counts = sp.csc_matrix((x, i, p), shape=tuple(dims))
# counts_df = pd.DataFrame.sparse.from_spmatrix(counts, index=rows, columns=cols)
# print(counts_df.shape)
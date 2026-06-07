
#### Setup ####

library(Seurat)
library(ggplot2)
library(dplyr)

setwd('C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/microglia_like_cells/analysis/single_nuc_genes')

hDRG <- readRDS('C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/microglia_like_cells/analysis/drg_single_nuc/srat_objs/hDRG.rds')
hDRG_meta <- readRDS("C:/Users/mazmsi/Box/gataca/personal_analysis/Khadijah/my_projects/microglia_like_cells/analysis/drg_single_nuc/srat_objs/hDRG_meta.rds")
hDRG <- AddMetaData( hDRG, hDRG_meta)

# Cell type order: bottom -> top on y-axis (matches image layout)
label_order <- c(
  'Adipocytes',
  'Endothelial',
  'Fibroblasts',
  'Granulocytes',
  'Lymphocytes',
  'MNPs',
  'MSchwann',
  'Mural',
  'Neurons',
  'NMSchwann',
  'SatGCs'
)

tol_gradient <- c(
  "#FFF9C4",
  "#F7EE55",
  "#F6C141",
  "#F1932D",
  "#E8601C",
  "#DC050C",
  "#72190E"
)



#### Plot 1 — MLC / Microglia Markers ####

genes1 <- c('ADGRG1', 'WNT5A', 'P2RY12', 'TMEM119', 'CX3CR1', 'CD83', 'CCL3')

plot1 <- DotPlot(hDRG, scale = F,
                 features = genes1,
                 group.by = 'HighName') +
  scale_color_gradientn(colours = tol_gradient, limits = c(0, 3), oob = scales::squish) +
  scale_size_continuous(limits = c(0, 60), range = c(0, 6)) +
  scale_y_discrete(limits = label_order) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
        axis.title.x = element_text(face = 'italic'),
        axis.text.y = element_text(hjust = 0)) +
  labs(x = NULL, y = NULL)

print(plot1)
ggsave('hdrg_mlc_markers.svg', plot = plot1, width = 6, height = 4)



#### Plot 2 — Pan-Immune / MNP Markers ####

genes2 <- c('ADGRG1', 'WNT5A', 'P2RY12', 'TMEM119', 'CX3CR1', 'CD83', 'CCL3', 
            'PTPRC', 'ITGAM', 'ITGAX', 'CD14', 'CD68', 'CD80', 'CD86', 'CD40', 'ICAM1', 'AIF1')

plot2 <- DotPlot(hDRG, scale = F,
                 features = genes2,
                 group.by = 'HighName') +
  scale_color_gradientn(colours = tol_gradient) +
  scale_y_discrete(limits = label_order) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
        axis.title.x = element_text(face = 'italic'),
        axis.text.y = element_text(hjust = 0)) +
  labs(x = NULL, y = NULL)

print(plot2)
ggsave('hdrg_immune_markers.svg', plot = plot2, width = 9, height = 3.5)

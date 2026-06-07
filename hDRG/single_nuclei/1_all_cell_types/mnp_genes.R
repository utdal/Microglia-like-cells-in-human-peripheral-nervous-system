
#### Setup ####

#Libraries
library(SingleCellExperiment)
library(Seurat)
library(slingshot)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(openxlsx)

hDRG <- readRDS("C:/Users/kxm130930/Box/gataca/personal_analysis/Khadijah/my_projects/microglia_like_cells/analysis/drg_single_nuc/srat_objs/hDRG/hDRG.rds")
hDRG <- NormalizeData( hDRG )

#Take a look
# Define the desired order for cell types
genes <- c( 'PTPRC', 'CD14', 'CD68', 'ITGAM', 'ITGAX','CD80', 'CD86', 'CD40', 
            'ICAM1', 'TLR1', 'TLR2', 'TLR3', 'TLR4', 'TLR9', 'AIF1',
            'WNT5A', 'KCNQ3', 'ADORA3', 'PADI2', 'ABCC4', 'ADGRG1', 
            'P2RY12', 'TMEM119', 'CX3CR1', 'TREM2', 'CD83', 'CCL3',
            'FABP7', 'S100B', 'THY1', 'ANPEP', 'PI16', 'CD3D', 'CD3E', 'CD3G',
            'PDGFRA', 'DPP4')


# Create the dot plot with reordered cell types and formatted gene names
your_plot <- DotPlot(hDRG, scale = F,
                     features = genes,
                     group.by = "Name_label")  +
  scale_color_gradientn(
    colours = c("#FFF9C4",
                "#F7EE55",  # yellow (Tol #8)
                "#F6C141",  # golden yellow (Tol #9)
                "#F1932D",  # orange (Tol #10)
                "#E8601C",  # dark orange (Tol #11)
                "#DC050C",  # red (Tol #12)
                "#72190E")  # dark red (Tol #13)
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),  # Rotate and italicize
        axis.title.x = element_text(face = "italic")) +  # Italicize x-axis title if needed
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_text(hjust = 0))

#Save
print(your_plot)
ggsave('mult_markers.svg', width = 12, height = 12)

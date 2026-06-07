
#Libraries
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(scales)
library(plyr)
library(harmony)
library(pheatmap)
library(patchwork)
library(ggalluvial)



# ============================================================
# PROCESSING MICROGLIA
# ============================================================

#### Load Data From Sun Et Al Paper ----

#Load full data / all MNPs
metaROSMAP <- readRDS(url("https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds?dl=0"))
countsROSMAP <- readRDS(url("https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.counts.rds?dl=0"))

#Load meta for their selected microglia after QC
#They had already done filtering in final_rosmap:
#QC filtering (counts>500 and mt<5)
#Filter for genes expressed in at least 50 cells
#DoubletFinder and manual contam cluster removal
metaMicroglia <- read.table(
  "https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.Microglia.6regions.seurat.harmony.selected.deidentified.metadata.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  stringsAsFactors = FALSE
)

#Make srat obj and clear mem
rosmap <- CreateSeuratObject( countsROSMAP, meta.data = metaROSMAP)
rm(countsROSMAP, metaROSMAP)
gc()

#Get the BAMs 
table(rosmap$seurat_clusters) # BAMS = clus 9, 4422 cells
Idents(rosmap) <- 'seurat_clusters'
bams <- subset( rosmap, idents = 9) # 4422 cells

#Get their final microglia
final_rosmap <- subset( rosmap, cells = rownames(metaMicroglia) )

#Merge and clear mem
brainMNPs <- merge( final_rosmap, bams )
rm( final_rosmap, bams, rosmap, metaMicroglia)
gc()



#### Cluster MNPs + UMAP ####

#Split obj into batch layers
brainMNPs <- JoinLayers(brainMNPs)
brainMNPs[['RNA']] <- split( brainMNPs[['RNA']] , f = brainMNPs$batch )

#Seurat pipeline w/Harmony
brainMNPs <- NormalizeData(brainMNPs)
brainMNPs <- FindVariableFeatures(brainMNPs)
brainMNPs <- ScaleData(brainMNPs)
brainMNPs <- RunPCA(brainMNPs, npcs = 100)
brainMNPs <- IntegrateLayers(brainMNPs, method = HarmonyIntegration,
                             orig.reduction = "pca", new.reduction = "harmony")
brainMNPs <- FindNeighbors( brainMNPs, reduction = 'harmony', dims = 1:80)
brainMNPs <- FindClusters( brainMNPs, resolution = 0.8)
brainMNPs <- RunUMAP(brainMNPs, reduction = "harmony", dims = 1:80, 
                     return.model = TRUE) #need model for mapping later
brainMNPs <- JoinLayers(brainMNPs)
DimPlot(brainMNPs, reduction = "umap", group.by = 'seurat_clusters', 
        label=T, raster=FALSE, alpha = 0.3)



#### Label Clusters ####

#Look from their markers in the paper/DEG list
genes <- unique( c(
  #BAM, Viral, Cycling, Stress
  "MRC1", "LYVE1", "IFI44L", "MX1", "EZH2", "BRIP1", "HSPB1", "HSPH1",
  #Surveillance, Homeostatic
  "P2RY12", "CX3CR1", "PRDM11", "TANC1",
  #Glycolytic, Lipid, Ribosome, Phagocytic
  "NAMPT", "SLC2A3", "GPNMB", "PTPRG", "RPL32", "RPL19", "CD163", "F13A1", 
  #Inflamm 1-3
  "TMEM163", "ERC2", "DUSP1", "SPON1", "CD83", "CCL3",
  #uG v BAM genes
  "ADGRG1", "WNT5A"
) )
DotPlot(brainMNPs, features = genes, scale = F, group.by = 'seurat_clusters') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Name the clusters
cluster_labels <- c(
  "0"  = "H.Surveillance",           # MG1 - P2RY12+, CX3CR1+ high
  "1"  = "H.Homeostatic",            # MG0
  "2"  = "D.InflammI",               # MG2 - TMEM163+, ERC2+
  "3"  = "H.Homeostatic",            # MG0
  "4"  = "D.LipidProcess",           # MG4 - GPNMB+
  "5"  = "D.InflammII",              # MG8 - SPON1+RGS1+
  "6"  = "D.Ribosome",               # MG3 - RPL+
  "7"  = "D.Phagocytic",             # MG5 - CD163+F13A1+
  "8"  = "BAM",                      # MG9 - MRC1+, LYVE1+
  "9"  = "D.InflammI",               # MG2 - TMEM163+, ERC2+, also HSP+
  "10" = "D.Ribosome",               # MG3 - RPL+
  "11" = "StressHS",                 # MG6 - HSP+
  "12" = "StressHS",                 # MG6 - HSP+
  "13" = "D.Glycolytic",             # MG7 - NAMPT+SLC2A3+
  "14" = "Antiviral",                # MG11 - IFI44L+
  "15" = "D.Glycolytic",             # MG7 - NAMPT+SLC2A3+
  "16" = "D.InflammIII",             # MG10 -  CCL3+CD83+IL1B+NFKB1+TNFAIP3+
  "17" = "D.Glycolytic",             # MG7 - NAMPT+SLC2A3+, also HSP+
  "18" = "D.InflammIII",             # MG10 - most similar to clus 16
  "19" = "Cycling",                  # MG12 - mitotic
  "20" = "H.Homeostatic",            # MG0
  "21" = "H.Homeostatic",            # MG0 - CTSS, CST3
  "22" = "D.InflammI"                # MG2 - TMEM163+, also HSP+
)
Idents(brainMNPs) <- 'seurat_clusters'
brainMNPs <- RenameIdents(brainMNPs, cluster_labels)
brainMNPs$Name <- Idents(brainMNPs)

#Order cell types
brainMNPs$Name <- factor(brainMNPs$Name, 
                          levels = c("BAM", "Antiviral", "Cycling", "StressHS",
                                     "H.Homeostatic", "H.Surveillance",
                                     "D.Glycolytic", "D.LipidProcess",
                                     "D.Ribosome", "D.Phagocytic", 
                                     "D.InflammI", "D.InflammII", "D.InflammIII"))



#### Clustering Plots ####

#Final umap
cbfColors <- c("#BBBBBB","#1965B0","#5289C7","#7BAFDE","#4EB265","#90C987","#CAE0AB",
  "#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C","#72190E")

DimPlot( brainMNPs, reduction = 'umap', group.by = 'Name', 
         alpha = 0.1, cols = cbfColors, raster = F, label = T) +
  labs(title = 'Human Brain MNPs', x = 'UMAP1', y = 'UMAP2')
ggsave( 'brain_umap.svg', width = 10, height = 8)

#Final dot plot
DotPlot(brainMNPs, features = genes, scale = F, group.by = 'Name') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('sun_markers_dot.svg', width = 10, height = 3.5)

#Save srat obj
saveRDS( brainMNPs, 'brainMNPs.rds' )



# ============================================================
# DAM GENES BUBBLE PLOTS
# ============================================================

#### Setup for Bubble Plots / DAM Genes Comparison ####

#Load hDRG obj
mlc <- readRDS("../hDRG/mlc.rds")

#DAM gene list
gene_order <- c(
  #Homeo
  "CX3CR1", "P2RY12", "TMEM119", "CSF1R", "CST3", "CTSS",
  "SPARC", "C1QB", "C1QA", "HEXB", "TMSB4X", "CTSD", "CST7",
  #DAM
  "AXL", "LILRB4", "TIMP2", "CLEC7A", "TREM2", "ITGAX",
  "LYZ", "APOE", "B2M", "TYROBP", "CTSB", "FTH1",
  "CTSL", "CD9", "LPL", "CSF1"
)

#List of ug clusters
ugClus <- c("D.InflammIII", "D.Glycolytic", "D.InflammI", 
            "D.Phagocytic", "D.InflammII", "D.LipidProcess", 
            "D.Ribosome", "H.Surveillance", "H.Homeostatic")

#Shared plot elements
color_scale <- scale_fill_gradient2(
  low      = "#2166AC",
  mid      = "lightyellow",
  high     = "#B2182B",
  midpoint = 0,
  limits   = c(-2, 2),
  oob      = scales::squish,
  name     = "avg log2FC",
  labels   = function(x) ifelse(x <= -2, "≤-2", ifelse(x >= 2, "≥2", x)),
  guide    = guide_colorbar(direction = "horizontal", title.position = "top")
)
size_scale <- scale_size_continuous(
  name   = "p.adj",
  breaks = c(1, 10, 100, 300),
  labels = c("1e-1", "1e-10", "1e-100", "<1e-100"),
  range  = c(0.5, 5),
  guide  = guide_legend(direction = "horizontal", title.position = "top")
)
base_theme <- theme_minimal() +
  theme(
    axis.text.x        = element_text(color = "black", face = "italic", angle = 30, hjust = 1),
    axis.text.y        = element_text(color = "black"),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray85"),
    panel.grid.minor   = element_blank(),
    text               = element_text(color = "black"),
    legend.position    = "bottom",
    legend.box         = "horizontal"
  )



#### Bubble Plot 1: DAM genes in brain MNPs (relative to H.Homeostatic) ####

#Set genes and clusters
genes <- gene_order[gene_order %in% rownames(brainMNPs)] #CST7, CTSD not in brainMNPs dataset
cluster_names_1 <- ugClus[-c(9)]

#Run DE for gene list
comparisons_1 <- lapply(cluster_names_1, function(cl) {
  list(ident.1 = cl, ident.2 = "H.Homeostatic", label = cl)
  })
all_results_1 <- lapply(comparisons_1, function(comp) {
    FindMarkers(
      brainMNPs,
      group.by        = 'Name',
      ident.1         = comp$ident.1,
      ident.2         = comp$ident.2,
      features        = genes,
      logfc.threshold = 0,
      min.pct         = 0,
      test.use        = "MAST",
      latent.vars     = "brainRegion"
    ) %>%
      tibble::rownames_to_column("gene") %>%
      mutate(comparison = comp$label)
  })
combined_df_1 <- bind_rows(all_results_1)

#Bubble plot
plot_df_1 <- combined_df_1 %>%
  mutate(
    neg_log10_padj = pmin(-log10(p_val_adj + 1e-300), 300),
    gene           = factor(gene, levels = gene_order),
    comparison     = factor(comparison, levels = cluster_names_1)
  )
ggplot(plot_df_1, aes(x = comparison, y = gene)) +
  geom_point(aes(size = neg_log10_padj, fill = avg_log2FC),
             shape = 21, color = "black", stroke = 0.3) +
  color_scale + size_scale + coord_flip() + base_theme +
  labs(x = NULL, y = NULL,
       title = "Expression of Microglial Genes Relative to H.Homeostatic Population")
ggsave('dam_genes_vsHomeo.svg', width = 8, height = 3.5)



#### Bubble Plot 2: DAM genes in brain MNPs (relative to H.Surveillance) ####

#Set genes and clusters
genes <- gene_order[gene_order %in% rownames(brainMNPs)] #CST7, CTSD not in brainMNPs dataset
cluster_names_2 <- ugClus[-c(8)]

#Run DE for gene list
comparisons_2 <- lapply(cluster_names_2, function(cl) {
  list(ident.1 = cl, ident.2 = "H.Surveillance", label = cl)
})
all_results_2 <- lapply(comparisons_2, function(comp) {
  FindMarkers(
    brainMNPs,
    group.by        = 'Name',
    ident.1         = comp$ident.1,
    ident.2         = comp$ident.2,
    features        = genes,
    logfc.threshold = 0,
    min.pct         = 0,
    test.use        = "MAST",
    latent.vars     = "brainRegion"
  ) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(comparison = comp$label)
})
combined_df_2 <- bind_rows(all_results_2)

#Bubble plot
plot_df_2 <- combined_df_2 %>%
  mutate(
    neg_log10_padj = pmin(-log10(p_val_adj + 1e-300), 300),
    gene           = factor(gene, levels = gene_order),
    comparison     = factor(comparison, levels = cluster_names_2)
  )
ggplot(plot_df_2, aes(x = comparison, y = gene)) +
  geom_point(aes(size = neg_log10_padj, fill = avg_log2FC),
             shape = 21, color = "black", stroke = 0.3) +
  color_scale + size_scale + coord_flip() + base_theme +
  labs(x = NULL, y = NULL,
       title = "Expression of Microglial Genes Relative to H.Surveillance Population")
ggsave('dam_genes_vsSurveil.svg', width = 8, height = 3.5)



#### Bubble Plot 3: DAM genes in MLCs ####

#Set genes and clusters 
cluster_names_mlc <- c("MLC(DA) vs MLC(NFKB2)", "MLC(H) vs MLC(NFKB2)",
                       "MLC(H) vs MLC(DA)", "MG(H) vs MG(D)")
genes <- gene_order[gene_order %in% rownames(mlc)]
  
#Run DE for gene list in MLCs
mlc_comparisons <- list(
    list(ident.2 = "MLC(H)",  ident.1 = "MLC(DA)",    label = "MLC(H) vs MLC(DA)"),
    list(ident.2 = "MLC(H)",  ident.1 = "MLC(NFKB2)", label = "MLC(H) vs MLC(NFKB2)"),
    list(ident.2 = "MLC(DA)", ident.1 = "MLC(NFKB2)", label = "MLC(DA) vs MLC(NFKB2)")
  )
mlc_results <- lapply(mlc_comparisons, function(comp) {
    FindMarkers(
      mlc,
      group.by        = 'Name_label',
      ident.1         = comp$ident.1,
      ident.2         = comp$ident.2,
      features        = genes,
      logfc.threshold = 0,
      min.pct         = 0,
      test.use        = "MAST",
      latent.vars     = "sample"
    ) %>%
      tibble::rownames_to_column("gene") %>%
      mutate(comparison = comp$label)
  })

#Run DE for gene list in MG
d_clusters <- c("D.InflammIII", "D.Glycolytic", "D.InflammI", "D.Phagocytic",
                  "D.InflammII", "D.LipidProcess", "D.Ribosome")
h_clusters <- c("H.Homeostatic", "H.Surveillance")
brainMNPs$Name2 <- ifelse( brainMNPs$Name %in% d_clusters, "MG(D)", brainMNPs$Name)
brainMNPs$Name2 <- ifelse( brainMNPs$Name %in% h_clusters, "MG(H)", brainMNPs$Name)
mg_result <- FindMarkers(
    brainMNPs,
    group.by        = 'Name2',
    ident.1         = "MG(D)",
    ident.2         = "MG(H)",
    features        = genes,
    logfc.threshold = 0,
    min.pct         = 0,
    test.use        = "MAST",
    latent.vars     = "brainRegion"
  ) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(comparison = "MG(H) vs MG(D)")
  
#Join Results
combined_df_3 <- bind_rows(mlc_results, list(mg_result))

#Plot
plot_df_3 <- combined_df_3 %>%
  mutate(
    neg_log10_padj = pmin(-log10(p_val_adj + 1e-300), 300),
    gene           = factor(gene, levels = genes),
    comparison     = factor(comparison, levels = cluster_names_mlc)
  )
ggplot(plot_df_3, aes(x = comparison, y = gene)) +
  geom_point(aes(size = neg_log10_padj, fill = avg_log2FC),
             shape = 21, color = "black", stroke = 0.3) +
  color_scale + size_scale + coord_flip() + base_theme +
  labs(x = NULL, y = NULL,
       title = "Expression of Microglial Genes in MLCs and MG")
ggsave('dam_genes_mlc.svg', width = 8, height = 2.5)



# ============================================================
# LABEL TRANSFER TO MLCS
# ============================================================

#### Seurat MapQuery pipeline ####

#Load hDRG MNPs
drgMNPs <- readRDS("../hDRG/immune.rds")
drgMNPs <- NormalizeData(drgMNPs)
# both ref and query are merged (not split in sample layers)

#Label transfer / map
anchors <- FindTransferAnchors(
  reference        = brainMNPs,
  query            = drgMNPs,
  reference.reduction  = "pca",
  dims             = 1:30
)
drgMNPs <- MapQuery(
  anchorset       = anchors,
  query           = drgMNPs,
  reference       = brainMNPs,
  refdata         = list(Name = "Name"),   # transfers Name labels
  reference.reduction = "pca",
  reduction.model = "umap"                 # your named UMAP in the ref
)
# Transferred labels land in:  drgMNPs$predicted.Name
# Prediction scores land in:   drgMNPs$predicted.Name.score
# Projected UMAP coords are in: drgMNPs[["ref.umap"]]

#Save meta
saveRDS( drgMNPs@meta.data, 'drgMNPs.meta.rds')
saveRDS( drgMNPs[['ref.umap']], 'ref.umap.rds')
saveRDS( drgMNPs[['ref.pca']], 'ref.pca.rds')


#### Mapped UMAP Plots ####

#Define colors for brain MNPs
cbfColors <- c("#BBBBBB","#1965B0","#5289C7","#7BAFDE","#4EB265","#90C987","#CAE0AB",
               "#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C","#72190E")
cell_types   <- levels(factor(brainMNPs$Name))
named_colors <- setNames(cbfColors, cell_types)

#Define colors for hDRG mNPs
name_label_colors <- c(
  'Mono(Early)'    = "#72190E",
  'Mono(Diff.)'    = "#DC050C",
  'MDM(GPNMB)'     = "#E8601C",
  'MDM(Early)'     = "#F6C141",
  'MDM(Mature)'    = "#F7EE55",
  'DendriticCells' = "#CAE0AB",
  'HeatShock'      = "#4EB265",
  'Mitotic'        = "#BBBBBB",
  'MLC(DA)'        = "#7BAFDE",
  'MLC(H)'         = "#1965B0",
  'MLC(NFKB2)'     = "#882E72"
)

# ── p1: Reference – Name ─────────────────────────────────────────────────
p1 <- DimPlot(brainMNPs, reduction = "umap", alpha = 0.3, 
              group.by = "Name", label = TRUE, repel = TRUE) +
  scale_color_manual(values = named_colors) +
  labs(title = NULL, x = "UMAP1", y = "UMAP2") +
  NoLegend()

# ── p2: Query – predicted Name ───────────────────────────────────────────
p2 <- DimPlot(drgMNPs, reduction = "ref.umap", alpha = 0.3,
              group.by = "predicted.Name") +
  scale_color_manual(values = named_colors) +
  labs(title = NULL, x = "UMAP1", y = "UMAP2") +
  NoLegend()

# ── p3: Query – original Name_label ───────────────────────────────────────
p3 <- DimPlot(drgMNPs, reduction = "ref.umap", alpha = 0.3,
              group.by = "Name_label") +
  scale_color_manual(values = name_label_colors) +
  labs(title = NULL, x = "UMAP1", y = "UMAP2")

# ── p4: Prediction score ───────────────────────────────────────────────────
p4 <- FeaturePlot(drgMNPs, reduction = "ref.umap", order = F, pt.size = 0.1,
                  features = "predicted.Name.score") +
  scale_color_gradientn(colors = c("#F7EE55", "#F6C141", "#DC050C", "#72190E")) +
  labs(title = NULL, x = "UMAP1", y = "UMAP2")

# ── Combine ────────────────────────────────────────────────────────────────
combined <- (p1 | p2) / (p4 | p3)
ggsave( 'mapquery_umap.svg', plot = combined, width = 17, height = 14)



#### Label Transfer Plots ####

# ── Shared data prep ───────────────────────────────────────────────────────
label_trans_df <- data.frame(
  Name_label = drgMNPs$Name_label,
  predicted  = ifelse(drgMNPs$predicted.Name.score >= 0.5,
                      as.character(drgMNPs$predicted.Name),
                      "Unassigned") #unassigned if score < 0.5
)
plot_df <- label_trans_df %>%
  dplyr::count(Name_label, predicted) %>%
  dplyr::group_by(Name_label) %>%
  dplyr::mutate(proportion = n / sum(n)) %>%
  dplyr::ungroup()

# ── Colors for predicted.Name (reuse named_colors + Unassigned) ──────────
pred_levels   <- unique(plot_df$predicted)
pred_colors   <- named_colors[pred_levels]
pred_colors["Unassigned"] <- "#DDDDDD"   # light grey for unassigned

# ── Define desired order ───────────────────────────────────────────────────
pred_order <- c("BAM", "Antiviral", "Cycling", "StressHS",
                "H.Homeostatic", "H.Surveillance",
                "D.Glycolytic", "D.LipidProcess",
                "D.Ribosome", "D.Phagocytic",
                "D.InflammI", "D.InflammII", "D.InflammIII", "Unassigned")
plot_df$predicted <- factor(plot_df$predicted, levels = pred_order)
pred_colors <- pred_colors[pred_order] # reorder pred_colors to match

# ── Plot 1: Heatmap ────────────────────────────────────────────────────────
p_heat <- ggplot(plot_df, aes(x = predicted, y = Name_label, fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(proportion > 0.05,
                               scales::percent(proportion, accuracy = 1), "")),
            size = 3, color = "white") +
  scale_x_discrete(limits = pred_order) +
  scale_fill_gradientn(colors = c("#F7F7F7", "#F7EE55", "#DC050C", "#72190E"),
                       labels = scales::percent) +
  labs(x = "CNS MNPs", y = "DRG MNPs", fill = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank())

# ── Plot 2: Stacked bar ────────────────────────────────────────────────────
p_bar <- ggplot(plot_df, aes(x = Name_label, y = proportion, fill = predicted)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pred_colors, breaks = pred_order) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion", fill = "CNS MNPs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())

# ── Plot 3: Bubble plot ────────────────────────────────────────────────────
score_df <- data.frame(
  Name_label = drgMNPs$Name_label,
  predicted  = ifelse(drgMNPs$predicted.Name.score >= 0.5,
                      as.character(drgMNPs$predicted.Name),
                      "Unassigned"),
  score      = drgMNPs$predicted.Name.score
) %>%
  dplyr::group_by(Name_label, predicted) %>%
  dplyr::summarise(mean_score = mean(score), .groups = "drop")

plot_df <- dplyr::left_join(plot_df, score_df, by = c("Name_label", "predicted"))
plot_df$predicted <- factor(plot_df$predicted, levels = pred_order)

p_bubble <- ggplot(plot_df, aes(x = predicted, y = Name_label,
                                size = proportion, color = mean_score)) +
  geom_point(alpha = 0.9) +
  scale_x_discrete(limits = pred_order) +
  scale_size_area(max_size = 14, labels = scales::percent) +
  scale_color_gradientn(colors = c("#F7F7F7", "#F7EE55", "#DC050C", "#72190E"),
                        limits = c(0, 1)) +
  labs(x = "CNS MNPs", y = "DRG MNPs",
       size = "Proportion", color = "Mean prediction\nscore") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_line(color = "grey92"))

# ── Save all three ─────────────────────────────────────────────────────────
ggsave("mapping_heatmap.svg",    plot = p_heat,   width = 8, height = 5.5)
ggsave("mapping_stackedbar.svg", plot = p_bar,    width = 4.5, height = 4)
ggsave("mapping_bubble.svg",     plot = p_bubble, width = 8, height = 5.5)

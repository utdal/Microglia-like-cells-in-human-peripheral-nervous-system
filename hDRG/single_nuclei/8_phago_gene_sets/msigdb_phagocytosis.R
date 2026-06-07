
#### Setup ####

#Libraries
library(SingleCellExperiment)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(readxl)
library(stringr)
library(msigdbr)
library(patchwork)
library(ggh4x)
library(tidyverse)

#Load full data
hDRG <- readRDS("../hDRG.rds")
hDRG <- NormalizeData(hDRG)



#### SEARCH FOR GENE SETS ####

# Load db
msig <- msigdbr(species = "Homo sapiens")

# Broad search terms to capture all relevant terms + save results
search_pattern <- paste(
  "phago", "efferocyt", "apoptotic.cell", "apoptotic.corp",
  "engulf", "clearance", "scaveng", "complement", "opsoniz",
  "fc_gamma", "fcgr", "fc_receptor", "opsoniz", "amyloid",
  "killing", "cytotox", "trail", "death_domain", "death_receptor",
  sep = "|"
)

#Do search
all_hits <- msig %>%
  filter(grepl(search_pattern, gs_name, ignore.case = TRUE)) %>%
  distinct(gs_name, gs_collection, gs_subcollection, gene_symbol) %>%
  group_by(gs_collection, gs_subcollection, gs_name) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  arrange(gs_collection, gs_subcollection, n_genes)

# Summary of libraries
all_hits %>%
  group_by(gs_collection, gs_subcollection) %>%
  summarise(n_sets = n(), .groups = "drop") %>%
  print()

#Take a look
print(all_hits, n = 200)

#Save
write.xlsx(all_hits, "msigdb_phagocytosis_search_all_hits.xlsx", row.names = FALSE)

#Manually annotated hits in the xlsx file and selected most relevant



#### DEFINE SELECTED GENE SETS ####

#Load selected hits
selected_hits <- read_excel("msigdb_phagocytosis_search_all_hits.xlsx", sheet = 2)[, c(3, 5)]

#Get genes
selected_genes <- msig %>%
  filter(gs_name %in% selected_hits$gs_name) %>%
  distinct(gs_name, gene_symbol) %>%
  left_join(selected_hits, by = "gs_name")

#Save to file
wb <- loadWorkbook("msigdb_phagocytosis_search_all_hits.xlsx")
addWorksheet(wb, "selected_genes")
writeData(wb, sheet = "selected_genes", selected_genes)
saveWorkbook(wb, "msigdb_phagocytosis_search_all_hits.xlsx", overwrite = TRUE)



#### ADD MODULE SCORES PER CATEGORY ####

#category names
categories <- unique(selected_hits$category)

for (cat in categories) {
  cat_genes <- selected_genes %>%
    filter(category == cat) %>%
    distinct(gene_symbol) %>%
    pull(gene_symbol)
  
  hDRG <- AddModuleScore(
    hDRG,
    features = list(cat_genes),
    name = cat,
    seed = 42
  )
  old_col <- paste0(cat, "1")
  hDRG@meta.data[[cat]] <- hDRG@meta.data[[old_col]]
  hDRG@meta.data[[old_col]] <- NULL
}
# Took note of genes not found, saved in missing_genes.txt

#Save
saveRDS(hDRG@meta.data, "hDRG_meta_scores.rds")



#### BARPLOTS PER CATEGORY - ALL CELLS ####

plot_list_all <- list()
for (cat in categories) {
  meta <- hDRG@meta.data %>%
    group_by(HighName) %>%
    summarise(
      mean_score = mean(.data[[cat]]),
      se         = sd(.data[[cat]]) / sqrt(n()),
      n_cells    = n()
    )
  
  plot_list_all[[cat]] <- ggplot(meta, aes(x = reorder(HighName, mean_score), y = mean_score)) +
    geom_col(aes(fill = mean_score)) +
    geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.3) +
    scale_fill_gradientn(
      colours = c("#FFF9C4", "#F7EE55", "#F6C141",
                  "#F1932D", "#E8601C", "#DC050C", "#72190E")
    ) +
    coord_flip() +
    labs(title = cat, x = NULL, y = "Mean Module Score") +
    theme_classic(base_size = 7) +
    theme(legend.position = "none")
}

combined_all <- wrap_plots(plot_list_all, ncol = 4)
ggsave("all_cells_highLevel.svg", plot = combined_all,
       width = 10, height = 9)



#### BARPLOTS PER CATEGORY - SELECTED CELLS ####

#Look at subtypes for MNPs and satellite cells
cells_of_interest <- c(
  'DendriticCells', 'Mono(Early)', 'Mono(Diff.)',
  'MDM(GPNMB)', 'MDM(Early)', 'MDM(Mature)',
  'MLC(H)', 'MLC(DA)', 'MLC(NFKB2)',
  'MNP.Mitotic', 'HeatShock', 'SGC.PILRB',
  'SGC.SERPINA5', 'NMSC.PXDN'
)
hDRG_sub <- subset(hDRG, cells = which(hDRG$Name_label %in% cells_of_interest))

plot_list_sub <- list()

for (cat in categories) {
  meta <- hDRG_sub@meta.data %>%
    group_by(Name_label) %>%
    summarise(
      mean_score = mean(.data[[cat]]),
      se         = sd(.data[[cat]]) / sqrt(n()),
      n_cells    = n()
    ) %>%
    arrange(desc(mean_score))
  
  plot_list_sub[[cat]] <- ggplot(meta, aes(x = reorder(Name_label, mean_score), y = mean_score)) +
    geom_col(aes(fill = mean_score)) +
    geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.3) +
    scale_fill_gradientn(
      colours = c("#FFF9C4", "#F7EE55", "#F6C141",
                  "#F1932D", "#E8601C", "#DC050C", "#72190E")
    ) +
    coord_flip() +
    labs(title = cat, x = NULL, y = "Mean Module Score") +
    theme_classic(base_size = 7) +
    theme(legend.position = "none")
}

combined_sub <- wrap_plots(plot_list_sub, ncol = 4)
ggsave("sub_highLevel.svg", plot = combined_sub, width = 10, height = 9)



#### DOTPLOT WITH Z-SCORES - ALL CELLS ####
z_threshold <- 1

dot_data_all <- purrr::map_dfr(categories, function(cat) {
  # cell-level z-scores across all cells for this category
  all_vals <- hDRG@meta.data[[cat]]
  all_z    <- scale(all_vals)[,1]
  
  hDRG@meta.data %>%
    mutate(cell_z = all_z) %>%
    group_by(HighName) %>%
    summarise(
      mean_score = mean(.data[[cat]]),
      pct_above  = mean(cell_z > z_threshold) * 100,
      n_cells    = n(),
      .groups    = "drop"
    ) %>%
    mutate(
      Category = cat,
      z_score  = scale(mean_score)[,1]  # z-score across cell types
    )
})

p_dot_all <- ggplot(dot_data_all, aes(x = Category, y = HighName, color = z_score, size = pct_above)) +
  geom_point() +
  scale_color_gradientn(
    colours = c("#FFF9C4", "#F7EE55", "#F6C141",
                "#F1932D", "#E8601C", "#DC050C", "#72190E"),
    name = "Z-score\n(mean per cell type)"
  ) +
  scale_size_continuous(
    name = paste0("% cells > ", z_threshold, " z"),
    range = c(0.5, 8)
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Module Scores by Cell Type and Category")

ggsave("dotplot_all.svg", plot = p_dot_all, width = 8, height = 5)



#### DOTPLOT WITH Z-SCORES - SELECTED CELLS ####
dot_data_sub <- purrr::map_dfr(categories, function(cat) {
  all_vals <- hDRG_sub@meta.data[[cat]]
  all_z    <- scale(all_vals)[,1]
  
  hDRG_sub@meta.data %>%
    mutate(cell_z = all_z) %>%
    group_by(Name_label) %>%
    summarise(
      mean_score = mean(.data[[cat]]),
      pct_above  = mean(cell_z > z_threshold) * 100,
      n_cells    = n(),
      .groups    = "drop"
    ) %>%
    mutate(
      Category = cat,
      z_score  = scale(mean_score)[,1]
    )
})

#order cell types
dot_data_sub <- dot_data_sub %>%
  mutate(Name_label = factor(Name_label, levels = rev(cells_of_interest)))

p_dot_sub <- ggplot(dot_data_sub, aes(x = Category, y = Name_label, color = z_score, size = pct_above)) +
  geom_point() +
  scale_color_gradientn(
    colours = c("#FFF9C4", "#F7EE55", "#F6C141",
                "#F1932D", "#E8601C", "#DC050C", "#72190E"),
    name = "Z-score\n(mean per cell type)"
  ) +
  scale_size_continuous(
    name = paste0("% cells > ", z_threshold, " z"),
    range = c(0.5, 8)
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Module Scores by Cell Type and Category")

ggsave("dotplot_sub.svg", plot = p_dot_sub, width = 10, height = 6)



#### GENE EXPRESSION TABLE (LONG) ####

# Define groups (needed for group lookup in table)
groups <- list(
  group1 = c("General Phagocytosis", "Phagocytosis Recognition",
             "Phagocytic Engulfing", "Phagosome", "Phagosome ROS"),
  group2 = c("Inhibition of Phagocytosis", "Inhibition of Phagosome"),
  group3 = c("Fc Phagocytosis", "Opsonization",
             "Pathogen Phagocytosis", "RBC Phagocytosis"),
  group4 = c("Cytotoxicity", "Apoptosis Induction",
             "Inhibition of Apoptotic Induction", "Apoptotic Cell Phagocytosis"),
  group5 = c("Scavenger Receptors", "Lipid Clearance",
             "Amyloid Phagocytosis", "Amyloid Response")
)

group_titles <- c(
  group1 = "Phagocytosis Mechanism",
  group2 = "Inhibition of Phagocytosis",
  group3 = "Other Phagocytosis",
  group4 = "Cytotoxicity and Apoptosis",
  group5 = "Clearance and Amyloid"
)

# Get normalized expression for all genes in selected sets
genes_to_plot <- unique(selected_genes$gene_symbol)
genes_present <- genes_to_plot[genes_to_plot %in% rownames(hDRG_sub)]

# Extract expression matrix once
expr_mat_sub <- hDRG_sub[['RNA']]$data[genes_present, ]

# Build full df for sub once
df_sub_full <- purrr::map_dfr(genes_present, function(g) {
  expr <- expr_mat_sub[g, ]
  data.frame(
    gene      = g,
    cell_type = hDRG_sub@meta.data$Name_label,
    expr      = as.numeric(expr)
  )
}) %>%
  group_by(gene, cell_type) %>%
  summarise(
    mean_exp = mean(expr),
    pct_exp  = mean(expr > 0) * 100,
    .groups  = "drop"
  ) %>%
  left_join(
    selected_genes %>% distinct(gene_symbol, category),
    by = c("gene" = "gene_symbol")
  )

# Build category -> group -> group_title lookup
group_lookup <- purrr::imap_dfr(groups, function(cats, grp) {
  data.frame(
    Category    = cats,
    Group       = grp,
    Group_Title = group_titles[[grp]],
    stringsAsFactors = FALSE
  )
})

# One row per gene x cell_type x gene_set
table_long <- df_sub_full %>%
  left_join(
    selected_genes %>% distinct(gene_symbol, gs_name, category),
    by = c("gene" = "gene_symbol", "category")
  ) %>%
  left_join(group_lookup, by = c("category" = "Category")) %>%  # fix case mismatch
  select(
    gene, gs_name, category, Group, Group_Title,
    cell_type, mean_exp, pct_exp
  ) %>%
  arrange(Group, category, gs_name, gene, cell_type)

wb <- loadWorkbook("msigdb_phagocytosis_search_all_hits.xlsx")
addWorksheet(wb, "gene_expression")
writeData(wb, "gene_expression", table_long)
saveWorkbook(wb, "msigdb_phagocytosis_search_all_hits.xlsx", overwrite = TRUE)

message("Saved: ", nrow(table_long), " rows, ", n_distinct(table_long$gene),
        " unique genes, ", n_distinct(table_long$gs_name), " gene sets")



#### COMPUTE AND SAVE PER-GENE-SET MODULE SCORES ####

gs_names <- unique(selected_hits$gs_name)

sink("gs_module_score_log.txt", split = TRUE)

min_genes <- 5
skipped_gs <- c()

for (gs in gs_names) {
  gs_genes <- selected_genes %>%
    filter(gs_name == gs) %>%
    pull(gene_symbol) %>%
    unique()
  
  gs_genes_present <- intersect(gs_genes, rownames(hDRG_sub))
  gs_genes_missing <- setdiff(gs_genes, rownames(hDRG_sub))
  
  if (length(gs_genes_present) < min_genes) {
    message("SKIPPED: ", gs, " — only ", length(gs_genes_present), "/", 
            length(gs_genes), " genes found")
    skipped_gs <- c(skipped_gs, gs)
    next
  }
  
  if (length(gs_genes_missing) > 0) {
    message("WARNING: ", gs, " — missing ", length(gs_genes_missing), "/", 
            length(gs_genes), " genes: ", 
            paste(gs_genes_missing, collapse = ", "))
  }
  
  gs_col <- make.names(gs)
  
  hDRG_sub <- AddModuleScore(
    hDRG_sub,
    features = list(gs_genes_present),
    name     = gs_col,
    seed     = 42
  )
  
  old_col <- paste0(gs_col, "1")
  hDRG_sub@meta.data[[gs_col]] <- hDRG_sub@meta.data[[old_col]]
  hDRG_sub@meta.data[[old_col]] <- NULL
  
  message("OK: ", gs, " — ", length(gs_genes_present), "/", 
          length(gs_genes), " genes used")
}

message("\n--- SUMMARY ---")
message("Completed: ", length(gs_names) - length(skipped_gs), " gene sets")
message("Skipped:   ", length(skipped_gs), " gene sets")
if (length(skipped_gs) > 0) {
  message("Skipped sets:\n", paste(" -", skipped_gs, collapse = "\n"))
}

sink()

saveRDS(hDRG_sub@meta.data, "hDRG_sub_meta_gs_scores.rds")



#### GENE SET DOTPLOT - 5 GROUPS, SELECTED CELLS ONLY ####

#Set thresh
z_threshold <- 1

# Only use gene sets that were successfully scored
gs_names_scored <- gs_names[purrr::map_lgl(gs_names, function(gs) {
  gs_col <- make.names(gs)
  !is.null(hDRG_sub@meta.data[[gs_col]])
})]
message("Using ", length(gs_names_scored), " / ", length(gs_names), " gene sets")

#Prep data
gs_dot_data <- purrr::map_dfr(gs_names_scored, function(gs) {
  gs_col     <- make.names(gs)
  score_vals <- hDRG_sub@meta.data[[gs_col]]
  z_vals     <- scale(score_vals)[, 1]
  
  data.frame(
    gs_name   = gs,
    cell_type = hDRG_sub@meta.data$Name_label,
    score     = score_vals,
    z_cell    = z_vals
  ) %>%
    group_by(gs_name, cell_type) %>%
    summarise(
      mean_score = mean(score),
      pct_above  = mean(z_cell > z_threshold) * 100,
      .groups    = "drop"
    )
}) %>%
  left_join(selected_hits %>% select(gs_name, category), by = "gs_name") %>%
  filter(!is.na(category))

# Z-score mean_score per gene set across cell types
gs_dot_data <- gs_dot_data %>%
  group_by(gs_name) %>%
  mutate(z_score = scale(mean_score)[, 1]) %>%
  ungroup()

# Order cell types
gs_dot_data <- gs_dot_data %>%
  mutate(cell_type = factor(cell_type, levels = cells_of_interest))

# Fix caps helper
fix_caps <- function(x) {
  x %>%
    gsub("^(Gobp|Gocc|Gomf|Kegg|Reactome|Wp|Hp|Hamai|Pid)\\b", "\\U\\1", ., perl = TRUE) %>%
    gsub("\\bRos\\b",   "ROS",   ., ignore.case = FALSE) %>%
    gsub("\\bRns\\b",   "RNS",   ., ignore.case = FALSE) %>%
    gsub("\\bTrail\\b", "TRAIL", ., ignore.case = FALSE) %>%
    gsub("\\bTp53\\b",  "TP53",  ., ignore.case = FALSE) %>%
    gsub("\\bVldl\\b",  "VLDL",  ., ignore.case = FALSE) %>%
    gsub("\\bLdl\\b",   "LDL",   ., ignore.case = FALSE) %>%
    gsub("\\bHdl\\b",   "HDL",   ., ignore.case = FALSE)
}

# Line height multiplier
lh <- 0.5

#Dotplot fxn
make_gs_dotplot <- function(df, title) {
  
  row_heights <- df %>%
    distinct(gs_name, category) %>%
    mutate(n_lines = str_count(as.character(gs_name), "\n") + 1) %>%
    group_by(category) %>%
    summarise(total_lines = sum(n_lines), .groups = "drop") %>%
    pull(total_lines)
  
  ggplot(df, aes(x = cell_type, y = gs_name)) +
    geom_point(aes(size = pct_above, color = z_score)) +
    scale_size_continuous(
      name   = paste0("% cells > ", z_threshold, " z"),
      range  = c(0.3, 5),
      breaks = c(10, 25, 50, 75, 100)
    ) +
    scale_color_gradientn(
      name    = "Z-score\n(mean per cell type)",
      colours = c("#FFF9C4", "#F7EE55", "#F6C141",
                  "#F1932D", "#E8601C", "#DC050C", "#72190E")
    ) +
    facet_grid(category ~ ., scales = "free_y", switch = "y") +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.x           = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y           = element_text(size = 6.5),
      strip.text.y.left     = element_text(angle = 0, hjust = 1, size = 7.5,
                                           face = "bold", lineheight = 1.2,
                                           margin = margin(t = 4, b = 4, r = 4, l = 4)),
      strip.placement       = "outside",
      strip.clip            = "off",
      strip.background      = element_rect(fill = "grey92", color = NA),
      strip.switch.pad.grid = unit(0, "cm"),
      panel.spacing         = unit(0.3, "lines"),
      legend.position       = "right",
      plot.title            = element_text(face = "bold", size = 11),
      plot.margin           = margin(t = 5, r = 5, b = 5, l = 0)
    ) +
    force_panelsizes(
      rows = unit(row_heights * lh, "cm"),
      cols = unit(0.55 * length(unique(df$cell_type)), "cm")
    )
}

#Run dotplots
dir.create("dotplots_genesets", showWarnings = FALSE)

for (grp in names(groups)) {
  cats  <- groups[[grp]]
  title <- group_titles[[grp]]
  
  df_grp <- gs_dot_data %>%
    filter(category %in% cats) %>%
    mutate(
      category = str_wrap(gsub("_", " ", as.character(category)) %>% str_to_title() %>% fix_caps(), width = 10) %>%
        factor(levels = unique(str_wrap(gsub("_", " ", cats) %>% str_to_title() %>% fix_caps(), width = 10)))
    )
  
  gs_order <- df_grp %>%
    distinct(gs_name, category) %>%
    arrange(category, gs_name) %>%
    pull(gs_name) %>%
    as.character() %>%
    unique()
  
  df_grp <- df_grp %>%
    mutate(
      gs_name = factor(
        str_wrap(gsub("_", " ", as.character(gs_name)) %>% str_to_title() %>% fix_caps(), width = 20),
        levels = str_wrap(gsub("_", " ", gs_order) %>% str_to_title() %>% fix_caps(), width = 20)
      )
    )
  
  p <- make_gs_dotplot(df_grp, title)
  
  total_lines <- df_grp %>%
    distinct(gs_name) %>%
    mutate(n_lines = str_count(as.character(gs_name), "\n") + 1) %>%
    pull(n_lines) %>%
    sum()
  
  n_cats <- length(unique(df_grp$category))
  
  ggsave(
    file.path("dotplots_genesets", paste0(grp, "_gs_dotplot.svg")),
    plot      = p,
    width     = length(unique(df_grp$cell_type)) * 0.03 + 10,
    height    = total_lines * lh + n_cats * 0.3 + 3
  )
  
  message("Saved ", grp, " — ", length(unique(df_grp$gs_name)), " gene sets")
}



#### DE OF SELECTED GENES IN MLC CELL TYPES ####

# MLC cell types
mlc_types <- c('MLC(H)', 'MLC(DA)', 'MLC(NFKB2)')

# MNP background
mnp_types <- c('DendriticCells', 'Mono(Early)', 'Mono(Diff.)',
               'MDM(GPNMB)', 'MDM(Early)', 'MDM(Mature)')

# Genes to test
genes_to_test <- unique(selected_genes$gene_symbol)
genes_present <- intersect(genes_to_test, rownames(hDRG_sub))

# Set identity to Name_label
Idents(hDRG_sub) <- "Name_label"

# Run DE for each MLC type vs MNP
de_results <- purrr::map_dfr(mlc_types, function(mlc) {
  
  message("Running DE for ", mlc, " vs MNP...")
  
  res <- FindMarkers(
    hDRG_sub,
    ident.1          = mlc,
    ident.2          = mnp_types,
    features         = genes_present,
    logfc.threshold  = 0,
    min.pct          = 0,
    test.use         = "wilcox",
    only.pos         = FALSE
  )
  
  res %>%
    rownames_to_column("gene") %>%
    mutate(MLC = mlc)
})

# Add gene set annotations
de_annotated <- de_results %>%
  left_join(
    selected_genes %>% distinct(gene_symbol, gs_name, category),
    by = c("gene" = "gene_symbol")
  ) %>%
  left_join(group_lookup, join_by("category" == "Category")) %>%
  rename(
    log2FC    = avg_log2FC,
    pct_MLC   = pct.1,
    pct_other = pct.2
  ) %>%
  select(
    gene, MLC, log2FC, p_val, p_val_adj,
    pct_MLC, pct_other,
    gs_name, category, Group, Group_Title
  ) %>%
  arrange(MLC, p_val_adj, desc(log2FC))

# Significant hits
de_sig <- de_annotated %>%
  filter(
    p_val_adj < 0.001,
    log2FC    > 1,
    pct_MLC   >= 0.25,
    pct_other <= 0.5
  )

# Print summary of all genes p<0.05 (590 genes)
de_results %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  distinct(gene) %>%
  nrow() %>%
  message("Unique genes with LFC > 0 and padj < 0.05: ", .)

# Print summary of all highly enriched genes (50 genes)
message("Significant hits: ", nrow(de_sig),
        " across ", n_distinct(de_sig$gene), " unique genes")
message("Per MLC type:")
de_sig %>% count(MLC) %>% print()

# Save to one xlsx
write.xlsx(
  list(
    all_genes     = de_results,
    all_annotated = de_annotated,
    sig_annotated = de_sig
  ),
  "MLC_DE_selected_genes.xlsx",
  row.names = FALSE
)



#### DOTPLOT - MLC DE SIGNIFICANT HITS ####

# Get sig hit genes in alpha order
sig_genes <- de_sig %>%
  distinct(gene) %>%
  arrange(gene) %>%
  pull(gene)

# Pull expression data for these genes from df_sub_full, all cell types
df_mlc_dot <- df_sub_full %>%
  filter(gene %in% sig_genes) %>%
  mutate(
    gene      = factor(gene, levels = rev(sig_genes)),
    cell_type = factor(cell_type, levels = cells_of_interest)
  )

# Gene -> category annotation (before cat_order so we can filter)
gene_cat_annot <- df_sub_full %>%
  filter(gene %in% sig_genes) %>%
  distinct(gene, category) %>%
  left_join(group_lookup, join_by("category" == "Category")) %>%
  mutate(gene = factor(gene, levels = rev(sig_genes)))

# Ordered categories (only those present in sig genes)
cat_order <- group_lookup %>%
  arrange(Group, Category) %>%
  pull(Category) %>%
  unique() %>%
  .[. %in% unique(gene_cat_annot$category)]

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

# Ordered colors (only used categories)
cat_colors <- setNames(
  colorRampPalette(cbfColors)(length(cat_order)),
  cat_order
)

# Annotation panel (left)
n_genes <- length(sig_genes)
n_cats  <- length(cat_order)
p_annot <- ggplot(gene_cat_annot, aes(x = category, y = gene, fill = category)) +
  geom_tile(width = 1, height = 1, color = NA) +
  geom_vline(xintercept = seq(0.5, n_cats + 0.5, by = 1), color = "grey80", linewidth = 0.2) +
  geom_hline(yintercept = seq(0.5, n_genes + 0.5, by = 1), color = "grey80", linewidth = 0.2) +
  scale_x_discrete(limits = cat_order, expand = c(0, 0)) +
  scale_y_discrete(limits = levels(gene_cat_annot$gene), expand = c(0, 0)) +
  scale_fill_manual(values = cat_colors, guide = "none", na.value = "white") +
  labs(x = NULL, y = NULL, title = "MLC Enriched Phagocytosis Genes") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y  = element_text(size = 6.5, face = "italic"),
    axis.ticks   = element_line(color = "black", linewidth = 0.3),
    axis.line    = element_line(color = "black", linewidth = 0.4),
    panel.ontop  = FALSE,
    plot.margin  = margin(t = 5, r = 0, b = 5, l = 5)
  )

# Main dotplot
p_main <- ggplot(df_mlc_dot, aes(x = cell_type, y = gene)) +
  geom_vline(xintercept = c(6.5, 9.5, 11.5), color = "grey60", linewidth = 0.4) +
  geom_point(aes(size = pct_exp, color = mean_exp)) +
  scale_size_continuous(
    name   = "% Expressed",
    range  = c(0.3, 5),
    breaks = c(10, 25, 50, 75, 100)
  ) +
  scale_color_gradientn(
    name    = "Mean Expr.",
    colours = c("#FFF9C4", "#F7EE55", "#F6C141",
                "#F1932D", "#E8601C", "#DC050C", "#72190E")
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y     = element_text(size = 6.5, face = "italic"),
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 11),
    plot.margin     = margin(t = 5, r = 5, b = 5, l = 0)
  )

# Combine
p_combined <- p_annot + p_main +
  plot_layout(widths = unit(c(length(cat_order) * 0.4, length(cells_of_interest) * 0.55), "cm"))

ggsave(
  "dotplot_MLC_DE_sig.svg",
  plot      = p_combined,
  width     = length(cells_of_interest) * 0.055 + length(cat_order) * 0.4 + 4,
  height    = length(sig_genes) * 0.1 + 3
)

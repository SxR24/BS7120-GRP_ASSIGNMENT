# =============================================================================
# Figure S1B
# =============================================================================

library(readr)
library(dplyr)
library(tibble)
library(FactoMineR)
library(gplots)
library(RColorBrewer)
library(dynamicTreeCut)


# Inputs
expr_matrix <- read_tsv("2026_UoL_prj_gene_matrix_f.tsv")
zone_meta   <- read_tsv("cellFAKE_assignments.txt")


# Prepare expression matrix
expr_df <- as.data.frame(expr_matrix)
rownames(expr_df) <- expr_df[, 1]          
expr_df <- expr_df[, -1]                   

expr_mat <- as.matrix(expr_df)
mode(expr_mat) <- "numeric"

cat("Total cells before filtering:", ncol(expr_mat), "\n")


# Filter samples using zone assignment file
samples_keep <- intersect(colnames(expr_mat), zone_meta$Sample)
expr_mat <- expr_mat[, samples_keep, drop = FALSE]

meta_sel <- zone_meta %>%
  filter(Sample %in% samples_keep) %>%
  select(Sample, Assigned_Zone, Assigned_CellType) %>%
  distinct() %>%
  column_to_rownames("Sample")

meta_sel <- meta_sel[colnames(expr_mat), , drop = FALSE]

cat("Samples retained after zone filtering:", ncol(expr_mat), "\n")
cat("Zones present:", paste(unique(meta_sel$Assigned_Zone), collapse = ", "), "\n")
cat("Cell types present:", paste(unique(meta_sel$Assigned_CellType), collapse = ", "), "\n")

# Define variable genes on filtered matrix
gene_vars <- apply(expr_mat, 1, var)
gene_expr_cells <- rowSums(expr_mat > 1)

variable_genes <- rownames(expr_mat)[
  gene_vars > 0.5 & gene_expr_cells >= 2
]

cat("Variable genes selected (all cells):", length(variable_genes), "\n")
expr_var <- expr_mat[variable_genes, ]

# PCA on ALL filtered samples
pca_res_all <- PCA(t(expr_var), scale.unit = TRUE, graph = FALSE)
loadings_all <- pca_res_all$var$coord

# Helper function for gene extraction
n_max <- 50   

extract_genes <- function(load_vec, threshold, decreasing, n_max) {
  idx <- if (decreasing) load_vec > threshold else load_vec < -threshold
  genes <- names(sort(load_vec[idx], decreasing = decreasing))
  if (length(genes) > n_max) genes <- genes[1:n_max]
  return(genes)
}

# Extract gene groups from all-cells PCA
groupA <- extract_genes(loadings_all[,1], 0.2, TRUE, n_max)   # PC1 pos
groupB <- extract_genes(loadings_all[,1], 0.2, FALSE, n_max)  # PC1 neg
groupC <- extract_genes(loadings_all[,2], 0.2, TRUE, n_max)   # PC2 pos
groupD <- extract_genes(loadings_all[,2], 0.2, FALSE, n_max)  # PC2 neg
groupE <- extract_genes(loadings_all[,3], 0.2, FALSE, n_max)  # PC3 neg
pc4_anti_all <- extract_genes(loadings_all[,4], 0.2, FALSE, n_max)  # PC4 neg

cat("\nAll-cells PCA genes:\n")
cat("  Group A (PC1 cor):", length(groupA), "\n")
cat("  Group B (PC1 anti):", length(groupB), "\n")
cat("  Group C (PC2 cor):", length(groupC), "\n")
cat("  Group D (PC2 anti):", length(groupD), "\n")
cat("  Group E (PC3 anti):", length(groupE), "\n")
cat("  PC4 anti (for Group F):", length(pc4_anti_all), "\n")

# Read NPC sample list and PCA on NPCs only
npc_file <- "Sample_lists_Neurons_vs_NPCs.tsv"
npc_df <- read_tsv(npc_file, col_names = TRUE)

npc_samples <- npc_df[[2]]
npc_samples <- npc_samples[!is.na(npc_samples)]
npc_samples <- intersect(npc_samples, colnames(expr_mat))

cat("\nNPC samples found:", length(npc_samples), "\n")

if (length(npc_samples) > 5) {
  expr_npc_full <- expr_mat[, npc_samples, drop = FALSE]
  
  npc_gene_vars <- apply(expr_npc_full, 1, var)
  npc_gene_expr_cells <- rowSums(expr_npc_full > 1)
  variable_genes_npc <- rownames(expr_npc_full)[
    npc_gene_vars > 0.5 & npc_gene_expr_cells >= 2
  ]
  cat("Variable genes in NPCs only:", length(variable_genes_npc), "\n")
  
  expr_npc_var <- expr_npc_full[variable_genes_npc, ]
  
  pca_res_npc <- PCA(t(expr_npc_var), scale.unit = TRUE, graph = FALSE)
  loadings_npc <- pca_res_npc$var$coord
  
  # NPC-derived groups
  pc1_cor_npc   <- extract_genes(loadings_npc[,1], 0.2, TRUE, n_max)   # PC1 cor NPC
  groupE        <- extract_genes(loadings_npc[,2], 0.2, TRUE, n_max)   # PC2 cor NPC
  groupG        <- extract_genes(loadings_npc[,2], 0.2, FALSE, n_max)  # PC2 anti NPC (Group G)
  
  # Combined Group F (PC1 cor NPC + PC4 anti all)
  groupF_combo <- unique(c(pc1_cor_npc, pc4_anti_all))
  
  cat("\nNPC-only PCA genes:\n")
  cat("  PC1 cor NPC:", length(pc1_cor_npc), "\n")
  cat("  PC2 cor NPC:", length(pc2_cor_npc), "\n")
  cat("  Group G (PC2 anti NPC):", length(groupG), "\n")
  cat("  Group F combo (PC1 cor NPC + PC4 anti all):", length(groupF_combo), "\n")
} else {
  stop("Not enough NPC samples for PCA. Found only ", length(npc_samples))
}

# Define final gene sections
gene_sections <- list(
  "Group A (PC1 cor all)"                = groupA,
  "Group B (PC1 anti all)"               = groupB,
  "Group C (PC2 cor all)"                = groupC,
  "Group D (PC2 anti all)"               = groupD,
  "Group E combo (PC3 anti all + PC2 cor NPC)" = groupE,
  "Group F combo (PC1 cor NPC + PC4 anti all)" = groupF_combo,
  "Group G (PC2 anti NPC)"               = groupG
)

# Within‑section hierarchical clustering of genes
cluster_section <- function(gene_list, expr_sub) {
  if (length(gene_list) <= 1) return(gene_list)
  gene_cor <- cor(t(expr_sub[gene_list, , drop = FALSE]), method = "pearson")
  gene_dist <- as.dist(1 - gene_cor)
  hc_genes <- hclust(gene_dist, method = "average")
  return(gene_list[hc_genes$order])
}

ordered_sections <- lapply(gene_sections, function(sec) {
  if (length(sec) > 0) cluster_section(sec, expr_mat) else character(0)
})

heat_genes_ordered <- unlist(ordered_sections)
cat("\nTotal heatmap rows (genes):", length(heat_genes_ordered), "\n")

# Write gene lists to CSV
 pc_lists <- list(
   "PC1 cor"      = groupA,
   "PC1 anti"     = groupB,
   "PC2 cor"      = groupC,
   "PC2 anti"     = groupD,
   "PC3 anti"     = groupE,
   "PC4 anti"     = pc4_anti_all,
   "PC1 rg cor"   = pc1_cor_npc,
   "PC2 rg cor"   = pc2_cor_npc,
   "PC2 rg anti"  = groupG
 )
 
 max_len <- max(sapply(pc_lists, length))
 pc_lists_padded <- lapply(pc_lists, function(x) {
   if (length(x) < max_len) c(x, rep(NA, max_len - length(x))) else x
 })
 pc_df <- as.data.frame(pc_lists_padded, stringsAsFactors = FALSE, check.names = FALSE)
 write.csv(pc_df, "D:/GroupA_University2026_Project/Figure_1C_gene_lists_by_PC.csv",
           row.names = FALSE, na = "")


# Prepare colour bars for heatmap
zone_colors <- c(
  "VZ"  = "#FFBF20",
  "iSVZ"= "#F77800",
  "oSVZ"= "#CD2624",
  "CP"  = "#7D007E"
)
all_zones <- unique(meta_sel$Assigned_Zone)
zone_colors <- zone_colors[names(zone_colors) %in% all_zones]
col_side <- zone_colors[meta_sel$Assigned_Zone]

section_names <- names(ordered_sections)
section_lengths <- sapply(ordered_sections, length)
row_section <- rep(section_names, times = section_lengths)

section_palette <- RColorBrewer::brewer.pal(min(8, length(section_names)), "Set1")
if (length(section_names) > 8) {
  section_palette <- colorRampPalette(section_palette)(length(section_names))
}
names(section_palette) <- section_names
row_side <- section_palette[row_section]

row_sep <- cumsum(section_lengths)
row_sep <- row_sep[-length(row_sep)]

# Hierarchical clustering of samples
expr_heat <- expr_mat[heat_genes_ordered, ]
cor_mat <- cor(expr_heat, method = "pearson")
cor_dist <- as.dist(1 - cor_mat)
hc_col <- hclust(cor_dist, method = "average")

# Dendrogram with dynamic clusters
dynamic_clusters <- cutreeDynamic(
  dendro = hc_col,
  distM = as.matrix(cor_dist),
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 5
)
names(dynamic_clusters) <- hc_col$labels
print(table(dynamic_clusters))

cluster_ids_unique <- sort(unique(dynamic_clusters))
n_clust <- length(cluster_ids_unique)
clust_cols <- RColorBrewer::brewer.pal(min(n_clust, 8), "Set1")
if (n_clust > 8) clust_cols <- rainbow(n_clust)
names(clust_cols) <- cluster_ids_unique

leaf_order <- hc_col$order
leaf_clusters <- dynamic_clusters[hc_col$labels[leaf_order]]
leaf_cols <- clust_cols[as.character(leaf_clusters)]

plot(hc_col, hang = -1, main = "Dendrogram with dynamic clusters (leaf colours)",
     labels = FALSE)
axis(1, at = 1:length(leaf_order), labels = FALSE)
text(x = 1:length(leaf_order),
     y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
     labels = hc_col$labels[leaf_order],
     col = leaf_cols,
     srt = 90, adj = 1, cex = 0.6, xpd = NA)
legend("topright",
       legend = paste("Cluster", names(clust_cols)),
       fill = clust_cols,
       border = NA,
       bty = "n",
       cex = 0.8,
       title = "Dynamic Clusters")

# Get column order from dendrogram
col_order <- hc_col$order
samples_ordered <- colnames(expr_heat)[col_order]

# Define sample groups (boundary SRR IDs)
group_definition <- list(
  "AP2" = c("SRR2967715", "SRR2967737"),
  "AP1" = c("SRR2967677", "SRR2967689"),
  "BP1" = c("SRR2967710", "SRR2967820"),
  "BP2" = c("SRR2967760", "SRR2967775"),
  "N1"  = c("SRR2967774", "SRR2967819"),
  "N2"  = c("SRR2967816", "SRR2967768"),
  "N3"  = c("SRR2967763", "SRR2967748")
)

# Find contiguous indices for each group in the clustered order
group_ranges <- list()
for (grp in names(group_definition)) {
  ids <- group_definition[[grp]]
  pos1 <- which(samples_ordered == ids[1])
  pos2 <- which(samples_ordered == ids[2])
  if (length(pos1) == 0 || length(pos2) == 0) next
  start_idx <- min(pos1, pos2)
  end_idx   <- max(pos1, pos2)
  group_ranges[[grp]] <- c(start_idx, end_idx)
  cat("Group", grp, ": indices", start_idx, "to", end_idx,
      "(", samples_ordered[start_idx], "to", samples_ordered[end_idx], ")\n")
}

# Assign custom colours
group_cols <- c(
  "AP2" = "#95F297",
  "AP1" = "#2C8334",
  "BP1" = "#15B5B5",
  "BP2" = "#9DFFFF",
  "N1"  = "#8EC9F5",
  "N2"  = "#4589EF",
  "N3"  = "#0545BB"
)
group_cols <- group_cols[names(group_ranges)]

# Define drawing function
draw_manual_bar <- function() {
  n_cols <- length(samples_ordered)
  usr <- par("usr")
  y_bottom <- usr[3] - 0.06 * (usr[4] - usr[3])
  y_top    <- usr[3] - 0.02 * (usr[4] - usr[3])
  
  for (grp in names(group_ranges)) {
    range_idx <- group_ranges[[grp]]
    rect(xleft   = range_idx[1] - 0.5,
         xright  = range_idx[2] + 0.5,
         ybottom = y_bottom,
         ytop    = y_top,
         col     = group_cols[grp],
         border  = "black",
         xpd     = NA)
  }
}

# Generate heatmap
hm_colors <- colorRampPalette(c("darkblue", "white", "red"))(100)

heatmap.2(expr_heat,
          Colv = as.dendrogram(hc_col),
          Rowv = FALSE,
          dendrogram = "column",
          col = hm_colors,
          trace = "none",
          density.info = "none",
          key = TRUE,
          keysize = 1.0,
          key.title = "log2(FPKM)",
          margins = c(12, 8),               
          labRow = "",
          labCol = "",
          ColSideColors = col_side,
          RowSideColors = row_side,
          rowsep = row_sep,
          sepcolor = "black",
          sepwidth = c(0.5, 0.5),
          main = "Figure S1B: Gene Expression Patterns in Human Fetal Neocortex",
          add.expr = draw_manual_bar() )


# Add legend for manual groups
legend("bottomleft",
       legend = names(group_cols),
       fill = group_cols,
       border = "black",
       bty = "n",
       cex = 0.8,
       title = "Sample Groups")

# Add text annotations inside row side color rectangles
n_rows <- length(heat_genes_ordered)
section_starts <- c(1, cumsum(section_lengths)[-length(section_lengths)] + 1)
section_ends <- cumsum(section_lengths)
section_midpoints <- (section_starts + section_ends) / 2

x_pos <- 0.04
y_positions <- c(0.8, 0.65, 0.5, 0.4, 0.3, 0.1, -0.08)

labels <- c(
  "Group A: GO: Cell cycle,\n  Forebrain Dev.",
  "Group B: GO: Neuron Diff.\n  and projection",
  "Group C: GO: Cell Signaling,\n  Vesicle transport",
  "Group D: GO: Neurogenesis,\n  Cell proliferation",
  "Group E: GO: Cell adhesion",
  "Group F: GO: Cell cycle,\n  mitosis",
  "Group G: GO: Neurogenesis,\n  Cell morph."
)

for (i in seq_along(y_positions)) {
  text(x = x_pos, y = y_positions[i],
       labels = labels[i],
       srt = 0, xpd = NA, col = "black", cex = 0.6, font = 2)
}

# Add legend for cortical zones
legend("topright",
       legend = names(zone_colors),
       fill = zone_colors,
       border = FALSE,
       bty = "n",
       cex = 0.8,
       title = "Cortical Zones")
# Second marker gene heatmap
target_genes <- c("PAX6", "GLI3", "SOX2", "VIM", "PROM1", "HES1",
                  "NEUROG2", "EOMES", "NEUROD4", "INSM1", "ASPM",
                  "CCNB1", "MKI67", "HES6", "NEUROD6", "BCL11B",
                  "MYT1L", "TBR1", "MEF2C")

genes_present <- intersect(target_genes, rownames(expr_mat))
genes_for_panel <- target_genes[target_genes %in% rownames(expr_mat)]

expr_panel <- expr_mat[genes_for_panel, hc_col$order, drop = FALSE]
col_side_panel <- col_side[hc_col$order]

draw_manual_bar_lines <- function(offset_lines = 3) {
  line_height <- strheight("M", units = "user")
  y_bottom <- par("usr")[3] - offset_lines * line_height
  y_top    <- par("usr")[3] - (offset_lines - 1) * line_height
  
  for (grp in names(group_ranges)) {
    range_idx <- group_ranges[[grp]]
    rect(xleft   = range_idx[1] - 0.5,
         xright  = range_idx[2] + 0.5,
         ybottom = y_bottom,
         ytop    = y_top,
         col     = group_cols[grp],
         border  = "black",
         xpd     = NA)
  }
}

heatmap.2(expr_panel,
          Colv = FALSE,
          Rowv = FALSE,
          dendrogram = "none",
          col = hm_colors,
          trace = "none",
          density.info = "none",
          key = TRUE,
          keysize = 1.0,
          key.title = "log2(FPKM)",
          margins = c(10, 8),            
          labRow = genes_for_panel,
          labCol = "",
          ColSideColors = col_side_panel,
          main = "Selected Marker Genes (same sample order as Fig S1B)",
          cexRow = 0.9,
          add.expr = draw_manual_bar_lines(offset_lines = 3))

# Assign samples to groups using dynamic cluster IDs and manual overrides
cluster_to_group <- c(
  "1" = "N3",
  "3" = "N2",
  "2" = "N1",
  "5" = "BP2",
  "4" = "BP1",
  "7" = "AP1",
  "6" = "AP2"
)

sample_clusters <- dynamic_clusters[colnames(expr_mat)]  # ensure order matches all samples
sample_groups <- cluster_to_group[as.character(sample_clusters)]

assignment_df <- data.frame(
  Sample = names(sample_groups),
  Group  = sample_groups,
  stringsAsFactors = FALSE
)

# 4. Reassign the four specific samples from N2 to N3
reassign_samples <- c("SRR2967763", "SRR2967681", "SRR2967748", "SRR2967680")
for (s in reassign_samples) {
  if (s %in% assignment_df$Sample) {
    old_group <- assignment_df$Group[assignment_df$Sample == s]
    assignment_df$Group[assignment_df$Sample == s] <- "N3"
    cat(sprintf("Reassigned sample %s from %s to N3\n", s, old_group))
  } else {
    warning(sprintf("Sample %s not found in expression matrix", s))
  }
}

# 5. Write TSV
output_tsv <- "D:/GroupA_University2026_Project/sample_group_assignments.tsv"
write_tsv(assignment_df, output_tsv)

cat("Sample-to-group assignment saved to:", output_tsv, "\n")
cat("Number of samples per group after reassignment:\n")
print(table(assignment_df$Group))

# Create sample-to-group assignment table

sample_groups <- rep(NA, length(samples_ordered))

for (grp in names(group_ranges)) {
  range_idx <- group_ranges[[grp]]
  sample_groups[range_idx[1]:range_idx[2]] <- grp
}

# Build a data frame
assignment_df <- data.frame(
  sample_id = samples_ordered,
  group = sample_groups,
  stringsAsFactors = FALSE
)

# Check for any unassigned samples
if (any(is.na(assignment_df$group))) {
  warning("Some samples were not assigned to any group.")
  print(assignment_df[is.na(assignment_df$group), ])
}

# Write to CSV
write.csv(assignment_df, "D:/GroupA_University2026_Project/OUTPUTS/Altsubpopulations.csv", row.names = FALSE, quote = FALSE)

cat("Sample group assignments saved to Altsubpopulations.csv\n")
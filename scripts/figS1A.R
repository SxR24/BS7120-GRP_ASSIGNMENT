dim(sc_fetal)
head(rownames(sc_fetal))
summary(as.numeric(as.matrix(sc_fetal[1:5,1:5])))
"F5_fetal_12wpc_c1" %in% rownames(sc_fetal)
library(FactoMineR)

expr_mat <- t(as.matrix(sc_fetal))
# Condition 1: expressed in > 2 cells
expressed_cells <- rowSums(expr_mat > 0)
keep_expr <- expressed_cells > 2

# Condition 2: variance > 0.5
gene_var <- apply(expr_mat, 1, var, na.rm = TRUE)
keep_var <- gene_var > 0.5
keep_var[is.na(keep_var)] <- FALSE

# Combine both conditions
variable_genes <- keep_expr & keep_var

# Create filtered matrix
expr_var <- expr_mat[variable_genes, ]

# Check number of genes retained
sum(variable_genes)
expr_for_pca <- t(expr_var)

library(FactoMineR)

pca_res <- PCA(expr_for_pca,
               scale.unit = TRUE,
               graph = FALSE)
pc_df <- as.data.frame(pca_res$ind$coord[, 1:2])
colnames(pc_df) <- c("PC1", "PC2")
pc_df$cell_id <- rownames(pc_df)
pc_df$experiment_group <- substr(pc_df$cell_id, 1, 1)
library(ggplot2)

p <- ggplot(pc_df, aes(x = PC1, y = PC2, color = experiment_group)) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "A" = "#66c2a5",
    "B" = "#41ae76",
    "C" = "#238b45",
    "D" = "#006d2c"
  )) +
  theme_classic() +
  labs(title = "PCA of Fetal Cerebral Cortex (12–13 wpc)",
       x = "PC1",
       y = "PC2")

print(p)
dim(pc_df)

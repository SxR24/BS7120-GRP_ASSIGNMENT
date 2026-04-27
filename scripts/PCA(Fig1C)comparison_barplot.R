# ===========================================================================
# PCA Gene-List Overlap Analysis & Plot
# ---------------------------------------------------------------------------
# This script compares gene lists derived from Principal Component PC
# across all pipelines (original, non-coding, protein-coding) and 
# quantifies the percentage overlap with the reference original dataset.
# ===========================================================================

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Inputs
test_files <- c(
  Original       = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/Oldpipeline_Figure_1C_gene_lists_by_PC.csv",
  Noncoding = "D:/GroupA_University2026_Project/AltnoncodingPCA.csv",
  Protein   = "D:/GroupA_University2026_Project/Figure_1C_gene_lists_by_PC_fetal.csv"
)

original_df <- read.csv("OriginalPCA.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Function to compute gene overlap for a single test file against the reference.
compute_overlap <- function(test_path, test_name, original_df) {
  test_df <- read.csv(test_path, stringsAsFactors = FALSE, check.names = FALSE)
  common_cols <- intersect(names(original_df), names(test_df))
  
  if (length(common_cols) == 0) {
    warning(paste("No common columns for", test_name))
    return(data.frame(
      TestDataset = character(), PC = character(),
      Original_Count = integer(), Test_Count = integer(),
      Overlap_Count = integer(), Percent_Overlap = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  results <- data.frame(TestDataset = character(), PC = character(),
                        Original_Count = integer(), Test_Count = integer(),
                        Overlap_Count = integer(), Percent_Overlap = numeric(),
                        stringsAsFactors = FALSE)
  
  for (col in common_cols) {
    orig_genes <- na.omit(original_df[[col]])
    orig_genes <- orig_genes[orig_genes != ""]
    test_genes <- na.omit(test_df[[col]])
    test_genes <- test_genes[test_genes != ""]
    overlap <- intersect(orig_genes, test_genes)
    pct <- if (length(test_genes) > 0) length(overlap) / length(test_genes) * 100 else NA
    
    results <- rbind(results, data.frame(
      TestDataset = test_name, PC = col,
      Original_Count = length(orig_genes), Test_Count = length(test_genes),
      Overlap_Count = length(overlap), Percent_Overlap = pct,
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}

# Apply the overlap function and combine results into one table
all_results <- bind_rows(lapply(names(test_files), function(nm) {
  compute_overlap(test_files[nm], nm, original_df)
}))
all_results <- all_results[!is.na(all_results$Percent_Overlap), ]

write.csv(all_results, "PCA_overlap_results_all.csv", row.names = FALSE)

# Optional: list genes that appear in the "Original" pipeline file but NOT in the reference original file,
# written out by PC column to a text file (not included in analysis hence optional).
non_overlap_file <- "non_overlap_genes.txt"
file_conn <- file(non_overlap_file, "w")
writeLines("Genes present in Original file but not in original pipeline file, by PC column", file_conn)
writeLines("==========================================================================\n", file_conn)

test_Original <- read.csv(test_files["Original"], stringsAsFactors = FALSE, check.names = FALSE)
common_cols <- intersect(names(original_df), names(test_Original))
for (col in common_cols) {
  orig_genes <- na.omit(original_df[[col]]); orig_genes <- orig_genes[orig_genes != ""]
  Original_genes  <- na.omit(test_Original[[col]]); Original_genes <- Original_genes[Original_genes != ""]
  Original_only <- setdiff(Original_genes, orig_genes)
  
  writeLines(paste0("PC column: ", col), file_conn)
  writeLines(paste0("Number of non‑overlapping genes: ", length(Original_only)), file_conn)
  if (length(Original_only) > 0) {
    writeLines("Genes:", file_conn)
    for (gene in sort(Original_only)) writeLines(paste0("  ", gene), file_conn)
  } else {
    writeLines("  (None – all genes in this column overlap)", file_conn)
  }
  writeLines("--------------------------------------------\n", file_conn)
}
close(file_conn)

# Tidy up PC labels for the plot (replace abbreviations with readable names),
label_map <- c("cor" = "correlating", "anti" = "anticorrelating", "rg" = "NPC")
replace_shorthand <- function(x, map) {
  for (pattern in names(map)) x <- gsub(pattern, map[[pattern]], x, fixed = TRUE)
  return(x)
}
all_results$PC_label <- replace_shorthand(all_results$PC, label_map)

pc_order_original <- names(original_df)
pc_order_original <- pc_order_original[pc_order_original %in% unique(all_results$PC)]
all_results$PC <- factor(all_results$PC, levels = pc_order_original)
all_results$PC_label <- factor(all_results$PC_label, 
                               levels = replace_shorthand(pc_order_original, label_map))

# Create a dodged bar chart showing percentage overlap per PC for each test dataset.
p <- ggplot(all_results, aes(x = PC_label, y = Percent_Overlap, fill = TestDataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percent_Overlap)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Original" = "#2CA02C", "Noncoding" = "#E69F00", "Protein" = "#56B4E9")) +
  labs(title = "Overlap of PC Gene Lists with Original",
       x = "Principal Component",
       y = "Percentage overlap (%)",
       fill = "Test dataset") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

print(p)
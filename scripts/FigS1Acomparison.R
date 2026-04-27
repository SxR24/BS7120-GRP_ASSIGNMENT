# ----------------------------------------------------------------------
# Bar plot of proportion overlap of cell types
# ----------------------------------------------------------------------

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# Inputs
ref_files <- list(
  sample    = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/Old_Sample_lists_Neurons_vs_NPCs.tsv",
  assign    = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/Old_cell_assignments.txt",
  subpop    = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/oldsubpopulations.csv"
)

noncoding_files <- list(
  sample    = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/noncoding_Sample_lists_Neurons_vs_NPCs.tsv",
  assign    = "altpipeline_noncoding_cell_assignments.txt",
  subpop    = "Altprotiensubpopulations.csv"
)

protein_files <- list(
  sample    = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/protein_Sample_lists_Neurons_vs_NPCs.tsv",
  assign    = "altpipeline_protien_cell_assignments.txt",
  subpop    = "Altnoncodingsubpopulations.csv"
)

# Extract unique sample IDs per category from a file
get_sample_categories <- function(file_path, type = "sample") {
  if (type == "sample") {
    df <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
    neurons <- unique(na.omit(df$Neurons))
    npcs    <- unique(na.omit(df$NPCs))
    return(list(Neurons = neurons, NPCs = npcs))
  }
  else if (type == "assign") {
    df <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
    zones <- split(df$Sample, df$Assigned_Zone)   # list of samples per zone
    cells <- split(df$Sample, df$Assigned_CellType)
    return(c(zones, cells))
  }
  else if (type == "subpop") {
    df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
    groups <- split(df$sample_id, df$group)
    return(groups)
  }
}

# Compute overlap for orignal dataset against reference
calc_overlap_for_categories <- function(test_categories, ref_categories) {
  categories <- intersect(names(test_categories), names(ref_categories))
  result <- data.frame()
  for (cat in categories) {
    test_set <- test_categories[[cat]]
    ref_set  <- ref_categories[[cat]]
    if (length(test_set) == 0) next
    overlap <- length(intersect(test_set, ref_set)) / length(test_set)
    result <- rbind(result, data.frame(Category = cat, Overlap = overlap))
  }
  return(result)
}

# Load reference categories
ref_sample <- get_sample_categories(ref_files$sample, "sample")
ref_assign <- get_sample_categories(ref_files$assign, "assign")
ref_subpop <- get_sample_categories(ref_files$subpop, "subpop")
ref_all <- c(ref_sample, ref_assign, ref_subpop)

# Process noncoding and protein datasets
process_dataset <- function(file_list, dataset_name) {
  result <- data.frame()
  if (!is.null(file_list$sample)) {
    test <- get_sample_categories(file_list$sample, "sample")
    ov <- calc_overlap_for_categories(test, ref_all)
    ov$BarType <- "Neurons/NPCs"
    result <- rbind(result, ov)
  }
  if (!is.null(file_list$assign)) {
    test <- get_sample_categories(file_list$assign, "assign")
    ov <- calc_overlap_for_categories(test, ref_all)
    zone_cats <- c("VZ", "iSVZ", "oSVZ", "CP")
    ov$BarType <- ifelse(ov$Category %in% zone_cats, "Assigned Zone", "Assigned Cell Type")
    result <- rbind(result, ov)
  }
  if (!is.null(file_list$subpop)) {
    test <- get_sample_categories(file_list$subpop, "subpop")
    ov <- calc_overlap_for_categories(test, ref_all)
    ov$BarType <- "Subpopulations"
    result <- rbind(result, ov)
  }
  result$Dataset <- dataset_name
  return(result)
}

noncoding_overlap <- process_dataset(noncoding_files, "Noncoding")
protein_overlap  <- process_dataset(protein_files, "Protein")

# Combine
plot_data <- bind_rows(noncoding_overlap, protein_overlap)

# Prepare for plotting: order categories as in original script
category_order <- c(
  "NPCs", "Neurons",
  "VZ", "iSVZ", "oSVZ", "CP",
  "aRG", "bRG", "neuron",
  "AP2", "AP1", "BP1", "BP2", "N1", "N2", "N3"
)
plot_data$Category <- factor(plot_data$Category, levels = category_order)

plot_data$BarType <- factor(plot_data$BarType,
                            levels = c("Neurons/NPCs", "Assigned Zone", "Assigned Cell Type", "Subpopulations"))

# Plot
ggplot(plot_data, aes(x = Category, y = Overlap, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(aes(label = scales::percent(Overlap, accuracy = 0.1)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Noncoding" = "#E69F00", "Protein" = "#56B4E9")) +
  facet_grid(~ BarType, scales = "free_x", space = "free_x") +
  labs(
    title = "Overlap of each population with the original pipeline",
    y = "Proportion overlapping with original pipeline",
    x = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )
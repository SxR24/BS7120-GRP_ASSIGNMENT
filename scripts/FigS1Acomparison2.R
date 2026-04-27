# ----------------------------------------------------------------------
# Comparison bar plot of proportions within all three pipelines
# ----------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# Inputs
sample_files <- c(
  Noncoding = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/Old_Sample_lists_Neurons_vs_NPCs.tsv",
  Protein   = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/protein_Sample_lists_Neurons_vs_NPCs.tsv",
  Original = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/noncoding_Sample_lists_Neurons_vs_NPCs.tsv"
)

calc_sample_props <- function(file_path, set_name) {
  df <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
  n_neurons <- length(unique(na.omit(df$Neurons)))
  n_npcs    <- length(unique(na.omit(df$NPCs)))
  total <- n_neurons + n_npcs
  data.frame(
    Set    = set_name,
    Bar    = "Neurons/NPCs",
    Type   = c("Neurons", "NPCs"),
    Proportion = c(n_neurons / total, n_npcs / total)
  )
}

sample_data <- map_dfr(names(sample_files), ~ calc_sample_props(sample_files[.x], .x))

assign_files <- c(
  Original = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/Old_cell_assignments.txt",
  Noncoding = "altpipeline_noncoding_cell_assignments.txt",
  Protein   = "altpipeline_protien_cell_assignments.txt"
)

calc_assign_props <- function(file_path, set_name, column_name, bar_label) {
  df <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
  counts <- table(df[[column_name]])
  total <- sum(counts)
  data.frame(
    Set    = set_name,
    Bar    = bar_label,
    Type   = names(counts),
    Proportion = as.numeric(counts) / total
  )
}

assign_data <- bind_rows(
  map_dfr(names(assign_files), function(set) {
    calc_assign_props(assign_files[set], set, "Assigned_Zone", "Assigned Zone")
  }),
  map_dfr(names(assign_files), function(set) {
    calc_assign_props(assign_files[set], set, "Assigned_CellType", "Assigned CellType")
  })
)

subpop_files <- c(
  Original = "D:/GroupA_University2026_Project/OUTPUTS/Metadata/oldsubpopulations.csv",
  Protein   = "Altprotiensubpopulations.csv",
  Noncoding   =  "Altnoncodingsubpopulations.csv"
)

calc_subpop_props <- function(file_path, set_name) {
  df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  df$group <- as.character(df$group)
  
  df$group[df$group == "4"] <- "Cluster 4"
  
  df$group[is.na(df$group)] <- "NA"
  
  counts <- table(df$group)
  total <- sum(counts)
  data.frame(
    Set    = set_name,
    Bar    = "Subpopulations",
    Type   = names(counts),
    Proportion = as.numeric(counts) / total
  )
}

subpop_data <- map_dfr(names(subpop_files), ~ calc_subpop_props(subpop_files[.x], .x))

# Combine all data
plot_data <- bind_rows(sample_data, assign_data, subpop_data)


plot_data$Bar <- factor(plot_data$Bar,
                        levels = c("Neurons/NPCs", "Assigned Zone", 
                                   "Assigned CellType", "Subpopulations"))

plot_data$Set <- factor(plot_data$Set, levels = c("Original", "Noncoding", "Protein"))


# 5. Set stacking order
type_order <- c(
  # Original bar
  "NPCs",
  "Neurons",
  # Zone bar
  "VZ",
  "iSVZ",
  "oSVZ",
  "CP",
  # CellType bar
  "aRG",
  "bRG",
  "neuron",
  # Subpopulations bar
  "AP2",
  "AP1",
  "BP1",
  "BP2",
  "N1",
  "N2",
  "N3",
  "Cluster 4"
)
plot_data$Type <- factor(plot_data$Type, levels = type_order)

# Custom colours
zone_colors <- c(
  "VZ"   = "#FFBF20",
  "iSVZ" = "#F77800",
  "oSVZ" = "#CD2624",
  "CP"   = "#7D007E"
)

celltype_colors <- c(
  "neuron" = "#2ca02c",
  "aRG"    = "#6bbf6b",
  "bRG"    = "#a5d9a5"
)

original_colors <- c(
  "Neurons" = "#1f77b4",
  "NPCs"    = "#75aad9"
)

subpop_colors <- c(
  "AP2" = "#95F297",
  "AP1" = "#2C8334",
  "BP1" = "#15B5B5",
  "BP2" = "#9DFFFF",
  "N1"  = "#8EC9F5",
  "N2"  = "#4589EF",
  "N3"  = "#0545BB",
  "Cluster 4"  = "#B0B0B0"
)

all_colors <- c(zone_colors, celltype_colors, original_colors, subpop_colors)

# Insert dummy types for legend gaps
dummy_gaps <- c("Assigned Zone", "Assigned Cell Type", "Subpopulations")
plot_data_with_gaps <- plot_data

gap_rows <- data.frame(
  Set = rep(levels(plot_data$Set)[1], 3),   
  Bar = factor("Neurons/NPCs", levels = levels(plot_data$Bar)),
  Type = dummy_gaps,
  Proportion = 0
)
plot_data_with_gaps <- bind_rows(plot_data_with_gaps, gap_rows)

type_order_with_gaps <- c(
  "NPCs", "Neurons",
  "Assigned Zone",
  "VZ", "iSVZ", "oSVZ", "CP",
  "Assigned Cell Type",
  "aRG", "bRG", "neuron",
  "Subpopulations",
  "AP2", "AP1", "BP1", "BP2", "N1", "N2", "N3", "Cluster 4"
)
plot_data_with_gaps$Type <- factor(plot_data_with_gaps$Type, levels = type_order_with_gaps)

# Add NA colours for gaps (invisible)
gap_colors <- c("Assigned Zone" = NA, "Assigned Cell Type" = NA, "Subpopulations" = NA)
all_colors_with_gaps <- c(all_colors, gap_colors)

# Plot
ggplot(plot_data_with_gaps, aes(x = Bar, y = Proportion, fill = Type)) +
  geom_col(width = 0.9, color = "white", linewidth = 0.2) +
  geom_text(data = subset(plot_data_with_gaps, Proportion > 0),  # skip dummies
            aes(label = scales::percent(Proportion, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), size = 3.5, color = "white") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = all_colors_with_gaps, name = "Category",
                    breaks = type_order_with_gaps) +  # ensure order
  facet_grid(~ Set, scales = "free_x", space = "free_x") +
  labs(title = "Proportions of identified cell types", x = NULL, y = "Proportion") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(3, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.major.x = element_blank(),
    legend.spacing.y = unit(0.1, "cm")   # optional: reduce space between entries
  )
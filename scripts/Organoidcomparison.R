# -------------------------------------------------------------------
# Compare cluster overlaps with descriptive labels and grouped bars
# -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)


# Read cluster files
read_cluster_file <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  clusters <- list()
  for (col in names(df)) {
    samples <- na.omit(df[[col]])
    samples <- samples[samples != ""]
    if (length(samples) > 0) clusters[[col]] <- samples
  }
  return(clusters)
}


# Define mapping table with labels (Discovered in Fig3D and S3H)

mappings <- tribble(
  ~Label,                      ~TestDataset, ~OrigClusters,              ~TestClusters,
  "Ventral forebrain NPCs",    "Noncoding",  list(c("C1","C4","C9")),    list(c("C1","C2")),
  "Ventral forebrain NPCs",    "Protein",    list(c("C1","C9","C4")),    list(c("C4","C5","C1")),
  "Dorsal forebrain NPCs",     "Noncoding",  list(c("C2","C3","C10")),   list(c("C3","C5")),
  "Dorsal forebrain NPCs",     "Protein",    list(c("C2","C10","C3")),   list(c("C6","C2")),
  "Ventral forebrain Neurons", "Noncoding", list(c("C7")),          list(c("C4")),
  "Ventral forebrain Neurons", "Protein",   list(c("C7")),          list(c("C3")),
  "Dorsal forebrain Neurons", "Noncoding", list(c("C8")),          list(c("C7")),
  "Dorsal forebrain Neurons", "Protein",   list(c("C8")),          list(c("C8")),
  "RSPO+",                     "Noncoding",  list(c("C5")),              list(c("C8")),
  "RSPO+",                     "Protein",    list(c("C5")),              list(c("C9")),
  "Mesenchymal cells",         "Noncoding",  list(c("C6")),              list(c("C6")),
  "Mesenchymal cells",         "Protein",    list(c("C6")),              list(c("C7"))
)


# Function to compute overlap for one row
compute_overlap_row <- function(orig_clusters, test_clusters_list, orig_vec, test_vec, test_name) {
  test_clusters <- if (test_name == "Noncoding") "ALtnoncodingOrganoidclusters.csv" else "ALtOrganoidclusters.csv"
  
  orig_combined <- unique(unlist(orig_clusters[orig_vec], use.names = FALSE))
  test_combined <- unique(unlist(test_clusters[test_vec], use.names = FALSE))
  
  if (length(test_combined) == 0) {
    return(data.frame(
      TestDataset = test_name,
      Orig_Size = length(orig_combined),
      Test_Size = 0,
      Overlap_Size = 0,
      Percent_Overlap = 0,
      stringsAsFactors = FALSE
    ))
  }
  
  overlap <- intersect(orig_combined, test_combined)
  data.frame(
    TestDataset = test_name,
    Orig_Size = length(orig_combined),
    Test_Size = length(test_combined),
    Overlap_Size = length(overlap),
    Percent_Overlap = length(overlap) / length(test_combined) * 100,
    stringsAsFactors = FALSE
  )
}

# Apply to each row using pmap_dfr
results <- pmap_dfr(
  list(
    Label = mappings$Label,
    TestDataset = mappings$TestDataset,
    OrigVec = mappings$OrigClusters,
    TestVec = mappings$TestClusters
  ),
  function(Label, TestDataset, OrigVec, TestVec) {
    compute_overlap_row("Organoidclusters.csv", NULL, unlist(OrigVec), unlist(TestVec), TestDataset) %>%
      mutate(Label = Label)
  }
)

results$Label <- factor(results$Label, levels = unique(mappings$Label))


# Assign numeric x positions for grouped bars
label_levels <- levels(results$Label)
results <- results %>%
  mutate(label_idx = as.numeric(Label),
         x_pos = ifelse(TestDataset == "Noncoding", label_idx * 3 - 1, label_idx * 3))


# Plot with custom x-axis breaks and labels
x_breaks <- seq_along(label_levels) * 3 - 0.5
x_labels <- label_levels

ggplot(results, aes(x = x_pos, y = Percent_Overlap, fill = TestDataset)) +
  geom_bar(stat = "identity", position = "identity", width = 0.8, color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", Percent_Overlap)),
            vjust = -0.5, size = 3.5) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.02, 0)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 105), expand = c(0, 0)) +
  scale_fill_manual(values = c("Noncoding" = "#E69F00", "Protein" = "#56B4E9")) +
  labs(title = "Overlap of cluster assignments with original pipeline",
       x = NULL,
       y = "Percentage overlap (%)",
       fill = "Test dataset") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom",
        panel.grid.major.x = element_blank())
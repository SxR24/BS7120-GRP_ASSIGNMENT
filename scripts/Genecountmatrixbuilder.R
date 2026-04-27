# Input
star_dir <- "D:/GroupA_University2026_Project/STAR_out/"

# Get a list of all sample folders (e.g., SRR*)
sample_folders <- list.dirs(star_dir, recursive = FALSE, full.names = TRUE)
sample_folders <- grep("SRR", sample_folders, value = TRUE)

# Initialize an empty list to hold count data
count_list <- list()

# Loop over each sample folder
for (folder in sample_folders) {
  file_path <- file.path(folder, "ReadsPerGene.out.tab")
  
  counts <- read.table(file_path, header = FALSE, row.names = 1, sep = "\t")[, 2, drop = FALSE]
  
  sample_name <- basename(folder)
  colnames(counts) <- sample_name
  
  count_list[[sample_name]] <- counts
}

# Merge all samples by gene IDs
count_matrix <- do.call(cbind, count_list)

head(count_matrix[, 1:5])

write.table(count_matrix, file = "counts.txt", sep = "\t", quote = FALSE, col.names = NA)
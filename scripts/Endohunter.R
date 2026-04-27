# --------------------------------------------------------------------------
# Identifies Cell discrepancies between matrices to find endothelial cell
# --------------------------------------------------------------------------

# Inputs
file_path <- "D:/GroupA_University2026_Project/PCA_Output/GSE75140_series_matrix.txt"
lines <- readLines(file_path, warn = FALSE)

title_row <- grep("^!Sample_title", lines, value = TRUE)
parts <- strsplit(title_row, "\t")[[1]]
sample_names <- gsub('^"|"$', '', parts[-1])   # remove surrounding quotes if any

# Keep only samples starting with "fetal" (case-insensitive)
fetal_samples <- sample_names[grepl("^fetal", sample_names, ignore.case = TRUE)]

# Remove the prefixes "fetal1_", "fetal2_", "fetal3_" from the beginning of each name
stripped_ids <- gsub("^fetal[123]_", "", fetal_samples, ignore.case = TRUE)

cat("Number of fetal samples found:", length(fetal_samples), "\n")
cat("First 5 original fetal sample names:\n")
print(head(fetal_samples, 5))
cat("\nFirst 5 stripped IDs:\n")
print(head(stripped_ids, 5))


master_file <- "D:/GroupA_University2026_Project/PCA_Output/GSE75140_hOrg.fetal.master.data.frame.txt"

master <- read.table(master_file, sep = "\t", header = TRUE, row.names = 1, 
                     check.names = FALSE, stringsAsFactors = FALSE)

# Now row names are the sample identifiers
master_rows <- rownames(master)
cat("\nNumber of rows (samples) in master data frame:", length(master_rows), "\n")
cat("First 5 master sample IDs (row names):\n")
print(head(master_rows, 5))

# Compare stripped IDs with original row names
stripped_lower <- tolower(stripped_ids)
master_lower <- tolower(master_rows)

missing_in_master <- stripped_ids[!stripped_lower %in% master_lower]

extra_in_master <- master_rows[!master_lower %in% stripped_lower]

# Print results
cat("\n========================================\n")
cat("COMPARISON RESULTS\n")
cat("========================================\n\n")

if (length(missing_in_master) > 0) {
  cat("Not found in papers matrix (", length(missing_in_master), " total):\n", sep = "")
  print(missing_in_master)
} else {
  cat("Endothelial cell found.\n")
}

cat("\n")

if (length(extra_in_master) > 0) {
  cat("Extra cell not found", length(extra_in_master), " total):\n", sep = "")
  print(extra_in_master)
} else {
  cat("extra cell found\n")
}
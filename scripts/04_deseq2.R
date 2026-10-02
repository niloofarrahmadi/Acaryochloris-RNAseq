# Differential expression analysis of Acaryochloris marina RNA-seq data
# Main comparison: Orange vs FarRed
# Analysis performed in RStudio with R 4.3.3

library(DESeq2)

# -----------------------------
# 1. Read count matrix
# -----------------------------

counts <- read.delim(
  "counts/gene_counts.txt",
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

count_matrix <- counts[, 7:12]

rownames(count_matrix) <- counts$Geneid

colnames(count_matrix) <- c(
  "FarRed_1",
  "Orange_1",
  "FarRed_2",
  "Orange_2",
  "FarRed_3",
  "Orange_3"
)

# -----------------------------
# 2. Define metadata
# -----------------------------

metadata <- data.frame(
  condition = factor(
    c(
      "FarRed",
      "Orange",
      "FarRed",
      "Orange",
      "FarRed",
      "Orange"
    ),
    levels = c("FarRed", "Orange")
  )
)

rownames(metadata) <- colnames(count_matrix)

# -----------------------------
# 3. Create DESeq2 dataset
# -----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

# -----------------------------
# 4. Low-count filtering
# -----------------------------

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# -----------------------------
# 5. Differential expression
# -----------------------------

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Orange", "FarRed")
)

# Save differential expression results
dir.create("results", showWarnings = FALSE)

write.csv(
  as.data.frame(res),
  "results/deseq2_results.csv"
)

# -----------------------------
# 6. Normalized counts
# -----------------------------

norm_counts <- counts(
  dds,
  normalized = TRUE
)

write.csv(
  as.data.frame(norm_counts),
  "results/normalized_counts.csv"
)

# -----------------------------
# 7. Variance-stabilizing transformation
# -----------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

# -----------------------------
# 8. PCA
# -----------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

write.csv(
  pca_data,
  "results/pca_data.csv"
)

# -----------------------------
# 9. Save size factors
# -----------------------------

write.csv(
  data.frame(
    sample = colnames(dds),
    size_factor = sizeFactors(dds)
  ),
  "results/size_factors.csv",
  row.names = FALSE
)

cat("Analysis completed successfully.\n")
cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering:", nrow(dds), "\n")

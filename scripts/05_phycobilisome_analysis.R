# Focused analysis of structural phycobilisome-associated genes
# Acaryochloris marina MBIC11017
# Analysis performed in RStudio with R 4.3.3

library(DESeq2)
library(pheatmap)

dir.create("results", showWarnings = FALSE)

# --------------------------------------------------
# 1. Load DESeq2 results and count data
# --------------------------------------------------

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

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Orange", "FarRed")
)

# --------------------------------------------------
# 2. Read current RefSeq annotation
# --------------------------------------------------

gff_lines <- readLines(
  gzfile("annotation/GCF_000018105.1_ASM1810v1_genomic.gff.gz")
)

cds_lines <- gff_lines[grepl("\tCDS\t", gff_lines)]

get_attr <- function(x, key) {
  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  hit <- parts[startsWith(parts, paste0(key, "="))]
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^", key, "="), "", hit[1])
}

cds_annotation <- data.frame(
  gene_id = sapply(cds_lines, get_attr, key = "Parent"),
  locus_tag = sapply(cds_lines, get_attr, key = "locus_tag"),
  gene = sapply(cds_lines, get_attr, key = "gene"),
  product = sapply(cds_lines, get_attr, key = "product"),
  protein_id = sapply(cds_lines, get_attr, key = "protein_id"),
  stringsAsFactors = FALSE
)

# Remove duplicated annotation records
cds_annotation <- unique(cds_annotation)

# Keep one annotation record per gene
cds_annotation <- cds_annotation[
  !duplicated(cds_annotation$gene_id),
]

# --------------------------------------------------
# 3. Merge DESeq2 results with annotation
# --------------------------------------------------

res_annotated <- merge(
  data.frame(
    gene_id = rownames(res),
    as.data.frame(res),
    row.names = NULL
  ),
  cds_annotation,
  by = "gene_id",
  all.x = TRUE
)

# --------------------------------------------------
# 4. Identify structural phycobilisome-associated genes
# --------------------------------------------------

structural_phyco <- res_annotated[
  grepl(
    "phycobilisome (rod-core )?linker polypeptide|phycocyanin subunit",
    res_annotated$product,
    ignore.case = TRUE
  ),
]

# Save focused gene results
write.csv(
  structural_phyco,
  "results/phycobilisome_structural_genes.csv",
  row.names = FALSE
)

# --------------------------------------------------
# 5. Summary of the 15-gene set
# --------------------------------------------------

summary_structural <- data.frame(
  n_genes = nrow(structural_phyco),
  n_positive = sum(
    structural_phyco$log2FoldChange > 0,
    na.rm = TRUE
  ),
  n_lfc_gt1 = sum(
    structural_phyco$log2FoldChange > 1,
    na.rm = TRUE
  ),
  n_sig = sum(
    structural_phyco$padj < 0.05,
    na.rm = TRUE
  ),
  n_sig_lfc_gt1 = sum(
    structural_phyco$padj < 0.05 &
      structural_phyco$log2FoldChange > 1,
    na.rm = TRUE
  ),
  mean_log2FC = mean(
    structural_phyco$log2FoldChange,
    na.rm = TRUE
  ),
  median_log2FC = median(
    structural_phyco$log2FoldChange,
    na.rm = TRUE
  )
)

write.csv(
  summary_structural,
  "results/phycobilisome_summary.csv",
  row.names = FALSE
)

# --------------------------------------------------
# 6. Normalized expression values for the 15 genes
# --------------------------------------------------

norm_counts <- counts(
  dds,
  normalized = TRUE
)

structural_ids <- structural_phyco$gene_id

structural_norm <- norm_counts[
  structural_ids,
  ,
  drop = FALSE
]

structural_norm <- data.frame(
  gene_id = rownames(structural_norm),
  structural_norm,
  row.names = NULL
)

write.csv(
  structural_norm,
  "results/phycobilisome_normalized_counts.csv",
  row.names = FALSE
)

# --------------------------------------------------
# 7. Variance-stabilized expression and heatmap
# --------------------------------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

mat <- assay(vsd)[structural_ids, ]

rownames(mat) <- structural_phyco$locus_tag

mat_z <- t(scale(t(mat)))

annotation_col <- data.frame(
  Condition = c(
    "FarRed",
    "Orange",
    "FarRed",
    "Orange",
    "FarRed",
    "Orange"
  )
)

rownames(annotation_col) <- colnames(mat_z)

png(
  "results/phycobilisome_heatmap.png",
  width = 1800,
  height = 1500,
  res = 200
)

pheatmap(
  mat_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Structural phycobilisome-associated genes",
  fontsize_row = 9,
  fontsize_col = 10,
  border_color = NA
)

dev.off()

# --------------------------------------------------
# 8. Sensitivity analysis
# --------------------------------------------------

dds_sens <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

keep_sens <- rowSums(
  counts(dds_sens) >= 10
) >= 4

dds_sens <- dds_sens[keep_sens, ]

dds_sens <- DESeq(dds_sens)

res_sens <- results(
  dds_sens,
  contrast = c("condition", "Orange", "FarRed")
)

sens_15 <- res_sens[
  structural_ids,
  c("log2FoldChange", "padj")
]

sens_15$locus_tag <- structural_phyco$locus_tag

sens_15 <- sens_15[
  ,
  c("locus_tag", "log2FoldChange", "padj")
]

write.csv(
  sens_15,
  "results/sensitivity_15_genes.csv",
  row.names = FALSE
)

sensitivity_summary <- data.frame(
  n_positive = sum(
    sens_15$log2FoldChange > 0,
    na.rm = TRUE
  ),
  n_lfc_gt1 = sum(
    sens_15$log2FoldChange > 1,
    na.rm = TRUE
  ),
  n_sig = sum(
    sens_15$padj < 0.05,
    na.rm = TRUE
  ),
  n_sig_lfc_gt1 = sum(
    sens_15$padj < 0.05 &
      sens_15$log2FoldChange > 1,
    na.rm = TRUE
  )
)

write.csv(
  sensitivity_summary,
  "results/sensitivity_summary.csv",
  row.names = FALSE
)

cat("Focused analysis completed successfully.\n")
cat("Structural genes:", nrow(structural_phyco), "\n")
cat("Sensitivity genes:", nrow(sens_15), "\n")

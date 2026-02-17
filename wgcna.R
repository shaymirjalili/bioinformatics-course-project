# ============================================================
# WGCNA Figure 4 Reproduction Script (FIXED)
# ============================================================
# Fixes applied:
# 1) Correctly aligns genes across datasets BEFORE cbind (prevents scrambled rows)
# 2) Collapses duplicate probes -> one gene (mean) instead of dropping randomly
# 3) Uses TOMsimilarityFromExpr (more stable, standard WGCNA)
# 4) Uses TOMplot (WGCNA-style heatmap with proper coloring)
# 5) Fixes dendrogram color labels (no mismatch in groupLabels)
# ============================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(GEOquery)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

set.seed(1234)
allowWGCNAThreads()

OUTDIR <- "geo_deg_out"
WGCNA_DIR <- file.path(OUTDIR, "wgcna_figure4")
dir.create(WGCNA_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
clean_symbols <- function(x) {
  x <- trimws(x)
  x <- sapply(strsplit(x, " /// "), `[`, 1)   # keep first symbol if "A /// B"
  x <- gsub("\\s+", "", x)
  x <- toupper(x)
  x[x == ""] <- NA_character_
  x
}

collapse_by_gene_mean <- function(expr_mat) {
  # rows=genes/probes, cols=samples, rownames are symbols (may repeat)
  rn <- rownames(expr_mat)
  keep <- !is.na(rn) & rn != ""
  expr_mat <- expr_mat[keep, , drop = FALSE]
  rn <- rn[keep]

  # mean across duplicates
  # rowsum() sums; divide by counts to get mean
  sums <- rowsum(expr_mat, group = rn, reorder = FALSE)
  counts <- as.numeric(table(rn))[match(rownames(sums), names(table(rn)))]
  means <- sweep(sums, 1, counts, "/")
  means
}

# ----------------------------
# Load overlapping DEGs (symbols)
# ----------------------------
common_up <- readLines(file.path(OUTDIR, "overlap_up_symbols.txt"))
common_down <- readLines(file.path(OUTDIR, "overlap_down_symbols.txt"))
all_overlapping_genes <- unique(clean_symbols(c(common_up, common_down)))
all_overlapping_genes <- all_overlapping_genes[!is.na(all_overlapping_genes)]

# ----------------------------
# Load expression data from all datasets
# ----------------------------
source("data_gathering.R")  # must define: get_eset(), maybe_log2(), extract_symbol()

gse_ids <- c("GSE3268", "GSE1987", "GSE31547", "GSE18842")
expr_list <- list()

for (gse_id in gse_ids) {
  message("Loading ", gse_id, " ...")
  eset <- get_eset(gse_id)

  expr <- exprs(eset)
  expr <- maybe_log2(expr)

  fdat <- fData(eset)
  symbols <- extract_symbol(fdat)  # vector aligned with rows(expr)

  if (is.null(symbols)) stop("No gene symbols extracted for ", gse_id)

  symbols <- clean_symbols(symbols)
  rownames(expr) <- symbols

  # collapse duplicate probes -> one gene
  expr_gene <- collapse_by_gene_mean(expr)

  # keep only overlapping genes (for speed), but do not force intersection yet
  expr_gene <- expr_gene[rownames(expr_gene) %in% all_overlapping_genes, , drop = FALSE]

  expr_list[[gse_id]] <- expr_gene
}

# ----------------------------
# CRITICAL: genes must exist in ALL datasets (intersection)
# Otherwise you get NAs / weird sparsity
# ----------------------------
genes_in_all <- Reduce(intersect, lapply(expr_list, rownames))
overlapping_in_data <- intersect(all_overlapping_genes, genes_in_all)

if (length(overlapping_in_data) < 30) {
  stop("Too few genes (", length(overlapping_in_data),
       ") remain after intersecting across datasets. WGCNA will be unstable.")
}

# Keep stable order
overlapping_in_data <- sort(overlapping_in_data)

# ----------------------------
# Combine expression matrices (ALIGNED rows!)
# ----------------------------
combined_expr <- NULL
for (gse_id in names(expr_list)) {
  expr <- expr_list[[gse_id]]

  # enforce identical gene order before cbind (THIS FIXES THE SCRAMBLING BUG)
  expr_subset <- expr[overlapping_in_data, , drop = FALSE]

  combined_expr <- if (is.null(combined_expr)) expr_subset else cbind(combined_expr, expr_subset)
}

# WGCNA needs: samples x genes
expr_wgcna <- t(combined_expr)

# Safety check
gsg <- goodSamplesGenes(expr_wgcna, verbose = 3)
if (!gsg$allOK) {
  expr_wgcna <- expr_wgcna[gsg$goodSamples, gsg$goodGenes]
}

message("WGCNA matrix: ", nrow(expr_wgcna), " samples x ", ncol(expr_wgcna), " genes")

# ----------------------------
# Figure 4A: Soft threshold selection
# ----------------------------
png(file.path(WGCNA_DIR, "figure_4a.png"), width = 1100, height = 550, res = 140)
par(mfrow = c(1, 2))

powers <- 1:30
sft <- pickSoftThreshold(expr_wgcna, powerVector = powers, networkType = "signed", verbose = 5)

# Scale independence
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale Independence",
     type = "n")
points(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       pch = 20, cex = 1.6)
abline(h = 0.5, col = "red", lwd = 2)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = sft$fitIndices[, 1], cex = 0.75, pos = 3)

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean Connectivity",
     type = "n")
points(sft$fitIndices[, 1], sft$fitIndices[, 5],
       pch = 20, cex = 1.6)
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = sft$fitIndices[, 1], cex = 0.75, pos = 3)

dev.off()

# Use the paper's power
soft_threshold <- 16

# ----------------------------
# Build TOM (standard, stable)
# ----------------------------
TOM <- TOMsimilarityFromExpr(expr_wgcna, power = soft_threshold, networkType = "signed")
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

minModuleSize <- 10
deepSplit <- 1

dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = deepSplit,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)
dynamicColors <- labels2colors(dynamicMods)

# (Optional but useful) merge similar modules
MEList <- moduleEigengenes(expr_wgcna, colors = dynamicColors)
MEs <- MEList$eigengenes
merge <- mergeCloseModules(expr_wgcna, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

module_summary <- table(mergedColors)

# ----------------------------
# Figure 4B: Dendrogram
# ----------------------------
png(file.path(WGCNA_DIR, "figure_4b.png"), width = 1100, height = 500, res = 140)
plotDendroAndColors(
  geneTree,
  cbind(dynamicColors, mergedColors),
  groupLabels = c("DynamicTreeCut", "Merged"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Cluster Dendrogram"
)
dev.off()

# ----------------------------
# Figure 4C: TOM heatmap (WGCNA style, colored correctly)
# ----------------------------
# This produces the "blocky" colored TOM plot like many papers.
png(file.path(WGCNA_DIR, "figure_4c.png"), width = 1100, height = 1100, res = 140)
TOMplot(dissTOM, geneTree, mergedColors, main = "Network heatmap plot, all genes")
dev.off()

# ----------------------------
# ZBTB16 module check
# ----------------------------
cat("\nWGCNA module summary (power=16, merged colors):\n")
print(module_summary)

if ("ZBTB16" %in% colnames(expr_wgcna)) {
  zbtb16_module <- mergedColors[which(colnames(expr_wgcna) == "ZBTB16")]
  cat("ZBTB16 is in module:", zbtb16_module, "\n")
} else {
  cat("ZBTB16 not found in expr_wgcna columns.\n")
}

# Save module assignments
module_df <- data.frame(
  gene = colnames(expr_wgcna),
  module_dynamic = dynamicColors,
  module_merged = mergedColors
)
write.csv(module_df, file.path(WGCNA_DIR, "gene_module_assignment.csv"), row.names = FALSE)

cat("\nAll figures and module assignments saved to:", normalizePath(WGCNA_DIR), "\n")

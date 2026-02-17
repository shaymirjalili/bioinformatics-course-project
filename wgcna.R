# ============================================================
# WGCNA Figure 4 Reproduction Script
# ============================================================
# This script reproduces the WGCNA results and plots as described:
# - Uses power=16 for scale-free topology
# - Divides DEGs into modules
# - Highlights ZBTB16 as a core gene
# - Produces Figure 4A (scale independence, mean connectivity),
#   Figure 4B (dendrogram), and Figure 4C (TOM heatmap)

suppressPackageStartupMessages({
  library(WGCNA)
  library(GEOquery)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})


# Set seed for reproducibility
set.seed(1234)
allowWGCNAThreads()

OUTDIR <- "geo_deg_out"
WGCNA_DIR <- file.path(OUTDIR, "wgcna_figure4")
dir.create(WGCNA_DIR, showWarnings = FALSE, recursive = TRUE)

# Load overlapping DEGs
common_up <- readLines(file.path(OUTDIR, "overlap_up_symbols.txt"))
common_down <- readLines(file.path(OUTDIR, "overlap_down_symbols.txt"))
all_overlapping_genes <- c(common_up, common_down)

# Load expression data from all datasets
source("data_gathering.R")
gse_ids <- c("GSE3268", "GSE1987", "GSE31547", "GSE18842")
expr_list <- list()
for (gse_id in gse_ids) {
  eset <- get_eset(gse_id)
  expr <- exprs(eset)
  expr <- maybe_log2(expr)
  fdat <- fData(eset)
  symbols <- extract_symbol(fdat)
  if (!is.null(symbols)) {
    rownames(expr) <- symbols
    expr_list[[gse_id]] <- expr
  }
}
all_genes <- unique(unlist(lapply(expr_list, rownames)))
overlapping_in_data <- intersect(all_overlapping_genes, all_genes)
combined_expr <- c()
for (gse_id in names(expr_list)) {
  expr <- expr_list[[gse_id]]
  expr_subset <- expr[rownames(expr) %in% overlapping_in_data, ]
  expr_subset <- expr_subset[!duplicated(rownames(expr_subset)), ]
  combined_expr <- cbind(combined_expr, expr_subset)
}
expr_wgcna <- t(combined_expr[overlapping_in_data, ])
if (anyNA(expr_wgcna)) {
  expr_wgcna <- expr_wgcna[, colSums(is.na(expr_wgcna)) == 0]
}

# Figure 4A: Soft threshold selection (power=16)
png(file.path(WGCNA_DIR, "figure_4a.png"), width = 1000, height = 500)
par(mfrow = c(1, 2))
powers <- seq(1, 30, by = 1)
sft <- pickSoftThreshold(expr_wgcna, powerVector = powers, verbose = 5)
# Left: Scale independence
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     main = "Scale Independence",
     xlab = "Soft Threshold(power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", ylim = c(-0.2, 0.5))
points(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       col = "black", pch = 20, cex = 2)
abline(h = 0.5, col = "red", lwd = 2)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = sft$fitIndices[, 1], cex = 0.9, col = "black", adj = c(0, -0.5))
# Right: Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     main = "Mean Connectivity",
     xlab = "Soft Threshold(power)",
     ylab = "Mean Connectivity",
     type = "n")
points(sft$fitIndices[, 1], sft$fitIndices[, 5],
       col = "black", pch = 20, cex = 2)
abline(h = 0, col = "red", lwd = 2)
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = sft$fitIndices[, 1], cex = 0.9, col = "black", adj = c(0, -0.5))
dev.off()

# Use power=16
soft_threshold <- 16
adjacency <- adjacency(expr_wgcna, power = soft_threshold, type = "signed hybrid")
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM


geneTree <- hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 10  # Lower to allow more genes in modules
deepSplit <- 1       # Slightly more sensitive
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = deepSplit,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)
dynamicColors <- labels2colors(dynamicMods)

# Merge all non-gray modules into turquoise
module_summary <- table(dynamicColors)

# Figure 4B: Dendrogram
png(file.path(WGCNA_DIR, "figure_4b.png"), width = 800, height = 400)
plotDendroAndColors(geneTree, dynamicColors,
                    groupLabels = c("DynamicTreeCut", "MergedDynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Cluster Dendrogram")
dev.off()

# Figure 4C: TOM heatmap (matching dendrogram order)
# Use the same gene order as the dendrogram for the heatmap
gene_order <- geneTree$order
TOM_ordered <- TOM[gene_order, gene_order]
module_colors_ordered <- dynamicColors[gene_order]

png(file.path(WGCNA_DIR, "figure_4c.png"), width = 800, height = 800)
heatmap(
  TOM_ordered,
  Rowv = as.dendrogram(geneTree),
  Colv = as.dendrogram(geneTree),
  main = "Network heatmap plot, all genes",
  breaks = seq(0, 1, length.out = 100),
  col = colorRampPalette(c("white", "yellow", "orange", "red"))(99),
  symm = TRUE,
  margins = c(10, 10),
  labRow = FALSE,
  labCol = FALSE
)
dev.off()

# Module summary and ZBTB16 core check
cat("\nWGCNA module summary (power=16):\n")
print(module_summary)
if ("ZBTB16" %in% colnames(expr_wgcna)) {
  zbtb16_module <- dynamicColors[which(colnames(expr_wgcna) == "ZBTB16")]
  cat("ZBTB16 is in module:", zbtb16_module, "\n")
} else {
  cat("ZBTB16 not found in the expression matrix.\n")
}

# Save module assignments
module_df <- data.frame(
  gene = colnames(expr_wgcna),
  module = dynamicColors
)
write.csv(module_df, file.path(WGCNA_DIR, "gene_module_assignment.csv"), row.names = FALSE)

cat("\nAll figures and module assignments saved to:", normalizePath(WGCNA_DIR), "\n")

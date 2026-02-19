# ============================================================
# Figure 6 Reproduction Script
# Immune infiltration analysis of ZBTB16 in lung cancer (ssGSEA)
# Panels:
#   A) Immune signatures by group (tumor vs normal)
#   B) Correlation lollipop (ZBTB16 vs immune signatures)
#   C) Scatter: ZBTB16 vs Mast cells
#   D) Scatter: ZBTB16 vs Eosinophils
# ============================================================

# ----------------------------
# User settings
# ----------------------------
FIG6_GSE_IDS <- c("GSE3268", "GSE1987", "GSE31547", "GSE18842")
FIG6_OUTDIR <- file.path("geo_deg_out", "figure6")
FIG6_IMMUNE_GMT <- file.path("resources", "immune_signatures_24cells.gmt")
FIG6_TARGET_GENE <- "ZBTB16"

# Optional explicit sample IDs (recommended if available)
FIG6_GROUPS_BY_GSE <- list(
  # "GSE31547" = list(
  #   normal = c("GSMxxxx", "GSMyyyy"),
  #   tumor  = c("GSMzzzz", "GSMwwww")
  # )
)

# ----------------------------
# Package setup
# ----------------------------
ensure_pkg <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))

  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Failed to install/load package: ", pkg)
  }
  invisible(TRUE)
}

ensure_pkg("GEOquery", bioc = TRUE)
ensure_pkg("GSVA", bioc = TRUE)
ensure_pkg("GSEABase", bioc = TRUE)
ensure_pkg("dplyr")
ensure_pkg("tidyr")
ensure_pkg("ggplot2")
ensure_pkg("stringr")
ensure_pkg("ggpubr")

suppressPackageStartupMessages({
  library(GEOquery)
  library(GSVA)
  library(GSEABase)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(ggpubr)
})

# ----------------------------
# Load helper functions
# ----------------------------
source("data_gathering.R")

# ----------------------------
# Local helpers
# ----------------------------
clean_symbols_fig6 <- function(x) {
  x <- trimws(x)
  x <- sapply(strsplit(x, " /// "), `[`, 1)
  x <- gsub("\\s+", "", x)
  x <- toupper(x)
  x[x == ""] <- NA_character_
  x
}

collapse_by_gene_mean <- function(expr_mat) {
  rn <- rownames(expr_mat)
  keep <- !is.na(rn) & rn != ""
  expr_mat <- expr_mat[keep, , drop = FALSE]
  rn <- rn[keep]

  sums <- rowsum(expr_mat, group = rn, reorder = FALSE)
  counts <- as.numeric(table(rn))[match(rownames(sums), names(table(rn)))]
  sweep(sums, 1, counts, "/")
}

infer_groups_local <- function(pdat) {
  fields <- intersect(c("title", "source_name_ch1", "characteristics_ch1"), colnames(pdat))
  if (length(fields) == 0) return(rep(NA_character_, nrow(pdat)))

  text <- apply(pdat[, fields, drop = FALSE], 1, paste, collapse = " | ")
  text_l <- tolower(text)

  is_tumor <- str_detect(text_l, "tumou?r|cancer|carcinoma|adenocarcinoma|nsclc|lusc|luad")
  is_normal <- str_detect(text_l, "normal|control|non[- ]?tumou?r|adjacent|healthy")

  group <- rep(NA_character_, nrow(pdat))
  group[is_tumor & !is_normal] <- "tumor"
  group[is_normal & !is_tumor] <- "normal"
  group
}

find_signature_name <- function(sig_names, candidates) {
  normalize_sig <- function(x) {
    gsub("[^a-z0-9]+", " ", tolower(x))
  }
  sig_l <- normalize_sig(sig_names)
  cand_l <- normalize_sig(candidates)

  idx <- match(cand_l, sig_l)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(sig_names[idx[1]])

  for (cnd in cand_l) {
    m <- grep(cnd, sig_l, fixed = TRUE)
    if (length(m) > 0) return(sig_names[m[1]])
  }

  NA_character_
}

make_scatter <- function(df, signature_name, panel_label, outfile) {
  x <- df$zbtb16_expr
  y <- df[[signature_name]]

  ctest <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))

  p <- ggplot(df, aes(x = zbtb16_expr, y = .data[[signature_name]])) +
    geom_point(color = "#66C2DC", alpha = 0.8, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "#4D7FB5", fill = "#9ABCE0") +
    theme_bw() +
    labs(
      title = panel_label,
      x = "The expression of ZBTB16\nLog2(TPM+1)",
      y = paste0("Enrichment of ", signature_name)
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      hjust = 1.05,
      vjust = -0.2,
      label = sprintf("Spearman\nR = %.3f\nP %s", ctest$estimate, ifelse(ctest$p.value < 0.001, "< 0.001", paste0("= ", signif(ctest$p.value, 3)))),
      size = 3.5
    )

  ggsave(outfile, plot = p, width = 5.2, height = 4.2, dpi = 300)
  p
}

# ----------------------------
# Input checks
# ----------------------------
if (!file.exists(FIG6_IMMUNE_GMT)) {
  stop(
    "Missing immune signature GMT file: ", FIG6_IMMUNE_GMT, "\n",
    "Add the GMT file and rerun."
  )
}

dir.create(FIG6_OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Load and merge expression data
# ----------------------------
expr_list <- list()
sample_meta_list <- list()

for (gse_id in FIG6_GSE_IDS) {
  message("Loading ", gse_id, " for Figure 6 analysis ...")
  eset <- get_eset(gse_id)

  expr <- exprs(eset)
  expr <- maybe_log2(expr)

  pdat <- pData(eset)
  fdat <- fData(eset)

  if (!all(rownames(expr) %in% rownames(fdat))) {
    stop("Feature IDs in expression matrix are not all present in feature data for ", gse_id)
  }
  fdat <- fdat[rownames(expr), , drop = FALSE]

  symbols <- extract_symbol(fdat)
  if (is.null(symbols)) stop("Could not extract gene symbols for ", gse_id)
  if (length(symbols) != nrow(expr)) {
    stop("Symbol length mismatch for ", gse_id, ": symbols=", length(symbols), ", expr rows=", nrow(expr))
  }

  rownames(expr) <- clean_symbols_fig6(symbols)
  expr_gene <- collapse_by_gene_mean(expr)

  group <- rep(NA_character_, nrow(pdat))
  if (gse_id %in% names(FIG6_GROUPS_BY_GSE) && length(FIG6_GROUPS_BY_GSE[[gse_id]]) > 0) {
    tumor_ids <- FIG6_GROUPS_BY_GSE[[gse_id]]$tumor
    normal_ids <- FIG6_GROUPS_BY_GSE[[gse_id]]$normal
    group[rownames(pdat) %in% tumor_ids] <- "tumor"
    group[rownames(pdat) %in% normal_ids] <- "normal"
    message(gse_id, ": groups assigned via FIG6_GROUPS_BY_GSE.")
  } else {
    group <- infer_groups_local(pdat)
    message(gse_id, ": groups inferred heuristically:")
    print(table(group, useNA = "ifany"))
  }

  # Prefix sample IDs with GSE to avoid collisions when merging studies.
  new_sample_ids <- paste(gse_id, colnames(expr_gene), sep = "::")
  colnames(expr_gene) <- new_sample_ids

  sample_meta_list[[gse_id]] <- data.frame(
    sample = new_sample_ids,
    group = group[match(sub("^.*::", "", new_sample_ids), rownames(pdat))],
    gse = gse_id,
    stringsAsFactors = FALSE
  )

  expr_list[[gse_id]] <- expr_gene
}

genes_in_all <- Reduce(intersect, lapply(expr_list, rownames))
if (length(genes_in_all) < 100) {
  stop(
    "Too few intersecting genes across datasets (", length(genes_in_all), "). ",
    "Check symbol mapping and datasets."
  )
}

genes_in_all <- sort(genes_in_all)
expr_gene <- do.call(
  cbind,
  lapply(expr_list, function(m) m[genes_in_all, , drop = FALSE])
)

sample_df <- do.call(rbind, sample_meta_list)
sample_df <- sample_df[match(colnames(expr_gene), sample_df$sample), , drop = FALSE]

if (!all(c("tumor", "normal") %in% sample_df$group)) {
  stop(
    "Could not assign both tumor and normal groups across merged datasets.\n",
    "Provide explicit GSM IDs in FIG6_GROUPS_BY_GSE for problematic datasets."
  )
}

# ----------------------------
# Run ssGSEA
# ----------------------------
gmt <- getGmt(FIG6_IMMUNE_GMT)
immune_sets <- geneIds(gmt)

message("Running ssGSEA on ", nrow(expr_gene), " genes x ", ncol(expr_gene), " samples ...")
expr_mat <- as.matrix(expr_gene)
gsva_formals <- names(formals(GSVA::gsva))

# GSVA >= 2.0 uses parameter objects; older versions use legacy arguments.
if ("param" %in% gsva_formals && exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
  ssgsea_param <- GSVA::ssgseaParam(
    exprData = expr_mat,
    geneSets = immune_sets,
    minSize = 2,
    normalize = TRUE
  )
  ssgsea_scores <- GSVA::gsva(ssgsea_param, verbose = FALSE)
} else {
  ssgsea_scores <- GSVA::gsva(
    expr = expr_mat,
    gset.idx.list = immune_sets,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = TRUE,
    min.sz = 2,
    ssgsea.norm = TRUE,
    verbose = FALSE
  )
}

# keep only non-NA group samples
keep_samples <- !is.na(sample_df$group)
sample_df <- sample_df[keep_samples, , drop = FALSE]
ssgsea_scores <- ssgsea_scores[, sample_df$sample, drop = FALSE]

sample_df$group <- factor(sample_df$group, levels = c("normal", "tumor"))

# ----------------------------
# Panel A: immune signatures by group
# ----------------------------
a_df <- as.data.frame(t(ssgsea_scores))
a_df$sample <- rownames(a_df)
a_df <- a_df %>%
  left_join(sample_df, by = "sample") %>%
  pivot_longer(cols = all_of(rownames(ssgsea_scores)), names_to = "signature", values_to = "score")

p_a <- ggplot(a_df, aes(x = signature, y = score, fill = group)) +
  geom_boxplot(width = 0.75, outlier.size = 0.4, alpha = 0.95) +
  stat_compare_means(aes(group = group), method = "wilcox.test", label = "p.signif", hide.ns = TRUE, size = 2.8) +
  scale_fill_manual(values = c("normal" = "#42B7D5", "tumor" = "#E64B35")) +
  theme_bw() +
  labs(title = "A", x = NULL, y = "ssGSEA score") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(FIG6_OUTDIR, "figure_6a.png"), plot = p_a, width = 12, height = 4.4, dpi = 300)

# ----------------------------
# Panel B: correlation lollipop
# ----------------------------
if (!(FIG6_TARGET_GENE %in% rownames(expr_gene))) {
  stop("Target gene ", FIG6_TARGET_GENE, " not found in expression matrix after symbol mapping.")
}

zbtb16_expr <- as.numeric(expr_gene[FIG6_TARGET_GENE, sample_df$sample])

b_df <- data.frame(
  signature = rownames(ssgsea_scores),
  rho = NA_real_,
  p_value = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(b_df))) {
  sig <- b_df$signature[i]
  y <- as.numeric(ssgsea_scores[sig, sample_df$sample])
  ct <- suppressWarnings(cor.test(zbtb16_expr, y, method = "spearman", exact = FALSE))
  b_df$rho[i] <- unname(ct$estimate)
  b_df$p_value[i] <- ct$p.value
}

b_df <- b_df %>%
  arrange(desc(rho)) %>%
  mutate(signature = factor(signature, levels = rev(signature)))

p_b <- ggplot(b_df, aes(y = signature, x = rho, color = p_value)) +
  geom_segment(aes(x = 0, xend = rho, yend = signature), linewidth = 0.7, color = "#A0A0A0") +
  geom_point(aes(size = abs(rho)), alpha = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_gradient(low = "#E64B35", high = "#4DBBD5") +
  theme_bw() +
  labs(title = "B", x = "Correlation", y = NULL, color = "P value", size = "|Cor|")

ggsave(file.path(FIG6_OUTDIR, "figure_6b.png"), plot = p_b, width = 6.2, height = 7.4, dpi = 300)

# ----------------------------
# Panels C and D: scatter plots
# ----------------------------
signature_names <- rownames(ssgsea_scores)
mast_name <- find_signature_name(signature_names, c("Mast cells", "Mast cell"))
eos_name <- find_signature_name(signature_names, c("Eosinophils", "Eosinophil"))

scatter_df <- as.data.frame(t(ssgsea_scores))
scatter_df$sample <- rownames(scatter_df)
scatter_df <- scatter_df[sample_df$sample, , drop = FALSE]
scatter_df$zbtb16_expr <- zbtb16_expr

if (!is.na(mast_name)) {
  make_scatter(
    df = scatter_df,
    signature_name = mast_name,
    panel_label = "C",
    outfile = file.path(FIG6_OUTDIR, "figure_6c.png")
  )
} else {
  warning("Mast cell signature not found in GMT; skipped panel C.")
}

if (!is.na(eos_name)) {
  make_scatter(
    df = scatter_df,
    signature_name = eos_name,
    panel_label = "D",
    outfile = file.path(FIG6_OUTDIR, "figure_6d.png")
  )
} else {
  warning("Eosinophil signature not found in GMT; skipped panel D.")
}

write.csv(b_df, file.path(FIG6_OUTDIR, "figure_6b_correlation_table.csv"), row.names = FALSE)
cat("\nFigure 6 outputs saved to:", normalizePath(FIG6_OUTDIR), "\n")

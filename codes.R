# ============================================================
# DEG analysis on GEO datasets (GSE3268, GSE1987, GSE31547, GSE18842)
# - downloads GEO
# - (optionally) log2-transforms if needed
# - runs limma differential expression
# - makes volcano plots
# - computes overlap of up/down DEGs across datasets
#
# NOTE: You MUST set the sample groups (tumor vs normal) for each dataset.
#       I give two options:
#         A) by providing GSM IDs (most reliable)
#         B) by regex matching in sample titles/characteristics (quick, may need tweaks)
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(VennDiagram)
  library(grid)
})

# ----------------------------
# User settings
# ----------------------------
GSES <- c("GSE3268", "GSE1987", "GSE31547", "GSE18842")

# DEG thresholds (match paper style)
LOGFC_CUTOFF <- 1
ADJ_P_CUTOFF <- 0.05

OUTDIR <- "geo_deg_out"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helper: robust gene mapping to SYMBOL (if available)
# ----------------------------
extract_symbol <- function(fdata) {
  # Common column names across platforms:
  candidates <- c("Gene.symbol", "Gene Symbol", "GENE_SYMBOL", "Symbol", "SYMBOL",
                  "gene_assignment", "GENE", "Gene")
  col <- intersect(candidates, colnames(fdata))
  if (length(col) == 0) return(NULL)
  sym <- as.character(fdata[[col[1]]])
  
  # Some annotations store "gene_assignment" like "xxx // SYMBOL // ..."
  if (col[1] == "gene_assignment") {
    # Try to parse out the symbol-like token after ' // '
    sym2 <- str_split_fixed(sym, " // ", 3)[,2]
    sym2[sym2 == ""] <- NA
    sym <- sym2
  }
  
  sym <- str_trim(sym)
  sym[sym == "" | sym == "---"] <- NA
  sym
}

# ----------------------------
# Helper: choose a GEO ExpressionSet (in case multiple platforms)
# ----------------------------
get_eset <- function(gse_id) {
  gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE)
  if (length(gse) > 1) {
    # pick the platform with most samples
    ns <- sapply(gse, ncol)
    eset <- gse[[which.max(ns)]]
    message(gse_id, ": multiple platforms found; using ", annotation(eset),
            " with n=", ncol(eset), " samples.")
  } else {
    eset <- gse[[1]]
  }
  eset
}

# ----------------------------
# Helper: log2 transform if needed
# ----------------------------
maybe_log2 <- function(expr) {
  qx <- quantile(expr, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  # heuristic: if high values typical of non-log microarray intensities
  if (qx[5] > 100 || (qx[6] - qx[1] > 50 && qx[2] > 0)) {
    message("Applying log2 transform (with offset if needed).")
    expr[expr <= 0] <- NA
    return(log2(expr))
  }
  expr
}

# ----------------------------
# Helper: build groups
# Option A (recommended): provide GSM IDs explicitly
# ----------------------------
# Fill these in if you want fully reproducible grouping:
GROUPS_BY_GSE <- list(
  "GSE3268" = list(
    normal = c("GSM73386", "GSM73388", "GSM73390", "GSM73392", "GSM73394"),
    tumor  = c("GSM73387", "GSM73389", "GSM73391", "GSM73393", "GSM73395")
  )
)


# ----------------------------
# Helper: build groups
# Option B (quick): regex match from phenoData columns
# You may need to tweak patterns per dataset.
# ----------------------------
infer_groups_quick <- function(pdat) {
  # Try several fields commonly present
  fields <- intersect(c("title", "source_name_ch1", "characteristics_ch1"), colnames(pdat))
  text <- apply(pdat[, fields, drop = FALSE], 1, paste, collapse = " | ")
  text_l <- tolower(text)
  
  # Heuristics: adjust these patterns if needed per dataset
  is_tumor  <- str_detect(text_l, "tumou?r|cancer|carcinoma|adenocarcinoma|nsclc|lusc|luad")
  is_normal <- str_detect(text_l, "normal|control|non[- ]?tumou?r|adjacent|healthy")
  
  group <- rep(NA_character_, nrow(pdat))
  group[is_tumor & !is_normal]  <- "tumor"
  group[is_normal & !is_tumor]  <- "normal"
  
  # If still ambiguous, return NA and let user fix
  group
}

# ----------------------------
# limma DEG runner
# ----------------------------
run_limma_deg <- function(expr, group, fdata, gse_id) {
  stopifnot(length(group) == ncol(expr))
  keep <- !is.na(group)
  expr <- expr[, keep, drop = FALSE]
  group <- factor(group[keep], levels = c("normal", "tumor"))
  
  if (any(table(group) < 3)) {
    warning(gse_id, ": very small group size after filtering: ", paste(table(group), collapse = ", "))
  }
  
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  fit <- lmFit(expr, design)
  cont <- makeContrasts(tumor_vs_normal = tumor - normal, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))
  
  tt <- topTable(fit2, coef = "tumor_vs_normal", number = Inf, sort.by = "P") %>%
    rownames_to_column("probe_id")
  
  # Add SYMBOL if available
  sym <- extract_symbol(fdata)
  if (!is.null(sym)) {
    ann <- tibble(probe_id = rownames(fdata), SYMBOL = sym)
    tt <- tt %>% left_join(ann, by = "probe_id")
  } else {
    tt$SYMBOL <- NA_character_
  }
  
  # If multiple probes map to same SYMBOL, keep the most significant per SYMBOL
  tt_sym <- tt %>%
    filter(!is.na(SYMBOL)) %>%
    group_by(SYMBOL) %>%
    arrange(adj.P.Val, P.Value) %>%
    slice(1) %>%
    ungroup()
  
  list(full = tt, by_symbol = tt_sym)
}

# ----------------------------
# Volcano plot
# ----------------------------
plot_volcano <- function(tt_sym, gse_id) {
  df <- tt_sym %>%
    mutate(
      neglog10p = -log10(adj.P.Val),
      status = case_when(
        adj.P.Val < ADJ_P_CUTOFF & logFC >=  LOGFC_CUTOFF ~ "Up",
        adj.P.Val < ADJ_P_CUTOFF & logFC <= -LOGFC_CUTOFF ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  p <- ggplot(df, aes(x = logFC, y = neglog10p)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed") +
    geom_hline(yintercept = -log10(ADJ_P_CUTOFF), linetype = "dashed") +
    labs(
      title = paste0(gse_id, " Volcano (symbol-level)"),
      x = "log2 Fold Change (Tumor vs Normal)",
      y = "-log10(adj.P.Val)"
    ) +
    theme_minimal(base_size = 12)
  
  p
}

# ----------------------------
# Main loop
# ----------------------------
deg_sets_up <- list()
deg_sets_down <- list()
deg_tables <- list()

for (gse_id in GSES) {
  message("\n=== Processing ", gse_id, " ===")
  eset <- get_eset(gse_id)
  
  expr <- exprs(eset)
  expr <- maybe_log2(expr)
  
  pdat <- pData(eset)
  fdat <- fData(eset)
  
  # Build group labels
  if (gse_id %in% names(GROUPS_BY_GSE) && length(GROUPS_BY_GSE[[gse_id]]) > 0) {
    tumor_ids  <- GROUPS_BY_GSE[[gse_id]]$tumor
    normal_ids <- GROUPS_BY_GSE[[gse_id]]$normal
    group <- rep(NA_character_, nrow(pdat))
    group[rownames(pdat) %in% tumor_ids]  <- "tumor"
    group[rownames(pdat) %in% normal_ids] <- "normal"
    message("Groups assigned via explicit GSM lists.")
  } else {
    group <- infer_groups_quick(pdat)
    message("Groups inferred via heuristics. Please verify!")
    message("Group counts (including NA):")
    print(table(group, useNA = "ifany"))
  }
  
  # IMPORTANT: ensure at least both groups exist
  if (!all(c("tumor", "normal") %in% group)) {
    stop(
      gse_id, ": Could not infer both tumor and normal groups.\n",
      "Fix by providing GSM IDs in GROUPS_BY_GSE[[\"", gse_id, "\"]]."
    )
  }
  
  res <- run_limma_deg(expr, group, fdat, gse_id)
  tt_sym <- res$by_symbol
  deg_tables[[gse_id]] <- tt_sym
  
  # Save tables
  write.csv(tt_sym, file = file.path(OUTDIR, paste0(gse_id, "_deg_symbol.csv")), row.names = FALSE)
  
  # Volcano
  vp <- plot_volcano(tt_sym, gse_id)
  ggsave(filename = file.path(OUTDIR, paste0(gse_id, "_volcano.png")),
         plot = vp, width = 7, height = 5, dpi = 200)
  
  # DEG sets
  up_syms <- tt_sym %>%
    filter(adj.P.Val < ADJ_P_CUTOFF, logFC >= LOGFC_CUTOFF) %>%
    pull(SYMBOL) %>% unique()
  
  down_syms <- tt_sym %>%
    filter(adj.P.Val < ADJ_P_CUTOFF, logFC <= -LOGFC_CUTOFF) %>%
    pull(SYMBOL) %>% unique()
  
  deg_sets_up[[gse_id]] <- up_syms
  deg_sets_down[[gse_id]] <- down_syms
  
  message(gse_id, ": Up=", length(up_syms), " Down=", length(down_syms))
}

# ----------------------------
# Overlap (intersection) across datasets
# ----------------------------
common_up <- Reduce(intersect, deg_sets_up)
common_down <- Reduce(intersect, deg_sets_down)

writeLines(common_up, con = file.path(OUTDIR, "overlap_up_symbols.txt"))
writeLines(common_down, con = file.path(OUTDIR, "overlap_down_symbols.txt"))

message("\n=== OVERLAP RESULTS ===")
message("Common UP genes (all datasets): ", length(common_up))
message("Common DOWN genes (all datasets): ", length(common_down))

# ----------------------------
# Venn diagrams (Up and Down)
# ----------------------------
venn_plot <- function(sets, title, filename) {
  # VennDiagram draws to a device; handle up to 4 sets nicely
  png(filename, width = 900, height = 700, res = 120)
  grid.newpage()
  vd <- venn.diagram(
    x = sets,
    filename = NULL,
    category.names = names(sets),
    main = title,
    main.cex = 1.4,
    cat.cex = 1.1,
    cex = 1.1,
    margin = 0.1
  )
  grid.draw(vd)
  dev.off()
}

venn_plot(deg_sets_up, "Overlap of UP DEGs (symbol-level)", file.path(OUTDIR, "venn_up.png"))
venn_plot(deg_sets_down, "Overlap of DOWN DEGs (symbol-level)", file.path(OUTDIR, "venn_down.png"))

# ----------------------------
# Quick printouts
# ----------------------------
cat("\nTop common UP genes:\n")
print(head(common_up, 30))

cat("\nTop common DOWN genes:\n")
print(head(common_down, 30))

# ============================================================
# HOW TO FIX GROUPING (important!)
# ============================================================
# If the heuristic grouping fails, do this:
# 1) Inspect sample metadata:
#    eset <- get_eset("GSE18842")
#    View(pData(eset)[, c("title","source_name_ch1","characteristics_ch1")])
#
# 2) Manually define GSM IDs:
#    GROUPS_BY_GSE[["GSE18842"]] <- list(
#      tumor  = c("GSM....", "GSM...."),
#      normal = c("GSM....", "GSM....")
#    )
# Then rerun.
# ============================================================
# ============================================================
# GO + KEGG enrichment for overlapping DEGs (common_up/common_down)
# Uses: clusterProfiler + org.Hs.eg.db
# Outputs:
#   - CSV tables for GO BP/CC/MF and KEGG
#   - dotplots saved as PNG
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ----------------------------
# Install/load required packages (run once if needed)
# ----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_bioc <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "pathview")
pkgs_cran <- c("readr")

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
}
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(readr)
})

# ----------------------------
# Inputs: your overlap gene symbols
# Assumes you already have:
#   common_up, common_down
# If you saved them to txt earlier, load like this:
# common_up   <- readLines(file.path(OUTDIR, "overlap_up_symbols.txt"))
# common_down <- readLines(file.path(OUTDIR, "overlap_down_symbols.txt"))
# ----------------------------

# Output folder (reuse yours or set new)
OUTDIR <- "geo_deg_out"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
ENRDIR <- file.path(OUTDIR, "enrichment")
dir.create(ENRDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helper: clean symbols like "AKAP2 /// PALM2-AKAP2"
# Keep only first symbol, drop empties
# ----------------------------
clean_symbols <- function(x) {
  x <- trimws(x)
  x <- sapply(strsplit(x, " /// "), `[`, 1)
  x <- trimws(x)
  x[x != "" & !is.na(x)]
}

common_up   <- unique(clean_symbols(common_up))
common_down <- unique(clean_symbols(common_down))
common_all  <- unique(c(common_up, common_down))

# ----------------------------
# Convert SYMBOL -> ENTREZID (needed for clusterProfiler)
# ----------------------------
sym2entrez <- function(symbols) {
  bitr(symbols,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db) %>%
    distinct(SYMBOL, ENTREZID)
}

map_up   <- sym2entrez(common_up)
map_down <- sym2entrez(common_down)
map_all  <- sym2entrez(common_all)

ent_up   <- unique(map_up$ENTREZID)
ent_down <- unique(map_down$ENTREZID)
ent_all  <- unique(map_all$ENTREZID)

# Optional: set a background/universe (recommended for microarrays)
# If you have a "tt_sym" from a dataset (symbol-level table), you can use:
# universe_symbols <- unique(tt_sym$SYMBOL)
# universe_entrez  <- unique(bitr(universe_symbols, "SYMBOL","ENTREZID", org.Hs.eg.db)$ENTREZID)
# If you don't, clusterProfiler will default to all annotated genes.
universe_entrez <- NULL

# ----------------------------
# Enrichment runners
# ----------------------------
run_go <- function(entrez, ont, universe = NULL) {
  enrichGO(
    gene          = entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,           # "BP" "CC" "MF"
    universe      = universe,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
}

run_kegg <- function(entrez, universe = NULL) {
  # KEGG uses organism code 'hsa' for human
  enrichKEGG(
    gene          = entrez,
    organism      = "hsa",
    keyType       = "kegg",
    universe      = universe,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  ) %>% setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

save_enrich <- function(enrich_obj, prefix) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    message(prefix, ": no significant terms.")
    return(invisible(NULL))
  }
  df <- as.data.frame(enrich_obj)
  write_csv(df, file.path(ENRDIR, paste0(prefix, ".csv")))
  return(df)
}

save_dotplot <- function(enrich_obj, prefix, show_n = 15) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) return(invisible(NULL))
  p <- dotplot(enrich_obj, showCategory = show_n) +
    ggtitle(prefix) +
    theme_minimal(base_size = 12)
  ggsave(file.path(ENRDIR, paste0(prefix, "_dotplot.png")), p, width = 9, height = 6, dpi = 200)
}

# ----------------------------
# Run enrichment (All overlap genes)
# ----------------------------
ego_bp_all <- run_go(ent_all, "BP", universe_entrez)
ego_cc_all <- run_go(ent_all, "CC", universe_entrez)
ego_mf_all <- run_go(ent_all, "MF", universe_entrez)
ekegg_all  <- run_kegg(ent_all, universe_entrez)

save_enrich(ego_bp_all, "GO_BP_all_overlap")
save_enrich(ego_cc_all, "GO_CC_all_overlap")
save_enrich(ego_mf_all, "GO_MF_all_overlap")
save_enrich(ekegg_all,  "KEGG_all_overlap")

save_dotplot(ego_bp_all, "GO_BP_all_overlap", show_n = 15)
save_dotplot(ego_cc_all, "GO_CC_all_overlap", show_n = 15)
save_dotplot(ego_mf_all, "GO_MF_all_overlap", show_n = 15)
save_dotplot(ekegg_all,  "KEGG_all_overlap",  show_n = 15)

# ----------------------------
# (Optional) Run enrichment separately for UP and DOWN
# ----------------------------
ego_bp_up <- run_go(ent_up, "BP", universe_entrez)
ego_bp_dn <- run_go(ent_down, "BP", universe_entrez)
ekegg_up  <- run_kegg(ent_up, universe_entrez)
ekegg_dn  <- run_kegg(ent_down, universe_entrez)

save_enrich(ego_bp_up, "GO_BP_up_overlap")
save_enrich(ego_bp_dn, "GO_BP_down_overlap")
save_enrich(ekegg_up,  "KEGG_up_overlap")
save_enrich(ekegg_dn,  "KEGG_down_overlap")

save_dotplot(ego_bp_up, "GO_BP_up_overlap", show_n = 15)
save_dotplot(ego_bp_dn, "GO_BP_down_overlap", show_n = 15)
save_dotplot(ekegg_up,  "KEGG_up_overlap",  show_n = 15)
save_dotplot(ekegg_dn,  "KEGG_down_overlap", show_n = 15)

# ----------------------------
# Print quick summary
# ----------------------------
cat("\nSaved enrichment results to:\n", normalizePath(ENRDIR), "\n")

cat("\nTop KEGG terms (all overlap):\n")
if (!is.null(ekegg_all) && nrow(as.data.frame(ekegg_all)) > 0) {
  print(as.data.frame(ekegg_all)[1:min(10, nrow(as.data.frame(ekegg_all))),
                                 c("ID","Description","p.adjust","Count")])
} else {
  cat("No significant KEGG terms.\n")
}
save_barplot <- function(enrich_obj, prefix, show_n = 15) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) return(invisible(NULL))
  
  p <- barplot(enrich_obj, showCategory = show_n) +
    ggtitle(prefix) +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(ENRDIR, paste0(prefix, "_barplot.png")),
         p, width = 9, height = 6, dpi = 200)
}
save_barplot(ego_bp_all, "GO_BP_all_overlap")
save_barplot(ego_cc_all, "GO_CC_all_overlap")
save_barplot(ego_mf_all, "GO_MF_all_overlap")
save_barplot(ekegg_all,  "KEGG_all_overlap")


# ============================================================
# Panels C & D (Chord + Circle plots like the paper)
# Uses GOplot (as in many papers) + your limma logFC + your enrichment
#
# What you need:
# 1) A DEG table with columns: SYMBOL, logFC, adj.P.Val
#    (Use one dataset e.g. GSE31547 like the paper, OR use your overlap genes
#     with logFC taken from a chosen dataset.)
# 2) An enrichment result data.frame (GO/KEGG) with:
#    Category (GO term / KEGG pathway), ID, Term, adj_pval, Genes
#
# Output:
#   - chord plot PNG (Panel C style)
#   - circle plot PNG (Panel D style)
# ============================================================
library(GOplot)
library(dplyr)
library(stringr)

enr_dir <- file.path("geo_deg_out", "enrichment")

# 1) Load ONE DEG table to provide logFC (choose any GEO you analyzed)
deg_df0 <- read.csv(file.path("geo_deg_out", "GSE31547_deg_symbol.csv"), stringsAsFactors = FALSE) %>%
  transmute(ID = SYMBOL, logFC = logFC) %>%
  filter(!is.na(ID), ID != "") %>%
  distinct(ID, .keep_all = TRUE)

# 2) Load overlap gene lists you saved earlier
common_up   <- readLines(file.path("geo_deg_out", "overlap_up_symbols.txt"))
common_down <- readLines(file.path("geo_deg_out", "overlap_down_symbols.txt"))

clean_symbols <- function(x) {
  x <- trimws(x)
  x <- sapply(strsplit(x, " /// "), `[`, 1)
  trimws(x)
}
overlap_all <- unique(c(clean_symbols(common_up), clean_symbols(common_down)))

# 3) Build DEG logFC table ONLY for overlap genes (this is key)
deg_df <- deg_df0 %>% filter(ID %in% overlap_all)

# If some overlap genes are missing from this dataset, that's okay, but we need enough:
if (nrow(deg_df) < 10) {
  stop("Too few overlap genes found in the selected DEG table. Try a different dataset (e.g., GSE18842).")
}

# 4) Load enrichment results (GO BP + KEGG)
go_bp <- read.csv(file.path(enr_dir, "GO_BP_all_overlap.csv"), stringsAsFactors = FALSE)
kegg  <- read.csv(file.path(enr_dir, "KEGG_all_overlap.csv"), stringsAsFactors = FALSE)

# Choose more terms so we have enough valid ones
TOP_GO <- 15
TOP_KEGG <- 5

go_bp_top <- go_bp %>% arrange(p.adjust) %>% slice_head(n = min(TOP_GO, nrow(go_bp)))
kegg_top  <- kegg  %>% arrange(p.adjust) %>% slice_head(n = min(TOP_KEGG, nrow(kegg)))

make_terms <- function(df) {
  df %>%
    transmute(
      category = ID,
      term     = Description,
      genes    = str_replace_all(geneID, "/", ","),
      adj_pval = p.adjust
    )
}

terms_goplot <- bind_rows(
  make_terms(go_bp_top),
  if (nrow(kegg_top) > 0) make_terms(kegg_top) else NULL
) %>%
  mutate(genes = sapply(strsplit(genes, ","), function(x) paste(trimws(x), collapse = ","))) %>%
  distinct(category, .keep_all = TRUE)

# 5) Keep only terms that have >= 2 genes in deg_df (overlap + logFC)
term_hits <- sapply(terms_goplot$genes, function(gs) {
  g <- strsplit(gs, ",")[[1]] |> trimws()
  sum(g %in% deg_df$ID)
})

terms_keep <- terms_goplot %>% mutate(hits = term_hits) %>% filter(hits >= 2)

if (nrow(terms_keep) < 2) {
  stop("Still <2 valid terms after filtering. Try: (1) TOP_GO larger, (2) switch DEG source dataset, (3) set hits>=1.")
}

# 6) Build circle data
circ <- circle_dat(terms_keep, deg_df)

# 7) PANEL C: Chord plot
cats <- unique(circ$process$category)
if (length(cats) < 2) stop("circ still has <2 categories; try lowering hits threshold to 1.")

png(file.path(enr_dir, "PANEL_C_chord.png"), width = 1400, height = 1000, res = 160)
ch <- chord_dat(circ, deg_df, cats)
GOChord(ch, space = 0.02, gene.order = "logFC", gene.space = 0.25, gene.size = 3)
dev.off()

# 8) PANEL D: Circle plot
safe_nsub <- min(10, min(terms_keep$hits))
png(file.path(enr_dir, "PANEL_D_circle.png"), width = 1400, height = 1000, res = 160)
GOCircle(circ, nsub = safe_nsub)
dev.off()

cat("Saved:\n",
    file.path(enr_dir, "PANEL_C_chord.png"), "\n",
    file.path(enr_dir, "PANEL_D_circle.png"), "\n")

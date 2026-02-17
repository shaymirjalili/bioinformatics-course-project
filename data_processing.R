# ============================================================
# Data Processing Module
# Functions for differential expression and enrichment analysis
# ============================================================

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
# Extract DEG gene sets
# ----------------------------
extract_deg_sets <- function(tt_sym, logfc_cutoff, adj_p_cutoff) {
  up_syms <- tt_sym %>%
    filter(adj.P.Val < adj_p_cutoff, logFC >= logfc_cutoff) %>%
    pull(SYMBOL) %>% unique()
  
  down_syms <- tt_sym %>%
    filter(adj.P.Val < adj_p_cutoff, logFC <= -logfc_cutoff) %>%
    pull(SYMBOL) %>% unique()
  
  list(up = up_syms, down = down_syms)
}

# ----------------------------
# Calculate overlap across datasets
# ----------------------------
calculate_overlap <- function(deg_sets_up, deg_sets_down) {
  common_up <- Reduce(intersect, deg_sets_up)
  common_down <- Reduce(intersect, deg_sets_down)
  
  list(up = common_up, down = common_down)
}

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

# ----------------------------
# Save enrichment results
# ----------------------------
save_enrich <- function(enrich_obj, prefix, output_dir) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    message(prefix, ": no significant terms.")
    return(invisible(NULL))
  }
  df <- as.data.frame(enrich_obj)
  write_csv(df, file.path(output_dir, paste0(prefix, ".csv")))
  return(df)
}

# ----------------------------
# Perform enrichment analysis
# ----------------------------
perform_enrichment <- function(common_up, common_down, enrdir, universe_entrez = NULL) {
  # Clean symbols
  common_up   <- unique(clean_symbols(common_up))
  common_down <- unique(clean_symbols(common_down))
  common_all  <- unique(c(common_up, common_down))
  
  # Convert to Entrez IDs
  map_up   <- sym2entrez(common_up)
  map_down <- sym2entrez(common_down)
  map_all  <- sym2entrez(common_all)
  
  ent_up   <- unique(map_up$ENTREZID)
  ent_down <- unique(map_down$ENTREZID)
  ent_all  <- unique(map_all$ENTREZID)
  
  # Run enrichment (All overlap genes)
  ego_bp_all <- run_go(ent_all, "BP", universe_entrez)
  ego_cc_all <- run_go(ent_all, "CC", universe_entrez)
  ego_mf_all <- run_go(ent_all, "MF", universe_entrez)
  ekegg_all  <- run_kegg(ent_all, universe_entrez)
  
  # Run enrichment separately for UP and DOWN
  ego_bp_up <- run_go(ent_up, "BP", universe_entrez)
  ego_bp_dn <- run_go(ent_down, "BP", universe_entrez)
  ekegg_up  <- run_kegg(ent_up, universe_entrez)
  ekegg_dn  <- run_kegg(ent_down, universe_entrez)
  
  # Return all results
  list(
    all = list(bp = ego_bp_all, cc = ego_cc_all, mf = ego_mf_all, kegg = ekegg_all),
    up = list(bp = ego_bp_up, kegg = ekegg_up),
    down = list(bp = ego_bp_dn, kegg = ekegg_dn)
  )
}

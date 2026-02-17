# ============================================================
# Plotting Module
# Functions for creating visualizations
# ============================================================

# ----------------------------
# Volcano plot
# ----------------------------
plot_volcano <- function(tt_sym, gse_id, logfc_cutoff, adj_p_cutoff) {
  df <- tt_sym %>%
    mutate(
      neglog10p = -log10(adj.P.Val),
      status = case_when(
        adj.P.Val < adj_p_cutoff & logFC >=  logfc_cutoff ~ "Up",
        adj.P.Val < adj_p_cutoff & logFC <= -logfc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  p <- ggplot(df, aes(x = logFC, y = neglog10p, color = status)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(
      values = c("Up" = "#E31A1C", "Down" = "#1F78B4", "NS" = "grey70"),
      labels = c("Up" = "Up-regulated", "Down" = "Down-regulated", "NS" = "Not significant")
    ) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(adj_p_cutoff), linetype = "dashed", color = "grey30") +
    labs(
      title = paste0(gse_id, " Volcano (symbol-level)"),
      x = "log2 Fold Change (Tumor vs Normal)",
      y = "-log10(adj.P.Val)",
      color = "Regulation"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  
  p
}

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

# ----------------------------
# Save dotplot
# ----------------------------
save_dotplot <- function(enrich_obj, prefix, output_dir, show_n = 15) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) return(invisible(NULL))
  p <- dotplot(enrich_obj, showCategory = show_n) +
    ggtitle(prefix) +
    theme_minimal(base_size = 12)
  ggsave(file.path(output_dir, paste0(prefix, "_dotplot.png")), p, width = 9, height = 6, dpi = 200)
}

# ----------------------------
# Save barplot
# ----------------------------
save_barplot <- function(enrich_obj, prefix, output_dir, show_n = 15) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) return(invisible(NULL))
  
  p <- barplot(enrich_obj, showCategory = show_n) +
    ggtitle(prefix) +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(output_dir, paste0(prefix, "_barplot.png")),
         p, width = 9, height = 6, dpi = 200)
}

# ----------------------------
# Create chord and circle plots
# ----------------------------
create_chord_circle_plots <- function(deg_file, overlap_up_file, overlap_down_file, enr_dir) {
  # Load DEG table for logFC values
  deg_df0 <- read.csv(deg_file, stringsAsFactors = FALSE) %>%
    transmute(ID = SYMBOL, logFC = logFC) %>%
    filter(!is.na(ID), ID != "") %>%
    distinct(ID, .keep_all = TRUE)
  
  # Load overlap gene lists
  common_up   <- readLines(overlap_up_file)
  common_down <- readLines(overlap_down_file)
  
  overlap_all <- unique(c(clean_symbols(common_up), clean_symbols(common_down)))
  
  # Build DEG logFC table ONLY for overlap genes
  deg_df <- deg_df0 %>% filter(ID %in% overlap_all)
  
  if (nrow(deg_df) < 10) {
    warning("Too few overlap genes found in the selected DEG table. Skipping chord/circle plots.")
    return(invisible(NULL))
  }
  
  # Load enrichment results (GO BP + KEGG)
  go_bp <- read.csv(file.path(enr_dir, "GO_BP_all_overlap.csv"), stringsAsFactors = FALSE)
  kegg  <- read.csv(file.path(enr_dir, "KEGG_all_overlap.csv"), stringsAsFactors = FALSE)
  
  # Choose terms
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
  
  # Keep only terms that have >= 2 genes in deg_df
  term_hits <- sapply(terms_goplot$genes, function(gs) {
    g <- strsplit(gs, ",")[[1]] |> trimws()
    sum(g %in% deg_df$ID)
  })
  
  terms_keep <- terms_goplot %>% mutate(hits = term_hits) %>% filter(hits >= 2)
  
  if (nrow(terms_keep) < 2) {
    warning("Not enough valid terms for chord/circle plots. Skipping.")
    return(invisible(NULL))
  }
  
  # Build circle data
  circ <- circle_dat(terms_keep, deg_df)
  
  # PANEL C: Chord plot
  cats <- unique(circ$process$category)
  if (length(cats) >= 2) {
    ch <- chord_dat(circ, deg_df, cats)
    if (!is.null(ch) && nrow(ch) > 0) {
      png(file.path(enr_dir, "PANEL_C_chord.png"), width = 1400, height = 1000, res = 160)
      GOChord(ch, space = 0.02, gene.order = "logFC", gene.space = 0.25, gene.size = 3)
      dev.off()
    } else {
      warning("Chord data is empty. Skipping PANEL_C_chord.png plot.")
    }
  } else {
    warning("Not enough categories for chord plot. Skipping PANEL_C_chord.png plot.")
  }

  # PANEL D: Circle plot
  safe_nsub <- min(10, min(terms_keep$hits))
  if (!is.null(circ) && !is.null(circ$process) && nrow(circ$process) > 0) {
    png(file.path(enr_dir, "PANEL_D_circle.png"), width = 1400, height = 1000, res = 160)
    tryCatch({
      GOCircle(circ, nsub = safe_nsub)
    }, error = function(e) {
      warning("GOCircle failed: ", e$message)
    })
    dev.off()
  } else {
    warning("Circle data is empty. Skipping PANEL_D_circle.png plot.")
  }

  cat("Saved:\n",
      file.path(enr_dir, "PANEL_C_chord.png"), "\n",
      file.path(enr_dir, "PANEL_D_circle.png"), "\n")
}

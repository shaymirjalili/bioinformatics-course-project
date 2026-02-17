# ============================================================
# Main Analysis Script
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

# ----------------------------
# Load required libraries
# ----------------------------
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
# Install/load enrichment packages (must be before sourcing data_processing.R)
# ----------------------------

# Only install BiocManager if not already available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

pkgs_bioc <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "pathview", "GOplot")
pkgs_cran <- c("readr")

# Try to install missing packages
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing bioconductor package: ", p)
    tryCatch({
      BiocManager::install(p, update = FALSE, ask = FALSE, force = TRUE)
    }, error = function(e) {
      message("Warning: Failed to install ", p)
    }, warning = function(w) {
      message("Warning during install of ", p, ": ", w$message)
    })
  }
}

for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing CRAN package: ", p)
    tryCatch({
      install.packages(p, repos = "http://cran.r-project.org")
    }, error = function(e) {
      message("Warning: Failed to install ", p)
    })
  }
}

# Check which packages are actually available now
enrichment_available <- requireNamespace("clusterProfiler", quietly = TRUE) &&
                        requireNamespace("org.Hs.eg.db", quietly = TRUE)

if (enrichment_available && requireNamespace("enrichplot", quietly = TRUE) && 
    requireNamespace("GOplot", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(GOplot)
  })
  # Load optional packages
  if (requireNamespace("readr", quietly = TRUE)) {
    suppressPackageStartupMessages(library(readr))
  }
} else {
  message("Note: Enrichment packages not fully available. Enrichment analysis will be skipped.")
  enrichment_available <- FALSE
}

# ----------------------------
# Source module files
# ----------------------------
source("data_gathering.R")
source("data_processing.R")
source("plotting.R")

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
# Main loop: Process each dataset
# ----------------------------
deg_sets_up <- list()
deg_sets_down <- list()
deg_tables <- list()

for (gse_id in GSES) {
  message("\n=== Processing ", gse_id, " ===")
  
  # Check if CSV file already exists
  csv_file <- file.path(OUTDIR, paste0(gse_id, "_deg_symbol.csv"))
  if (file.exists(csv_file)) {
    message("CSV file already exists for ", gse_id, ". Skipping download and processing.")
    tt_sym <- read.csv(csv_file)
    deg_tables[[gse_id]] <- tt_sym
    
    # Extract DEG sets from existing data
    deg_sets <- extract_deg_sets(tt_sym, LOGFC_CUTOFF, ADJ_P_CUTOFF)
    deg_sets_up[[gse_id]] <- deg_sets$up
    deg_sets_down[[gse_id]] <- deg_sets$down
    
    message(gse_id, ": Up=", length(deg_sets$up), " Down=", length(deg_sets$down))
    next
  }
  
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
  
  # Volcano plot
  vp <- plot_volcano(tt_sym, gse_id, LOGFC_CUTOFF, ADJ_P_CUTOFF)
  ggsave(filename = file.path(OUTDIR, paste0(gse_id, "_volcano.png")),
         plot = vp, width = 7, height = 5, dpi = 200)
  
  # DEG sets
  deg_sets <- extract_deg_sets(tt_sym, LOGFC_CUTOFF, ADJ_P_CUTOFF)
  deg_sets_up[[gse_id]] <- deg_sets$up
  deg_sets_down[[gse_id]] <- deg_sets$down
  
  message(gse_id, ": Up=", length(deg_sets$up), " Down=", length(deg_sets$down))
}

# ----------------------------
# Overlap (intersection) across datasets
# ----------------------------
overlap <- calculate_overlap(deg_sets_up, deg_sets_down)
common_up <- overlap$up
common_down <- overlap$down

writeLines(common_up, con = file.path(OUTDIR, "overlap_up_symbols.txt"))
writeLines(common_down, con = file.path(OUTDIR, "overlap_down_symbols.txt"))

message("\n=== OVERLAP RESULTS ===")
message("Common UP genes (all datasets): ", length(common_up))
message("Common DOWN genes (all datasets): ", length(common_down))

# Venn diagrams
venn_plot(deg_sets_up, "Overlap of UP DEGs (symbol-level)", file.path(OUTDIR, "venn_up.png"))
venn_plot(deg_sets_down, "Overlap of DOWN DEGs (symbol-level)", file.path(OUTDIR, "venn_down.png"))

# Quick printouts
cat("\nTop common UP genes:\n")
print(head(common_up, 30))

cat("\nTop common DOWN genes:\n")
print(head(common_down, 30))

# ============================================================
# GO + KEGG enrichment for overlapping DEGs
# ============================================================

if (enrichment_available) {
  message("\nPerforming enrichment analysis...")
  
  # Create enrichment output directory
  ENRDIR <- file.path(OUTDIR, "enrichment")
  dir.create(ENRDIR, showWarnings = FALSE, recursive = TRUE)
  
  # Perform enrichment analysis
  enrich_results <- perform_enrichment(common_up, common_down, ENRDIR)
  
  # Save enrichment results
  save_enrich(enrich_results$all$bp, "GO_BP_all_overlap", ENRDIR)
  save_enrich(enrich_results$all$cc, "GO_CC_all_overlap", ENRDIR)
  save_enrich(enrich_results$all$mf, "GO_MF_all_overlap", ENRDIR)
  save_enrich(enrich_results$all$kegg, "KEGG_all_overlap", ENRDIR)
  
  save_enrich(enrich_results$up$bp, "GO_BP_up_overlap", ENRDIR)
  save_enrich(enrich_results$down$bp, "GO_BP_down_overlap", ENRDIR)
  save_enrich(enrich_results$up$kegg, "KEGG_up_overlap", ENRDIR)
  save_enrich(enrich_results$down$kegg, "KEGG_down_overlap", ENRDIR)
  
  # Save dotplots
  save_dotplot(enrich_results$all$bp, "GO_BP_all_overlap", ENRDIR)
  save_dotplot(enrich_results$all$cc, "GO_CC_all_overlap", ENRDIR)
  save_dotplot(enrich_results$all$mf, "GO_MF_all_overlap", ENRDIR)
  save_dotplot(enrich_results$all$kegg, "KEGG_all_overlap", ENRDIR)
  
  save_dotplot(enrich_results$up$bp, "GO_BP_up_overlap", ENRDIR)
  save_dotplot(enrich_results$down$bp, "GO_BP_down_overlap", ENRDIR)
  save_dotplot(enrich_results$up$kegg, "KEGG_up_overlap", ENRDIR)
  save_dotplot(enrich_results$down$kegg, "KEGG_down_overlap", ENRDIR)
  
  # Save barplots
  save_barplot(enrich_results$all$bp, "GO_BP_all_overlap", ENRDIR)
  save_barplot(enrich_results$all$cc, "GO_CC_all_overlap", ENRDIR)
  save_barplot(enrich_results$all$mf, "GO_MF_all_overlap", ENRDIR)
  save_barplot(enrich_results$all$kegg, "KEGG_all_overlap", ENRDIR)
  
  # Print summary
  cat("\nSaved enrichment results to:\n", normalizePath(ENRDIR), "\n")
  
  cat("\nTop KEGG terms (all overlap):\n")
  if (!is.null(enrich_results$all$kegg) && nrow(as.data.frame(enrich_results$all$kegg)) > 0) {
    print(as.data.frame(enrich_results$all$kegg)[1:min(10, nrow(as.data.frame(enrich_results$all$kegg))),
                                                  c("ID","Description","p.adjust","Count")])
  } else {
    cat("No significant KEGG terms.\n")
  }
  
  # Panels C & D (Chord + Circle plots)
  create_chord_circle_plots(
    deg_file = file.path(OUTDIR, "GSE31547_deg_symbol.csv"),
    overlap_up_file = file.path(OUTDIR, "overlap_up_symbols.txt"),
    overlap_down_file = file.path(OUTDIR, "overlap_down_symbols.txt"),
    enr_dir = ENRDIR
  )
} else {
  message("\nSkipping enrichment analysis (required packages not available).")
}

message("\n=== INITIAL ANALYSIS COMPLETE ===")
message("All results saved to: ", normalizePath(OUTDIR))

# ============================================================
# Run WGCNA Figure 4 Reproduction Analysis
# ============================================================
message("\n=== RUNNING WGCNA FIGURE 4 REPRODUCTION ===")
source("wgcna.R")
message("\n=== WGCNA FIGURE 4 REPRODUCTION COMPLETE ===")

message("\n=== FULL ANALYSIS PIPELINE COMPLETE ===")


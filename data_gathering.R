# ============================================================
# Data Gathering Module
# Functions for downloading and extracting data from GEO
# ============================================================

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
# Helper: clean symbols like "AKAP2 /// PALM2-AKAP2"
# Keep only first symbol, drop empties
# ----------------------------
clean_symbols <- function(x) {
  x <- trimws(x)
  x <- sapply(strsplit(x, " /// "), `[`, 1)
  x <- trimws(x)
  x[x != "" & !is.na(x)]
}

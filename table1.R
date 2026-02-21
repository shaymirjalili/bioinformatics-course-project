# ============================================================
# Consolidated Analysis Script: Table 1 + DEG Analysis (Fixed)
# ============================================================

# ----------------------------
# 1. Load Required Libraries
# ----------------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# ----------------------------
# 2. Internal Helper Functions
# ----------------------------

get_eset <- function(gse_id) {
  gse <- getGEO(gse_id, destdir = ".", getGPL = TRUE)
  return(if (length(gse) > 1) gse[[1]] else gse[[1]])
}

maybe_log2 <- function(expr) {
  qx <- as.numeric(quantile(expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (LogC) {
    expr[expr <= 0] <- NaN
    expr <- log2(expr)
    message("Applied log2 transformation.")
  }
  return(expr)
}

infer_groups_quick <- function(pdat) {
  combined_text <- tolower(apply(pdat, 1, paste, collapse = " "))
  group <- rep(NA, nrow(pdat))
  group[grepl("normal|control|healthy", combined_text)] <- "normal"
  group[grepl("tumor|cancer|adenocarcinoma|carcinoma", combined_text)] <- "tumor"
  return(group)
}

run_limma_deg <- function(expr, group, fdat, gse_id) {
  valid_idx <- !is.na(group)
  expr <- expr[, valid_idx]
  group <- group[valid_idx]
  
  design <- model.matrix(~0 + factor(group))
  colnames(design) <- c("normal", "tumor")
  
  fit <- lmFit(expr, design)
  cont.matrix <- makeContrasts(tumor-normal, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  tt <- topTable(fit2, coef=1, number=Inf)
  
  
  symbol_col <- grep("SYMBOL|Gene Symbol", colnames(fdat), ignore.case = TRUE, value = TRUE)[1]
  if(is.na(symbol_col)) symbol_col <- "ID"
  
  tt$Symbol <- as.character(fdat[rownames(tt), symbol_col])
  

  tt$Symbol <- sapply(strsplit(tt$Symbol, " /// "), `[`, 1)
  
  
  tt <- tt %>% filter(!is.na(Symbol) & Symbol != "" & !grepl("---", Symbol))
  
  
  tt_sym <- tt %>%
    group_by(Symbol) %>%
    summarise(
      logFC = logFC[which.max(abs(logFC))],
      adj.P.Val = min(adj.P.Val)
    ) %>% ungroup()
  
  return(list(by_symbol = tt_sym))
}

extract_deg_sets <- function(tt_sym, logfc, adjp) {
  up <- tt_sym %>% filter(logFC >= logfc & adj.P.Val < adjp) %>% pull(Symbol)
  down <- tt_sym %>% filter(logFC <= -logfc & adj.P.Val < adjp) %>% pull(Symbol)
  return(list(up = unique(up), down = unique(down)))
}

# ----------------------------
# 3. User Settings & Main Loop
# ----------------------------
GSES <- c("GSE3268", "GSE1987", "GSE31547", "GSE18842")

LOGFC_CUTOFF <- 0.7 
ADJ_P_CUTOFF <- 0.05
OUTDIR <- "geo_deg_out"
dir.create(OUTDIR, showWarnings = FALSE)

GROUPS_BY_GSE <- list(
  "GSE3268" = list(
    normal = c("GSM73386", "GSM73388", "GSM73390", "GSM73392", "GSM73394"),
    tumor  = c("GSM73387", "GSM73389", "GSM73391", "GSM73393", "GSM73395")
  )
)

deg_sets_up <- list()
deg_sets_down <- list()
table1_list <- list()

for (gse_id in GSES) {
  message("\n=== Processing ", gse_id, " ===")
  eset <- get_eset(gse_id)
  expr <- maybe_log2(exprs(eset))
  pdat <- pData(eset)
  fdat <- fData(eset)
  
  if (gse_id %in% names(GROUPS_BY_GSE)) {
    group <- rep(NA, nrow(pdat))
    group[rownames(pdat) %in% GROUPS_BY_GSE[[gse_id]]$tumor]  <- "tumor"
    group[rownames(pdat) %in% GROUPS_BY_GSE[[gse_id]]$normal] <- "normal"
  } else {
    group <- infer_groups_quick(pdat)
  }
  
  # Statistics for Table 1
  counts <- table(group)
  table1_list[[gse_id]] <- data.frame(
    Accession = gse_id,
    Normal = as.numeric(counts["normal"]),
    Tumor = as.numeric(counts["tumor"]),
    Total = sum(counts)
  )
  
  res <- run_limma_deg(expr, group, fdat, gse_id)
  deg_sets <- extract_deg_sets(res$by_symbol, LOGFC_CUTOFF, ADJ_P_CUTOFF)
  
  deg_sets_up[[gse_id]] <- deg_sets$up
  deg_sets_down[[gse_id]] <- deg_sets$down
  message(gse_id, " Done. Up: ", length(deg_sets$up), " Down: ", length(deg_sets$down))
}

# ----------------------------
# 4. Final Overlap & Output
# ----------------------------
common_up <- sort(Reduce(intersect, deg_sets_up))
common_down <- sort(Reduce(intersect, deg_sets_down))

message("\n=== FINAL OVERLAP RESULTS ===")
message("Common UP: ", length(common_up))
message("Common DOWN: ", length(common_down))

# ساخت جدول نهایی با پر کردن فضاهای خالی (Padding)
max_len <- max(length(common_up), length(common_down))
final_table <- data.frame(
  `Up_Regulated_DEGs` = c(common_up, rep("", max_len - length(common_up))),
  `Down_Regulated_DEGs` = c(common_down, rep("", max_len - length(common_down))),
  stringsAsFactors = FALSE
)


print(final_table)


write.csv(final_table, file = file.path(OUTDIR, "Table1_DEGs_Final.csv"), row.names = FALSE)
message("\nSuccess! Full list saved to: ", file.path(OUTDIR, "Table1_DEGs_Final.csv"))
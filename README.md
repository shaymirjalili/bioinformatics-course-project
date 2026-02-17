# Bioinformatics Course Project

## Project Structure

This project is organized into modular files for better maintainability and clarity:

### Main Files

- **`main.R`** - Main orchestration script that runs the entire analysis pipeline
- **`data_gathering.R`** - Functions for downloading and extracting data from GEO
- **`data_processing.R`** - Functions for differential expression and enrichment analysis
- **`plotting.R`** - Functions for creating all visualizations
- **`codes.R`** - Original monolithic script (kept for reference)

## File Descriptions

### 1. data_gathering.R
Contains functions for:
- Downloading GEO datasets (`get_eset`)
- Log2 transformation (`maybe_log2`)
- Gene symbol extraction (`extract_symbol`)
- Group inference from metadata (`infer_groups_quick`)
- Symbol cleaning utilities (`clean_symbols`)

### 2. data_processing.R
Contains functions for:
- Running limma differential expression analysis (`run_limma_deg`)
- Extracting DEG gene sets (`extract_deg_sets`)
- Calculating overlap across datasets (`calculate_overlap`)
- Converting symbols to Entrez IDs (`sym2entrez`)
- GO and KEGG enrichment analysis (`run_go`, `run_kegg`)
- Performing complete enrichment workflow (`perform_enrichment`)

### 3. plotting.R
Contains functions for:
- Volcano plots with color coding (`plot_volcano`)
- Venn diagrams (`venn_plot`)
- Enrichment dotplots (`save_dotplot`)
- Enrichment barplots (`save_barplot`)
- Chord and circle plots (`create_chord_circle_plots`)

### 4. main.R
Main workflow that:
1. Loads required libraries
2. Sources all module files
3. Sets up configuration and settings
4. Processes each GEO dataset
5. Calculates gene overlaps
6. Performs enrichment analysis
7. Creates all visualizations

## How to Run

### Option 1: Run the modular version (recommended)
```r
Rscript main.R
```

### Option 2: Run the original single file
```r
Rscript codes.R
```

### Option 3: Interactive in RStudio
```r
source("main.R")
```

## Requirements

### CRAN Packages
- dplyr
- ggplot2
- stringr
- tibble
- readr

### Bioconductor Packages
- GEOquery
- limma
- VennDiagram
- clusterProfiler
- org.Hs.eg.db
- enrichplot
- pathview
- GOplot

## Configuration

Edit the settings in `main.R`:

```r
GSES <- c("GSE3268", "GSE1987", "GSE31547", "GSE18842")
LOGFC_CUTOFF <- 1
ADJ_P_CUTOFF <- 0.05
OUTDIR <- "geo_deg_out"
```

### Defining Sample Groups

Option A (Recommended - Explicit GSM IDs):
```r
GROUPS_BY_GSE <- list(
  "GSE3268" = list(
    normal = c("GSM73386", "GSM73388", "GSM73390"),
    tumor  = c("GSM73387", "GSM73389", "GSM73391")
  )
)
```

Option B: Automatic inference from metadata (may need verification)

## Output Structure

```
geo_deg_out/
├── GSE3268_deg_symbol.csv       # DEG results per dataset
├── GSE3268_volcano.png          # Volcano plot per dataset
├── GSE1987_deg_symbol.csv
├── GSE1987_volcano.png
├── ...
├── overlap_up_symbols.txt       # Common up-regulated genes
├── overlap_down_symbols.txt     # Common down-regulated genes
├── venn_up.png                  # Venn diagram for up genes
├── venn_down.png                # Venn diagram for down genes
└── enrichment/
    ├── GO_BP_all_overlap.csv
    ├── GO_BP_all_overlap_dotplot.png
    ├── GO_BP_all_overlap_barplot.png
    ├── KEGG_all_overlap.csv
    ├── PANEL_C_chord.png
    └── PANEL_D_circle.png
```

## Benefits of Modular Structure

1. **Maintainability** - Easy to modify specific parts without affecting others
2. **Reusability** - Functions can be imported and used in other projects
3. **Testing** - Individual modules can be tested separately
4. **Collaboration** - Multiple people can work on different modules
5. **Clarity** - Clear separation of concerns makes code easier to understand
6. **Debugging** - Easier to isolate and fix issues

## Analysis Workflow

1. **Data Gathering**
   - Download datasets from GEO
   - Apply log2 transformation if needed
   - Extract gene symbols
   - Assign tumor/normal groups

2. **Data Processing**
   - Run limma differential expression analysis
   - Extract significant DEGs
   - Calculate overlap across datasets
   - Perform GO/KEGG enrichment

3. **Visualization**
   - Generate volcano plots (blue for down, red for up)
   - Create Venn diagrams
   - Generate enrichment plots
   - Create chord and circle plots

## Troubleshooting

If group inference fails, manually define GSM IDs:
```r
eset <- get_eset("GSE18842")
View(pData(eset)[, c("title","source_name_ch1","characteristics_ch1")])
```

Then add to GROUPS_BY_GSE in main.R.

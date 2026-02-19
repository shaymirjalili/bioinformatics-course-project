# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# ==============================================================================
# 1. GENERATE SYNTHETIC TCGA PAN-CANCER DATASET
# ==============================================================================

set.seed(123) # Ensure reproducibility

# List of cancer types observed in Figure 7A
cancer_types <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", 
                  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", 
                  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", 
                  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

# Function to generate synthetic expression data for a cancer type
generate_cancer_data <- function(type, n_tumor=50, n_normal=20, mean_diff=0) {
  # Base expression
  base_expr <- runif(1, 1, 4) 
  
  # Tumor samples
  tumor_expr <- rnorm(n_tumor, mean = base_expr + mean_diff, sd = 1.0)
  tumor_expr <- pmax(0, tumor_expr) # Ensure non-negative
  
  df_tumor <- data.frame(
    CancerType = type,
    SampleType = "Tumor",
    Expression = tumor_expr,
    PatientID = paste0(type, "_T_", 1:n_tumor)
  )
  
  # Normal samples
  normal_expr <- rnorm(n_normal, mean = base_expr, sd = 0.8)
  normal_expr <- pmax(0, normal_expr)
  
  df_normal <- data.frame(
    CancerType = type,
    SampleType = "Normal",
    Expression = normal_expr,
    PatientID = paste0(type, "_N_", 1:n_normal)
  )
  
  return(rbind(df_tumor, df_normal))
}

# Generate data for all cancer types with varying differences to mimic the plot
# Some have higher Normal (blue), some higher Tumor (red), some similar
all_data_list <- list()

for (ct in cancer_types) {
  # Randomly decide trend: 
  # -1 = Tumor lower (common in ZBTB16 according to plot, e.g. BLCA, BRCA, LUAD)
  # 0 = No difference
  # 1 = Tumor higher
  trend <- sample(c(-2, -1.5, -0.5, 0, 0.5), 1)
  
  # Generate matched-like counts for some, just group for others
  # Figure 7A shows varying boxes. 
  df <- generate_cancer_data(ct, n_tumor=sample(40:100, 1), n_normal=sample(10:50, 1), mean_diff=trend)
  all_data_list[[ct]] <- df
}

pan_cancer_data <- do.call(rbind, all_data_list)

# Make SampleType a factor with Normal first so it appears blue/left
pan_cancer_data$SampleType <- factor(pan_cancer_data$SampleType, levels = c("Normal", "Tumor"))

# Save the synthetic dataset
write.csv(pan_cancer_data, "tcga_pan_cancer_synthetic_data.csv", row.names = FALSE)
cat("Synthetic dataset created: tcga_pan_cancer_synthetic_data.csv\n")


# ==============================================================================
# 2. GENERATE FIGURE 7A: GROUPED COMPARISON
# ==============================================================================

# Define colors roughly matching the plot: Normal (Light Blue), Tumor (Red/Orange)
colors_groups <- c("Normal" = "#00BFC4", "Tumor" = "#F8766D")

p1 <- ggplot(pan_cancer_data, aes(x = CancerType, y = Expression, fill = SampleType)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = colors_groups) +
  labs(y = "The expression of ZBTB16\nLog2(TPM+1)", x = NULL, title = "A") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top",
        panel.grid.major.x = element_blank()) +
  stat_compare_means(aes(group = SampleType), label = "p.signif", method = "t.test", hide.ns = TRUE)

ggsave("Figure7A_PanCancer_Boxplot.png", plot = p1, width = 15, height = 6)
cat("Figure 7A created: Figure7A_PanCancer_Boxplot.png\n")


# ==============================================================================
# 3. GENERATE FIGURE 7B: PAIRED SAMPLES COMPARISON
# ==============================================================================

# Create Paired Data
# Figure 7B subsets specific cancers and shows lines connecting Normal -> Tumor
paired_cancers <- c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KICH", 
                    "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PCPG", "PRAD", 
                    "READ", "SARC", "SKCM", "STAD", "THCA", "THYM", "UCEC")

paired_data_list <- list()

for (ct in paired_cancers) {
  n_pairs <- sample(15:40, 1)
  
  # Base for this cancer
  base <- runif(1, 2, 5)
  diff <- sample(c(-2, -1, 0, 0.5), 1) # Generally downregulation
  
  # Generate pairs
  normal_vals <- rnorm(n_pairs, mean = base, sd = 0.5)
  tumor_vals <- normal_vals + rnorm(n_pairs, mean = diff, sd = 0.3)
  
  # Create a dataframe for this cancer's pairs
  df_pairs <- data.frame(
    PatientID = rep(paste0(ct, "_P", 1:n_pairs), 2),
    CancerType = ct,
    SampleType = rep(c("Normal", "Tumor"), each = n_pairs),
    Expression = c(normal_vals, tumor_vals)
  )
  paired_data_list[[ct]] <- df_pairs
}

paired_data <- do.call(rbind, paired_data_list)
paired_data$SampleType <- factor(paired_data$SampleType, levels = c("Normal", "Tumor"))

# Plotting Paired Data
# We use geom_point and geom_line grouped by patient
p2 <- ggplot(paired_data, aes(x = SampleType, y = Expression)) +
  # Draw lines first so points are on top
  geom_line(aes(group = PatientID), color = "gray50", alpha = 0.5) +
  geom_point(aes(color = SampleType), size = 2) +
  scale_color_manual(values = c("Normal" = "#00BFC4", "Tumor" = "#ED1C24")) +
  facet_grid(~CancerType, switch = "x") + # Split by cancer type
  theme_bw() +
  labs(y = "The expression of ZBTB16\nLog2(TPM+1)", x = NULL, title = "B") +
  theme(
    axis.text.x = element_blank(), # Hide x text inside facets, use strip text
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(0, "lines"), # Reduce space between facets
    legend.position = "top"
  ) +
  stat_compare_means(paired = TRUE, label = "p.signif", method = "t.test", 
                     comparisons = list(c("Normal", "Tumor")), 
                     label.y.npc = 0.9) 

# Note: The 'stat_compare_means' in facet_grid might be tricky to align exactly like the original 
# which acts like a categorical x-axis. 
# Alternative approach to match visual closer: Plot x=CancerType, and dodge points.
# But "Paired" requires connecting lines. The faceted approach is cleaner for paired lines in ggplot.

ggsave("Figure7B_PanCancer_Paired.png", plot = p2, width = 15, height = 6)
cat("Figure 7B created: Figure7B_PanCancer_Paired.png\n")

#-----------------
# LOAD LIBRARIES
#-----------------

library(ggpubr)
library(dplyr)
library(readr)
library(ggrepel)
library(ggplot2)
library(vctrs)
library(purrr)

#-----------------
# FUNCTIONS
#-----------------

# Load TF table
load_tf_table <- function(filepath) {
  TF_all <- read_delim(filepath, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  TF_all$TumorFraction <- as.numeric(gsub(",", ".", TF_all$TumorFraction))
  TF_all$SampleName <- gsub("_R$", "", TF_all$SampleName)
  return(TF_all)
}

# Prepare data subsets
filter_tf_data <- function(df, file_type, tech) {
  df %>% filter(file == file_type & Technology == tech)
}

# Join and arrange TF values for comparison
prepare_tf_comparison <- function(df1, df2) {
  inner_join(
    df1 %>% select(SampleName, illu_tf = TumorFraction),
    df2 %>% select(SampleName, ont_tf = TumorFraction),
    by = "SampleName"
  ) %>% arrange(SampleName)
}

# Plot correlation
plot_tf_correlation <- function(data, xvar, yvar, labelvar, title, out_file, label_size = 4, max_overlaps = 1) {
  p <- ggplot(data, aes_string(x = xvar, y = yvar, label = labelvar)) +
    geom_point(size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed", fill = "grey88", color = "grey", size = 0.5) +
    geom_text_repel(size = label_size, max.overlaps = max_overlaps, na.rm = TRUE) +
    stat_cor(method = "pearson", size = 6, label.x = min(data[[xvar]]), label.y = max(data[[yvar]])) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)
    ) +
    labs(x = "Tumor Fraction - ichorCNA (Illumina)", 
         y = "Tumor Fraction - ichorCNA (Nanopore)", 
         title = title)
  
  ggsave(p, filename = out_file, dpi = 600)
  print(p)
}

# Fit model and compute residual labels
add_residual_labels <- function(data, xvar, yvar, threshold_quantile = 0.85) {
  model <- lm(as.formula(paste(yvar, "~", xvar)), data = data)
  data$residual <- abs(resid(model))
  threshold <- quantile(data$residual, threshold_quantile)
  data$label <- ifelse(data$residual > threshold, data$SampleName, NA)
  return(data)
}

# Load .seg files and process column names
load_and_clean_seg_files <- function(path, reference_sample) {
  files <- list.files(path = path, pattern = "\\.seg$", full.names = TRUE)
  data_list <- lapply(files, read_delim, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  sample_names <- gsub(".cna.seg$", "", basename(files))
  sample_names <- gsub("_R$", "", sample_names)
  sample_names <- gsub("_downsampled$", "", sample_names)
  
  names(data_list) <- sample_names
  
  if (reference_sample %in% names(data_list)) {
    cleaned_colnames <- gsub(paste0("^", reference_sample, "_R\\."), "", colnames(data_list[[reference_sample]]))
    cleaned_colnames <- gsub(paste0("^", reference_sample, "\\."), "", cleaned_colnames)
    cleaned_colnames <- gsub(paste0("^", reference_sample, "_downsampled\\."), "", cleaned_colnames)
    
    data_list <- lapply(data_list, function(df) {
      colnames(df) <- cleaned_colnames
      return(df)
    })
  } else {
    warning("Reference sample not found")
  }
  
  combined <- bind_rows(data_list, .id = "sample_id")
  combined$chr <- as.factor(paste0("chr", combined$chr))
  return(combined)
}

# Compute per-bin correlations per sample
compute_bin_correlations <- function(merged_data) {
  merged_data %>%
    group_by(sample_id) %>%
    filter(n() >= 2) %>%
    summarise(
      R = cor.test(illu_logR, ont_logR, method = "pearson")$estimate,
      p_value = cor.test(illu_logR, ont_logR, method = "pearson")$p.value,
      .groups = "drop"
    )
}

#-----------------
# MAIN ANALYSIS
#-----------------

# Load TF data
TF_all <- load_tf_table("samples info - TF_info_all.tsv")

# Filter subsets
illu_orig <- filter_tf_data(TF_all, "original_size", "Illumina")
ont_orig  <- filter_tf_data(TF_all, "original_size", "Nanopore")
illu_down <- filter_tf_data(TF_all, "downsampled", "Illumina")
ont_down  <- filter_tf_data(TF_all, "ont_downsampled", "Nanopore")

# Update downsampled sets
illu_down <- rbind(illu_down, illu_orig %>% filter(SampleName %in% c("IMN029", "IMN041", "IMN072", "IMN3110", "IMN3179")))
ont_down  <- rbind(ont_down, ont_orig %>% filter(!SampleName %in% c("IMN029", "IMN041", "IMN072", "IMN3110", "IMN3179")))

# Prepare comparisons
my_data <- prepare_tf_comparison(illu_orig, ont_orig)

illu_down$SampleName <- gsub("_downsampled$", "", illu_down$SampleName)
ont_down$SampleName <- gsub("_downsampled$", "", ont_down$SampleName)
data_downsampled <- prepare_tf_comparison(illu_down, ont_down)

# Plot original size correlation
my_data_labeled <- add_residual_labels(my_data, "illu_tf", "ont_tf")
plot_tf_correlation(my_data, "illu_tf", "ont_tf", "SampleName", "Original Size Correlation", "plots_TFG/original_cov_more_labels.jpeg")
plot_tf_correlation(my_data_labeled, "illu_tf", "ont_tf", "label", "Original Correlation (Outliers labeled)", "plots_TFG/original_cov_less_labels.jpeg")

# Plot downsampled correlation
data_down_labeled <- add_residual_labels(data_downsampled, "illu_tf", "ont_tf")
plot_tf_correlation(data_downsampled, "illu_tf", "ont_tf", "SampleName", "Downsampled Correlation", "plots_TFG/downsampled_cov_more_labels.jpeg")
plot_tf_correlation(data_down_labeled, "illu_tf", "ont_tf", "label", "Downsampled Correlation (Outliers labeled)", "plots_TFG/downsampled_cov_less_labels.jpeg")

###-------------------------------------
# BINS CORRELATION ORIGINAL COVERAGE
###-------------------------------------

# Load and process .seg files
illumina_combined_seg <- load_and_clean_seg_files("/home/alessandra/Projects/FDP/all_cna_segs/all_illumina_original_size", "IMN029")
nanopore_combined_seg <- load_and_clean_seg_files("/home/alessandra/Projects/FDP/all_cna_segs/all_nanopore_original_size", "IMN041")

# Merge bins and compute correlation
merged_bins <- inner_join(
  illumina_combined_seg %>% select(sample_id, chr, start, end, logR, logR_Copy_Number) %>% dplyr::rename(illu_logR = logR),
  nanopore_combined_seg %>% select(sample_id, chr, start, end, logR, logR_Copy_Number) %>% dplyr::rename(ont_logR = logR),
  by = c("sample_id", "chr", "start", "end")
)

# Compute correlations per sample
correlations_per_sample <- compute_bin_correlations(merged_bins)

# Plot example correlation
ggplot(merged_bins, aes(x = illu_logR, y = ont_logR)) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "grey", size = 1) +
  geom_point(alpha = 0.2, color = "grey30", size = 2) +
  stat_cor(method = "pearson") +
  theme_classic() +
  labs(x = "Log2 ratio Illumina", y = "Log2 ratio Nanopore")

# Test for specific samples
for (sample in c("IMN4425", "IMN3179")) {
  sample_data <- merged_bins %>% filter(sample_id == sample)
  test_result <- cor.test(sample_data$illu_logR, sample_data$ont_logR, method = "pearson")
  print(paste("Sample:", sample))
  print(test_result)
}


###-------------------------------------
# BINS CORRELATION MATCHED COVERAGE
###-------------------------------------

illumina_combined_seg_down <- load_and_clean_seg_files("/home/alessandra/Projects/FDP/all_cna_segs/samples_with_similar_cov/illumina", "IMN029")
nanopore_combined_seg_down <- load_and_clean_seg_files("/home/alessandra/Projects/FDP/all_cna_segs/samples_with_similar_cov/nanopore", "IMN029")

merged_bins_down <- inner_join(
  illumina_combined_seg_down %>% select(sample_id, chr, start, end, logR, logR_Copy_Number) %>% dplyr::rename(illu_logR = logR),
  nanopore_combined_seg_down %>% select(sample_id, chr, start, end, logR, logR_Copy_Number) %>% dplyr::rename(ont_logR = logR),
  by = c("sample_id", "chr", "start", "end")
)

correlations_per_sample_down <- compute_bin_correlations(merged_bins)
View(correlations_per_sample_down)

ggplot(merged_bins, aes(x = illu_logR, y = ont_logR)) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "grey", size = 1) +
  geom_point(alpha = 0.2, color = "grey30", size = 2) +
  stat_cor(method = "pearson") +
  theme_classic() +
  labs(x = "Log2 ratio Illumina", y = "Log2 ratio Nanopore")

for (sample in c("IMN4425", "IMN3179")) {
  sample_data <- merged_bins %>% filter(sample_id == sample)
  test_result <- cor.test(sample_data$illu_logR, sample_data$ont_logR, method = "pearson")
  print(paste("Sample:", sample))
  print(test_result)
}


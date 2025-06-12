#============================================
# LOAD LIBRARIES
#============================================

library(data.table)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(paletteer)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(stringr)
library(magick)
library(patchwork)
library(ggrepel)

#============================================
# PARAMETERS & SAMPLE INFO
#============================================

bin_size <- 100000
chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

sample_info <- tibble::tibble(
  sample_id = c("IMN259", "IMN041", "IMN072", "IMN4048", "IMN3952", "IMN3852", 
                "IMN3110", "IMN4004", "IMN3179", "IMN4521", "IMN1339", 
                "IMN3549", "IMN4425", "IMN4426", "IMN4405", "IMN029", "IMN3573", 
                "IMN3469", "IMN508", "IMN4478", "IMN4541", "IMN1269", "IMN083"),
  tumor_type = c("Colorectal", "Colorectal", "Neuroendocrine", "Prostate", "Breast",
                 "Prostate", "Colorectal", "Prostate", "Colorectal", "Vater Ampoule", "Cervix",
                 "Urinary Tract", "Colorectal", "Gastric", "Gastric", "Non-Small Cell Lung Cancer",
                 "Liver", "Melanoma", "Colorectal", "Vater Ampoule", "Non-Small Cell Lung Cancer",
                 "Colorectal", "Neuroendocrine"),
  timepoint = c("EOT", "Baseline","Baseline","EOT","EOT","Baseline","Baseline","Post1","Baseline","EOT",
                "EOT","EOT", "EOT","EOT","Post1w","Baseline","EOT","Post5","Post3","Post1w","Post1w","EOT","Post1w"))



#============================================
# FUNCTIONS
#============================================

# Process file individually
process_methylation_file <- function(file, replace_name) {
  message("Processing file: ", file)
  data <- fread(file)
  sample_name <- gsub(replace_name, "", basename(file))
  data[, chr := paste0("chr", chr)]
  data <- data[chr %in% chromosomes]
  data[, sample_id := sample_name]
  data[, bin := floor(pos / bin_size) * bin_size]
  
  summary <- data[, .(methylation = sum(X) / sum(N)), by = .(sample_id, chr, bin)]
  return(summary)
}

# Bind all files summaries
summarise_all_methylation <- function(path, pattern, replace_name) {
  files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  all_summaries <- lapply(files, process_methylation_file, replace_name = replace_name)
  final_summary <- rbindlist(all_summaries)
  final_summary <- merge(final_summary, sample_info, by = "sample_id")
  return(final_summary)
}

create_methylation_matrix <- function(summary) {
  meth_matrix <- summary %>%
    unite(region, chr, bin, sep = "_") %>%
    pivot_wider(names_from = region, values_from = methylation)
  
  meth_mat <- as.matrix(meth_matrix[,-c(1,2,3)])
  rownames(meth_mat) <- meth_matrix$sample_id
  return(meth_mat)
}


generate_heatmap <- function(meth_mat, sample_info, output_file) {
  sample_info <- sample_info %>% filter(sample_id %in% rownames(meth_mat))
  sample_info <- sample_info[match(rownames(meth_mat), sample_info$sample_id), ]
  
  tumor_types <- unique(sample_info$tumor_type)
  tumor_colors <- setNames(
    colorRampPalette(paletteer::paletteer_d("RColorBrewer::Spectral"))(length(tumor_types)),
    tumor_types
  )
  
  time_points <- unique(sample_info$timepoint)
  timepoint_colors <- c(
    "Baseline" = "#AC6A9F",
    "Post1" = "#E3732A",
    "Post3" = "#F69541",
    "Post5" = "#EFB27E",
    "Post1w" = "#EC867E",
    "EOT" = "#B81840"
  )
  
  ha <- rowAnnotation(
    "Tumor Type" = sample_info$tumor_type,
    "Time Point" = sample_info$timepoint,
    col = list("Tumor Type" = tumor_colors,
               "Time Point" = timepoint_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  
  column_chr <- factor(sapply(strsplit(colnames(meth_mat), "_"), `[`, 1),
                       levels = chromosomes)
  row_split <- sample_info$tumor_type
  
  col_fun <- colorRamp2(c(0, 0.5, 1), c("#235789", "white", "#8B0000"))
  
  hm <- Heatmap(meth_mat,
                name = "6hmC",
                col = col_fun,
                left_annotation = ha,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                column_title_rot = 90,
                show_column_names = FALSE,
                show_row_names = TRUE,
                column_split = column_chr,
                row_split = row_split,
                show_row_dend = FALSE,
                heatmap_legend_param = list(
                  title = "Methylation beta value",
                  title_gp = gpar(fontsize = 14, fontface = "bold"),
                  labels_gp = gpar(fontsize = 13),
                  legend_direction = "horizontal",
                  legend_width = unit(5, "cm")
                ),
                row_names_gp = gpar(fontsize = 14),
                column_title_gp = gpar(fontsize = 14),
                row_title = NULL,
                use_raster = TRUE)
  
  lgd_time <- Legend(
    title = "Time Point",
    labels = names(timepoint_colors),
    legend_gp = gpar(fill = timepoint_colors),
    ncol = 3,
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 13)
  )
  

  lgd_tumor <- Legend(
    title = "Tumor Type",
    labels = names(tumor_colors),
    legend_gp = gpar(fill = tumor_colors),
    ncol = 4,
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 13)
  )
  
  png(output_file, width = 20, height = 9, units = "in", res = 600)
  
  draw(hm,
       heatmap_legend_side = "bottom",
       annotation_legend_list = c(lgd_tumor,lgd_time),
       annotation_legend_side = "bottom",
       merge_legends = TRUE)
  
  dev.off()
}


#============================================
# PIPELINE FUNCTION
#============================================

run_heatmap_pipeline <- function(path, output_file, pattern, replace_name) {
  message("Summarising methylation data from: ", path)
  meth_summary <- summarise_all_methylation(path, pattern, replace_name)
  
  message("Creating methylation matrix...")
  meth_mat <- create_methylation_matrix(meth_summary)
  
  message("Generating heatmap and saving to: ", output_file)
  generate_heatmap(meth_mat, sample_info, output_file)
  
  message("Done! Heatmap saved as: ", output_file)
  
  return(list(meth_summary = meth_summary, meth_mat = meth_mat))
}

#============================================
# GENERATE HEATMAP FOR 3 TYPES OF METHYLATION
#============================================

# 6mC (cg) - m
mod_6mc <- run_heatmap_pipeline("wgmeth_nanopore",
                     "heatmap_cg_m.png", "\\_wgmeth$", "\\.?_m_wgmeth$")


mod_6mc$meth_summary$beta_category <- cut(
  mod_6mc$meth_summary$methylation,
  breaks = c(-Inf, 0, 0.1, 0.5, 1),
  labels = c("0", "0-0.1", "0.1-0.5", "0.5-1")
)


table_6mc_beta_counts <- mod_6mc$meth_summary %>%
  group_by(beta_category) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100,
         label = paste0(round(percentage, 1), "%"))


#---------------------------------
# 6hmC (c) - h
mod_6hmc <- run_heatmap_pipeline("wgmeth_h_nanopore",
                      "heatmap_c_h.png", "\\_h_wgmeth$", "\\.?_h_wgmeth$")

mod_6hmc$meth_summary$beta_category <- cut(
  mod_6hmc$meth_summary$methylation,
  breaks = c(-Inf, 0, 0.1, 0.5, 1),
  labels = c("0", "0-0.1", "0.1-0.5", "0.5-1")
)

table_6hmC_beta_counts <- mod_6hmc$meth_summary %>% group_by(beta_category) %>% summarise(count = n())

#### genes 

sub_6hmc_higher = mod_6hmc$meth_summary  %>% filter(methylation >=0.5 )

bins_gr_6hmc <- GRanges(
  seqnames = sub_6hmc_higher$chr,
  ranges = IRanges(
    start = sub_6hmc_higher$bin,
    end = sub_6hmc_higher$bin + 1000000  
  ),
  sample_id = sub_6hmc_higher$sample_id,
  methylation = sub_6hmc_higher$methylation,
  tumor_type = sub_6hmc_higher$tumor_type
)

# equal names in seqlevels (chr)
seqlevels(bins_gr_6hmc) <- gsub("^chr", "", seqlevels(bins_gr))

overlaps_6hmc <- findOverlaps(bins_gr_6hmc, genes_gr)

result_df_6hmc <- data.frame(
  sample_id = bins_gr_6hmc[queryHits(overlaps_6hmc)]$sample_id,
  tumor_type = bins_gr_6hmc[queryHits(overlaps_6hmc)]$tumor_type,
  methylation = bins_gr_6hmc[queryHits(overlaps_6hmc)]$methylation,
  gene_name = genes_gr[subjectHits(overlaps_6hmc)]$gene_name
)


result_df2_6hmc <- data.frame(
  sample_id = bins_gr_6hmc[queryHits(overlaps_6hmc)]$sample_id,
  tumor_type = bins_gr_6hmc[queryHits(overlaps_6hmc)]$tumor_type,
  chr = as.character(seqnames(bins_gr_6hmc[queryHits(overlaps_6hmc)])),
  bin_start = start(bins_gr_6hmc[queryHits(overlaps_6hmc)]),
  bin_end = end(bins_gr_6hmc[queryHits(overlaps_6hmc)]),
  methylation = bins_gr_6hmc[queryHits(overlaps_6hmc)]$methylation,
  gene_name = genes_gr[subjectHits(overlaps_6hmc)]$gene_name
)
result_df2_6hmc <- na.omit(result_df2_6hmc)

write.csv(result_df2_6hmc, file = "genes_6hmc_result_df.csv", row.names = FALSE)

table(result_df2_6hmc$gene_name)

#------------------------------
# 6mA (a) - a
mod_6ma <- run_heatmap_pipeline("wgmeth_a_nanopore",
                      "heatmap_a_a.png", "\\_h_wgmeth$", "\\.?_h_wgmeth$")
mod_6ma$meth_summary$beta_category <- cut(
  mod_6ma$meth_summary$methylation,
  breaks = c(-Inf, 0, 0.1, 0.5, 1),
  labels = c("0", "0-0.1", "0.1-0.5", "0.5-1")
)

View(mod_6ma$meth_summary)

table_6ma_beta_counts <- mod_6ma$meth_summary %>% group_by(beta_category) %>% summarise(count = n())

library(rtracklayer)
library(GenomicRanges)
sub_6ma_higher = mod_6ma$meth_summary  %>% filter(methylation >=0.5 )

bins_gr <- GRanges(
  seqnames = sub_6ma_higher$chr,
  ranges = IRanges(
    start = sub_6ma_higher$bin,
    end = sub_6ma_higher$bin + 1000000  
  ),
  sample_id = sub_6ma_higher$sample_id,
  methylation = sub_6ma_higher$methylation,
  tumor_type = sub_6ma_higher$tumor_type
)


genes_gr <- import("/home/alessandra/Reference/bwa_GRCh38/Homo_sapiens.GRCh38.113.chr.gtf")
genes_gr <- genes_gr[genes_gr$type == "gene"]

# equal names in seqlevels (chr)
seqlevels(bins_gr) <- gsub("^chr", "", seqlevels(bins_gr))

overlaps <- findOverlaps(bins_gr, genes_gr)

result_df <- data.frame(
  sample_id = bins_gr[queryHits(overlaps)]$sample_id,
  tumor_type = bins_gr[queryHits(overlaps)]$tumor_type,
  methylation = bins_gr[queryHits(overlaps)]$methylation,
  gene_name = genes_gr[subjectHits(overlaps)]$gene_name
)

result_df2 <- data.frame(
  sample_id = bins_gr[queryHits(overlaps)]$sample_id,
  tumor_type = bins_gr[queryHits(overlaps)]$tumor_type,
  chr = as.character(seqnames(bins_gr[queryHits(overlaps)])),
  bin_start = start(bins_gr[queryHits(overlaps)]),
  bin_end = end(bins_gr[queryHits(overlaps)]),
  methylation = bins_gr[queryHits(overlaps)]$methylation,
  gene_name = genes_gr[subjectHits(overlaps)]$gene_name
)

write.csv(result_df, file = "genes_6ma_result_df.csv", row.names = FALSE)
write.csv(result_df2, file = "genes_6ma_result_df2.csv", row.names = FALSE)


saveRDS(mod_6ma, file = "mod_6ma.rds")
saveRDS(mod_6hmc, file = "mod_6hmc.rds")
saveRDS(mod_6mc, file = "mod_6mc.rds")


#------------------------------
## PLOT DISTRIBUTION
# labels
mod_6mc$meth_summary$methylation_type <- "6mC"
mod_6hmc$meth_summary$methylation_type <- "6hmC"
mod_6ma$meth_summary$methylation_type <- "6mA"

meth_summary_combined <- bind_rows(mod_6mc$meth_summary, mod_6hmc$meth_summary, mod_6ma$meth_summary)

custom_colors <- c("6mC" = "#D53E4FFF", "6hmC" = "#5E4FA2FF", "6mA" = "#3288BDFF")
order_mod <-  c("6mC", "6hmC", "6mA")

meth_summary_combined$methylation_type = as.factor(meth_summary_combined$methylation_type)
meth_summary_combined$methylation_type <- factor(
  meth_summary_combined$methylation_type,
  levels = order_mod
)

dist_meth <- ggplot(meth_summary_combined %>% filter(methylation_type == '6mA'), aes(x = methylation, fill = methylation_type)) +
  geom_density(aes(color = after_scale(fill)), alpha = 0.3, size = 0.5, show.legend = FALSE) +
  facet_wrap(~methylation_type, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Distribution - Beta values",
       x = "Beta value",
       y = "Density",
       fill = "Methylation mod:") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 12, face = "bold"),         
    axis.title = element_text(size = 14, face = "bold"),        
    strip.text = element_text(size = 13, face = "bold"),        
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5) 
  )

ggsave(dist_meth, filename = 'distribution_6mA.png', height = 6, width = 6, dpi = 600)





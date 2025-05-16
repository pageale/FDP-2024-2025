
#============================================
# LIBRARIES
#============================================
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(circlize)
library(ComplexHeatmap)
library(forcats)
library(karyoploteR)
library(cowplot)
library(magick)
library(grid)
library(textshape)
library(RColorBrewer)
library(paletteer)

#============================================
# READING SEGS FILES
#============================================

read_and_clean_seg <- function(seg_files_path, ref_sample, suffix_pattern) {
  seg_files <- list.files(path = seg_files_path, pattern = "\\.seg$", full.names = TRUE)
  seg_data_list <- lapply(seg_files, read_delim, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                          show_col_types = FALSE)
  
  sample_names <- gsub(".cna.seg$", "", basename(seg_files))
  sample_names <- gsub("_R$", "", sample_names)
  sample_names <- gsub("_downsampled$", "", sample_names)
  
  names(seg_data_list) <- sample_names
  
  if (ref_sample %in% names(seg_data_list)) {
    cleaned_colnames <- gsub(suffix_pattern, "", colnames(seg_data_list[[ref_sample]]))
    seg_data_list <- lapply(seg_data_list, function(df) {
      colnames(df) <- cleaned_colnames
      df
    })
  } else {
    warning(paste("Reference sample not found in", seg_files_path))
  }
  
  combined_seg <- bind_rows(seg_data_list, .id = "sample_id")
  combined_seg <- combined_seg %>%
    mutate(chr = factor(paste0("chr", chr))) %>%
    arrange(sample_id)
  
  return(combined_seg)
}

illumina_combined_seg <- read_and_clean_seg(
  "/home/alessandra/Projects/FDP/all_cna_segs/all_illumina_original_size",
  "IMN029","^IMN029_R\\.")

nanopore_combined_seg <- read_and_clean_seg(
  "/home/alessandra/Projects/FDP/all_cna_segs/all_nanopore_original_size",
  "IMN029","^IMN029\\.")

nanopore_similar_coverage <- read_and_clean_seg(
  "/home/alessandra/Projects/FDP/all_cna_segs/samples_with_similar_cov/nanopore",
  "IMN083","^IMN083\\.")

illumina_similar_coverage <- read_and_clean_seg(
  "/home/alessandra/Projects/FDP/all_cna_segs/samples_with_similar_cov/illumina",
  "IMN029","^IMN029_R\\.")

#============================================
# BINDING SEGMENTS BY BINS
#============================================

merged_bins_original_coverage <- inner_join(
  illumina_combined_seg %>% select(sample_id, chr, start, end, logR) %>% dplyr::rename(illu_logR = logR),
  nanopore_combined_seg %>% select(sample_id, chr, start, end, logR) %>% dplyr::rename(ont_logR = logR),
  by = c("sample_id", "chr", "start", "end")
)

merged_bins_downsampled <- inner_join(
  nanopore_similar_coverage %>% select(sample_id, chr, start, end, logR) %>% dplyr::rename(ont_down_logR = logR),
  illumina_similar_coverage %>% select(sample_id, chr, start, end, logR) %>% dplyr::rename(illu_down_logR = logR),
  by = c("sample_id", "chr", "start", "end")
)

# Squish logR values
merged_bins_scaled_ori <- merged_bins_original_coverage %>%
  mutate(
    illu_logR_squished = squish(illu_logR, c(-2, 2)),
    ont_logR_squished  = squish(ont_logR, c(-2, 2))
  )

merged_bins_scaled_downsampled<- merged_bins_downsampled %>%
  mutate(
    illu_logR_squished = squish(illu_down_logR, c(-2, 2)),
    ont_logR_squished  = squish(ont_down_logR, c(-2, 2))
  )

#============================================
# GENERATE HEATMAP
#============================================
crear_heatmap_cna <- function(merged_bins_scaled, output_name, show_row_names = FALSE) {
  
  ## create heatmap matrix
  heatmap_data <- merged_bins_scaled %>%
    select(sample_id, chr, start, illu_logR_squished, ont_logR_squished) %>%
    pivot_longer(cols = c(illu_logR_squished, ont_logR_squished),
                 names_to = "sequencing", values_to = "logR") %>%
    mutate(
      sequencing = ifelse(sequencing == "illu_logR_squished", "Illumina", "Nanopore"),
      chr_pos = paste0(chr, ":", start),
      sample_seq = paste0(sample_id, "_", sequencing)
    )
  
  heatmap_matrix <- heatmap_data %>%
    select(sample_seq, chr_pos, logR) %>%
    pivot_wider(names_from = chr_pos, values_from = logR) %>%
    column_to_rownames("sample_seq") %>%
    as.matrix()
  
  ## annotations
  anno_df <- data.frame(sample_seq = rownames(heatmap_matrix)) %>%
    mutate(
      sequencing = factor(gsub("^.*_", "", sample_seq), levels = c("Illumina", "Nanopore"))
    )
  
  tech_colors <- c("Illumina" = "#eb9307", "Nanopore" = "#5BB0BA")
  
  row_anno <- rowAnnotation(
    df = anno_df[, "sequencing", drop = FALSE],
    col = list(sequencing = tech_colors),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      sequencing = list(title = NULL)  # acá quitamos el título de la leyenda
    )
  )
  
  ## split columns by chr
  chr_labels <- gsub(":.*$", "", colnames(heatmap_matrix))
  col_split <- factor(chr_labels, levels = unique(chr_labels))
  
  ## split rows by sample ID
  samples <- gsub("_(Illumina|Nanopore)$", "", rownames(heatmap_matrix))
  
  ## heatmap
  ht <- Heatmap(
    heatmap_matrix,
    name = "LogR",
    col = colorRamp2(c(-2, 0, 1, 2), c("#0000FF", "white", "#FF0000", "#8B0000")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = samples,
    column_split = col_split,
    column_title_rot = 90,
    show_column_names = FALSE,
    show_row_names = show_row_names,
    row_names_gp = gpar(fontsize = 7),
    column_title_gp = gpar(fontsize = 7),
    row_names_side = "left",
    left_annotation = row_anno,
    row_title = NULL,
    heatmap_legend_param = list(direction = "horizontal")
  )
  
  ## karyoplot
  png("karyo.png", width=1550, height=150)
  kp <- plotKaryotype(plot.type=5, chromosomes=c("chr1", "chr2", "chr3",
                                                 "chr4", "chr5", "chr6", "chr7", "chr8", 
                                                 "chr9", "chr10", "chr11", "chr12", "chr13", 
                                                 "chr14", "chr15", "chr16", "chr17", "chr18", 
                                                 "chr19", "chr20", "chr21", "chr22", "chrX"))
  dev.off()
  
  karyo_image <- image_read("karyo.png")
  ht_grob <- grid.grabExpr(draw(ht, 
                                heatmap_legend_side = "bottom", 
                                annotation_legend_side = "bottom"))
  
  final_plot <- plot_grid(
    ggdraw() + draw_image(karyo_image),
    ggdraw() + draw_grob(ht_grob),
    ncol = 1,
    rel_heights = c(0.12, 1)
  )
  
  print(final_plot)
  
  ggsave(paste0(output_name, ".pdf"), final_plot, width = 16, height = 9, dpi = 300)
  ggsave(paste0(output_name, ".png"), final_plot, width = 16, height = 9, dpi = 300)
  ggsave(paste0(output_name, ".jpeg"), final_plot, width = 16, height = 9, dpi = 300)
}


#============================================
# HEATMAPS ORIGINAL AND DOWNSAMPLED
#============================================
crear_heatmap_cna(merged_bins_scaled_ori, "original_coverage_horizontal")

crear_heatmap_cna(merged_bins_scaled_downsampled, "downsampled_coverage_horizontal")


#============================================
# HEATMAPS ORIGINAL AND DOWNSAMPLED BY TUMOR TYPE
#============================================
# Colorectal	IMN259
# Colorectal	IMN041
# Neuroendrocrine	IMN072
# Colorectal	IMN603
# prostate	IMN4048
# breast	IMN3952
# prostate	IMN3852
# Colorectal	IMN3110
# prostate	IMN4004
# Colorectal	IMN3179
# vater_ampoule	IMN4521
# cervix	IMN1339
# urinary_tract	IMN3549
# Colorectal	IMN4425
# Gastric	IMN4426
# Gastric	IMN4405
# Non-Small Cell Lung Cancer	IMN029
# liver	IMN3573
# Melanoma	IMN3469
# Colorectal	IMN508
# vater_ampoule	IMN4478
# Non-Small Cell Lung Cancer	IMN4541
# Colorectal	IMN1269
# Neuroendrocrine	IMN083


crear_heatmap_cna_tumor_type <- function(merged_bins_scaled, output_name, show_row_names = FALSE) {
  
  ## create heatmap matrix
  heatmap_data <- merged_bins_scaled %>%
    select(sample_id, chr, start, illu_logR_squished, ont_logR_squished) %>%
    pivot_longer(cols = c(illu_logR_squished, ont_logR_squished),
                 names_to = "sequencing", values_to = "logR") %>%
    mutate(
      sequencing = ifelse(sequencing == "illu_logR_squished", "Illumina", "Nanopore"),
      chr_pos = paste0(chr, ":", start),
      sample_seq = paste0(sample_id, "_", sequencing)
    )
  
  heatmap_matrix <- heatmap_data %>%
    select(sample_seq, chr_pos, logR) %>%
    pivot_wider(names_from = chr_pos, values_from = logR) %>%
    column_to_rownames("sample_seq") %>%
    as.matrix()
  
  ## annotations
  tumor_info <- c(
    "IMN259" = "Colorectal", "IMN041" = "Colorectal", "IMN072" = "Neuroendrocrine", "IMN603" = "Colorectal",
    "IMN4048" = "Prostate", "IMN3952" = "Breast", "IMN3852" = "Prostate", "IMN3110" = "Colorectal",
    "IMN4004" = "Prostate", "IMN3179" = "Colorectal", "IMN4521" = "Vater Ampoule", "IMN1339" = "Cervix",
    "IMN3549" = "Urinary Tract", "IMN4425" = "Colorectal", "IMN4426" = "Gastric", "IMN4405" = "Gastric",
    "IMN029" = "Non-Small Cell Lung Cancer", "IMN3573" = "Liver", "IMN3469" = "Melanoma", "IMN508" = "Colorectal",
    "IMN4478" = "Vater Ampoule", "IMN4541" = "Non-Small Cell Lung Cancer", "IMN1269" = "Colorectal",
    "IMN083" = "Neuroendrocrine"
  )
  
  anno_df <- data.frame(sample_seq = rownames(heatmap_matrix)) %>%
    mutate(
      sequencing = factor(gsub("^.*_", "", sample_seq), levels = c("Illumina", "Nanopore")),
      sample_id = gsub("_(Illumina|Nanopore)$", "", sample_seq),
      tumor_type = tumor_info[sample_id]
    )
  
  tech_colors <- c("Illumina" = "#eb9307", "Nanopore" = "#5BB0BA")
  
  # tumor_colors <- c(
  #   "Colorectal" = "#A0C878",
  #   "Neuroendrocrine" = "#E6F5D0",
  #   "Prostate" = "#8073AC",
  #   "Breast" = "#F38C79",
  #   "Vater Ampoule" = "#276419",
  #   "Cervix" = "#C2A5CF",
  #   "Urinary Tract" = "#E3C89B",
  #   "Gastric" = "#A53860",
  #   "Non-Small Cell Lung Cancer" = "#FDE0EF",
  #   "Liver" = "#F1B6DA",
  #   "Melanoma" = "#DE77AE"
  # )
  
  tumor_types <- unique(anno_df$tumor_type)

  tumor_colors <- setNames(
    paletteer::paletteer_d("RColorBrewer::Spectral", n = length(tumor_types)),
    tumor_types
  )
  
  # Anotación de filas
  row_anno <- rowAnnotation(
    Sequencing = anno_df$sequencing,
    TumorType = anno_df$tumor_type,
    col = list(
      Sequencing = tech_colors,
      TumorType = tumor_colors
    ),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Sequencing = list(title = NULL),
      TumorType = list(title = "Tumor Type")
    )
  )

  
  ## split columns by chr
  chr_labels <- gsub(":.*$", "", colnames(heatmap_matrix))
  col_split <- factor(chr_labels, levels = unique(chr_labels))
  
  ## split rows by sample ID
  samples <- gsub("_(Illumina|Nanopore)$", "", rownames(heatmap_matrix))
  
  ## heatmap
  ht <- Heatmap(
    heatmap_matrix,
    name = "LogR",
    col = colorRamp2(c(-2,-1, 0, 1, 2), c("#235789","#3674B5", "white", "#FF0000", "#8B0000")),
    cluster_rows = T,
    show_row_dend = F,
    cluster_columns = FALSE,
    row_split = samples,
    column_split = col_split,
    column_title_rot = 90,
    show_column_names = FALSE,
    show_row_names = show_row_names,
    row_names_gp = gpar(fontsize = 7),
    column_title_gp = gpar(fontsize = 7),
    row_names_side = "left",
    left_annotation = row_anno,
    row_title = NULL,
    heatmap_legend_param = list(direction = "horizontal")
  )
  
  ## karyoplot
  png("karyo.png", width=1550, height=150)
  kp <- plotKaryotype(plot.type=5, chromosomes=c("chr1", "chr2", "chr3",
                                                 "chr4", "chr5", "chr6", "chr7", "chr8", 
                                                 "chr9", "chr10", "chr11", "chr12", "chr13", 
                                                 "chr14", "chr15", "chr16", "chr17", "chr18", 
                                                 "chr19", "chr20", "chr21", "chr22", "chrX"))
  dev.off()
  
  karyo_image <- image_read("karyo.png")
  ht_grob <- grid.grabExpr(draw(ht, 
                                heatmap_legend_side = "right", 
                                annotation_legend_side = "right"))
  
  final_plot <- plot_grid(
    ggdraw() + draw_image(karyo_image),
    ggdraw() + draw_grob(ht_grob),
    ncol = 1,
    rel_heights = c(0.12, 1)
  )
  
  print(final_plot)
  
  ggsave(paste0(output_name, ".pdf"), final_plot, width = 16, height = 9, dpi = 300)
  ggsave(paste0(output_name, ".png"), final_plot, width = 16, height = 9, dpi = 300)
  ggsave(paste0(output_name, ".jpeg"), final_plot, width = 16, height = 9, dpi = 300)
}

crear_heatmap_cna_tumor_type(merged_bins_scaled_ori, "original_coverage_horizontal_tumor")

crear_heatmap_cna_tumor_type(merged_bins_scaled_downsampled, "downsampled_coverage_horizontal_tumor")


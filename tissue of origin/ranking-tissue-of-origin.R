#-----------------------------
# LOAD LIBRARIES
#----------------------------

library(ggplot2)
library(dplyr)
library(readr)
library(RColorBrewer)
library(paletteer)

#-----------------------------
# LOAD DATA
#----------------------------

p160_correlation_illumina <- read_delim("illumina/tables/p160/p160_correlation_results.tsv", 
                                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)
p160_correlation_nanopore <- read_delim("nanopore/tables/p160/p160_correlation_results.tsv", 
                                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)

p780_correlation_illumina <- read_delim("illumina/tables/p780/p780_correlation_results.tsv", 
                                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)
p780_correlation_nanopore <- read_delim("nanopore/tables/p780/p780_correlation_results.tsv", 
                                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)

sample_order <- c("IMN072", "IMN083")
sample_order_p780 <- c("IMN3852", "IMN4004", "IMN4048")

#-----------------------------
# FUNCTION FOR RANKING THE TISSUES AND EXTRACT THE TOP 10
#----------------------------

rank_and_top10 <- function(df, sample_order, palette_name, n_top = 10) {
  ranked <- df %>%
    filter(sample_id %in% sample_order) %>%
    group_by(sample_id) %>%
    arrange(desc(abs(correlation))) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  top10 <- ranked %>%
    filter(sample_id == sample_order[1]) %>%
    arrange(rank) %>%
    dplyr::slice(1:n_top) %>%
    pull(tissue)
  
  color_palette <- paletteer::paletteer_d(palette_name, n = n_top)
  all_tissues <- unique(ranked$tissue)
  
  my_colors <- setNames(rep("#D3D3D3", length(all_tissues)), all_tissues)
  my_colors[top10] <- color_palette
  
  ranked <- ranked %>%
    mutate(tissue_label = ifelse(tissue %in% top10, tissue, ""))
  
  list(ranked = ranked, top10 = top10, colors = my_colors)
}


#-----------------------------
# FUCNTION FOR PLOTTING
#-----------------------------


plot_ranking <- function(ranked_df, top10, my_colors, title_text, sample_order) {
  
  n_tissues <- n_distinct(ranked_df$tissue)
  
  p <- ggplot(ranked_df, aes(x = sample_id, y = rank, group = tissue, color = tissue)) +
    geom_line(aes(size = tissue %in% top10), show.legend = FALSE) +
    geom_point(size = 6, show.legend = FALSE) +
    geom_text(data = subset(ranked_df, tissue %in% top10),
              aes(label = rank), color = "white", size = 2.8, vjust = 0.5) +
    scale_color_manual(values = my_colors) +
    scale_size_manual(values = c("TRUE" = 0.9, "FALSE" = 0.3)) +
    scale_y_reverse(
      breaks = 1:n_tissues,
      labels = ranked_df %>% 
        filter(sample_id == sample_order[1]) %>% 
        arrange(rank) %>% 
        pull(tissue)
    ) +
    scale_x_discrete(expand = c(0, 0.1)) +
    labs(title = title_text, x = NULL, y = "Rank") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 8),
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.border = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Add secondary axis only if there’s at least 2 samples
  if (length(sample_order) >= 2) {
    p <- p + scale_y_reverse(
      breaks = 1:n_tissues,
      labels = ranked_df %>% 
        filter(sample_id == sample_order[1]) %>% 
        arrange(rank) %>% 
        pull(tissue),
      sec.axis = sec_axis(
        trans = ~.,
        breaks = 1:n_tissues,
        labels = ranked_df %>%
          filter(sample_id == sample_order[length(sample_order)]) %>%
          arrange(rank) %>%
          pull(tissue)
      )
    ) +
      theme(
        axis.text.y = element_text(size = 12),
        axis.text.y.right = element_text(size = 12)
      )
  }
  
  return(p)
}


#-------------------------------------------------------
# P160 - Illumina
illumina_results <- rank_and_top10(p160_correlation_illumina, sample_order, "ggthemes::Classic_10")
p160_illumina_plot <- plot_ranking(illumina_results$ranked, illumina_results$top10, illumina_results$colors, 
                                   "P160 - Illumina", sample_order)

ggsave("plots/P160_Illumina.png", p160_illumina_plot, width = 6, height = 9, dpi = 600)

#-------------------------------------------------------
# P160 - Nanopore
nanopore_results <- rank_and_top10(p160_correlation_nanopore, sample_order, "ggthemes::Classic_10")
p160_nanopore_plot <- plot_ranking(nanopore_results$ranked, nanopore_results$top10, nanopore_results$colors, 
                                   "P160 - Nanopore", sample_order)

ggsave("plots/P160_Nanopore.png", p160_nanopore_plot, width = 6, height = 9, dpi = 600)

#-------------------------------------------------------
# P780 - Illumina
illumina_results_p780 <- rank_and_top10(p780_correlation_illumina, sample_order_p780, "ggthemes::Classic_10")
p780_illumina_plot <- plot_ranking(illumina_results_p780$ranked, illumina_results_p780$top10, illumina_results_p780$colors, 
                                   "P780 - Illumina", sample_order_p780)

ggsave("plots/P780_Illumina.png", p780_illumina_plot, width = 6, height = 9, dpi = 600)

#-------------------------------------------------------
# P780 - Nanopore
nanopore_results_p780 <- rank_and_top10(p780_correlation_nanopore, sample_order_p780, "ggthemes::Classic_10")
p780_nanopore_plot <- plot_ranking(nanopore_results_p780$ranked, nanopore_results_p780$top10, nanopore_results_p780$colors, 
                                   "P780 - Nanopore", sample_order_p780)

ggsave("plots/P780_Nanopore.png", p780_nanopore_plot, width = 6, height = 9, dpi = 600)

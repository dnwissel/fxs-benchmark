plot_figure_s01_j <- function(output_pdf, output_svg) {
  sample_names_6 <- c("E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3")
  config <- list(cpm_threshold = 1, cpm_sample_threshold = 3)
  anno_data <- data.frame(Condition = sample_names_6)
  rownames(anno_data) <- sample_names_6
  anno_colours <- list(`Condition` = c(
    "E3-1" = "#994f00ff", "E3-2" = "#994f00b2", "E3-3" = "#994f0066",
    "isoB11-1" = "#006CD1", "isoB11-2" = "#006cd1b2", "isoB11-3" = "#006cd166"
  ))
  n_colors <- 100
  color_palette <- colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "PuBuGn"))(n_colors)
  breaks_list <- seq(from = 0.9, to = 1.0, length.out = n_colors + 1)

  get_heatmap_row <- function(file_list, method_names) {
    plot_list <- lapply(seq_along(file_list), function(ix) {
      counts_df <- vroom::vroom(file_list[ix], show_col_types = FALSE)
      if ("transcript_id" %in% colnames(counts_df)) {
        counts_data <- counts_df %>% dplyr::select(-transcript_id)
      } else {
        counts_data <- counts_df
      }
      if ("gene_id" %in% colnames(counts_data)) {
        counts_data <- counts_data %>% dplyr::select(-gene_id)
      }
      colnames(counts_data) <- sample_names_6
      cpm_val <- edgeR::cpm(counts_data)
      mask <- rowSums(cpm_val > config$cpm_threshold) >= config$cpm_sample_threshold
      cpm_filtered <- cpm_val[mask, ]
      cormat <- round(cor(cpm_filtered, method = "pearson"), 3)
      p <- pheatmap::pheatmap(
        cormat,
        color = color_palette,
        breaks = breaks_list,
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        show_colnames = FALSE,
        show_rownames = FALSE,
        annotation_col = anno_data,
        annotation_row = anno_data,
        annotation_colors = anno_colours,
        display_numbers = TRUE,
        number_color = "black",
        main = method_names[ix],
        legend = TRUE,
        annotation_legend = TRUE,
        fontsize = 10,
        silent = TRUE
      )
      ggplotify::as.ggplot(p)
    })

    cowplot::plot_grid(plotlist = plot_list, nrow = 1, scale = 0.95)
  }

  names_bulk_heatmap <- c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish")
  pb_bulk_tx_heatmaps <- c("results/format/run_bambu_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv", "results/format/run_isoquant_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv", "results/format/run_isosceles_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv", "results/format/run_kallisto_long_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv", "results/format/run_miniquant_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv", "results/format/run_oarfish_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
  p_hm_pb_bulk <- get_heatmap_row(pb_bulk_tx_heatmaps, names_bulk_heatmap)

  panel_plot <- p_hm_pb_bulk
  ggplot2::ggsave(output_pdf, panel_plot, dpi = 600, width = 12, height = 8)
  ggplot2::ggsave(output_svg, panel_plot, dpi = 600, width = 12, height = 8)
  return(0)
}

if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

suppressPackageStartupMessages({
  library(cowplot)
  library(ggpubfigs)
  library(RColorBrewer)
  library(ggpubr)
  library(ggplot2)
  library(dplyr)
  library(ggplotify)
  library(forcats)
  library(vroom)
  library(scales)
  library(biomaRt)
  library(edgeR)
  library(limma)
  library(pheatmap)
  library(ggrepel)
  library(stringr)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_S01/panel_J.pdf"
  output_svg <- "results/figure_S01/panel_J.svg"
}

plot_figure_s01_j(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

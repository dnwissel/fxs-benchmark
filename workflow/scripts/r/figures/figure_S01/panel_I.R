plot_figure_s01_i <- function(output_pdf, output_svg) {

  names_sc <- c("Bambu", "Isosceles", "Kallisto")
  ont_sc_tx <- c("results/format_quantify_subsampled/gencode/ont/1_60000000.0/bambu/pseudobulk/transcript_counts_formatted.tsv", "results/format_quantify_subsampled/gencode/ont/1_60000000.0/isosceles/pseudobulk/transcript_counts_formatted.tsv", "results/format_quantify_subsampled/gencode/ont/1_60000000.0/kallisto/pseudobulk/transcript_counts_formatted.tsv")
  p_mds_ont_sc <- get_mds_plot(ont_sc_tx, names_sc)

  if (!exists("p_mds_ont_sc")) {
    stop("Panel object not found: p_mds_ont_sc", call. = FALSE)
  }

  panel_plot <- get("p_mds_ont_sc")
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
  output_pdf <- "results/figure_S01/panel_I.pdf"
  output_svg <- "results/figure_S01/panel_I.svg"
}

plot_figure_s01_i(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


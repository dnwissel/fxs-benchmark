plot_figure_s01_f <- function(output_pdf, output_svg) {

  names_bulk <- c("Bambu", "Isoquant", "Kallisto", "Miniquant", "Oarfish")
  pb_bulk_tx <- c("results/format/run_bambu_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv", "results/format/run_isoquant_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv", "results/format/run_kallisto_long_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv", "results/format/run_miniquant_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv", "results/format/run_oarfish_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
  p_mds_pb_bulk <- get_mds_plot(pb_bulk_tx, names_bulk)

  if (!exists("p_mds_pb_bulk")) {
    stop("Panel object not found: p_mds_pb_bulk", call. = FALSE)
  }

  panel_plot <- get("p_mds_pb_bulk")
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
  output_pdf <- "results/figure_S01/panel_F.pdf"
  output_svg <- "results/figure_S01/panel_F.svg"
}

plot_figure_s01_f(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


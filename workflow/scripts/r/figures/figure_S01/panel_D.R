plot_figure_s01_d <- function(output_pdf, output_svg) {

  names_bulk <- c("Bambu", "Isoquant", "Kallisto", "Miniquant", "Oarfish")
  pb_bulk_gene <- c("results/format/run_bambu_lr/1_15000000.0/pb/gencode/gene_counts_formatted.tsv", "results/format/run_isoquant_lr/1_15000000.0/pb/gencode/gene_counts_formatted.tsv", "results/format/run_kallisto_long_lr/1_15000000.0/pb/gencode/gene_counts_formatted.tsv", "results/format/run_miniquant_lr/1_15000000.0/pb/gencode/gene_counts_formatted.tsv", "results/format/run_oarfish_lr/1_15000000.0/pb/gencode/gene_counts_formatted.tsv")
  ont_bulk_gene <- c("results/format/run_bambu_lr/1_15000000.0/ont/gencode/gene_counts_formatted.tsv", "results/format/run_isoquant_lr/1_15000000.0/ont/gencode/gene_counts_formatted.tsv", "results/format/run_kallisto_long_lr/1_15000000.0/ont/gencode/gene_counts_formatted.tsv", "results/format/run_miniquant_lr/1_15000000.0/ont/gencode/gene_counts_formatted.tsv", "results/format/run_oarfish_lr/1_15000000.0/ont/gencode/gene_counts_formatted.tsv")
  ref_kinnex <- vroom::vroom(pb_bulk_gene[1], show_col_types = FALSE)
  mart <- useEnsembl("ensembl", "hsapiens_gene_ensembl")
  z <- getBM(c("ensembl_gene_id", "hgnc_symbol"), "ensembl_gene_id", sapply(strsplit(ref_kinnex$gene_id, "\\."), function(x) x[[1]]), mart)
  z$hgnc_symbol[which(z$hgnc_symbol == "")] <- z$ensembl_gene_id[which(z$hgnc_symbol == "")]
  gene_name_map <- data.frame(gene_id = sapply(strsplit(ref_kinnex$gene_id, "\\."), function(x) x[[1]])) %>% left_join(data.frame(gene_id = z$ensembl_gene_id, gene_name = z$hgnc_symbol), multiple = "first") %>% pull(gene_name)
  p_vol_ont_bulk <- get_volcano_plot(ont_bulk_gene, names_bulk, gene_name_map)

  if (!exists("p_vol_ont_bulk")) {
    stop("Panel object not found: p_vol_ont_bulk", call. = FALSE)
  }

  panel_plot <- get("p_vol_ont_bulk")
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
  output_pdf <- "results/figure_S01/panel_D.pdf"
  output_svg <- "results/figure_S01/panel_D.svg"
}

plot_figure_s01_d(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


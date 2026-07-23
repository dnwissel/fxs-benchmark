plot_figure_06_a <- function(output_pdf, output_svg) {

  kinnex <- vroom::vroom("results/format/run_isosceles_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  ont <- vroom::vroom("results/format/run_isosceles_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  kinnex_pseudobulk <- vroom::vroom("results/format_quantify_subsampled/gencode/pb/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  ont_pseudobulk <- vroom::vroom("results/format_quantify_subsampled/gencode/ont/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  relevant_transcripts <- c("ENST00000335327.6", "ENST00000361042.8", "ENST00000671038.1", "ENST00000496788.1")
  plt_frame <- data.frame(transcript_id = c(substr(kinnex$transcript_id[which(kinnex$transcript_id %in% relevant_transcripts)], 1, nchar(kinnex$transcript_id[which(kinnex$transcript_id %in% relevant_transcripts)]) - 2)), tech = c(rep("PB (Bulk)", 12), rep("ONT (Bulk)", 12), rep("PB (SC)", 12), rep("ONT (SC)", 12)), replicate = c(rep(c(rep("E3-1", 4), rep("E3-2", 4), rep("E3-3", 4)), 4)), cpm = c(unlist(lapply(1:3, function(x) unlist(edgeR::cpm(kinnex[, -(1:2)])[which(kinnex$transcript_id %in% relevant_transcripts), 
      ])[, x])), unlist(lapply(1:3, function(x) unlist(edgeR::cpm(ont[, -(1:2)])[which(kinnex$transcript_id %in% relevant_transcripts), ])[, x])), unlist(lapply(1:3, function(x) unlist(edgeR::cpm(kinnex_pseudobulk[, -(1:2)])[which(kinnex$transcript_id %in% relevant_transcripts), ])[, x])), unlist(lapply(1:3, function(x) unlist(edgeR::cpm(ont_pseudobulk[, -(1:2)])[which(kinnex$transcript_id %in% relevant_transcripts), ])[, x]))))
  cols <- c("#01799bff", "#e21b92ff", "#01799b80", "#e21b9280")
  plt_frame$tech <- factor(plt_frame$tech, levels = c("ONT (Bulk)", "PB (Bulk)", "ONT (SC)", "PB (SC)"))
  cpm_plot <- ggplot(plt_frame, aes(x = transcript_id, y = cpm, fill = tech)) + geom_bar(stat = "identity", position = "dodge") + facet_wrap(~replicate) + scale_fill_manual(values = cols) + theme_big_simple() + scale_y_continuous(trans = "pseudo_log") + labs(y = "CPM", x = "", fill = "", alpha = "") + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

  if (!exists("cpm_plot")) {
    stop("Panel object not found: cpm_plot", call. = FALSE)
  }

  panel_plot <- get("cpm_plot")
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
  library(magick)
  library(pdftools)
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(ggpubfigs)
  library(ggpubr)
  library(glue)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_06/panel_A.pdf"
  output_svg <- "results/figure_06/panel_A.svg"
}

plot_figure_06_a(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


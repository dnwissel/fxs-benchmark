plot_figure_04_e <- function(output_pdf, output_svg) {

  cols <- carto_pal(7, "Safe")
  sheets <- excel_sheets("results/reads_progression_updated.xlsx")
  read_numbers <- setNames(lapply(sheets, function(x) read_excel("results/reads_progression_updated.xlsx", sheet = x)), sheets)
  read_numbers <- read_numbers[["60M"]]
  quant_methods <- c("bambu", "isosceles", "kallisto", "oarfish")
  num_quant_methods_bulk <- 7
  ont_gencode_cts_list <- lapply(quant_methods, function(m) {
      path <- paste0("results/format_quantify_subsampled/gencode/ont/1_60000000.0/", m, "/pseudobulk/transcript_counts_formatted.tsv")
      gencode_cts <- read.csv(path, sep = "\t", header = TRUE)
      gencode_cts[is.na(gencode_cts)] <- 0
      return(gencode_cts)
  })
  pb_gencode_cts_list <- lapply(quant_methods, function(m) {
      path <- paste0("results/format_quantify_subsampled/gencode/pb/1_60000000.0/", m, "/pseudobulk/transcript_counts_formatted.tsv")
      gencode_cts <- read.csv(path, sep = "\t", header = TRUE)
      gencode_cts[is.na(gencode_cts)] <- 0
      return(gencode_cts)
  })
  aligned_reads_ont <- read_numbers %>% filter(Stage == "Deduplication") %>% pull(ONT)
  ont_quant_eff_df <- data.frame(Efficiency = sapply(quant_methods, function(m) {
      sum(ont_gencode_cts_list[[m]][, c(-1)])/sum(aligned_reads_ont)
  })) %>% rownames_to_column("Tool")
  ont_quant_eff_df$Tech <- "ONT"
  aligned_reads_pb <- read_numbers %>% filter(Stage == "Transcriptome alignment") %>% pull(PB)
  pb_quant_eff_df <- data.frame(Efficiency = sapply(quant_methods, function(m) {
      sum(pb_gencode_cts_list[[m]][, c(-1)])/sum(aligned_reads_pb)
  })) %>% rownames_to_column("Tool")
  pb_quant_eff_df$Tech <- "PB"
  quant_eff_df <- rbind(ont_quant_eff_df, pb_quant_eff_df)
  quant_eff_df$Tool <- factor(quant_eff_df$Tool, levels = quant_methods)
  cols <- carto_pal(num_quant_methods_bulk, "Safe")[c(1, 3, 4, 6)]
  panel_a <- ggplot(quant_eff_df, aes(x = Tool, fill = Tool, y = Efficiency)) + geom_bar(stat = "identity") + facet_wrap(~Tech) + theme_big_simple() + scale_fill_manual(values = cols) + labs(y = "Efficiency", x = "", fill = "") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

  if (!exists("panel_a")) {
    stop("Panel object not found: panel_a", call. = FALSE)
  }

  panel_plot <- get("panel_a")
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
  library(rcartocolor)
  library(readxl)
  library(scran)
  library(scuttle)
  library(spgs)
  library(tibble)
  library(tidyr)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_04/panel_E.pdf"
  output_svg <- "results/figure_04/panel_E.svg"
}

plot_figure_04_e(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


plot_figure_03_e <- function(output_pdf, output_svg) {

  path_master <- "results/format/run_bambu_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv"
  path_master <- "results/format/run_miniquant_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv"
  gencode_master <- read.csv(path_master, sep = "\t", header = TRUE)
  master_transcripts <- gencode_master$transcript_id
  quant_methods <- c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish")
  ont_gencode_cts_list <- lapply(tolower(quant_methods), function(m) {
      if (m == "kallisto") {
          path <- paste0("results/format/run_", m, "_long_lr//1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")
      }
      else {
          path <- paste0("results/format/run_", m, "_lr//1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")
      }
      gencode_cts <- read.csv(path, sep = "\t", header = TRUE)
      if (!all(master_transcripts %in% gencode_cts$transcript_id)) {
          concat_frame <- gencode_master[which(!(master_transcripts %in% gencode_cts$transcript_id)), ]
          concat_frame[, -(1:2)] <- 0
          gencode_cts <- rbind(gencode_cts, concat_frame)
          gencode_cts <- gencode_cts[match(gencode_master$transcript_id, gencode_cts$transcript_id), ] %>% tibble::remove_rownames()
      }
      gencode_cts[is.na(gencode_cts)] <- 0
      return(gencode_cts)
  })
  pb_gencode_cts_list <- lapply(tolower(quant_methods), function(m) {
      if (m == "kallisto") {
          path <- paste0("results/format/run_", m, "_long_lr//1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
      }
      else {
          path <- paste0("results/format/run_", m, "_lr//1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
      }
      gencode_cts <- read.csv(path, sep = "\t", header = TRUE)
      if (!all(master_transcripts %in% gencode_cts$transcript_id)) {
          concat_frame <- gencode_master[which(!(master_transcripts %in% gencode_cts$transcript_id)), ]
          concat_frame[, -(1:2)] <- 0
          gencode_cts <- rbind(gencode_cts, concat_frame)
          gencode_cts <- gencode_cts[match(gencode_master$transcript_id, gencode_cts$transcript_id), ] %>% tibble::remove_rownames()
      }
      gencode_cts[is.na(gencode_cts)] <- 0
      return(gencode_cts)
  })
  ont_quant_eff_df <- data.frame(Efficiency = sapply(quant_methods, function(m) {
      sum(ont_gencode_cts_list[[m]][, c(-(1:2))])/(1.5e+07 * 6)
  })) %>% rownames_to_column("Tool")
  ont_quant_eff_df$Tech <- "ONT"
  pb_quant_eff_df <- data.frame(Efficiency = sapply(quant_methods, function(m) {
      sum(pb_gencode_cts_list[[m]][, c(-(1:2))])/(1.5e+07 * 6)
  })) %>% rownames_to_column("Tool")
  pb_quant_eff_df$Tech <- "PB"
  quant_eff_df <- rbind(ont_quant_eff_df, pb_quant_eff_df)
  quant_eff_df$Tool <- factor(quant_eff_df$Tool, levels = c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish"))
  cols <- carto_pal(7, "Safe")
  quant_eff_df$Tech <- factor(quant_eff_df$Tech, levels = c("ONT", "PB"))
  panel_f <- ggplot(quant_eff_df, aes(x = Tool, fill = Tool, y = Efficiency)) + geom_bar(stat = "identity") + facet_grid(~Tech, scales = "free", space = "free_x") + theme_big_simple() + scale_fill_manual(values = cols) + labs(y = "Efficiency", x = "", fill = "") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if (!exists("panel_f")) {
    stop("Panel object not found: panel_f", call. = FALSE)
  }

  panel_plot <- get("panel_f")
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
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(ggpubfigs)
  library(ggpubr)
  library(glue)
  library(rcartocolor)
  library(tibble)
  library(tidyr)
  library(tximeta)
  library(DRIMSeq)
  library(stageR)
  library(stringr)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplotify)
  library(readxl)
  library(rtracklayer)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_03/panel_E.pdf"
  output_svg <- "results/figure_03/panel_E.svg"
}

plot_figure_03_e(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


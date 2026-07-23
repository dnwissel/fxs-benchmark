plot_figure_01_j <- function(output_pdf, output_svg) {
  read_numbers <- read_excel("results/reads_progression_updated.xlsx", sheet = "60M") %>%
    tidyr::pivot_longer(cols = c("ONT", "PB"), names_to = "Tech")

  raw_numbers <- data.frame(
    Sample = rep(unique(read_numbers$Sample), 2),
    Stage = "Raw",
    Tech = rep(c("ONT", "PB"), each = 6),
    value = 60000000
  )

  path_pb <- "results/format_quantify_subsampled/gencode/pb/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv"
  path_ont <- "results/format_quantify_subsampled/gencode/ont/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv"

  ont_cts <- colSums(read.csv(path_ont, sep = "\t", header = TRUE)[, 2:7])
  pb_cts <- colSums(read.csv(path_pb, sep = "\t", header = TRUE)[, 2:7])

  cts_numbers <- data.frame(
    Sample = rep(unique(read_numbers$Sample), 2),
    Stage = "Quantification",
    Tech = rep(c("ONT", "PB"), each = 6),
    value = c(ont_cts, pb_cts)
  )

  read_numbers_all <- rbind(read_numbers, raw_numbers, cts_numbers) %>%
    mutate(
      Stage_ordered = case_when(
        Tech == "ONT" & Stage == "Transcriptome alignment" ~ "3_Transcriptome alignment",
        Tech == "ONT" & Stage == "Deduplication" ~ "4_Deduplication",
        Tech == "PB" & Stage == "Deduplication" ~ "3_Deduplication",
        Tech == "PB" & Stage == "Transcriptome alignment" ~ "4_Transcriptome alignment",
        Stage == "Raw" ~ "1_Raw",
        Stage == "Demultiplexing" ~ "2_Demultiplexing",
        Stage == "Quantification" ~ "5_Quantification",
        TRUE ~ Stage
      ),
      Stage = factor(Stage, levels = c("Raw", "Demultiplexing", "Deduplication", "Transcriptome alignment", "Quantification")),
      Stage_ordered = factor(
        Stage_ordered,
        levels = c(
          "1_Raw", "2_Demultiplexing",
          "3_Deduplication", "3_Transcriptome alignment",
          "4_Deduplication", "4_Transcriptome alignment",
          "5_Quantification"
        )
      )
    )

  cols <- c(
    "Raw" = "#4575b4",
    "Demultiplexing" = "#91bfdb",
    "Deduplication" = "#e0f3f8",
    "Transcriptome alignment" = "#fee090",
    "Quantification" = "#fdae61"
  )

  panel_j <- ggplot(read_numbers_all, aes(x = Sample, y = value, fill = Stage)) +
    geom_col(position = position_dodge(width = 0.9), color = "black", aes(group = Stage_ordered)) +
    facet_wrap(~ Tech, scales = "free_x", strip.position = "top") +
    theme_big_simple() +
    scale_fill_manual(values = cols) +
    labs(x = NULL, y = "Number of reads", fill = "")

  ggsave(output_pdf, panel_j, dpi = 600, width = 12, height = 6)
  ggsave(output_svg, panel_j, dpi = 600, width = 12, height = 6)
  return(0)
}

if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggpubfigs)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_01/panel_J.pdf"
  output_svg <- "results/figure_01/panel_J.svg"
}

status <- plot_figure_01_j(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

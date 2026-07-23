plot_figure_01_i <- function(output_pdf, output_svg) {
  genome_df <- vroom::vroom("results/qc/genome_alignment_efficiency_benchmark.csv")
  trans_df <- vroom::vroom("results/qc/transcriptome_alignment_efficiency_benchmark.csv")

  plot_data <- genome_df %>%
    dplyr::select(tech, type, sample, Raw_Total = total_reads, Genome_Mapped = primary_mapped) %>%
    inner_join(
      trans_df %>% dplyr::select(tech, type, sample, Transcriptome_Mapped = primary_mapped),
      by = c("tech", "type", "sample")
    ) %>%
    mutate(
      Raw = 1000000,
      Genome = (Genome_Mapped / Raw_Total) * 1000000,
      Transcriptome = (Transcriptome_Mapped / Raw_Total) * 1000000
    ) %>%
    pivot_longer(cols = c(Raw, Genome, Transcriptome), names_to = "Stage", values_to = "ReadCount") %>%
    mutate(
      Stage = factor(Stage, levels = c("Raw", "Genome", "Transcriptome")),
      group_id = paste(tech, type, sample, sep = "_")
    )

  plot_data$label <- case_when(
    plot_data$tech == "ont" & plot_data$type == "single_cell" ~ "ONT (SC)",
    plot_data$tech == "ont" & plot_data$type != "single_cell" ~ "ONT (Bulk)",
    plot_data$tech == "pb" & plot_data$type == "single_cell" ~ "PB (SC)",
    plot_data$tech == "pb" & plot_data$type != "single_cell" ~ "PB (Bulk)"
  )

  plot_data$label <- factor(
    plot_data$label,
    levels = c("ONT (Bulk)", "PB (Bulk)", "ONT (SC)", "PB (SC)")
  )

  panel_i <- ggplot(plot_data, aes(x = ReadCount, y = sample, color = Stage)) +
    geom_line(aes(group = group_id), color = "grey70", size = 0.8) +
    geom_point(size = 3) +
    facet_wrap(~ label, nrow = 1) +
    scale_x_reverse() +
    scale_x_continuous(labels = scales::label_number(scale = 1e-6)) +
    scale_color_manual(values = c("Raw" = "#2A9D8F", "Genome" = "#E9C46A", "Transcriptome" = "#E76F51")) +
    labs(x = "Primary alignments (M)", y = "Replicate", color = "Alignment stage") +
    theme_big_simple()

  ggsave(output_pdf, panel_i, dpi = 600, width = 14, height = 6)
  ggsave(output_svg, panel_i, dpi = 600, width = 14, height = 6)
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
  library(vroom)
  library(ggpubfigs)
  library(scales)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_01/panel_I.pdf"
  output_svg <- "results/figure_01/panel_I.svg"
}

status <- plot_figure_01_i(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

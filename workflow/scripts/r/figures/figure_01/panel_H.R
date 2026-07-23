plot_figure_01_h <- function(output_pdf, output_svg) {
  priming_data <- vroom::vroom("results/qc/internal_priming_benchmark.csv")
  priming_data$label <- case_when(
    priming_data$tech == "ont" & priming_data$type == "single_cell" ~ "ONT (SC)",
    priming_data$tech == "ont" & priming_data$type != "single_cell" ~ "ONT (Bulk)",
    priming_data$tech == "pb" & priming_data$type == "single_cell" ~ "PB (SC)",
    priming_data$tech == "pb" & priming_data$type != "single_cell" ~ "PB (Bulk)"
  )

  tech_colors <- c(
    "ONT (Bulk)" = "#01799bff",
    "PB (Bulk)" = "#e21b92ff",
    "ONT (SC)" = "#01799b80",
    "PB (SC)" = "#e21b9280"
  )

  panel_h <- ggplot(priming_data, aes(x = label, y = proportion, fill = label)) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), width = 0.7) +
    stat_summary(
      fun.data = mean_sdl,
      fun.args = list(mult = 1),
      geom = "errorbar",
      position = position_dodge(width = 0.8),
      width = 0.2
    ) +
    geom_point(
      position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1),
      color = "black",
      alpha = 0.6
    ) +
    theme_big_simple() +
    scale_fill_manual(values = tech_colors) +
    labs(y = "Internally primed reads", x = "")

  ggsave(output_pdf, panel_h, dpi = 600, width = 8, height = 6)
  ggsave(output_svg, panel_h, dpi = 600, width = 8, height = 6)
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
  library(vroom)
  library(ggpubfigs)
  library(ggpubr)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_01/panel_H.pdf"
  output_svg <- "results/figure_01/panel_H.svg"
}

status <- plot_figure_01_h(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

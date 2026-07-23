plot_figure_06_c <- function(output_pdf, output_svg) {

  pb_track_plot <- ggdraw() + draw_image(magick::image_read_pdf("results/track-division-pb-colored.pdf", density = 300))

  if (!exists("pb_track_plot")) {
    stop("Panel object not found: pb_track_plot", call. = FALSE)
  }

  panel_plot <- get("pb_track_plot")
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
  output_pdf <- "results/figure_06/panel_C.pdf"
  output_svg <- "results/figure_06/panel_C.svg"
}

plot_figure_06_c(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


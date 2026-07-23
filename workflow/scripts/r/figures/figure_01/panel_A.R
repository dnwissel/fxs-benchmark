plot_figure_01_a <- function(output_pdf, output_svg, design_pdf) {
  panel_a <- ggdraw() +
    draw_image(magick::image_read_pdf(design_pdf, density = 300))

  ggsave(output_pdf, panel_a, dpi = 600, width = 10, height = 10)
  ggsave(output_svg, panel_a, dpi = 600, width = 10, height = 10)
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
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  design_pdf <- if ("design_pdf" %in% names(snakemake@input)) snakemake@input[["design_pdf"]] else "resources/roche_fxs_experimental_design.pdf"
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  design_pdf <- "resources/roche_fxs_experimental_design.pdf"
  output_pdf <- "results/figure_01/panel_A.pdf"
  output_svg <- "results/figure_01/panel_A.svg"
}

status <- plot_figure_01_a(output_pdf = output_pdf, output_svg = output_svg, design_pdf = design_pdf)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

safe_cbind <- function(data_list, type, tech) {
  df <- data.frame(V1 = unlist(data_list))
  df$type <- type
  df$tech <- tech
  df
}

plot_figure_01_b <- function(output_pdf, output_svg, depths_json) {
  params <- rjson::fromJSON(file = depths_json)
  sample_replicates <- c("E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3")

  plt_frame_depth <- data.frame(
    rbind(
      safe_cbind(params$bulk$quant_gencode$illumina, "GENCODE (Counts)", "Illmn (Bulk)"),
      safe_cbind(params$bulk$quant_gencode$ont, "GENCODE (Counts)", "ONT (Bulk)"),
      safe_cbind(params$bulk$quant_gencode$pb, "GENCODE (Counts)", "PB (Bulk)"),
      safe_cbind(params$bulk$quant_sirv$illumina, "SIRV (Counts)", "Illmn (Bulk)"),
      safe_cbind(params$bulk$quant_sirv$ont, "SIRV (Counts)", "ONT (Bulk)"),
      safe_cbind(params$bulk$quant_sirv$pb, "SIRV (Counts)", "PB (Bulk)"),
      safe_cbind(params$bulk$quant_ercc$illumina, "ERCC (Counts)", "Illmn (Bulk)"),
      safe_cbind(params$bulk$quant_ercc$ont, "ERCC (Counts)", "ONT (Bulk)"),
      safe_cbind(params$bulk$quant_ercc$pb, "ERCC (Counts)", "PB (Bulk)"),
      safe_cbind(params$bulk$raw$illumina, "Raw reads", "Illmn (Bulk)"),
      safe_cbind(params$bulk$raw$ont, "Raw reads", "ONT (Bulk)"),
      safe_cbind(params$bulk$raw$pb, "Raw reads", "PB (Bulk)"),
      safe_cbind(params$single_cell$raw$illumina, "Raw reads", "Illmn (SC)"),
      safe_cbind(params$single_cell$raw$ont, "Raw reads", "ONT (SC)"),
      safe_cbind(params$single_cell$raw$pb, "Raw reads", "PB (SC)"),
      safe_cbind(params$single_cell$quant$illumina, "GENCODE (Counts)", "Illmn (SC)"),
      safe_cbind(params$single_cell$quant$ont, "GENCODE (Counts)", "ONT (SC)"),
      safe_cbind(params$single_cell$quant$pb, "GENCODE (Counts)", "PB (SC)")
    )
  )

  plt_frame_depth$rep <- rep(sample_replicates, nrow(plt_frame_depth) / length(sample_replicates))
  plt_frame_depth$V1 <- as.numeric(plt_frame_depth$V1)
  plt_frame_depth$type <- factor(
    plt_frame_depth$type,
    levels = c("Raw reads", "GENCODE (Counts)", "SIRV (Counts)", "ERCC (Counts)")
  )

  tech_colors <- c(
    "Illmn (Bulk)" = "#f8a21fff",
    "ONT (Bulk)" = "#01799bff",
    "PB (Bulk)" = "#e21b92ff",
    "Illmn (SC)" = "#f8a21f80",
    "ONT (SC)" = "#01799b80",
    "PB (SC)" = "#e21b9280"
  )

  panel_b <- ggplot(plt_frame_depth, aes(x = V1, y = fct_rev(rep))) +
    geom_line(aes(group = rep), color = "grey", linewidth = 1.5, alpha = 0.7) +
    geom_point(aes(color = tech), size = 4) +
    facet_wrap(~ type, scales = "free_x", nrow = 1) +
    scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
    scale_color_manual(values = tech_colors) +
    labs(x = "Number of reads (M)", y = "", color = "Platform") +
    theme_big_simple() +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())

  ggsave(output_pdf, panel_b, dpi = 600, width = 16, height = 6)
  ggsave(output_svg, panel_b, dpi = 600, width = 16, height = 6)
  return(0)
}

if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubfigs)
  library(forcats)
  library(scales)
  library(rjson)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  depths_json <- if ("depths_json" %in% names(snakemake@input)) snakemake@input[["depths_json"]] else "results/qc/roche_depths.json"
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  depths_json <- "results/qc/roche_depths.json"
  output_pdf <- "results/figure_01/panel_B.pdf"
  output_svg <- "results/figure_01/panel_B.svg"
}

status <- plot_figure_01_b(output_pdf = output_pdf, output_svg = output_svg, depths_json = depths_json)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

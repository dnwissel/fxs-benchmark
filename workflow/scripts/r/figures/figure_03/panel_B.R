plot_figure_03_b <- function(output_pdf, output_svg) {

  cols <- carto_pal(7, "Safe")
  num_quant_methods_bulk <- 7
  cols <- carto_pal(num_quant_methods_bulk, "Safe")
  cols <- carto_pal(7, "Safe")
  plt_frame <- data.frame()
  for (depth in c("1000000.0", "2500000.0", "5000000.0", "10000000.0", "15000000.0")) {
      for (method in c("run_oarfish_lr", "run_miniquant_lr", "run_isoquant_lr", "run_isosceles_lr", "run_kallisto_long_lr", "run_bambu_lr")) {
          for (type in c("ont", "pb")) {
              for (rep in c("E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3")) {
                  plt_frame <- rbind(plt_frame, data.frame(method = method, type = type, rep = rep, depth = as.numeric(depth), time = vroom::vroom(paste0("results/benchmarks/quantify_bulk_downsampled/", method, "/", "1_", depth, "/gencode/", rep, "/", type, "/", rep, ".log"))$s, mem = vroom::vroom(paste0("results/benchmarks/quantify_bulk_downsampled/", method, "/", "1_", depth, "/gencode/", rep, "/", type, "/", rep, ".log"))$max_uss))
              }
          }
      }
  }
  plt_frame$type <- ifelse(plt_frame$type == "ont", "ONT", "PB")
  plt_frame$method <- ifelse(plt_frame$method == "run_bambu_lr", "Bambu", ifelse(plt_frame$method == "run_isoquant_lr", "Isoquant", ifelse(plt_frame$method == "run_kallisto_long_lr", "Kallisto", ifelse(plt_frame$method == "run_oarfish_lr", "Oarfish", ifelse(plt_frame$method == "run_miniquant_lr", "Miniquant", ifelse(plt_frame$method == "run_isosceles_lr", "Isosceles", "Salmon"))))))
  plt_frame$method <- factor(plt_frame$method, levels = c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish", "Salmon"))
  plt_frame$type <- factor(plt_frame$type, levels = c("ONT", "PB"))
  plt_frame$depth <- factor(sprintf("%.1e", plt_frame$depth), levels = c("1.0e+06", "2.5e+06", "5.0e+06", "1.0e+07", "1.5e+07"))
  time <- plt_frame %>% ggplot(aes(x = depth, y = time/60, color = method, group = method)) + geom_point(stat = "summary", fun = mean, size = 4) + stat_summary(fun.y = mean, geom = "line", linewidth = 2, lty = 2) + stat_summary(fun.data = mean_sdl, geom = "pointrange", linewidth = 1) + theme_big_simple() + labs(fill = "", y = "Runtime (m)", x = "Subsampled reads", color = "") + scale_color_manual(values = cols) + facet_wrap(~type) + scale_y_log10()

  if (!exists("time")) {
    stop("Panel object not found: time", call. = FALSE)
  }

  panel_plot <- get("time")
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
  output_pdf <- "results/figure_03/panel_B.pdf"
  output_svg <- "results/figure_03/panel_B.svg"
}

plot_figure_03_b(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


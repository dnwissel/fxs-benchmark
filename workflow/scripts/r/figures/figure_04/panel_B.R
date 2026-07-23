plot_figure_04_b <- function(output_pdf, output_svg) {

  cols <- carto_pal(7, "Safe")
  plt_frame <- data.frame()
  for (depth in c("10000000.0", "20000000.0", "30000000.0", "60000000.0")) {
      for (method in c("quantify_sc_run_bambu", "quantify_sc_run_isosceles", "quantify_kallisto_sc_lr", "quantify_oarfish_sc_lr")) {
          for (type in c("ont", "pb")) {
              for (rep in c("E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3")) {
                  plt_frame <- rbind(plt_frame, data.frame(method = method, type = type, rep = rep, depth = as.numeric(depth), time = vroom::vroom(paste0("results/timing_single_cell/", method, "/", type, "/", "1_", depth, "/", rep, ".log"))$s, mem = vroom::vroom(paste0("results/timing_single_cell/", method, "/", type, "/", "1_", depth, "/", rep, ".log"))$max_uss))
              }
          }
      }
  }
  plt_frame$type <- ifelse(plt_frame$type == "ont", "ONT", "PB")
  plt_frame$method <- ifelse(plt_frame$method == "quantify_kallisto_sc_lr", "Kallisto", ifelse(plt_frame$method == "quantify_oarfish_sc_lr", "Oarfish", ifelse(plt_frame$method == "quantify_sc_run_bambu", "Bambu", "Isosceles")))
  plt_frame$method <- factor(plt_frame$method, levels = c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish", "Salmon"))
  plt_frame$type <- factor(plt_frame$type, levels = c("ONT", "PB"))
  plt_frame$depth <- factor(sprintf("%.1e", plt_frame$depth), levels = c("1.0e+07", "2.0e+07", "3.0e+07", "6.0e+07"))
  memory <- plt_frame %>% ggplot(aes(x = depth, y = mem/1000, color = method, group = method)) + geom_point(stat = "summary", fun = mean, size = 4) + stat_summary(fun.y = mean, geom = "line", linewidth = 2, lty = 2) + stat_summary(fun.data = mean_sdl, geom = "pointrange", linewidth = 1) + theme_big_simple() + labs(fill = "", y = "Memory (USS GB)", x = "Subsampled reads", color = "") + scale_color_manual(values = cols[c(1, 3, 4, 6)]) + facet_wrap(~type) + scale_y_log10()

  if (!exists("memory")) {
    stop("Panel object not found: memory", call. = FALSE)
  }

  panel_plot <- get("memory")
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
  output_pdf <- "results/figure_04/panel_B.pdf"
  output_svg <- "results/figure_04/panel_B.svg"
}

plot_figure_04_b(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


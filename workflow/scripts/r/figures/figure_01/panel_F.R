plot_figure_01_f <- function(output_pdf, output_svg) {
  sample_replicates <- c("E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3")
  plt_frame_qc <- data.frame()

  for (type in c("bulk", "single_cell")) {
    for (day in sample_replicates) {
      for (tech in c("ont", "pb")) {
        tmp_frame <- vroom::vroom(
          paste0("results/qc/qc_prepare_sampled_read_quality_frame/", tech, "/1000000.0/", type, "/", day, ".tsv")
        )
        plt_frame_qc <- rbind(
          plt_frame_qc,
          data.frame(
            tech = paste0(tech, " (", type, ")"),
            sample = day,
            snv_rate = tmp_frame$snv_rate,
            indel_rate = tmp_frame$indel_rate
          )
        )
      }
    }
  }

  for (day in sample_replicates) {
    tmp_bulk <- vroom::vroom(
      paste0("results/qc/qc_prepare_sampled_read_quality_frame_short_read_bulk/1000000.0/bulk/", day, ".tsv")
    )
    plt_frame_qc <- rbind(
      plt_frame_qc,
      data.frame(
        tech = "Illumina (Bulk)",
        sample = day,
        snv_rate = tmp_bulk$snv_rate,
        indel_rate = tmp_bulk$indel_rate
      )
    )

    tmp_sc <- vroom::vroom(
      paste0("results/qc/qc_prepare_sampled_read_quality_frame_short_read_single_cell/1000000.0/single_cell/", day, ".tsv")
    )
    plt_frame_qc <- rbind(
      plt_frame_qc,
      data.frame(
        tech = "Illumina (Single-cell)",
        sample = day,
        snv_rate = tmp_sc$snv_rate,
        indel_rate = tmp_sc$indel_rate
      )
    )
  }

  plt_frame_qc <- plt_frame_qc %>%
    mutate(
      tech = case_when(
        tech == "ont (bulk)" ~ "ONT (Bulk)",
        tech == "ont (single_cell)" ~ "ONT (SC)",
        tech == "pb (bulk)" ~ "PB (Bulk)",
        tech == "pb (single_cell)" ~ "PB (SC)",
        tech == "Illumina (Bulk)" ~ "Illmn (Bulk)",
        tech == "Illumina (Single-cell)" ~ "Illmn (SC)"
      )
    )

  tech_colors <- c(
    "Illmn (Bulk)" = "#f8a21fff",
    "ONT (Bulk)" = "#01799bff",
    "PB (Bulk)" = "#e21b92ff",
    "Illmn (SC)" = "#f8a21f80",
    "ONT (SC)" = "#01799b80",
    "PB (SC)" = "#e21b9280"
  )

  panel_f <- plt_frame_qc %>%
    pivot_longer(cols = c(snv_rate, indel_rate), names_to = "error_type", values_to = "value") %>%
    mutate(
      error_type = recode(
        error_type,
        `snv_rate` = "SNV Rate",
        `indel_rate` = "Indel Rate"
      )
    ) %>%
    ggplot(aes(x = value, y = sample, color = tech)) +
    geom_violin(scale = "width") +
    facet_wrap(~ error_type) +
    coord_cartesian(xlim = c(0, 0.05)) +
    theme_big_simple() +
    labs(y = "", x = "Relative Rate", color = "") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.text = element_text(face = "bold")
    ) +
    scale_color_manual(values = tech_colors)

  ggsave(output_pdf, panel_f, dpi = 600, width = 12, height = 8)
  ggsave(output_svg, panel_f, dpi = 600, width = 12, height = 8)
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
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_01/panel_F.pdf"
  output_svg <- "results/figure_01/panel_F.svg"
}

status <- plot_figure_01_f(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

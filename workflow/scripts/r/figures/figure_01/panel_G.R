extract_and_normalize_coverage <- function(file_path) {
  lines <- readLines(file_path)

  if (length(lines) < 2) {
    stop(paste("Input file", file_path, "is empty or missing the data line."), call. = FALSE)
  }

  data_line <- lines[2]
  parts <- strsplit(data_line, "\t")[[1]]
  raw_values <- as.numeric(parts[2:length(parts)])

  if (any(is.na(raw_values))) {
    stop(paste("Error: Non-numeric values found in data line of", file_path), call. = FALSE)
  }

  max_val <- max(raw_values, na.rm = TRUE)
  if (max_val == 0) {
    return(raw_values)
  }

  raw_values / max_val
}

load_bias_data <- function(tech, platform, path_prefix, sample_names, bins) {
  lapply(bins, function(bin) {
    df <- data.frame(lapply(
      paste0(path_prefix, tech, "/1000000.0/", platform, "/", sample_names, "/", bin, "/three_prime_bias.geneBodyCoverage.txt"),
      extract_and_normalize_coverage
    ))
    colnames(df) <- sample_names
    df
  })
}

plot_figure_01_g <- function(output_pdf, output_svg) {
  sample_names <- c("E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3")
  bins <- c("lt1kb", "1-2kb", "2-3kb", "3-6kb", "gt6kb")

  pb_bulk <- load_bias_data("pb", "bulk", "results/qc/qc_calculate_three_prime_bias_binned/", sample_names, bins)
  ont_bulk <- load_bias_data("ont", "bulk", "results/qc/qc_calculate_three_prime_bias_binned/", sample_names, bins)
  pb_sc <- load_bias_data("pb", "single_cell", "results/qc/qc_calculate_three_prime_bias_binned/", sample_names, bins)
  ont_sc <- load_bias_data("ont", "single_cell", "results/qc/qc_calculate_three_prime_bias_binned/", sample_names, bins)
  ill_bulk <- load_bias_data("", "bulk", "results/qc/qc_calculate_three_prime_bias_binned_short_reads/", sample_names, bins)
  ill_sc <- load_bias_data("", "single_cell", "results/qc/qc_calculate_three_prime_bias_binned_short_reads/", sample_names, bins)

  plt_frame_bias <- rbind(
    cbind(pb_bulk[[1]], type = "PB (Bulk)", length = "<= 1 kb"),
    cbind(pb_bulk[[2]], type = "PB (Bulk)", length = "(1 kb, 2 kb]"),
    cbind(pb_bulk[[3]], type = "PB (Bulk)", length = "(2 kb, 3 kb]"),
    cbind(pb_bulk[[4]], type = "PB (Bulk)", length = "(3 kb, 6 kb]"),
    cbind(pb_bulk[[5]], type = "PB (Bulk)", length = "> 6 kb"),
    cbind(ont_bulk[[1]], type = "ONT (Bulk)", length = "<= 1 kb"),
    cbind(ont_bulk[[2]], type = "ONT (Bulk)", length = "(1 kb, 2 kb]"),
    cbind(ont_bulk[[3]], type = "ONT (Bulk)", length = "(2 kb, 3 kb]"),
    cbind(ont_bulk[[4]], type = "ONT (Bulk)", length = "(3 kb, 6 kb]"),
    cbind(ont_bulk[[5]], type = "ONT (Bulk)", length = "> 6 kb"),
    cbind(pb_sc[[1]], type = "PB (SC)", length = "<= 1 kb"),
    cbind(pb_sc[[2]], type = "PB (SC)", length = "(1 kb, 2 kb]"),
    cbind(pb_sc[[3]], type = "PB (SC)", length = "(2 kb, 3 kb]"),
    cbind(pb_sc[[4]], type = "PB (SC)", length = "(3 kb, 6 kb]"),
    cbind(pb_sc[[5]], type = "PB (SC)", length = "> 6 kb"),
    cbind(ont_sc[[1]], type = "ONT (SC)", length = "<= 1 kb"),
    cbind(ont_sc[[2]], type = "ONT (SC)", length = "(1 kb, 2 kb]"),
    cbind(ont_sc[[3]], type = "ONT (SC)", length = "(2 kb, 3 kb]"),
    cbind(ont_sc[[4]], type = "ONT (SC)", length = "(3 kb, 6 kb]"),
    cbind(ont_sc[[5]], type = "ONT (SC)", length = "> 6 kb"),
    cbind(ill_bulk[[1]], type = "Illmn (Bulk)", length = "<= 1 kb"),
    cbind(ill_bulk[[2]], type = "Illmn (Bulk)", length = "(1 kb, 2 kb]"),
    cbind(ill_bulk[[3]], type = "Illmn (Bulk)", length = "(2 kb, 3 kb]"),
    cbind(ill_bulk[[4]], type = "Illmn (Bulk)", length = "(3 kb, 6 kb]"),
    cbind(ill_bulk[[5]], type = "Illmn (Bulk)", length = "> 6 kb"),
    cbind(ill_sc[[1]], type = "Illmn (SC)", length = "<= 1 kb"),
    cbind(ill_sc[[2]], type = "Illmn (SC)", length = "(1 kb, 2 kb]"),
    cbind(ill_sc[[3]], type = "Illmn (SC)", length = "(2 kb, 3 kb]"),
    cbind(ill_sc[[4]], type = "Illmn (SC)", length = "(3 kb, 6 kb]"),
    cbind(ill_sc[[5]], type = "Illmn (SC)", length = "> 6 kb")
  )

  plt_frame_bias$x <- 1:100
  plt_frame_bias <- plt_frame_bias %>%
    pivot_longer(cols = starts_with(c("E3", "iso")), names_to = "day", values_to = "val")

  panel_g <- plt_frame_bias %>%
    filter(length != "> 6 kb") %>%
    ggplot(aes(x = as.numeric(x), y = as.numeric(val), col = day)) +
    geom_path(linewidth = 2) +
    facet_grid(vars(type), vars(length)) +
    theme_big_simple() +
    labs(x = "Transcript percentile (5' to 3')", y = "Normalized coverage", color = "")

  ggsave(output_pdf, panel_g, dpi = 600, width = 18, height = 10)
  ggsave(output_svg, panel_g, dpi = 600, width = 18, height = 10)
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
  library(ggpubfigs)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_01/panel_G.pdf"
  output_svg <- "results/figure_01/panel_G.svg"
}

status <- plot_figure_01_g(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

plot_figure_s02_a <- function(output_pdf, output_svg) {

  oarfish_kinnex_sc <- vroom::vroom("results/format_quantify_subsampled/gencode/pb/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  oarfish_ont_sc <- vroom::vroom("results/format_quantify_subsampled/gencode/ont/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  detected_pb_sc <- oarfish_kinnex_sc$transcript_id[apply(edgeR::cpm(oarfish_kinnex_sc[, 3:5]), 1, function(x) all(x >= 1))]
  detected_ont_sc <- oarfish_ont_sc$transcript_id[apply(edgeR::cpm(oarfish_ont_sc[, 3:5]), 1, function(x) all(x >= 1))]
  coverage_thresholds <- c(0.25, 0.5, 0.75, 0.9)
  k_values <- c(1, 3, 5)
  oarfish_cleanup_df <- bind_rows(summarize_sc_cleanup(detected_pb_sc, "PB (SC)"), summarize_sc_cleanup(detected_ont_sc, "ONT (SC)")) %>% mutate(facet_label = paste0("Coverage >= ", coverage_filter, ", k >= ", k_filter), facet_label = factor(facet_label, levels = as.vector(outer(coverage_thresholds, k_values, function(cov, k) paste0("Coverage >= ", cov, ", k >= ", k)))), Detection = factor(Detection, levels = c("Shared with >=1 bulk", "Unique to SC")))
  sc_filter_cleanup_plot <- ggplot(oarfish_cleanup_df, aes(x = Technology, y = Count, fill = Detection)) + geom_col(position = "stack") + facet_wrap(~facet_label, nrow = 2) + theme_big_simple() + labs(x = "", y = "Detected transcripts", fill = "Detection class") + scale_fill_manual(values = ggpubfigs::friendly_pals$bright_seven[1:2])

  if (!exists("sc_filter_cleanup_plot")) {
    stop("Panel object not found: sc_filter_cleanup_plot", call. = FALSE)
  }

  panel_plot <- get("sc_filter_cleanup_plot")
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
  library(ggplot2)
  library(Biostrings)
  library(tidyr)
  library(rtracklayer)
  library(dplyr)
  library(fishpond)
  library(SummarizedExperiment)
  library(ggpubfigs)
  library(tximeta)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_S02/panel_A.pdf"
  output_svg <- "results/figure_S02/panel_A.svg"
}

plot_figure_s02_a(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


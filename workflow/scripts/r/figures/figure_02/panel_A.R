plot_figure_02_a <- function(output_pdf, output_svg) {

  kinnex <- vroom::vroom("results/format/run_isosceles_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  ont <- vroom::vroom("results/format/run_isosceles_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  illumina_tpm_path <- "results/format/run_salmon_illumina_tpm/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv"
  illumina_tpm <- vroom::vroom(illumina_tpm_path) %>% arrange(desc(transcript_id))
  kinnex_pseudobulk <- vroom::vroom("results/format_quantify_subsampled/gencode/pb/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  ont_pseudobulk <- vroom::vroom("results/format_quantify_subsampled/gencode/ont/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  detected_transcripts <- list(illumina_tpm$transcript_id[apply(as.matrix(illumina_tpm[, 3:5]), 1, function(x) all(x >= 1))], ont$transcript_id[apply(edgeR::cpm(ont[, 3:5]), 1, function(x) all(x >= 1))], ont_pseudobulk$transcript_id[apply(edgeR::cpm(ont_pseudobulk[, 3:5]), 1, function(x) all(x >= 1))], kinnex$transcript_id[apply(edgeR::cpm(kinnex[, 3:5]), 1, function(x) all(x >= 1))], kinnex_pseudobulk$transcript_id[apply(edgeR::cpm(kinnex_pseudobulk[, 3:5]), 1, function(x) all(x >= 1))])
  all_detected_transcripts <- Reduce(union, detected_transcripts)
  upset_plt_frame <- data.frame(`Illumina (Bulk)` = as.integer(all_detected_transcripts %in% detected_transcripts[[1]]), `ONT (Bulk)` = as.integer(all_detected_transcripts %in% detected_transcripts[[2]]), `ONT (SC)` = as.integer(all_detected_transcripts %in% detected_transcripts[[3]]), `PB (Bulk)` = as.integer(all_detected_transcripts %in% detected_transcripts[[4]]), `PB (SC)` = as.integer(all_detected_transcripts %in% detected_transcripts[[5]]), check.names = FALSE)
  upset_data <- upset_plt_frame %>% mutate(transcript_id = row_number()) %>% pivot_longer(-transcript_id, names_to = "Tech", values_to = "Detected") %>% filter(Detected == 1) %>% group_by(transcript_id) %>% summarize(Technologies = list(Tech)) %>% mutate(Technologies = map(Technologies, sort)) %>% filter(map_int(Technologies, length) > 0)
  upset <- upset_data %>% ggplot(aes(x = Technologies)) + geom_bar(fill = "grey30") + scale_x_upset(n_intersections = 15, reverse = FALSE) + geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) + theme_big_simple() + labs(title = "Transcript detection across technologies (Illumina: TPM >= 1; ONT/PB: CPM >= 1; each in all three E3 replicates).", y = "Intersection size") + scale_y_continuous(limits = c(0, 25000), expand = c(0, 0))

  if (!exists("upset")) {
    stop("Panel object not found: upset", call. = FALSE)
  }

  panel_plot <- get("upset")
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
  library(tidyverse)
  library(SummarizedExperiment)
  library(ggpubfigs)
  library(ggupset)
  library(fishpond)
  library(tximeta)
  library(tximport)
  library(rtracklayer)
  library(Biostrings)
  library(dplyr)
  library(cowplot)
  library(scales)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_02/panel_A.pdf"
  output_svg <- "results/figure_02/panel_A.svg"
}

plot_figure_02_a(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


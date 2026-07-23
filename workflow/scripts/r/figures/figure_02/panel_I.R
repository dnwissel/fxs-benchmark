plot_figure_02_i <- function(output_pdf, output_svg) {

  config <- list(log_expr_pseudo = 0.1)
  transcriptome <- readDNAStringSet("results/gencode_transcriptome.fa")
  transcriptome_gtf <- import("results/gencode.v45.primary_assembly.annotation.gtf.gz")
  kinnex <- vroom::vroom("results/format/run_isosceles_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  ont <- vroom::vroom("results/format/run_isosceles_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  illumina_tpm_path <- "results/format/run_salmon_illumina_tpm/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv"
  illumina_tpm <- vroom::vroom(illumina_tpm_path) %>% arrange(desc(transcript_id))
  kinnex_pseudobulk <- vroom::vroom("results/format_quantify_subsampled/gencode/pb/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  ont_pseudobulk <- vroom::vroom("results/format_quantify_subsampled/gencode/ont/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv") %>% arrange(desc(transcript_id))
  detected_transcripts <- list(illumina_tpm$transcript_id[apply(as.matrix(illumina_tpm[, 3:5]), 1, function(x) all(x >= 1))], ont$transcript_id[apply(edgeR::cpm(ont[, 3:5]), 1, function(x) all(x >= 1))], ont_pseudobulk$transcript_id[apply(edgeR::cpm(ont_pseudobulk[, 3:5]), 1, function(x) all(x >= 1))], kinnex$transcript_id[apply(edgeR::cpm(kinnex[, 3:5]), 1, function(x) all(x >= 1))], kinnex_pseudobulk$transcript_id[apply(edgeR::cpm(kinnex_pseudobulk[, 3:5]), 1, function(x) all(x >= 1))])
  all_detected_transcripts <- Reduce(union, detected_transcripts)
  detail_frame <- data.frame(transcript_id = all_detected_transcripts)
  detail_frame <- detail_frame %>% left_join(data.frame(transcript_id = sapply(strsplit(names(transcriptome), "\\ "), function(x) x[[1]]), gc_content = letterFrequency(transcriptome, letters = "GC", as.prob = TRUE)[, 1], length = nchar(transcriptome)))
  detail_frame$pb_bulk_detection <- ifelse(detail_frame$transcript_id %in% detected_transcripts[[1]] & detail_frame$transcript_id %in% detected_transcripts[[2]] & detail_frame$transcript_id %in% detected_transcripts[[3]] & detail_frame$transcript_id %in% detected_transcripts[[5]] & !detail_frame$transcript_id %in% detected_transcripts[[4]], "All other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[4]], "PB (Bulk)", "Other"))
  detail_frame$ont_bulk_detection <- ifelse(detail_frame$transcript_id %in% detected_transcripts[[1]] & detail_frame$transcript_id %in% detected_transcripts[[3]] & detail_frame$transcript_id %in% detected_transcripts[[4]] & detail_frame$transcript_id %in% detected_transcripts[[5]] & !detail_frame$transcript_id %in% detected_transcripts[[2]], "All other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[2]], "ONT (Bulk)", "Other"))
  detail_frame <- data.frame(transcript_id = all_detected_transcripts)
  detail_frame <- detail_frame %>% left_join(data.frame(transcript_id = sapply(strsplit(names(transcriptome), "\\ "), function(x) x[[1]]), gc_content = letterFrequency(transcriptome, letters = "GC", as.prob = TRUE)[, 1], length = nchar(transcriptome)))
  gtf <- transcriptome_gtf
  gene_bounds <- gtf[gtf$type == "gene"] %>% as.data.frame() %>% dplyr::select(gene_id, gene_start = start, gene_end = end, gene_strand = strand) %>% distinct()
  tx_bounds <- gtf[gtf$type == "transcript"] %>% as.data.frame() %>% dplyr::select(transcript_id, gene_id, tx_start = start, tx_end = end)
  detail_frame <- detail_frame %>% inner_join(tx_bounds, by = "transcript_id") %>% inner_join(gene_bounds, by = "gene_id") %>% mutate(gene_length = gene_end - gene_start, rel_pos_3p = if_else(gene_strand == "+", (tx_end - gene_start)/gene_length, (gene_end - tx_start)/gene_length)) %>% mutate(rel_start = if_else(gene_strand == "+", (tx_start - gene_start)/gene_length, (gene_end - tx_end)/gene_length), rel_end = if_else(gene_strand == "+", (tx_end - gene_start)/gene_length, (gene_end - tx_start)/gene_length))
  detail_frame$detection <- ifelse(detail_frame$transcript_id %in% detected_transcripts[[3]] & detail_frame$transcript_id %in% detected_transcripts[[5]] & (detail_frame$transcript_id %in% detected_transcripts[[2]] & detail_frame$transcript_id %in% detected_transcripts[[4]]), "Bulk and SC", ifelse(detail_frame$transcript_id %in% detected_transcripts[[3]] & detail_frame$transcript_id %in% detected_transcripts[[5]], "SC", "Other"))
  all_detected_transcripts <- Reduce(union, detected_transcripts)
  detail_frame <- data.frame(transcript_id = all_detected_transcripts)
  detail_frame <- detail_frame %>% left_join(data.frame(transcript_id = sapply(strsplit(names(transcriptome), "\\ "), function(x) x[[1]]), gc_content = letterFrequency(transcriptome, letters = "GC", as.prob = TRUE)[, 1], length = nchar(transcriptome))) %>% left_join(transcriptome_gtf %>% as.data.frame %>% filter(type == "transcript") %>% dplyr::select(seqnames, start, end, gene_id, transcript_id, strand) %>% left_join(transcriptome_gtf %>% as.data.frame %>% filter(type == "gene") %>% dplyr::select(seqnames, 
      start, end, gene_id, width), by = "gene_id") %>% mutate(relative_five_prime_position = ifelse(strand == "+", (start.x - start.y)/(width), (end.y - end.x)/(width)), relative_three_prime_position = ifelse(strand == "+", (end.x - start.y)/(width), (end.y - start.x)/(width))))
  mask_bulk <- Reduce(intersect, lapply(list(kinnex, ont), function(counts) {
      which(apply(edgeR::cpm(counts[, 3:5]), 1, function(x) all(x >= 1)))
  }))
  mask_bulk <- Reduce(intersect, list(mask_bulk, which(apply(as.matrix(illumina_tpm[, 3:5]), 1, function(x) all(x >= 1)))))
  mask_pseudobulk <- Reduce(intersect, lapply(list(kinnex_pseudobulk, ont_pseudobulk), function(counts) {
      which(apply(edgeR::cpm(counts[, 3:5]), 1, function(x) all(x >= 1)))
  }))
  mask <- mask_bulk[mask_bulk %in% mask_pseudobulk]
  master_wide <- data.frame(transcript_id = illumina_tpm[mask, ]$transcript_id, ill_bulk = log2(illumina_tpm[mask, ]$`E3-1` + config$log_expr_pseudo), pb_bulk = log2(edgeR::cpm(kinnex[mask, ]$`E3-1`, log = FALSE) + config$log_expr_pseudo), ont_bulk = log2(edgeR::cpm(ont[mask, ]$`E3-1`, log = FALSE) + config$log_expr_pseudo), pb_pb = log2(edgeR::cpm(kinnex_pseudobulk[mask, ]$`E3-1`, log = FALSE) + config$log_expr_pseudo), ont_pb = log2(edgeR::cpm(ont_pseudobulk[mask, ]$`E3-1`, log = FALSE) + config$log_expr_pseudo)) %>% 
      mutate(set_ont_bulk = calc_set_multi(ont_bulk, list(ill_bulk, pb_bulk, pb_pb, ont_pb)), set_ont_pb = calc_set_multi(ont_pb, list(ill_bulk, pb_bulk, ont_bulk)), set_pb_bulk = calc_set_multi(pb_bulk, list(ill_bulk, ont_bulk, pb_pb, ont_pb)), set_pb_pb = calc_set_multi(pb_pb, list(ill_bulk, pb_bulk, ont_bulk)))
  master_prepared <- master_wide %>% dplyr::rename(val_ontBulk = ont_bulk, val_pbBulk = pb_bulk, val_ontPB = ont_pb, val_pbPB = pb_pb, set_ontBulk = set_ont_bulk, set_pbBulk = set_pb_bulk, set_ontPB = set_ont_pb, set_pbPB = set_pb_pb)
  facet_frame <- master_prepared %>% pivot_longer(cols = starts_with("val_") | starts_with("set_"), names_to = c(".value", "tech_key"), names_sep = "_") %>% mutate(facet_label = case_match(tech_key, "ontBulk" ~ "ONT Bulk", "pbBulk" ~ "PB Bulk", "ontPB" ~ "ONT SC", "pbPB" ~ "PB SC"))
  plot_data <- facet_frame %>% left_join(detail_frame %>% dplyr::select(transcript_id, length), by = "transcript_id")
  quant_difference_gene_body_length <- plot_data %>% filter(facet_label %in% c("ONT Bulk", "PB Bulk")) %>% ggplot(aes(x = length, fill = set)) + geom_histogram(alpha = 0.5, position = "identity", bins = 50) + scale_x_log10() + scale_fill_manual(values = c(`Query higher` = "#e41a1c", `Query lower` = "#377eb8", Similar = "grey80")) + geom_vline(xintercept = 1750, lty = 2) + facet_wrap(~facet_label, nrow = 1) + theme_big_simple() + labs(x = "Transcript length (nts)", y = "Count", fill = "Bias") + labs(title = "Under- and overquantification due to length biases in PB (Bulk) and ONT (Bulk)")

  if (!exists("quant_difference_gene_body_length")) {
    stop("Panel object not found: quant_difference_gene_body_length", call. = FALSE)
  }

  panel_plot <- get("quant_difference_gene_body_length")
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
  output_pdf <- "results/figure_02/panel_I.pdf"
  output_svg <- "results/figure_02/panel_I.svg"
}

plot_figure_02_i(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


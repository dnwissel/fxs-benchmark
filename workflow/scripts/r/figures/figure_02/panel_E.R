plot_figure_02_e <- function(output_pdf, output_svg) {

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
  backup_frame <- detail_frame
  backup_frame$coverage <- backup_frame$relative_three_prime_position - backup_frame$relative_five_prime_position
  detail_frame <- data.frame(transcript_id = all_detected_transcripts)
  detail_frame <- detail_frame %>% left_join(data.frame(transcript_id = sapply(strsplit(names(transcriptome), "\\ "), function(x) x[[1]]), gc_content = letterFrequency(transcriptome, letters = "GC", as.prob = TRUE)[, 1], length = nchar(transcriptome)))
  coldata_ill <- data.frame(files = c("results/run_salmon_illumina/1_15000000.0/gencode/E3-1/illumina/quant.sf", "results/run_salmon_illumina/1_15000000.0/gencode/E3-2/illumina/quant.sf", "results/run_salmon_illumina/1_15000000.0/gencode/E3-3/illumina/quant.sf"), names = c("E3-1", "E3-2", "E3-3"))
  se <- tximeta(coldata_ill)
  inf_ill <- fishpond::computeInfRV(se)
  inf_var_frame <- data.frame(transcript_id = rownames(rowData(inf_ill)), mean_inf_rv_bulk = unname(rowData(inf_ill)[, 1]), mean_inf_rv_pseudobulk = NA, type = "Illmn")
  detail_frame <- detail_frame %>% left_join(inf_var_frame %>% filter(type == "Illmn") %>% dplyr::select(transcript_id, mean_inf_rv_bulk)) %>% left_join(backup_frame %>% dplyr::select(transcript_id, coverage))
  detail_frame$pb_bulk_detection <- ifelse((detail_frame$transcript_id %in% detected_transcripts[[1]] | detail_frame$transcript_id %in% detected_transcripts[[5]] | detail_frame$transcript_id %in% detected_transcripts[[2]] | detail_frame$transcript_id %in% detected_transcripts[[3]]) & detail_frame$transcript_id %in% detected_transcripts[[4]], "PB (Bulk) and other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[4]], "PB (Bulk)", "Other"))
  detail_frame$ont_bulk_detection <- ifelse((detail_frame$transcript_id %in% detected_transcripts[[1]] | detail_frame$transcript_id %in% detected_transcripts[[4]] | detail_frame$transcript_id %in% detected_transcripts[[3]] | detail_frame$transcript_id %in% detected_transcripts[[5]]) & detail_frame$transcript_id %in% detected_transcripts[[2]], "ONT (Bulk) and other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[2]], "ONT (Bulk)", "Other"))
  detail_frame$ont_pseudobulk_detection <- ifelse((detail_frame$transcript_id %in% detected_transcripts[[1]] | detail_frame$transcript_id %in% detected_transcripts[[4]] | detail_frame$transcript_id %in% detected_transcripts[[2]] | detail_frame$transcript_id %in% detected_transcripts[[5]]) & detail_frame$transcript_id %in% detected_transcripts[[3]], "ONT (SC) and other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[3]], "ONT (SC)", "Other"))
  detail_frame$pb_pseudobulk_detection <- ifelse((detail_frame$transcript_id %in% detected_transcripts[[1]] | detail_frame$transcript_id %in% detected_transcripts[[4]] | detail_frame$transcript_id %in% detected_transcripts[[2]] | detail_frame$transcript_id %in% detected_transcripts[[3]]) & detail_frame$transcript_id %in% detected_transcripts[[5]], "PB (SC) and other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[5]], "PB (SC)", "Other"))
  ont_bulk_unique <- detail_frame %>% filter(ont_bulk_detection != "Other") %>% ggplot(aes(x = ont_bulk_detection, y = length)) + geom_violin(linewidth = 1, scale = "width") + theme_big_simple() + labs(y = "Length (nts)", x = "Detected by", title = "ONT (Bulk) unique detections due to\n focus n on shorter transcripts") + scale_y_log10()

  if (!exists("ont_bulk_unique")) {
    stop("Panel object not found: ont_bulk_unique", call. = FALSE)
  }

  panel_plot <- get("ont_bulk_unique")
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
  output_pdf <- "results/figure_02/panel_E.pdf"
  output_svg <- "results/figure_02/panel_E.svg"
}

plot_figure_02_e(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


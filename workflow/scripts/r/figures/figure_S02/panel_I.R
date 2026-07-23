plot_figure_s02_i <- function(output_pdf, output_svg) {

  config <- list(log_expr_pseudo = 0.1)
  transcriptome <- readDNAStringSet("results/gencode_transcriptome.fa")
  transcriptome_gtf <- import("results/gencode.v45.primary_assembly.annotation.gtf.gz")
  mapping <- data.frame(transcript_id = transcriptome_gtf$transcript_id, gene_id = transcriptome_gtf$gene_id) %>% distinct %>% drop_na
  transcript_annotation_frame <- data.frame(transcript_id = sapply(strsplit(names(transcriptome), "\\ "), function(x) x[[1]]), transcript_length = nchar(as.character(transcriptome))) %>% left_join(mapping)
  coldata_ill <- data.frame(files = c("results/run_oarfish_lr_bootstrapped/1_15000000.0/gencode/E3-1/pb/E3-1.quant", "results/run_oarfish_lr_bootstrapped/1_15000000.0/gencode/E3-2/pb/E3-2.quant", "results/run_oarfish_lr_bootstrapped/1_15000000.0/gencode/E3-3/pb/E3-3.quant"), names = c("E3-1", "E3-2", "E3-3"))
  se <- tximeta(coldata_ill, type = "oarfish", skipMeta = TRUE)
  inf_rv_bulk <- fishpond::computeInfRV(se)
  transcript_annotation_frame <- transcript_annotation_frame %>% left_join(data.frame(transcript_id = names(rowData(inf_rv_bulk)$meanInfRV), transcript_illumina_mean_inf_rv = unname(rowData(inf_rv_bulk)$meanInfRV))) %>% left_join(transcriptome_gtf %>% as.data.frame %>% filter(type == "transcript") %>% dplyr::select(seqnames, start, end, gene_id, transcript_id, strand) %>% left_join(transcriptome_gtf %>% as.data.frame %>% filter(type == "gene") %>% dplyr::select(seqnames, start, end, gene_id, width), 
      by = "gene_id") %>% mutate(relative_five_prime_position = ifelse(strand == "+", (start.x - start.y)/(width), (end.y - end.x)/(width)), relative_three_prime_position = ifelse(strand == "+", (end.x - start.y)/(width), (end.y - start.x)/(width)))) %>% drop_na() %>% mutate(coverage = relative_three_prime_position - relative_five_prime_position) %>% dplyr::select(transcript_id, transcript_length, transcript_illumina_mean_inf_rv, coverage, gene_id)
  for (method in c("isosceles", "oarfish", "bambu", "kallisto")) {
      if (method != "kallisto") {
          kinnex <- vroom::vroom(paste0("results/format/run_", method, "_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")) %>% arrange(desc(transcript_id))
          ont <- vroom::vroom(paste0("results/format/run_", method, "_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")) %>% arrange(desc(transcript_id))
      }
      else {
          kinnex <- vroom::vroom(paste0("results/format/run_", method, "_long_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")) %>% arrange(desc(transcript_id))
          ont <- vroom::vroom(paste0("results/format/run_", method, "_long_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")) %>% arrange(desc(transcript_id))
      }
      illumina <- vroom::vroom(paste0("results/format/run_salmon_illumina_corrected/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv")) %>% arrange(desc(transcript_id))
      illumina_tpm_path <- "results/format/run_salmon_illumina_tpm/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv"
      illumina_tpm <- vroom::vroom(illumina_tpm_path) %>% arrange(desc(transcript_id))
      kinnex_pseudobulk <- vroom::vroom(paste0("results/format_quantify_subsampled/gencode/pb/1_60000000.0/", method, "/pseudobulk/transcript_counts_formatted.tsv")) %>% arrange(desc(transcript_id))
      ont_pseudobulk <- vroom::vroom(paste0("results/format_quantify_subsampled/gencode/ont/1_60000000.0/", method, "/pseudobulk/transcript_counts_formatted.tsv")) %>% arrange(desc(transcript_id))
      detected_transcripts <- list(illumina_tpm$transcript_id[apply(as.matrix(illumina_tpm[, 3:5]), 1, function(x) all(x >= 1))], ont$transcript_id[apply(edgeR::cpm(ont[, 3:5]), 1, function(x) all(x >= 1))], ont_pseudobulk$transcript_id[apply(edgeR::cpm(ont_pseudobulk[, 3:5]), 1, function(x) all(x >= 1))], kinnex$transcript_id[apply(edgeR::cpm(kinnex[, 3:5]), 1, function(x) all(x >= 1))], kinnex_pseudobulk$transcript_id[apply(edgeR::cpm(kinnex_pseudobulk[, 3:5]), 1, function(x) all(x >= 1))])
      all_detected_transcripts <- Reduce(union, detected_transcripts)
      detail_frame <- data.frame(transcript_id = all_detected_transcripts)
      detail_frame$pb_bulk_detection <- ifelse(detail_frame$transcript_id %in% detected_transcripts[[1]] & detail_frame$transcript_id %in% detected_transcripts[[2]] & detail_frame$transcript_id %in% detected_transcripts[[3]] & detail_frame$transcript_id %in% detected_transcripts[[5]] & !detail_frame$transcript_id %in% detected_transcripts[[4]], "All other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[4]], "PB (Bulk)", "Other"))
      detail_frame$ont_bulk_detection <- ifelse(detail_frame$transcript_id %in% detected_transcripts[[1]] & detail_frame$transcript_id %in% detected_transcripts[[3]] & detail_frame$transcript_id %in% detected_transcripts[[4]] & detail_frame$transcript_id %in% detected_transcripts[[5]] & !detail_frame$transcript_id %in% detected_transcripts[[2]], "All other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[2]], "ONT (Bulk)", "Other"))
      detail_frame$pseudobulk_detection <- ifelse(detail_frame$transcript_id %in% detected_transcripts[[3]] & detail_frame$transcript_id %in% detected_transcripts[[5]] & (detail_frame$transcript_id %in% detected_transcripts[[2]] & detail_frame$transcript_id %in% detected_transcripts[[4]]), "Bulk and SC", ifelse(detail_frame$transcript_id %in% detected_transcripts[[3]] & detail_frame$transcript_id %in% detected_transcripts[[5]], "SC", "Other"))
      detail_frame$pb_exclusive_bulk_detection <- ifelse((detail_frame$transcript_id %in% detected_transcripts[[1]] | detail_frame$transcript_id %in% detected_transcripts[[5]] | detail_frame$transcript_id %in% detected_transcripts[[2]] | detail_frame$transcript_id %in% detected_transcripts[[3]]) & detail_frame$transcript_id %in% detected_transcripts[[4]], "PB (Bulk) and other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[4]], "PB (Bulk)", "Other"))
      detail_frame$ont_exclusive_bulk_detection <- ifelse((detail_frame$transcript_id %in% detected_transcripts[[1]] | detail_frame$transcript_id %in% detected_transcripts[[4]] | detail_frame$transcript_id %in% detected_transcripts[[3]] | detail_frame$transcript_id %in% detected_transcripts[[5]]) & detail_frame$transcript_id %in% detected_transcripts[[2]], "ONT (Bulk) and other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[2]], "ONT (Bulk)", "Other"))
      detail_frame$ont_exclusive_pseudobulk_detection <- ifelse((detail_frame$transcript_id %in% detected_transcripts[[1]] | detail_frame$transcript_id %in% detected_transcripts[[4]] | detail_frame$transcript_id %in% detected_transcripts[[2]] | detail_frame$transcript_id %in% detected_transcripts[[5]]) & detail_frame$transcript_id %in% detected_transcripts[[3]], "ONT (SC) and other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[3]], "ONT (SC)", "Other"))
      detail_frame$pb_exclusive_pseudobulk_detection <- ifelse((detail_frame$transcript_id %in% detected_transcripts[[1]] | detail_frame$transcript_id %in% detected_transcripts[[4]] | detail_frame$transcript_id %in% detected_transcripts[[2]] | detail_frame$transcript_id %in% detected_transcripts[[3]]) & detail_frame$transcript_id %in% detected_transcripts[[5]], "PB (SC) and other techs", ifelse(detail_frame$transcript_id %in% detected_transcripts[[5]], "PB (SC)", "Other"))
      detail_frame$quant_method <- method
      combined_frame <- rbind(combined_frame, detail_frame)
  }
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
  plot_data <- facet_frame %>% left_join(transcript_annotation_frame)
  coverage_difference_gene_body_length <- ggplot(plot_data, aes(x = coverage, fill = set)) + geom_histogram(alpha = 0.5, position = "identity", bins = 50) + scale_fill_manual(values = c(`Query higher` = "#e41a1c", `Query lower` = "#377eb8", Similar = "grey80")) + facet_wrap(~facet_label, nrow = 1) + theme_big_simple() + scale_y_continuous(trans = "pseudo_log") + labs(x = "Gene coverage", y = "Count", fill = "Bias") + labs(title = "")

  if (!exists("coverage_difference_gene_body_length")) {
    stop("Panel object not found: coverage_difference_gene_body_length", call. = FALSE)
  }

  panel_plot <- get("coverage_difference_gene_body_length")
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
  output_pdf <- "results/figure_S02/panel_I.pdf"
  output_svg <- "results/figure_S02/panel_I.svg"
}

plot_figure_s02_i(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


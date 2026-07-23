plot_figure_s02_f <- function(output_pdf, output_svg) {

  transcriptome <- readDNAStringSet("results/gencode_transcriptome.fa")
  transcriptome_gtf <- import("results/gencode.v45.primary_assembly.annotation.gtf.gz")
  mapping <- data.frame(transcript_id = transcriptome_gtf$transcript_id, gene_id = transcriptome_gtf$gene_id) %>% distinct %>% drop_na
  transcript_annotation_frame <- data.frame(transcript_id = sapply(strsplit(names(transcriptome), "\\ "), function(x) x[[1]]), transcript_length = nchar(as.character(transcriptome))) %>% left_join(mapping)
  coldata_ill <- data.frame(files = c("results/run_oarfish_lr_bootstrapped/1_15000000.0/gencode/E3-1/pb/E3-1.quant", "results/run_oarfish_lr_bootstrapped/1_15000000.0/gencode/E3-2/pb/E3-2.quant", "results/run_oarfish_lr_bootstrapped/1_15000000.0/gencode/E3-3/pb/E3-3.quant"), names = c("E3-1", "E3-2", "E3-3"))
  se <- tximeta(coldata_ill, type = "oarfish", skipMeta = TRUE)
  inf_rv_bulk <- fishpond::computeInfRV(se)
  transcript_annotation_frame <- transcript_annotation_frame %>% left_join(data.frame(transcript_id = names(rowData(inf_rv_bulk)$meanInfRV), transcript_illumina_mean_inf_rv = unname(rowData(inf_rv_bulk)$meanInfRV))) %>% left_join(transcriptome_gtf %>% as.data.frame %>% filter(type == "transcript") %>% dplyr::select(seqnames, start, end, gene_id, transcript_id, strand) %>% left_join(transcriptome_gtf %>% as.data.frame %>% filter(type == "gene") %>% dplyr::select(seqnames, start, end, gene_id, width), 
      by = "gene_id") %>% mutate(relative_five_prime_position = ifelse(strand == "+", (start.x - start.y)/(width), (end.y - end.x)/(width)), relative_three_prime_position = ifelse(strand == "+", (end.x - start.y)/(width), (end.y - start.x)/(width)))) %>% drop_na() %>% mutate(coverage = relative_three_prime_position - relative_five_prime_position) %>% dplyr::select(transcript_id, transcript_length, transcript_illumina_mean_inf_rv, coverage, gene_id)
  combined_frame <- data.frame()
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
  combined_frame$facet_label <- case_when(combined_frame$quant_method == "bambu" ~ "Bambu", combined_frame$quant_method == "oarfish" ~ "Oarfish", combined_frame$quant_method == "kallisto" ~ "Kallisto", combined_frame$quant_method == "isosceles" ~ "Isosceles")
  combined_frame$facet_label <- factor(combined_frame$facet_label, levels = c("Bambu", "Isosceles", "Kallisto", "Oarfish"))
  ont_bulk_unique <- combined_frame %>% left_join(transcript_annotation_frame) %>% filter(ont_exclusive_bulk_detection != "Other") %>% ggplot(aes(x = ont_exclusive_bulk_detection, y = transcript_length)) + geom_violin(scale = "width") + theme_big_simple() + facet_wrap(~facet_label, nrow = 1) + labs(y = "Length (nts)", x = "Detected by", title = "ONT (Bulk) unique detections due to focus n on shorter transcripts") + scale_y_log10()

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
  output_pdf <- "results/figure_S02/panel_F.pdf"
  output_svg <- "results/figure_S02/panel_F.svg"
}

plot_figure_s02_f(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


plot_figure_04_c <- function(output_pdf, output_svg) {

  config <- list(log_expr_pseudo = 0.1)
  quant_methods <- c("bambu", "isosceles", "kallisto", "oarfish")
  master_mapping <- gencode_master %>% dplyr::select(transcript_id, gene_id) %>% distinct %>% drop_na
  quant_methods <- c("bambu", "isosceles", "kallisto", "oarfish")
  ont_tx_cts_list <- vector(mode = "list", length = length(quant_methods))
  for (m in quant_methods) {
      ont_path <- paste0("results/format_quantify_subsampled/gencode/ont/1_60000000.0/", m, "/pseudobulk/transcript_counts_formatted.tsv")
      ont_tx_cts <- vroom::vroom(ont_path) %>% left_join(master_mapping) %>% dplyr::select(gene_id, transcript_id, "E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3") %>% tibble::column_to_rownames("transcript_id") %>% dplyr::select(-gene_id)
      ont_tx_cts[is.na(ont_tx_cts)] <- 0
      ont_tx_cts_list[[m]] <- ont_tx_cts
  }
  keep_ont_tx <- Reduce(union, lapply(ont_tx_cts_list, function(counts) {
      names(which(apply(edgeR::cpm(counts[, 1:3], normalized.lib.sizes = TRUE), 1, function(x) all(x >= 1))))
  }))
  pb_tx_cts_list <- vector(mode = "list", length = length(quant_methods))
  for (m in quant_methods) {
      pb_path <- paste0("results/format_quantify_subsampled/gencode/pb/1_60000000.0/", m, "/pseudobulk/transcript_counts_formatted.tsv")
      pb_tx_cts <- vroom::vroom(pb_path) %>% left_join(master_mapping) %>% dplyr::select(gene_id, transcript_id, "E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3") %>% tibble::column_to_rownames("transcript_id") %>% dplyr::select(-gene_id)
      pb_tx_cts[is.na(pb_tx_cts)] <- 0
      pb_tx_cts_list[[m]] <- pb_tx_cts
  }
  keep_pb_tx <- Reduce(union, lapply(pb_tx_cts_list, function(counts) {
      names(which(apply(edgeR::cpm(counts[, 1:3], normalized.lib.sizes = TRUE), 1, function(x) all(x >= 1))))
  }))
  ill_tx_tpm_ref_path <- "results/format/run_salmon_illumina_tpm/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv"
  ill_tx_tpm_ref <- read.csv(ill_tx_tpm_ref_path, sep = "\t", header = TRUE) %>% tibble::column_to_rownames("transcript_id") %>% dplyr::select(-gene_id)
  keep_tx_bulk <- Reduce(union, lapply(list(ill_tx_tpm_ref), function(counts) {
      names(which(apply(counts[, 1:3], 1, function(x) all(x >= 1))))
  }))
  keep_tx <- Reduce(intersect, list(keep_tx_bulk, keep_ont_tx, keep_pb_tx))
  ont_tx_cts_list <- lapply(ont_tx_cts_list, function(m) m[keep_tx, ])
  pb_tx_cts_list <- lapply(pb_tx_cts_list, function(m) m[keep_tx, ])
  ill_tx_tpm_ref <- ill_tx_tpm_ref[keep_tx, ]
  ont_tx_cts_list <- lapply(names(ont_tx_cts_list), function(m) {
      colnames(ont_tx_cts_list[[m]]) <- paste("ONT", m, colnames(ont_tx_cts_list[[m]]), sep = "_")
      return(ont_tx_cts_list[[m]])
  })
  pb_tx_cts_list <- lapply(names(pb_tx_cts_list), function(m) {
      colnames(pb_tx_cts_list[[m]]) <- paste("PB", m, colnames(pb_tx_cts_list[[m]]), sep = "_")
      return(pb_tx_cts_list[[m]])
  })
  ont_tx_cts <- do.call(cbind, ont_tx_cts_list)
  ont_tx_cts <- ont_tx_cts %>% dplyr::select(which(grepl("E3-1", colnames(ont_tx_cts))))
  pb_tx_cts <- do.call(cbind, pb_tx_cts_list)
  pb_tx_cts <- pb_tx_cts %>% dplyr::select(which(grepl("E3.1", colnames(pb_tx_cts))))
  tx_tpm_ref <- ill_tx_tpm_ref %>% dplyr::select(which(grepl("E3.1", colnames(ill_tx_tpm_ref))))
  tx_cts <- cbind(ont_tx_cts, pb_tx_cts)
  tx_cpm <- log2(edgeR::cpm(tx_cts, log = FALSE, normalized.lib.sizes = TRUE) + config$log_expr_pseudo)
  tx_abs_quant_df <- tx_cpm %>% as.data.frame() %>% tibble::rownames_to_column("transcript_id") %>% tidyr::pivot_longer(cols = -transcript_id, names_to = "Sample", values_to = "CPM") %>% dplyr::mutate(Technology = sapply(strsplit(Sample, "_"), function(x) x[[1]]), Method = sapply(strsplit(Sample, "_"), function(x) x[[2]]), Replicate = sapply(Sample, function(x) substr(x, nchar(x), nchar(x)))) %>% dplyr::select(-Sample) %>% tidyr::pivot_wider(id_cols = c(transcript_id, Method, Replicate), names_from = Technology, 
      values_from = CPM) %>% arrange(transcript_id)
  tx_log_tpm_ref <- log2(tx_tpm_ref + config$log_expr_pseudo)
  tx_abs_quant_ref_df <- tx_log_tpm_ref %>% as.data.frame() %>% tibble::rownames_to_column("transcript_id") %>% tidyr::pivot_longer(cols = -transcript_id, names_to = "Sample", values_to = "CPM") %>% dplyr::mutate(Technology = sapply(strsplit(Sample, "_"), function(x) paste(x[1:2], collapse = "_")), Method = sapply(strsplit(Sample, "_"), function(x) x[[3]]), Replicate = sapply(Sample, function(x) substr(x, nchar(x), nchar(x)))) %>% dplyr::select(-Sample) %>% tidyr::pivot_wider(id_cols = c(transcript_id, 
      Method, Replicate), names_from = Technology, values_from = CPM) %>% group_by(transcript_id) %>% dplyr::slice(rep(1:n(), times = length(quant_methods))) %>% ungroup() %>% dplyr::select(transcript_id, Illumina_bulk) %>% arrange(transcript_id)
  tx_abs_quant_df <- cbind(tx_abs_quant_df, tx_abs_quant_ref_df)
  tx_abs_quant_df <- tx_abs_quant_df[, !duplicated(names(tx_abs_quant_df))]
  tx_abs_quant_df$Method <- stringr::str_to_title(tx_abs_quant_df$Method)
  tx_abs_quant_df$Method <- factor(tx_abs_quant_df$Method, levels = c("Bambu", "Isosceles", "Kallisto", "Oarfish"))
  panel_e3 <- ggplot(tx_abs_quant_df, aes(x = Illumina_bulk, y = ONT)) + geom_density_2d_filled(contour_var = "ndensity", alpha = 0.9) + stat_cor(aes(label = after_stat(r.label), group = 1), show.legend = FALSE, method = "pearson", size = 6, label.x.npc = 0.025, label.y.npc = 1, geom = "label", fill = "white") + scale_fill_viridis_d(option = "A", name = "Normalized density") + labs(x = "Illumina bulk transcript abundance (logTPM)", y = "ONT transcript counts (logCPM)") + ggpubfigs::theme_big_simple() + 
      theme(legend.position = "bottom") + facet_wrap(~Method, nrow = 1)

  if (!exists("panel_e3")) {
    stop("Panel object not found: panel_e3", call. = FALSE)
  }

  panel_plot <- get("panel_e3")
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
  output_pdf <- "results/figure_04/panel_C.pdf"
  output_svg <- "results/figure_04/panel_C.svg"
}

plot_figure_04_c(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


plot_figure_03_d <- function(output_pdf, output_svg) {

  config <- list(log_expr_pseudo = 0.1)
  path_master <- "results/format/run_bambu_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv"
  path_master <- "results/format/run_miniquant_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv"
  gencode_master <- read.csv(path_master, sep = "\t", header = TRUE)
  path_master <- "results/format/run_bambu_lr/1_15000000.0/ont/gencode/gene_counts_formatted.tsv"
  gencode_master_gene <- read.csv(path_master, sep = "\t", header = TRUE)
  master_genes <- gencode_master_gene$gene_id
  master_transcripts <- gencode_master$transcript_id
  quant_methods <- c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish")
  cols <- carto_pal(7, "Safe")
  quant_methods <- c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish")
  ill_gene_tpm_path <- "results/format/run_salmon_illumina_tpm/1_15000000.0/illumina/gencode/gene_counts_formatted.tsv"
  ill_gene_tpm <- read.csv(ill_gene_tpm_path, sep = "\t", header = TRUE) %>% tibble::column_to_rownames("gene_id")
  keep_ill_gene <- names(which(apply(ill_gene_tpm[, 1:3], 1, function(x) all(x >= 1))))
  ont_gene_cts_list <- vector(mode = "list", length = length(quant_methods))
  for (m in quant_methods) {
      if (m == "Kallisto") {
          ont_path <- paste0("results/format/run_", tolower(m), "_long_lr//1_15000000.0/ont/gencode/gene_counts_formatted.tsv")
      }
      else {
          ont_path <- paste0("results/format/run_", tolower(m), "_lr//1_15000000.0/ont/gencode/gene_counts_formatted.tsv")
      }
      ont_gene_cts <- read.csv(ont_path, sep = "\t", header = TRUE)
      if (!all(master_genes %in% ont_gene_cts$gene_id)) {
          concat_frame <- gencode_master_gene[which(!(master_genes %in% ont_gene_cts$gene_id)), ]
          concat_frame[, -(1)] <- 0
          gencode_cts <- rbind(ont_gene_cts, concat_frame)
          ont_gene_cts <- gencode_cts[match(gencode_master_gene$gene_id, gencode_cts$gene_id), ] %>% tibble::remove_rownames()
      }
      ont_gene_cts[is.na(ont_gene_cts)] <- 0
      ont_gene_cts_list[[m]] <- ont_gene_cts %>% tibble::column_to_rownames("gene_id")
  }
  keep_ont_gene <- Reduce(union, lapply(ont_gene_cts_list, function(counts) {
      names(which(apply(edgeR::cpm(counts[, 1:3], normalized.lib.sizes = TRUE), 1, function(x) all(x >= 1))))
  }))
  pb_gene_cts_list <- vector(mode = "list", length = length(quant_methods))
  for (m in quant_methods) {
      print(m)
      if (m == "Kallisto") {
          pb_path <- paste0("results/format/run_", tolower(m), "_long_lr//1_15000000.0/pb/gencode/gene_counts_formatted.tsv")
      }
      else {
          pb_path <- paste0("results/format/run_", tolower(m), "_lr//1_15000000.0/pb/gencode/gene_counts_formatted.tsv")
      }
      pb_gene_cts <- read.csv(pb_path, sep = "\t", header = TRUE)
      if (!all(master_genes %in% pb_gene_cts$gene_id)) {
          concat_frame <- gencode_master_gene[which(!(master_genes %in% pb_gene_cts$gene_id)), ]
          concat_frame[, -(1)] <- 0
          gencode_cts <- rbind(pb_gene_cts, concat_frame)
          pb_gene_cts <- gencode_cts[match(gencode_master_gene$gene_id, gencode_cts$gene_id), ] %>% tibble::remove_rownames()
      }
      pb_gene_cts[is.na(pb_gene_cts)] <- 0
      pb_gene_cts_list[[m]] <- pb_gene_cts %>% tibble::column_to_rownames("gene_id")
  }
  keep_pb_gene <- Reduce(union, lapply(pb_gene_cts_list, function(counts) {
      names(which(apply(edgeR::cpm(counts[, 1:3], normalized.lib.sizes = TRUE), 1, function(x) all(x >= 1))))
  }))
  keep_gene <- intersect(keep_ill_gene, intersect(keep_ont_gene, keep_pb_gene))
  ill_gene_tpm <- ill_gene_tpm[keep_gene, ]
  ont_gene_cts_list <- lapply(ont_gene_cts_list, function(m) m[keep_gene, ])
  pb_gene_cts_list <- lapply(pb_gene_cts_list, function(m) m[keep_gene, ])
  ont_gene_cts_list <- lapply(names(ont_gene_cts_list), function(m) {
      colnames(ont_gene_cts_list[[m]]) <- paste("ONT", m, colnames(ont_gene_cts_list[[m]]), sep = "_")
      return(ont_gene_cts_list[[m]])
  })
  pb_gene_cts_list <- lapply(names(pb_gene_cts_list), function(m) {
      colnames(pb_gene_cts_list[[m]]) <- paste("PB", m, colnames(pb_gene_cts_list[[m]]), sep = "_")
      return(pb_gene_cts_list[[m]])
  })
  ill_gene_tpm <- ill_gene_tpm %>% dplyr::select(which(grepl("E3", colnames(ill_gene_tpm))))
  ill_gene_tpm <- do.call(cbind, replicate(length(quant_methods), ill_gene_tpm, simplify = FALSE))
  ont_gene_cts <- do.call(cbind, ont_gene_cts_list)
  ont_gene_cts <- ont_gene_cts %>% dplyr::select(which(grepl("E3", colnames(ont_gene_cts))))
  pb_gene_cts <- do.call(cbind, pb_gene_cts_list)
  pb_gene_cts <- pb_gene_cts %>% dplyr::select(which(grepl("E3", colnames(pb_gene_cts))))
  gene_log_expr <- cbind(log2(ill_gene_tpm + config$log_expr_pseudo), log2(edgeR::cpm(cbind(ont_gene_cts, pb_gene_cts), log = FALSE, normalized.lib.sizes = TRUE) + config$log_expr_pseudo))
  gene_abs_quant_df <- gene_log_expr %>% as.data.frame() %>% tibble::rownames_to_column("gene_id") %>% tidyr::pivot_longer(cols = -gene_id, names_to = "Sample", values_to = "log_expr") %>% dplyr::mutate(Technology = sapply(strsplit(Sample, "_"), function(x) x[[1]]), Method = sapply(strsplit(Sample, "_"), function(x) x[[2]]), Replicate = sapply(Sample, function(x) substr(x, nchar(x), nchar(x)))) %>% dplyr::select(-Sample) %>% tidyr::pivot_wider(id_cols = c(gene_id, Method, Replicate), names_from = Technology, 
      values_from = log_expr)
  gene_abs_quant_df$Method <- factor(gene_abs_quant_df$Method, levels = c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish"))
  coldata_ill <- data.frame(files = c("results/run_salmon_illumina/1_15000000.0/gencode/E3-1/illumina/quant.sf", "results/run_salmon_illumina/1_15000000.0/gencode/E3-2/illumina/quant.sf", "results/run_salmon_illumina/1_15000000.0/gencode/E3-3/illumina/quant.sf"), names = c("E3-1", "E3-2", "E3-3"))
  se <- tximeta(coldata_ill)
  se <- se[grepl("^ENST", rownames(se))]
  illumina_se <- fishpond::computeInfRV(se)
  keep_tx_inf_var <- names(illumina_se@elementMetadata@listData$meanInfRV[illumina_se@elementMetadata@listData$meanInfRV < 1])
  ill_tx_tpm_path <- "results/format/run_salmon_illumina_tpm/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv"
  ill_tx_tpm <- read.csv(ill_tx_tpm_path, sep = "\t", header = TRUE) %>% tibble::column_to_rownames("transcript_id") %>% dplyr::select(-gene_id)
  ill_tx_tpm[is.na(ill_tx_tpm)] <- 0
  keep_ill_tx_cts <- names(which(apply(ill_tx_tpm[, 1:3], 1, function(x) all(x >= 1))))
  ont_tx_cts_list <- vector(mode = "list", length = length(quant_methods))
  for (m in quant_methods) {
      if (m == "Kallisto") {
          ont_path <- paste0("results/format/run_", tolower(m), "_long_lr//1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")
      }
      else {
          ont_path <- paste0("results/format/run_", tolower(m), "_lr//1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")
      }
      ont_tx_cts <- read.csv(ont_path, header = TRUE, sep = "\t")
      if (!all(master_transcripts %in% ont_tx_cts$transcript_id)) {
          concat_frame <- gencode_master[which(!(master_transcripts %in% ont_tx_cts$transcript_id)), ]
          concat_frame[, -(1:2)] <- 0
          ont_tx_cts <- rbind(ont_tx_cts, concat_frame)
          ont_tx_cts <- ont_tx_cts[match(gencode_master$transcript_id, ont_tx_cts$transcript_id), ] %>% tibble::remove_rownames()
      }
      ont_tx_cts[is.na(ont_tx_cts)] <- 0
      ont_tx_cts_list[[m]] <- ont_tx_cts %>% tibble::column_to_rownames("transcript_id") %>% dplyr::select(-gene_id)
  }
  keep_ont_tx_cts <- Reduce(union, lapply(ont_tx_cts_list, function(counts) {
      names(which(apply(edgeR::cpm(counts[, 1:3], normalized.lib.sizes = TRUE), 1, function(x) all(x >= 1))))
  }))
  pb_tx_cts_list <- vector(mode = "list", length = length(quant_methods))
  for (m in quant_methods) {
      if (m == "Kallisto") {
          pb_path <- paste0("results/format/run_", tolower(m), "_long_lr//1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
      }
      else {
          pb_path <- paste0("results/format/run_", tolower(m), "_lr//1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
      }
      pb_tx_cts <- read.csv(pb_path, sep = "\t", header = TRUE)
      if (!all(master_transcripts %in% pb_tx_cts$transcript_id)) {
          concat_frame <- gencode_master[which(!(master_transcripts %in% pb_tx_cts$transcript_id)), ]
          concat_frame[, -(1:2)] <- 0
          pb_tx_cts <- rbind(pb_tx_cts, concat_frame)
          pb_tx_cts <- pb_tx_cts[match(gencode_master$transcript_id, pb_tx_cts$transcript_id), ] %>% tibble::remove_rownames()
      }
      pb_tx_cts[is.na(pb_tx_cts)] <- 0
      pb_tx_cts_list[[m]] <- pb_tx_cts %>% tibble::column_to_rownames("transcript_id") %>% dplyr::select(-gene_id)
  }
  keep_pb_tx_cts <- Reduce(union, lapply(pb_tx_cts_list, function(counts) {
      names(which(apply(edgeR::cpm(counts[, 1:3], normalized.lib.sizes = TRUE), 1, function(x) all(x >= 1))))
  }))
  keep_tx <- intersect(intersect(keep_ill_tx_cts, keep_tx_inf_var), intersect(keep_ont_tx_cts, keep_pb_tx_cts))
  ill_tx_tpm <- ill_tx_tpm[keep_tx, ]
  ont_tx_cts_list <- lapply(ont_tx_cts_list, function(m) m[keep_tx, ])
  pb_tx_cts_list <- lapply(pb_tx_cts_list, function(m) m[keep_tx, ])
  ont_tx_cts_list <- lapply(names(ont_tx_cts_list), function(m) {
      colnames(ont_tx_cts_list[[m]]) <- paste("ONT", m, colnames(ont_tx_cts_list[[m]]), sep = "_")
      return(ont_tx_cts_list[[m]])
  })
  pb_tx_cts_list <- lapply(names(pb_tx_cts_list), function(m) {
      colnames(pb_tx_cts_list[[m]]) <- paste("PB", m, colnames(pb_tx_cts_list[[m]]), sep = "_")
      return(pb_tx_cts_list[[m]])
  })
  ill_tx_tpm <- ill_tx_tpm %>% dplyr::select(which(grepl("E3", colnames(ill_tx_tpm))))
  ill_tx_tpm <- do.call(cbind, replicate(length(quant_methods), ill_tx_tpm, simplify = FALSE))
  ont_tx_cts <- do.call(cbind, ont_tx_cts_list)
  ont_tx_cts <- ont_tx_cts %>% dplyr::select(which(grepl("E3", colnames(ont_tx_cts))))
  pb_tx_cts <- do.call(cbind, pb_tx_cts_list)
  pb_tx_cts <- pb_tx_cts %>% dplyr::select(which(grepl("E3", colnames(pb_tx_cts))))
  tx_log_expr <- cbind(log2(ill_tx_tpm + config$log_expr_pseudo), log2(edgeR::cpm(cbind(ont_tx_cts, pb_tx_cts), log = FALSE, normalized.lib.sizes = TRUE) + config$log_expr_pseudo))
  tx_abs_quant_df <- tx_log_expr %>% as.data.frame() %>% tibble::rownames_to_column("transcript_id") %>% tidyr::pivot_longer(cols = -transcript_id, names_to = "Sample", values_to = "log_expr") %>% dplyr::mutate(Technology = sapply(strsplit(Sample, "_"), function(x) x[[1]]), Method = sapply(strsplit(Sample, "_"), function(x) x[[2]]), Replicate = sapply(Sample, function(x) substr(x, nchar(x), nchar(x)))) %>% dplyr::select(-Sample) %>% tidyr::pivot_wider(id_cols = c(transcript_id, Method, Replicate), 
      names_from = Technology, values_from = log_expr)
  gene_abs_quant_df$Method <- factor(gene_abs_quant_df$Method, levels = c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish"))
  num_quant_methods_bulk <- 7
  cols <- carto_pal(num_quant_methods_bulk, "Safe")
  plt_frame_transcript <- tx_abs_quant_df
  plt_frame_gene <- gene_abs_quant_df
  calculate_correlations <- function(df) {
      df %>% group_by(Method) %>% summarise(ONT_vs_Illumina = cor(ONT, Illumina, method = "pearson"), PB_vs_Illumina = cor(PB, Illumina, method = "pearson"), .groups = "drop")
  }
  data_list <- list(Gene = plt_frame_gene, Transcript = plt_frame_transcript)
  combined_corr <- bind_rows(lapply(data_list, calculate_correlations), .id = "Level")
  corr_long <- combined_corr %>% pivot_longer(cols = c(ONT_vs_Illumina, PB_vs_Illumina), names_to = "Platform", values_to = "Correlation") %>% mutate(Platform = str_remove(Platform, "_vs_Illumina"), Heatmap_Column = str_c(Method, " (", Platform, ")")) %>% dplyr::select(Level, Method, Platform, Heatmap_Column, Correlation)
  heatmap_matrix_df <- corr_long %>% pivot_wider(names_from = Level, values_from = Correlation)
  heatmap_matrix <- heatmap_matrix_df %>% column_to_rownames("Heatmap_Column")
  methods_order <- c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish")
  col_order_ont <- paste0(methods_order, " (ONT)")
  col_order_pb <- paste0(methods_order, " (PB)")
  final_col_order <- c(col_order_ont, col_order_pb)
  heatmap_matrix <- heatmap_matrix[final_col_order, ]
  heatmap_matrix <- heatmap_matrix[, c("Gene", "Transcript")]
  annotation_df <- data.frame(Platform = str_extract(rownames(heatmap_matrix), "ONT|PB"), Method = str_remove(rownames(heatmap_matrix), " \\(ONT\\)| \\(PB\\)"))
  my_annotation_colors <- list(Platform = c(ONT = "#01799bff", PB = "#e21b92ff"), Method = c(Bambu = cols[1], Isoquant = cols[2], Isosceles = cols[3], Kallisto = cols[4], Miniquant = cols[5], Oarfish = cols[6]))
  my_colors <- colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)
  pheatmap_plot <- pheatmap(heatmap_matrix, main = "", display_numbers = TRUE, number_format = "%.2f", fontsize_number = 20, cluster_rows = FALSE, cluster_cols = FALSE, color = my_colors, border_color = "white", show_colnames = TRUE, show_rownames = FALSE, angle_col = 0, fontsize = 24, fontsize_col = 24, annotation_row = annotation_df, annotation_colors = my_annotation_colors, filename = NA, silent = TRUE)
  ggplot_heatmap <- as.ggplot(pheatmap_plot$gtable)

  if (!exists("ggplot_heatmap")) {
    stop("Panel object not found: ggplot_heatmap", call. = FALSE)
  }

  panel_plot <- get("ggplot_heatmap")
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
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(ggpubfigs)
  library(ggpubr)
  library(glue)
  library(rcartocolor)
  library(tibble)
  library(tidyr)
  library(tximeta)
  library(DRIMSeq)
  library(stageR)
  library(stringr)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplotify)
  library(readxl)
  library(rtracklayer)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_03/panel_D.pdf"
  output_svg <- "results/figure_03/panel_D.svg"
}

plot_figure_03_d(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


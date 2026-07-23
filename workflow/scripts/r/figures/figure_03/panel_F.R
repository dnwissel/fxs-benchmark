plot_figure_03_f <- function(output_pdf, output_svg) {

  ill_gencode_cts <- vroom::vroom("results/format/run_salmon_illumina_corrected/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv")
  ill_gencode_tpm_path <- "results/format/run_salmon_illumina_tpm/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv"
  ill_gencode_tpm <- vroom::vroom(ill_gencode_tpm_path)
  path_master <- "results/format/run_bambu_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv"
  path_master <- "results/format/run_miniquant_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv"
  gencode_master <- read.csv(path_master, sep = "\t", header = TRUE)
  master_transcripts <- gencode_master$transcript_id
  quant_methods <- c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish")
  ont_gencode_cts_list <- lapply(tolower(quant_methods), function(m) {
      if (m == "kallisto") {
          path <- paste0("results/format/run_", m, "_long_lr//1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")
      }
      else {
          path <- paste0("results/format/run_", m, "_lr//1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")
      }
      gencode_cts <- read.csv(path, sep = "\t", header = TRUE)
      if (!all(master_transcripts %in% gencode_cts$transcript_id)) {
          concat_frame <- gencode_master[which(!(master_transcripts %in% gencode_cts$transcript_id)), ]
          concat_frame[, -(1:2)] <- 0
          gencode_cts <- rbind(gencode_cts, concat_frame)
          gencode_cts <- gencode_cts[match(gencode_master$transcript_id, gencode_cts$transcript_id), ] %>% tibble::remove_rownames()
      }
      gencode_cts[is.na(gencode_cts)] <- 0
      return(gencode_cts)
  })
  pb_gencode_cts_list <- lapply(tolower(quant_methods), function(m) {
      if (m == "kallisto") {
          path <- paste0("results/format/run_", m, "_long_lr//1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
      }
      else {
          path <- paste0("results/format/run_", m, "_lr//1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
      }
      gencode_cts <- read.csv(path, sep = "\t", header = TRUE)
      if (!all(master_transcripts %in% gencode_cts$transcript_id)) {
          concat_frame <- gencode_master[which(!(master_transcripts %in% gencode_cts$transcript_id)), ]
          concat_frame[, -(1:2)] <- 0
          gencode_cts <- rbind(gencode_cts, concat_frame)
          gencode_cts <- gencode_cts[match(gencode_master$transcript_id, gencode_cts$transcript_id), ] %>% tibble::remove_rownames()
      }
      gencode_cts[is.na(gencode_cts)] <- 0
      return(gencode_cts)
  })
  cols <- carto_pal(7, "Safe")
  design_vec <- c(0, 0, 0, 1, 1, 1)
  cts_ill <- ill_gencode_cts %>% dplyr::select(-c(gene_id)) %>% tibble::column_to_rownames("transcript_id")
  keep_ill_dte <- ill_gencode_tpm$transcript_id[apply(as.matrix(ill_gencode_tpm[, 3:5]), 1, function(x) all(x >= 1))]
  ill_stats <- get_edgeR_table(cts_ill, design_vec, keep_ids = keep_ill_dte) %>% dplyr::mutate(Tech = "Illumina", Method = "Salmon")
  ont_stats_list <- lapply(names(ont_gencode_cts_list), function(m) {
      cts <- ont_gencode_cts_list[[m]] %>% dplyr::select(-c(gene_id)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("transcript_id")
      get_edgeR_table(cts, design_vec) %>% dplyr::mutate(Tech = "ONT", Method = m)
  })
  ont_all_stats <- bind_rows(ont_stats_list)
  pb_stats_list <- lapply(names(pb_gencode_cts_list), function(m) {
      cts <- pb_gencode_cts_list[[m]] %>% dplyr::select(-c(gene_id)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("transcript_id")
      get_edgeR_table(cts, design_vec) %>% dplyr::mutate(Tech = "PB", Method = m)
  })
  pb_all_stats <- bind_rows(pb_stats_list)
  all_dte_data <- bind_rows(ill_stats, ont_all_stats, pb_all_stats)
  fdr_cutoff <- 0.01
  logfc_cutoff <- 1
  all_dte_data <- all_dte_data %>% dplyr::mutate(is_sig = FDR < fdr_cutoff & abs(logFC) >= logfc_cutoff)
  truth_table <- all_dte_data %>% group_by(transcript_id) %>% summarise(data = list(data.frame(Tech, Method, logFC, is_sig))) %>% ungroup() %>% mutate(classification = purrr::map_chr(data, function(df) {
      detected_techs <- unique(df$Tech)
      n_detected <- length(detected_techs)
      if (n_detected <= 1) 
          return("Unknown")
      best_class <- if (n_detected == 3) 
          "Silver-Negative"
      else "Unknown"
      drivers <- df[df$is_sig, ]
      if (nrow(drivers) == 0) 
          return(best_class)
      check_consistency <- function(target_data, driver_lfc, diff_cutoff) {
          if (nrow(target_data) == 0) 
              return(FALSE)
          valid_support <- target_data %>% dplyr::filter(abs(logFC) > 1, sign(logFC) == sign(driver_lfc), abs(logFC - driver_lfc) < diff_cutoff)
          return(nrow(valid_support) > 0)
      }
      for (r in 1:nrow(drivers)) {
          d_tech <- drivers$Tech[r]
          d_lfc <- drivers$logFC[r]
          d_meth <- drivers$Method[r]
          all_techs <- c("Illumina", "ONT", "PB")
          other_techs <- setdiff(all_techs, d_tech)
          t1_data <- df[df$Tech == other_techs[1] & df$Method != d_meth, ]
          t2_data <- df[df$Tech == other_techs[2] & df$Method != d_meth, ]
          t1_sig <- if (nrow(t1_data) > 0) 
              any(t1_data$is_sig)
          else FALSE
          t2_sig <- if (nrow(t2_data) > 0) 
              any(t2_data$is_sig)
          else FALSE
          if (t1_sig || t2_sig) {
              best_class <- "High Confidence"
              break
          }
          t1_med <- check_consistency(t1_data, d_lfc, diff_cutoff = 2)
          t2_med <- check_consistency(t2_data, d_lfc, diff_cutoff = 2)
          if (t1_med || t2_med) {
              if (best_class != "High Confidence") 
                  best_class <- "Medium Confidence"
              next
          }
          t1_low <- check_consistency(t1_data, d_lfc, diff_cutoff = 4)
          t2_low <- check_consistency(t2_data, d_lfc, diff_cutoff = 4)
          if (t1_low || t2_low) {
              if (best_class %in% c("Silver-Negative", "Unknown")) 
                  best_class <- "Low Confidence"
          }
      }
      return(best_class)
  }))
  evaluated_calls <- all_dte_data %>% dplyr::filter(is_sig) %>% dplyr::left_join(truth_table %>% dplyr::select(transcript_id, classification), by = "transcript_id")
  fdr_stats <- evaluated_calls %>% filter(Tech != "Illumina") %>% dplyr::group_by(Tech, Method) %>% dplyr::summarise(TP_Total = sum(grepl("Confidence", classification)), FP_Total = sum(classification == "Silver-Negative"), .groups = "drop") %>% dplyr::mutate(Total_Calls = TP_Total + FP_Total, FDR = ifelse(Total_Calls > 0, FP_Total/Total_Calls, 0), FDR_Label = sprintf("%.2f", FDR)) %>% dplyr::rename(Tool = Method)
  tp_plot_df <- evaluated_calls %>% dplyr::filter(grepl("Confidence", classification)) %>% dplyr::group_by(Tech, Method, classification) %>% dplyr::summarise(Count = n(), .groups = "drop") %>% dplyr::rename(Tool = Method, Confidence = classification)
  tool_levels <- c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish")
  tp_plot_df$Tool <- factor(stringr::str_to_title(tp_plot_df$Tool), levels = tool_levels)
  fdr_stats$Tool <- factor(stringr::str_to_title(fdr_stats$Tool), levels = tool_levels)
  tp_plot_df$Confidence <- factor(tp_plot_df$Confidence, levels = c("Low Confidence", "Medium Confidence", "High Confidence"))
  tp_plot_df <- tp_plot_df %>% filter(Tech != "Illumina")
  panel_g <- ggplot(tp_plot_df, aes(x = Tool, y = Count, fill = Tool)) + geom_bar(aes(alpha = Confidence), stat = "identity", position = "stack") + facet_wrap(~Tech) + geom_text(data = fdr_stats, aes(x = Tool, y = TP_Total, label = FDR_Label), vjust = -0.5, size = 5, inherit.aes = FALSE) + theme_big_simple() + scale_fill_manual(values = cols) + scale_alpha_manual(values = c(`High Confidence` = 1, `Medium Confidence` = 0.7, `Low Confidence` = 0.2), breaks = c("High Confidence", "Medium Confidence", 
      "Low Confidence")) + labs(y = "Replicable DTE calls", x = "", fill = "", alpha = "", subtitle = "Numbers above bars represent irreproducibility rates") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom", legend.box = "horizontal") + guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2, nrow = 2)) + coord_cartesian(ylim = c(0, 45))

  if (!exists("panel_g")) {
    stop("Panel object not found: panel_g", call. = FALSE)
  }

  panel_plot <- get("panel_g")
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
  output_pdf <- "results/figure_03/panel_F.pdf"
  output_svg <- "results/figure_03/panel_F.svg"
}

plot_figure_03_f(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}


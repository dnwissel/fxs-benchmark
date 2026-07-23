plot_figure_s04_a <- function(output_pdf, output_svg) {

  num_quant_methods_bulk <- 7
  quant_methods_bulk <- c("bambu", "isoquant", "isosceles", "kallisto", "miniquant", "oarfish")
  ont_bulk_gencode_cts_list <- lapply(quant_methods_bulk, function(m) {
      if (m == "kallisto") {
          path <- paste0("results/format/run_", m, "_long_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")
      }
      else {
          path <- paste0("results/format/run_", m, "_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")
      }
      gencode_cts <- read.csv(path, sep = "\t", header = TRUE)
      gencode_cts[is.na(gencode_cts)] <- 0
      return(gencode_cts)
  })
  pb_bulk_gencode_cts_list <- lapply(quant_methods_bulk, function(m) {
      if (m == "kallisto") {
          path <- paste0("results/format/run_", m, "_long_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
      }
      else {
          path <- paste0("results/format/run_", m, "_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
      }
      gencode_cts <- read.csv(path, sep = "\t", header = TRUE)
      gencode_cts[is.na(gencode_cts)] <- 0
      return(gencode_cts)
  })
  ill_bulk_gencode_cts <- vroom::vroom("results/format/run_salmon_illumina_corrected/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv")
  ill_bulk_gencode_tpm_path <- "results/format/run_salmon_illumina_tpm/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv"
  ill_bulk_gencode_tpm <- vroom::vroom(ill_bulk_gencode_tpm_path)
  quant_methods_pseudobulk <- c("bambu", "isosceles", "kallisto", "oarfish")
  meta <- read.csv("results/FXS_metadata.csv")
  ont_bulk_stats_list <- lapply(quant_methods_bulk, function(m) {
      meta_lib <- meta %>% filter(method == "ONT" & type == "bulk")
      cond <- factor(meta_lib$condition)
      cond <- relevel(cond, ref = "E3")
      ts <- get_dtu_table(ont_bulk_gencode_cts_list[[m]], cond) %>% dplyr::mutate(Tech = "ONT", Type = "Bulk", Method = m)
      return(ts)
  })
  ont_bulk_stats <- Reduce(rbind, ont_bulk_stats_list)
  pb_bulk_stats_list <- lapply(quant_methods_bulk, function(m) {
      meta_lib <- meta %>% filter(method == "PacBio" & type == "bulk")
      cond <- factor(meta_lib$condition)
      cond <- relevel(cond, ref = "E3")
      ts <- get_dtu_table(pb_bulk_gencode_cts_list[[m]], cond) %>% dplyr::mutate(Tech = "PB", Type = "Bulk", Method = m)
      return(ts)
  })
  pb_bulk_stats <- Reduce(rbind, pb_bulk_stats_list)
  meta_lib <- meta %>% filter(method == "Illumina" & type == "bulk")
  cond <- factor(meta_lib$condition)
  cond <- relevel(cond, ref = "E3")
  keep_ill_bulk <- ill_bulk_gencode_tpm$transcript_id[apply(as.matrix(ill_bulk_gencode_tpm[, 3:5]), 1, function(x) all(x >= 1))]
  ill_bulk_stats <- get_dtu_table(ill_bulk_gencode_cts, cond, keep_ids = keep_ill_bulk) %>% dplyr::mutate(Tech = "Illumina", Type = "Bulk", Method = "Salmon")
  ont_pseudobulk_stats_list <- lapply(quant_methods_pseudobulk, function(m) {
      meta_lib <- meta %>% filter(method == "ONT" & type == "sc")
      cond <- factor(meta_lib$condition)
      cond <- relevel(cond, ref = "E3")
      if (m == "isosceles") {
          ont_pseudobulk_gencode_cts_list[[m]] <- ont_pseudobulk_gencode_cts_list[[m]][-which(is.na(ont_pseudobulk_gencode_cts_list[[m]]$gene_id)), ]
      }
      ts <- get_dtu_table(ont_pseudobulk_gencode_cts_list[[m]], cond) %>% dplyr::mutate(Tech = "ONT", Type = "Pseudobulk", Method = m)
      return(ts)
  })
  ont_pseudobulk_stats <- Reduce(rbind, ont_pseudobulk_stats_list)
  pb_pseudobulk_stats_list <- lapply(quant_methods_pseudobulk, function(m) {
      meta_lib <- meta %>% filter(method == "PacBio" & type == "sc")
      cond <- factor(meta_lib$condition)
      cond <- relevel(cond, ref = "E3")
      if (m == "isosceles") {
          pb_pseudobulk_gencode_cts_list[[m]] <- pb_pseudobulk_gencode_cts_list[[m]][-which(is.na(pb_pseudobulk_gencode_cts_list[[m]]$gene_id)), ]
      }
      ts <- get_dtu_table(pb_pseudobulk_gencode_cts_list[[m]], cond) %>% dplyr::mutate(Tech = "PB", Type = "Pseudobulk", Method = m)
      return(ts)
  })
  pb_pseudobulk_stats <- Reduce(rbind, pb_pseudobulk_stats_list)
  all_dtu_data <- Reduce(rbind, list(ont_bulk_stats, pb_bulk_stats, ill_bulk_stats, ont_pseudobulk_stats, pb_pseudobulk_stats))
  fdr_cutoff <- 0.01
  logfc_cutoff <- 1
  all_dtu_data <- all_dtu_data %>% dplyr::mutate(is_sig = FDR < fdr_cutoff & abs(logFC) >= logfc_cutoff)
  truth_table <- all_dtu_data %>% group_by(transcript_id) %>% summarise(data = list(data.frame(Tech, Type, Method, logFC, is_sig))) %>% ungroup() %>% mutate(classification = purrr::map_chr(data, function(df) {
      detected_techs <- unique(paste(df$Tech, df$Type, sep = "_"))
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
          return(nrow(valid_support) > 1)
      }
      for (r in 1:nrow(drivers)) {
          d_tech <- drivers$Tech[r]
          d_lfc <- drivers$logFC[r]
          d_meth <- drivers$Method[r]
          d_type <- drivers$Type[r]
          all_techs <- c("Illumina", "ONT", "PB")
          all_types <- c("Bulk", "Pseudobulk")
          other_techs <- setdiff(all_techs, d_tech)
          t_data <- df[(df$Tech %in% other_techs | df$Type != d_type) & df$Method != d_meth, ]
          t_sig <- if (nrow(t_data) > 1) 
              any(t_data$is_sig)
          else FALSE
          if (t_sig) {
              best_class <- "High Confidence"
              break
          }
          t_med <- check_consistency(t_data, d_lfc, diff_cutoff = 2)
          if (t_med) {
              if (best_class != "High Confidence") 
                  best_class <- "Medium Confidence"
              next
          }
          t_low <- check_consistency(t_data, d_lfc, diff_cutoff = 4)
          if (t_low) {
              if (best_class %in% c("Silver-Negative", "Unknown")) 
                  best_class <- "Low Confidence"
          }
      }
      return(best_class)
  }))
  evaluated_calls <- all_dtu_data %>% dplyr::filter(is_sig) %>% dplyr::left_join(truth_table %>% dplyr::select(transcript_id, classification), by = "transcript_id")
  fdr_stats <- evaluated_calls %>% filter(Tech != "Illumina" & Type == "Pseudobulk") %>% dplyr::group_by(Tech, Method) %>% dplyr::summarise(TP_Total = sum(grepl("Confidence", classification)), FP_Total = sum(classification == "Silver-Negative"), .groups = "drop") %>% dplyr::mutate(Total_Calls = TP_Total + FP_Total, FDR = ifelse(Total_Calls > 0, FP_Total/Total_Calls, 0), FDR_Label = sprintf("%.2f", FDR)) %>% dplyr::rename(Tool = Method)
  tp_plot_df <- evaluated_calls %>% dplyr::filter(grepl("Confidence", classification)) %>% dplyr::group_by(Tech, Type, Method, classification) %>% dplyr::summarise(Count = n(), .groups = "drop") %>% dplyr::rename(Tool = Method, Confidence = classification)
  tool_levels <- c("Bambu", "Isosceles", "Kallisto", "Oarfish")
  tp_plot_df$Tool <- factor(stringr::str_to_title(tp_plot_df$Tool), levels = tool_levels)
  fdr_stats$Tool <- factor(stringr::str_to_title(fdr_stats$Tool), levels = tool_levels)
  tp_plot_df$Confidence <- factor(tp_plot_df$Confidence, levels = c("Low Confidence", "Medium Confidence", "High Confidence"))
  tp_plot_df <- tp_plot_df %>% filter(Tech != "Illumina" & Type == "Pseudobulk")
  cols <- carto_pal(num_quant_methods_bulk, "Safe")[c(1, 3, 4, 6)]
  dtu_plot <- ggplot(tp_plot_df, aes(x = Tool, y = Count, fill = Tool)) + geom_bar(aes(alpha = Confidence), stat = "identity", position = "stack") + facet_wrap(~Tech) + geom_text(data = fdr_stats, aes(x = Tool, y = TP_Total, label = FDR_Label), vjust = -0.5, size = 5, inherit.aes = FALSE) + theme_big_simple() + scale_fill_manual(values = cols) + scale_alpha_manual(values = c(`High Confidence` = 1, `Medium Confidence` = 0.7, `Low Confidence` = 0.2), breaks = c("High Confidence", "Medium Confidence", 
      "Low Confidence")) + labs(y = "Replicable DTU calls", x = "", fill = "", alpha = "", subtitle = "Numbers above bars represent irreproducibility rates") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom", legend.box = "horizontal") + guides(fill = guide_legend(order = 1, nrow = 2), alpha = guide_legend(order = 2, nrow = 2))
  dtu_plot <- dtu_plot + coord_cartesian(ylim = c(0, max(tp_plot_df$Count) + 2)) + scale_y_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 2))

  if (!exists("dtu_plot")) {
    stop("Panel object not found: dtu_plot", call. = FALSE)
  }

  panel_plot <- get("dtu_plot")
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
  library(ComplexUpset)
  library(rcartocolor)
  library(ggplot2)
  library(ggpubfigs)
  library(ggpubr)
  library(glue)
  library(readxl)
  library(dplyr)
  library(ggupset)
  library(scran)
  library(stageR)
  library(tidyr)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_S04/panel_A.pdf"
  output_svg <- "results/figure_S04/panel_A.svg"
}

plot_figure_s04_a(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

